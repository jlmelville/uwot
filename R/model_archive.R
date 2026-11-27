# The ustar size field holds eleven octal digits per member.
ustar_size_limit <- 2^33

find_pax_tar <- function() {
  windows <- .Platform$OS.type == "windows"
  system_tar <- if (windows) {
    file.path(Sys.getenv("SystemRoot"), "System32", "tar.exe")
  }
  paths <- c(
    Sys.getenv(c("tar", "TAR")),
    system_tar,
    Sys.which(c("tar", "gtar", "bsdtar"))
  )
  for (path in unique(paths[file.exists(paths)])) {
    version <- tryCatch(
      suppressWarnings(system2(
        path,
        "--version",
        stdout = TRUE,
        stderr = TRUE
      )),
      error = function(e) character()
    )
    gnu <- any(grepl("GNU tar", version, fixed = TRUE))
    if (gnu || any(grepl("bsdtar", version, fixed = TRUE))) {
      # GNU tar interprets Windows drive letters as remote hosts by default.
      flags <- if (windows && gnu) "--force-local"
      return(list(command = shQuote(path), flags = flags))
    }
  }
  NULL
}

windows_unicode_path <- function(...) {
  .Platform$OS.type == "windows" &&
    any(grepl("[^[:ascii:]]", enc2utf8(c(...)), perl = TRUE))
}

# On Windows, tempfile(), gzfile() and system2() translate paths through
# LC_CTYPE. Switch only for paths the caller's locale cannot represent.
with_archive_locale <- function(file, code) {
  if (
    .Platform$OS.type == "windows" &&
      !l10n_info()[["UTF-8"]] &&
      anyNA(iconv(
        enc2utf8(c(file, getwd(), tempdir())),
        from = "UTF-8",
        to = ""
      ))
  ) {
    old <- Sys.getlocale("LC_CTYPE")
    on.exit(Sys.setlocale("LC_CTYPE", old))
    suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
    if (!l10n_info()[["UTF-8"]]) {
      stop(
        "This Windows R installation cannot use UTF-8 for the archive path. ",
        "Use a path representable in ",
        old,
        ".",
        call. = FALSE
      )
    }
  }
  force(code)
}

# Windows tar programs can misread Unicode arguments. Let R open the archive
# and pass the file handle directly to tar. R's tar setting can include a
# quoted executable path followed by flags, which system2 takes separately.
windows_model_tar <- function(tar, args, stdin = "", stdout = "") {
  tar <- trimws(tar)
  program <- regmatches(tar, regexpr('^("[^"]+"|[^[:space:]]+)', tar))
  options <- substring(tar, nchar(program) + 1L)
  system2(
    sub('^"(.*)"$', "\\1", program),
    c(options, args),
    stdin = stdin,
    stdout = stdout
  )
}

# GNU tar requires an explicit decompressor when the archive is on stdin.
tar_compression_flag <- function(file) {
  magic <- paste(as.character(readBin(file, "raw", 6L)), collapse = "")
  flags <- c(
    "1f8b" = "-z",
    "1f9d" = "-z",
    "425a68" = "-j",
    "fd377a585a00" = "-J",
    "28b52ffd" = "--zstd"
  )
  unname(flags[startsWith(magic, names(flags))])
}

write_model_archive <- function(archive) {
  members <- list.files("uwot", full.names = TRUE)
  large <- which(file.size(members) >= ustar_size_limit)
  tar <- Sys.getenv("tar")
  flags <- if (is_win7()) "--force-local" else ""
  if (length(large)) {
    pax <- find_pax_tar()
    if (is.null(pax)) {
      stop(
        "Model file ",
        members[large[1]],
        " is at least 8 GiB. ",
        "Install GNU tar or bsdtar to save it.",
        call. = FALSE
      )
    }
    tar <- pax$command
    flags <- paste(c("--format=pax", pax$flags), collapse = " ")
  }
  if (windows_unicode_path(archive) && nzchar(tar) && tar != "internal") {
    return(windows_model_tar(
      tar,
      c(flags, "-cf", "-", "uwot"),
      stdout = archive
    ))
  }
  utils::tar(archive, files = "uwot", tar = tar, extra_flags = flags)
}

extract_model_archive <- function(file, exdir, verbose) {
  # uwot writes uncompressed archives. The total size is an upper bound
  # on the largest member, so a large legacy archive may still be plain ustar.
  large <- isTRUE(file.size(file) >= ustar_size_limit)
  pax <- if (large) find_pax_tar()
  tar <- if (is.null(pax)) Sys.getenv("TAR") else pax$command
  flags <- if (is.null(pax) && is_win7()) "--force-local" else pax$flags
  tryCatch(
    {
      unicode <- windows_unicode_path(file, exdir)
      if (unicode && !nzchar(tar)) {
        tar <- if (nzchar(Sys.which("tar.exe"))) "tar.exe" else "internal"
      }
      status <- if (unicode && tar != "internal") {
        wd <- setwd(exdir)
        on.exit(setwd(wd))
        args <- c(flags, tar_compression_flag(file), "-xf", "-")
        if (verbose) {
          message("untar: using ", tar, " ", paste(args, collapse = " "))
        }
        windows_model_tar(tar, args, stdin = file)
      } else {
        utils::untar(
          file,
          exdir = exdir,
          verbose = verbose,
          tar = tar,
          extras = flags
        )
      }
      if (status != 0L) {
        stop("tar exited with status ", status)
      }
    },
    error = function(e) {
      hint <- if (large && is.null(pax)) {
        " Install GNU tar or bsdtar to read members of 8 GiB or more."
      }
      stop(
        "Could not extract model archive ",
        file,
        ": ",
        conditionMessage(e),
        hint,
        call. = FALSE
      )
    }
  )
}
