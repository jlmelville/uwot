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

write_model_archive <- function(archive) {
  members <- list.files("uwot", full.names = TRUE)
  large <- which(file.size(members) >= ustar_size_limit)
  if (!length(large)) {
    flags <- if (is_win7()) "--force-local" else ""
    return(utils::tar(archive, files = "uwot", extra_flags = flags))
  }

  tar <- find_pax_tar()
  if (is.null(tar)) {
    stop(
      "Model file ",
      members[large[1]],
      " is at least 8 GiB. ",
      "Install GNU tar or bsdtar to save it.",
      call. = FALSE
    )
  }
  utils::tar(
    archive,
    files = "uwot",
    tar = tar$command,
    extra_flags = paste(c("--format=pax", tar$flags), collapse = " ")
  )
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
      status <- utils::untar(
        file,
        exdir = exdir,
        verbose = verbose,
        tar = tar,
        extras = flags
      )
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
