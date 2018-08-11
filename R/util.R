stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
# appears only if called from an environment where a logical verbose = TRUE
# OR force = TRUE
tsmessage <- function(..., domain = NULL, appendLF = TRUE, force = FALSE) {
  verbose <- get0("verbose", envir = sys.parent())

  if (force || (!is.null(verbose) && verbose)) {
    message(stime(), " ", ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# log vector information
summarize <- function(X, msg = "") {
  summary_X <- summary(X, digits = max(3, getOption("digits") - 3))
  tsmessage(msg, ": ", paste(names(summary_X), ":", summary_X, "|",
    collapse = ""
  ),
  force = get0("verbose", envir = sys.parent())
  )
}

# http://rorynolan.rbind.io/2018/05/08/rcsetseed/
get_seed <- function() {
  sample.int(.Machine$integer.max, 1)
}

# pluralize("thread", 1) => "1 thread"
# pluralize("thread", 2) => "2 threads"
pluralize <- function(str, n, prefix = NULL, inc_num = TRUE) {
  if (n == 0) {
    return("")
  }
  ret <- paste0(str, ifelse(n != 1, "s", ""))
  if (inc_num) {
    ret <- paste0(n, " ", ret)
  }
  if (!is.null(prefix)) {
    ret <- paste0(prefix, " ", ret)
  }
  ret
}
