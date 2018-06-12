stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE) {
  message(stime(), " ", ..., domain = domain, appendLF = appendLF)
  utils::flush.console()
}

vtsmessage <- function(..., domain = NULL, appendLF = TRUE, verbose = FALSE) {
  if (verbose) {
    tsmessage(..., domain = domain, appendLF = appendLF)
  }
}

# log vector information
summarize <- function(X, msg = "") {
  summary_X <- summary(X, digits = max(3, getOption("digits") - 3))
  tsmessage(msg, ": ", paste(names(summary_X), ":", summary_X, "|",
                             collapse = ""))
}
