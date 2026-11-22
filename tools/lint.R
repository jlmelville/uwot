if (as.character(utils::packageVersion("lintr")) != "3.4.0") {
  stop("Install lintr 3.4.0 to match CI.")
}
lints <- lintr::lint_package()
print(lints)
quit(status = as.integer(length(lints) > 0L))
