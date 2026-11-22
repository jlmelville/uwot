library(uwot)

# Exercise automatic selection in two separate dependency environments.
mode <- match.arg(commandArgs(trailingOnly = TRUE)[1L], c("full", "minimal"))
backends <- c("RcppHNSW", "rnndescent")
available <- vapply(backends, requireNamespace, logical(1), quietly = TRUE)
if (mode == "full") {
  stopifnot(all(available))
} else {
  stopifnot(!any(available))
}

fit <- function(method = NULL) {
  umap2(
    iris[, 1:4],
    nn_method = method,
    init = "random",
    n_neighbors = 5,
    n_epochs = 3,
    n_threads = 1,
    n_sgd_threads = 1,
    ret_model = TRUE,
    verbose = FALSE,
    seed = 42
  )
}
model <- fit()
expected_index <- if (mode == "full") "hnswv1" else "annoyv1"
stopifnot(
  identical(model$nn_index$type, expected_index),
  all(is.finite(model$embedding))
)
message("Default nearest neighbor index: ", expected_index)
if (mode == "full") {
  model <- fit("nndescent")
  stopifnot(
    identical(model$nn_method, "nndescent"),
    all(is.finite(model$embedding))
  )
  message("Explicit nearest neighbor backend: nndescent")
}
