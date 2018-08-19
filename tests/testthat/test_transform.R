library(uwot)
context("Transform")

graph <- V_asymm + diag(1, nrow(V_asymm), ncol(V_asymm))
dV <- as.matrix(graph)
vdV <- as.vector(t(dV))
dgraph <- matrix(vdV[vdV > 0], byrow = TRUE, nrow = 10)
dgraph <- t(apply(dgraph, 1, function(x) {
  sort(x, decreasing = TRUE)
}))

graph <- Matrix::t(graph)


train_embedding <- matrix(1:20, nrow = 10)

av <- matrix(c(
  6.00, 16.00,
  4.75, 14.75,
  4.00, 14.00,
  6.50, 16.50,
  5.25, 15.25,
  5.00, 15.00,
  5.50, 15.50,
  6.00, 16.00,
  4.50, 14.50,
  4.75, 14.75
), nrow = 10, byrow = TRUE)
embedding <- init_new_embedding(train_embedding, nn,
  graph = NULL,
  weighted = FALSE,
  n_threads = 0,
  verbose = FALSE
)
expect_equal(embedding, av, check.attributes = FALSE)

wav <- matrix(c(
  4.774600, 14.77460,
  5.153800, 15.15380,
  4.120000, 14.12000,
  5.485100, 15.48510,
  4.573100, 14.57310,
  4.000000, 14.00000,
  5.138362, 15.13836,
  5.184333, 15.18433,
  5.191600, 15.19160,
  5.166667, 15.16667
), nrow = 10, byrow = TRUE)
embedding <- init_new_embedding(train_embedding, nn,
  graph = dgraph,
  weighted = TRUE,
  n_threads = 0,
  verbose = FALSE
)
expect_equal(embedding, wav, check.attributes = FALSE, tol = 1e-5)


# Check threaded code
RcppParallel::setThreadOptions(numThreads = 1)
embedding <- init_new_embedding(train_embedding, nn,
  graph = NULL,
  weighted = FALSE,
  n_threads = 1,
  verbose = FALSE
)
expect_equal(embedding, av, check.attributes = FALSE)

embedding <- init_new_embedding(train_embedding, nn,
  graph = dgraph,
  weighted = TRUE,
  n_threads = 1,
  verbose = FALSE
)
expect_equal(embedding, wav, check.attributes = FALSE, tol = 1e-5)


iris10_range <- scale_input(iris10, scale_type = "range", ret_model = TRUE)
iris10_rtrans <- apply_scaling(iris10, attr_to_scale_info(iris10_range))
expect_equal(iris10_range, iris10_rtrans, check.attributes = FALSE)

iris10_maxabs <- scale_input(iris10, scale_type = "maxabs", ret_model = TRUE)
iris10_matrans <- apply_scaling(iris10, attr_to_scale_info(iris10_maxabs))
expect_equal(iris10_maxabs, iris10_matrans, check.attributes = FALSE)

iris10_scale <- scale_input(iris10, scale_type = "scale", ret_model = TRUE)
iris10_strans <- apply_scaling(iris10, attr_to_scale_info(iris10_scale))
expect_equal(iris10_scale, iris10_strans, check.attributes = FALSE)

iris10_zv_col <- iris10
iris10_zv_col[, 3] <- 10
iris10zvc_scale <- scale_input(iris10_zv_col,
  scale_type = "scale",
  ret_model = TRUE
)
# scale the original iris10 here on purpose to check that full-variance column
# is correctly removed
iris10_zvstrans <- apply_scaling(iris10, attr_to_scale_info(iris10zvc_scale))
expect_equal(iris10zvc_scale, iris10_zvstrans, check.attributes = FALSE)

iris10_none <- scale_input(iris10, scale_type = FALSE, ret_model = TRUE)
expect_null(attr_to_scale_info(iris10_none))
