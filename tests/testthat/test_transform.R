library(uwot)
library(RSpectra)
context("Transform")

diagonal1s <- as(Matrix::drop0(diag(1, nrow(V_asymm), ncol(V_asymm))), "generalMatrix")
graph <- V_asymm + diagonal1s
dV <- as.matrix(graph)
vdV <- as.vector(t(dV))
dgraph <- matrix(vdV[vdV > 0], byrow = TRUE, nrow = 10)
dgraph <- apply(dgraph, 1, function(x) {
  sort(x, decreasing = TRUE)
})

graph <- Matrix::t(graph)
nnt <- nn_graph_t(nn)

train_embedding <- t(matrix(1:20, nrow = 10))

av <- t(matrix(c(
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
), nrow = 10, byrow = TRUE))
embedding <- init_new_embedding(train_embedding, as.vector(nnt$idx), ncol(nnt$idx),
  graph = NULL,
  weighted = FALSE,
  n_threads = 0,
  verbose = FALSE
)
expect_equal(embedding, av, check.attributes = FALSE)

wav <- t(matrix(c(
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
), nrow = 10, byrow = TRUE))
embedding <- init_new_embedding(train_embedding, as.vector(nnt$idx), ncol(nnt$idx),
  graph = dgraph,
  weighted = TRUE,
  n_threads = 0,
  verbose = FALSE
)
expect_equal(embedding, wav, check.attributes = FALSE, tol = 1e-5)


# Check threaded code
embedding <- init_new_embedding(train_embedding, as.vector(nnt$idx), ncol(nnt$idx),
  graph = NULL,
  weighted = FALSE,
  n_threads = 1,
  verbose = FALSE
)
expect_equal(embedding, av, check.attributes = FALSE)

embedding <- init_new_embedding(train_embedding, as.vector(nnt$idx), ncol(nnt$idx),
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


iris10_colrange <- scale_input(iris10, scale_type = "colrange", ret_model = TRUE)
iris10_crtrans <- apply_scaling(iris10, attr_to_scale_info(iris10_colrange))
expect_equal(iris10_colrange, iris10_crtrans, check.attributes = FALSE)

# test pca transform works
iris10pca <- pca_init(iris10, ndim = 2, ret_extra = TRUE)
iris10pcat <- apply_pca(iris10, iris10pca)
expect_equal(iris10pca$scores, iris10pcat, check.attributes = FALSE)

# #64 (and some #81)
test_that("can use pre-calculated neighbors in transform", {
  set.seed(1337)
  X_train <- as.matrix(iris[c(1:10, 51:60), -5])
  X_test <- as.matrix(iris[101:110, -5])
  iris_train_nn <- annoy_nn(
    X = X_train, k = 4,
    metric = "euclidean", n_threads = 0,
    ret_index = TRUE
  )
  # (81) test row names are found if it's just the dist matrix of the NN graph
  row.names(iris_train_nn$dist) <- row.names(X_train)
  iris_umap_train <- umap(
    X = NULL, nn_method = iris_train_nn, ret_model = TRUE,
    n_neighbors = 4
  )
  expect_equal(row.names(iris_umap_train$embedding), row.names(X_train))

  query_ref_nn <- annoy_search(
    X = X_test, k = 4,
    ann = iris_train_nn$index, n_threads = 0
  )
  # (81) test row names are found if it's just the index matrix of the NN graph
  row.names(query_ref_nn$dist) <- row.names(X_test)

  iris_umap_test <- umap_transform(
    X = NULL, model = iris_umap_train,
    nn_method = query_ref_nn, ret_extra = c("nn")
  )
  expect_ok_matrix(iris_umap_test$embedding)
  expect_equal(row.names(iris_umap_test$embedding), row.names(X_test))
  expect_equal(iris_umap_test$nn$precomputed$idx, query_ref_nn$idx)
  expect_equal(iris_umap_test$nn$precomputed$dist, query_ref_nn$dist)

  # also test that we can provide our own input and it's unchanged with 0 epochs
  nr <- nrow(query_ref_nn$idx)
  nc <- ncol(iris_umap_train$embedding)
  test_init <- matrix(rnorm(nr * nc), nrow = nr, ncol = nc)
  # set init row name and then set the NN dist names back to NULL to test
  # we can get row names from init matrix if needed
  row.names(test_init) <- row.names(X_test)
  row.names(query_ref_nn$dist) <- NULL

  iris_umap_test_rand0 <- umap_transform(
    X = NULL, model = iris_umap_train,
    nn_method = query_ref_nn,
    init = test_init, n_epochs = 0
  )
  expect_equal(iris_umap_test_rand0, test_init)
})

test_that("equivalent results with nn graph or sparse distance matrix", {
  set.seed(42)
  iris_even <- iris[seq(2, 75, 2), ]
  iris_odd <- iris[seq(1, 25, 2), ]

  iris_even_nn <- uwot:::annoy_nn(
    X = uwot:::x2m(iris_even),
    k = 10,
    metric = "euclidean",
    ret_index = TRUE
  )
  row.names(iris_even_nn$idx) <- row.names(iris_even)
  row.names(iris_even_nn$dist) <- row.names(iris_even)

  iris_odd_nn <- annoy_search(
    X = uwot:::x2m(iris_odd),
    k = 10,
    ann = iris_even_nn$index
  )
  row.names(iris_odd_nn$idx) <- row.names(iris_odd)
  row.names(iris_odd_nn$dist) <- row.names(iris_odd)

  iris_even_nn$index <- NULL

  iris_even_umap <-
    umap(
      X = NULL,
      nn_method = iris_even_nn,
      ret_model = TRUE
    )

  set.seed(42)
  iris_odd_transform_nn_graph <-
    umap_transform(X = NULL, iris_even_umap, nn_method = iris_odd_nn)
  expect_ok_matrix(iris_odd_transform_nn_graph, nrow(iris_odd), 2)
  expect_equal(row.names(iris_odd_transform_nn_graph), row.names(iris_odd))

  iris_odd_nn_sp <-
    t(uwot:::nng_to_sparse(iris_odd_nn$idx, as.vector(iris_odd_nn$dist),
      self_nbr = FALSE,
      max_nbr_id = nrow(iris_even)
    ))
  row.names(iris_odd_nn_sp) <- row.names(iris_even_umap$embedding)
  colnames(iris_odd_nn_sp) <- row.names(iris_odd)

  set.seed(42)
  iris_odd_transform_sp <- umap_transform(
    X = NULL, iris_even_umap,
    nn_method = iris_odd_nn_sp,
    ret_extra = c("nn")
  )
  expect_ok_matrix(iris_odd_transform_sp$embedding, nrow(iris_odd), 2)
  expect_equal(row.names(iris_odd_transform_sp$embedding), row.names(iris_odd))
  expect_equal(iris_odd_transform_sp$embedding, iris_odd_transform_nn_graph)

  expect_equal(iris_odd_transform_sp$nn$precomputed, iris_odd_nn_sp)
})

test_that("n_components can be > n_neighbors (#102)", {
  train <- iris[1:20, ]
  test <- iris[101:110, ]
  set.seed(42)
  train_umap <-
    umap(
      train,
      n_components = 4,
      ret_model = TRUE,
      y = train$Petal.Length,
      init = "rand",
      n_neighbors = 3
    )
  set.seed(42)
  test_umap <- umap_transform(test, train_umap)

  expect_equal(dim(test_umap), c(10, 4))
})

test_that("return transform fgraph (#104)", {
  train <- iris[1:20, ]
  test <- iris[101:110, ]
  set.seed(42)
  train_umap <-
    umap(
      train,
      ret_model = TRUE,
      n_neighbors = 3
    )
  set.seed(42)
  test_umap <- umap_transform(test, train_umap,
    ret_extra = c("fgraph", "localr", "sigma", "nn")
  )
  expect_is(test_umap, "list")
  expect_ok_matrix(test_umap$embedding)
  expect_equal(dim(test_umap$embedding), c(10, 2))
  expect_is(test_umap$fgraph, "Matrix")
  expect_equal(dim(test_umap$fgraph), c(10, 20))
  expect_is(test_umap$localr, "numeric")
  expect_is(test_umap$sigma, "numeric")
  expect_is(test_umap$rho, "numeric")
  expect_equal(length(test_umap$localr), 10)
  expect_equal(length(test_umap$sigma), 10)
  expect_equal(length(test_umap$rho), 10)
  expect_equal(dim(test_umap$nn$euclidean$idx), c(10, 3))
  expect_equal(dim(test_umap$nn$euclidean$dist), c(10, 3))
})

# regression tests the bug reported in #103 where ai and aj were transposed and
# also the data being transformed is larger than the original data
# this at best leads to a wide ring structure being formed for some of the
# transformed data (but may also lead to NaN or a seg fault). This test checks
# that the transformed data doesn't cover a large range, which would be
# diagnostic of the ring forming: with the error present, the range of
# coordinates is around c(-40, 40) vs c(-6, 6) otherwise
test_that("leopold transform (#103)", {
  iris_species_12 <- iris[1:100, ]
  iris_species_3 <- iris[101:150, ]

  set.seed(42)
  iris_s3_leopold <- umap(iris_species_3, dens_scale = 1, ret_model = TRUE)

  set.seed(42)
  iris_s12_transform <- umap_transform(iris_species_12, iris_s3_leopold)

  transform_range <- range(iris_s12_transform)
  # as long as the coordinates aren't in the c(-40, 40) range then we are
  # probably ok
  expect_gt(transform_range[1], -10.0)
  expect_lt(transform_range[2], 10.0)
})


test_that("can transform with binary edge weights", {
  iris_species_12 <- iris[1:100, ]
  iris_species_3 <- iris[101:150, ]

  set.seed(42)
  iris_s3 <- umap(iris_species_3, binary_edge_weights = TRUE, ret_model = TRUE)
  expect_true(iris_s3$binary_edge_weights)

  iris_s12_transform <- umap_transform(iris_species_12, iris_s3,
    ret_extra = c("fgraph")
  )
  expect_true(all(iris_s12_transform$fgraph@x == 1))
})

test_that("transform can set or inherit model seed", {
  iris_species_12 <- iris[1:100, ]
  iris_species_3 <- iris[101:150, ]

  # transform inherits seed from model by default
  iris_model <-
    umap(
      iris_species_12,
      seed = 42,
      ret_model = TRUE,
      n_sgd_threads = 1
    )
  iris_transform1 <-
    umap_transform(iris_species_3, iris_model, n_sgd_threads = 1)
  iris_transform2 <-
    umap_transform(iris_species_3, iris_model, n_sgd_threads = 1)
  expect_equal(iris_transform1, iris_transform2)

  iris_transform3 <-
    umap_transform(iris_species_3,
      iris_model,
      seed = 42,
      n_sgd_threads = 1
    )
  expect_equal(iris_transform1, iris_transform3)

  # external seed setting should be same as internal
  set.seed(42)
  iris_model2 <-
    umap(iris_species_12,
      ret_model = TRUE,
      n_sgd_threads = 1
    )
  set.seed(42)
  iris_transform4 <-
    umap_transform(iris_species_3, iris_model2, n_sgd_threads = 1)
  expect_equal(iris_transform1, iris_transform4)

  # setting seed explicitly overrides model seed, gives different results
  iris_transform5 <-
    umap_transform(iris_species_3,
      iris_model,
      seed = 123,
      n_sgd_threads = 1
    )
  diff15 <- iris_transform1 - iris_transform5
  expect_gt(sqrt(sum(diff15 * diff15) / length(diff15)), 0.01)

  # force model seed setting off
  iris_transform6 <-
    umap_transform(iris_species_3,
      iris_model,
      seed = FALSE,
      n_sgd_threads = 1
    )
  diff16 <- iris_transform1 - iris_transform6
  expect_gt(sqrt(sum(diff16 * diff16) / length(diff16)), 0.01)
  iris_transform7 <-
    umap_transform(iris_species_3,
      iris_model,
      seed = FALSE,
      n_sgd_threads = 1
    )
  diff17 <- iris_transform1 - iris_transform7
  expect_gt(sqrt(sum(diff17 * diff17) / length(diff17)), 0.01)
  # and transforms are different from each other
  diff67 <- iris_transform6 - iris_transform7
  expect_gt(sqrt(sum(diff67 * diff67) / length(diff67)), 0.01)
})

#118 fgraph must be transposed even if n_epochs = 0
test_that("graph dim is consistent when n_epochs = 0", {
  iris_species_12 <- iris[1:100, ]
  iris_species_3 <- iris[101:150, ]

  iris_model <-
    umap(
      iris_species_12,
      ret_model = TRUE,
      n_epochs = 0,
      batch = TRUE
    )

  iris_transform_10 <- umap_transform(iris_species_3,
    iris_model,
    n_epochs = 10,
    ret_extra = "fgraph"
  )

  iris_transform_0 <- umap_transform(iris_species_3,
    iris_model,
    n_epochs = 0,
    ret_extra = "fgraph"
  )

  expect_equal(
    dim(iris_transform_10$fgraph),
    dim(iris_transform_0$fgraph)
  )
  expect_equal(dim(iris_transform_10$fgraph), c(50, 100))


  #118/129 and also without batch
  iris_model_no_batch <-
    umap(
      iris_species_12,
      ret_model = TRUE,
      n_epochs = 0,
      batch = FALSE
    )

  iris_transform_10_no_batch <- umap_transform(iris_species_3,
                                      iris_model_no_batch,
                                      n_epochs = 10,
                                      ret_extra = "fgraph"
  )

  iris_transform_0_no_batch <- umap_transform(iris_species_3,
                                     iris_model_no_batch,
                                     n_epochs = 0,
                                     ret_extra = "fgraph"
  )

  expect_equal(
    dim(iris_transform_10_no_batch$fgraph),
    dim(iris_transform_0_no_batch$fgraph)
  )
  expect_equal(dim(iris_transform_10_no_batch$fgraph), c(50, 100))
})
