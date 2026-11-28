library(uwot)
library(RSpectra)

test_that("new coordinates use weighted or unweighted neighbor averages", {
  diagonal1s <- as(
    Matrix::drop0(diag(1, nrow(V_asymm), ncol(V_asymm))),
    "generalMatrix"
  )
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

  av <- t(matrix(
    c(
      6.00,
      16.00,
      4.75,
      14.75,
      4.00,
      14.00,
      6.50,
      16.50,
      5.25,
      15.25,
      5.00,
      15.00,
      5.50,
      15.50,
      6.00,
      16.00,
      4.50,
      14.50,
      4.75,
      14.75
    ),
    nrow = 10,
    byrow = TRUE
  ))
  embedding <- init_new_embedding(
    train_embedding,
    as.vector(nnt$idx),
    ncol(nnt$idx),
    graph = NULL,
    weighted = FALSE,
    n_threads = 0,
    verbose = FALSE
  )
  expect_equal(embedding, av)

  wav <- t(matrix(
    c(
      4.774600,
      14.77460,
      5.153800,
      15.15380,
      4.120000,
      14.12000,
      5.485100,
      15.48510,
      4.573100,
      14.57310,
      4.000000,
      14.00000,
      5.138362,
      15.13836,
      5.184333,
      15.18433,
      5.191600,
      15.19160,
      5.166667,
      15.16667
    ),
    nrow = 10,
    byrow = TRUE
  ))
  embedding <- init_new_embedding(
    train_embedding,
    as.vector(nnt$idx),
    ncol(nnt$idx),
    graph = dgraph,
    weighted = TRUE,
    n_threads = 0,
    verbose = FALSE
  )
  expect_equal(embedding, wav, tolerance = 1e-5)

  # Check threaded code
  embedding <- init_new_embedding(
    train_embedding,
    as.vector(nnt$idx),
    ncol(nnt$idx),
    graph = NULL,
    weighted = FALSE,
    n_threads = 1,
    verbose = FALSE
  )
  expect_equal(embedding, av)

  embedding <- init_new_embedding(
    train_embedding,
    as.vector(nnt$idx),
    ncol(nnt$idx),
    graph = dgraph,
    weighted = TRUE,
    n_threads = 1,
    verbose = FALSE
  )
  expect_equal(embedding, wav, tolerance = 1e-5)
})


test_that("stored scaling parameters reproduce the training transformation", {
  iris10_range <- scale_input(iris10, scale_type = "range", ret_model = TRUE)
  iris10_rtrans <- apply_scaling(iris10, attr_to_scale_info(iris10_range))
  expect_equal(
    iris10_range,
    iris10_rtrans,
    ignore_attr = c("scaled:range:max", "scaled:range:min")
  )

  iris10_maxabs <- scale_input(iris10, scale_type = "maxabs", ret_model = TRUE)
  iris10_matrans <- apply_scaling(iris10, attr_to_scale_info(iris10_maxabs))
  expect_equal(iris10_maxabs, iris10_matrans, ignore_attr = "scaled:maxabs")

  iris10_scale <- scale_input(iris10, scale_type = "scale", ret_model = TRUE)
  iris10_strans <- apply_scaling(iris10, attr_to_scale_info(iris10_scale))
  expect_equal(iris10_scale, iris10_strans, ignore_attr = "scaled:nzvcols")

  iris10_zv_col <- iris10
  iris10_zv_col[, 3] <- 10
  iris10zvc_scale <- scale_input(
    iris10_zv_col,
    scale_type = "scale",
    ret_model = TRUE
  )
  # scale the original iris10 here on purpose to check that full-variance column
  # is correctly removed
  iris10_zvstrans <- apply_scaling(iris10, attr_to_scale_info(iris10zvc_scale))
  expect_equal(iris10zvc_scale, iris10_zvstrans, ignore_attr = "scaled:nzvcols")

  iris10_none <- scale_input(iris10, scale_type = FALSE, ret_model = TRUE)
  expect_null(attr_to_scale_info(iris10_none))

  iris10_colrange <- scale_input(
    iris10,
    scale_type = "colrange",
    ret_model = TRUE
  )
  iris10_crtrans <- apply_scaling(iris10, attr_to_scale_info(iris10_colrange))
  expect_equal(
    iris10_colrange,
    iris10_crtrans,
    ignore_attr = c("scaled:colrange:max", "scaled:colrange:min")
  )
})

test_that("stored PCA reproduces the training scores", {
  # test pca transform works
  iris10pca <- pca_init(iris10, ndim = 2, ret_extra = TRUE)
  iris10pcat <- apply_pca(iris10, iris10pca)
  expect_equal(iris10pca$scores, iris10pcat, ignore_attr = "dimnames")
})

# #64 (and some #81)
test_that("can use pre-calculated neighbors in transform", {
  set.seed(1337)
  X_train <- as.matrix(iris[c(1:10, 51:60), -5])
  X_test <- as.matrix(iris[101:110, -5])
  iris_train_nn <- annoy_nn(
    X = X_train,
    k = 4,
    metric = "euclidean",
    n_threads = 0,
    ret_index = TRUE
  )
  # (81) test row names are found if it's just the dist matrix of the NN graph
  row.names(iris_train_nn$dist) <- row.names(X_train)
  iris_umap_train <- umap(
    X = NULL,
    nn_method = iris_train_nn,
    ret_model = TRUE,
    n_neighbors = 4
  )
  expect_equal(row.names(iris_umap_train$embedding), row.names(X_train))

  query_ref_nn <- annoy_search(
    X = X_test,
    k = 4,
    ann = iris_train_nn$index,
    n_threads = 0
  )
  # (81) test row names are found if it's just the index matrix of the NN graph
  row.names(query_ref_nn$dist) <- row.names(X_test)

  iris_umap_test <- umap_transform(
    X = NULL,
    model = iris_umap_train,
    nn_method = query_ref_nn,
    ret_extra = c("nn")
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
    X = NULL,
    model = iris_umap_train,
    nn_method = query_ref_nn,
    init = test_init,
    n_epochs = 0
  )
  expect_equal(iris_umap_test_rand0, test_init)
})

test_that("equivalent results with nn graph or sparse distance matrix", {
  set.seed(42)
  iris_even <- as.matrix(iris[seq(2, 75, 2), 1:4])
  iris_odd <- as.matrix(iris[seq(1, 25, 2), 1:4])

  iris_even_nn <- uwot:::annoy_nn(
    X = iris_even,
    k = 10,
    metric = "euclidean",
    ret_index = TRUE
  )
  row.names(iris_even_nn$idx) <- row.names(iris_even)
  row.names(iris_even_nn$dist) <- row.names(iris_even)

  iris_odd_nn <- annoy_search(
    X = iris_odd,
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
    t(uwot:::nng_to_sparse(
      iris_odd_nn$idx,
      as.vector(iris_odd_nn$dist),
      self_nbr = FALSE,
      max_nbr_id = nrow(iris_even)
    ))
  row.names(iris_odd_nn_sp) <- row.names(iris_even_umap$embedding)
  colnames(iris_odd_nn_sp) <- row.names(iris_odd)

  set.seed(42)
  iris_odd_transform_sp <- umap_transform(
    X = NULL,
    iris_even_umap,
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
  test_umap <- umap_transform(
    test,
    train_umap,
    ret_extra = c("fgraph", "localr", "sigma", "nn")
  )
  expect_type(test_umap, "list")
  expect_ok_matrix(test_umap$embedding)
  expect_equal(dim(test_umap$embedding), c(10, 2))
  expect_s4_class(test_umap$fgraph, "Matrix")
  expect_equal(dim(test_umap$fgraph), c(10, 20))
  expect_true(is.numeric(test_umap$localr))
  expect_true(is.numeric(test_umap$sigma))
  expect_true(is.numeric(test_umap$rho))
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

  iris_s12_transform <- umap_transform(
    iris_species_12,
    iris_s3,
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
    umap_transform(iris_species_3, iris_model, seed = 42, n_sgd_threads = 1)
  expect_equal(iris_transform1, iris_transform3)

  # external seed setting should be same as internal
  set.seed(42)
  iris_model2 <-
    umap(iris_species_12, ret_model = TRUE, n_sgd_threads = 1)
  set.seed(42)
  iris_transform4 <-
    umap_transform(iris_species_3, iris_model2, n_sgd_threads = 1)
  expect_equal(iris_transform1, iris_transform4)

  # setting seed explicitly overrides model seed, gives different results
  iris_transform5 <-
    umap_transform(iris_species_3, iris_model, seed = 123, n_sgd_threads = 1)
  diff15 <- iris_transform1 - iris_transform5
  expect_gt(sqrt(sum(diff15 * diff15) / length(diff15)), 0.01)

  # force model seed setting off
  iris_transform6 <-
    umap_transform(iris_species_3, iris_model, seed = FALSE, n_sgd_threads = 1)
  diff16 <- iris_transform1 - iris_transform6
  expect_gt(sqrt(sum(diff16 * diff16) / length(diff16)), 0.01)
  iris_transform7 <-
    umap_transform(iris_species_3, iris_model, seed = FALSE, n_sgd_threads = 1)
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

  iris_transform_10 <- umap_transform(
    iris_species_3,
    iris_model,
    n_epochs = 10,
    ret_extra = "fgraph"
  )

  iris_transform_0 <- umap_transform(
    iris_species_3,
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

  iris_transform_10_no_batch <- umap_transform(
    iris_species_3,
    iris_model_no_batch,
    n_epochs = 10,
    ret_extra = "fgraph"
  )

  iris_transform_0_no_batch <- umap_transform(
    iris_species_3,
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

test_that("scaled PCA models project original-format training rows", {
  m <- fit_small_model(scale = TRUE, pca = 2, pca_method = "svd")
  q <- query_neighbors(m)
  expect_equal(unname(q$idx[, 1]), seq_len(nrow(iris10)))
  expect_equal(unname(q$dist[, 1]), rep(0, nrow(iris10)), tolerance = 1e-6)
})

test_that("HNSW projection defaults its search budget when nn_args is absent", {
  skip_if_not_installed("RcppHNSW")
  model <- fit_small_model(nn_method = "hnsw", nn_args = NULL)
  rows <- c(7L, 3L, 5L)
  query <- iris10[rows, , drop = FALSE]
  neighbors <- query_neighbors(model, query)
  expect_equal(unname(neighbors$idx[, 1]), rows)
  model$nn_args <- list(ef = 10)
  expect_equal(query_neighbors(model, query), neighbors)
})

test_that("each backend queries the fitted PCA space", {
  for (backend in c("annoy", "hnsw", "nndescent")) {
    if (backend == "hnsw" && !is_installed("RcppHNSW")) {
      next
    }
    if (backend == "nndescent" && !is_installed("rnndescent")) {
      next
    }
    m <- fit_small_model(nn_method = backend, pca = 2, pca_method = "svd")
    q <- query_neighbors(m)
    expect_equal(unname(q$idx[, 1]), seq_len(nrow(iris10)))
    expect_equal(unname(q$dist[, 1]), rep(0, nrow(iris10)), tolerance = 1e-5)
  }
})

test_that("a single feature block precedes PCA in fit and projection", {
  cols <- 2:4
  for (pca in list(NULL, 2L)) {
    m <- fit_small_model(
      metric = list(euclidean = cols),
      pca = pca,
      pca_method = "svd",
      ret_nn = TRUE
    )
    ref <- fit_small_model(
      iris10[, cols],
      pca = pca,
      pca_method = "svd",
      ret_nn = TRUE
    )
    expect_equal(m$nn[[1]], ref$nn[[1]])
    expected <- query_neighbors(ref, iris10[, cols])
    expect_equal(query_neighbors(m), expected)
    archive <- tempfile()
    save_uwot(m, archive, unload = TRUE)
    m <- load_uwot(archive)
    expect_equal(query_neighbors(m), expected)
    unload_uwot(m)
    unlink(archive)
  }
})

test_that("models fitted from a sparse graph project dense query neighbors", {
  model <- fit_small_model(NULL, nn_method = dmiris10z, seed = 42)
  idx <- rbind(first = c(1L, 4L), second = c(3L, 7L))
  query <- list(idx = idx, dist = rbind(c(1, 2), c(2, 3)))
  coords <- umap_transform(
    NULL,
    model,
    nn_method = query,
    init = "average",
    n_epochs = 0,
    n_threads = 0
  )
  expected <- rbind(
    first = colMeans(model$embedding[idx[1, ], , drop = FALSE]),
    second = colMeans(model$embedding[idx[2, ], , drop = FALSE])
  )
  expect_identical(rownames(coords), rownames(idx))
  expect_equal(unname(coords), unname(expected), tolerance = 1e-6)
})

test_that("sparse projection averages each query's own neighbors", {
  m <- fit_small_model(iris10[1:4, ], n_components = 1)
  m$embedding <- matrix(c(0, 10, 20, 30))
  s <- Matrix::sparseMatrix(i = 1:4, j = c(1, 2, 2, 2), x = c(1, 1, 2, 3))
  project <- function(...) {
    umap_transform(
      model = m,
      nn_method = s,
      ...,
      n_epochs = 0,
      n_threads = 0,
      ret_extra = "fgraph",
      verbose = FALSE
    )
  }
  average <- project(init = "average")$embedding
  expect_equal(as.vector(average), c(0, 20))
  expect_equal(project(init_weighted = FALSE)$embedding, average)

  q <- project(init = "weighted")
  expect_equal(project()$embedding, q$embedding)
  expect_equal(
    project(init = "weighted", init_weighted = FALSE)$embedding,
    q$embedding
  )
  expected <- as.vector(q$fgraph %*% m$embedding) / Matrix::rowSums(q$fgraph)
  expect_equal(as.vector(q$embedding), expected, tolerance = 1e-6)
})

test_that("batch transform learns at the explicitly requested rate", {
  m <- fit_small_model(batch = TRUE)
  initial <- matrix(0, 3, 2)
  project <- function(...) {
    umap_transform(
      iris10[1:3, ],
      m,
      init = initial,
      n_epochs = 4,
      n_threads = 0,
      learning_rate = 0,
      verbose = FALSE,
      ...
    )
  }
  expect_equal(unname(project()), initial)
  expect_gt(max(abs(project(opt_args = list(alpha = 0.1)))), 0)
})

test_that("saved models retain their RNG and calibration bandwidth", {
  m <- fit_small_model(rng_type = "deterministic", bandwidth = 0.5)
  archive <- tempfile()
  on.exit(unlink(archive))
  save_uwot(m, archive, unload = TRUE)
  m <- load_uwot(archive)
  on.exit(unload_uwot(m), add = TRUE, after = FALSE)
  expect_identical(m$rng_type, "deterministic")
  expect_equal(m$bandwidth, 0.5)
  x <- iris10[1:3, ] + 0.15
  project <- function(model, seed, ...) {
    umap_transform(
      x,
      model,
      seed = seed,
      n_epochs = 4,
      n_threads = 0,
      verbose = FALSE,
      ...
    )
  }
  expect_equal(project(m, 41), project(m, 42))
  for (bandwidth in c(0.5, 1)) {
    if (bandwidth == 1) {
      m$bandwidth <- NULL
      m$rng_type <- NULL
      expect_false(isTRUE(all.equal(project(m, 41), project(m, 42))))
    }
    q <- project(m, 42, ret_extra = c("sigma", "nn"))
    d <- q$nn[[1]]$dist
    expected <- smooth_knn(
      as.vector(t(d)),
      ncol(d),
      TRUE,
      log2(ncol(d)) * bandwidth,
      local_connectivity = 0,
      n_threads = 0,
      grain_size = 1,
      verbose = FALSE,
      ret_sigma = TRUE
    )
    expect_equal(q$sigma, expected$sigma)
  }
})

test_that("only a shared embedding excludes a matching negative index", {
  update <- function(row, batch, shared) {
    result <- optimize_layout_r(
      matrix(c(1, 1), 1),
      if (!shared) matrix(c(0, 2), 1),
      row,
      1L,
      c(0L, as.integer(row == 0), 1L),
      2L,
      2L,
      2L,
      1,
      "tumap",
      list(gamma = 1),
      0.1,
      list(method = "sgd", alpha = 0.1),
      NULL,
      2,
      "deterministic",
      batch,
      0L,
      1L,
      shared,
      FALSE
    )
    result[1, row + 1L]
  }
  for (batch in c(FALSE, TRUE)) {
    expect_equal(update(0L, batch, FALSE), update(1L, batch, FALSE))
    expect_equal(update(0L, batch, TRUE), 1)
  }
})

test_that("models without an epoch count use the size of the query", {
  model <- fit_small_model()
  model$n_epochs <- NULL
  neighbors <- query_neighbors(model, iris10[1:3, ])

  project <- function(n, expected_epochs) {
    rows <- rep_len(seq_len(3), n)
    graph <- lapply(neighbors, function(x) x[rows, , drop = FALSE])
    last_epoch <- 0L
    embedding <- umap_transform(
      model = model,
      nn_method = graph,
      n_threads = 0,
      epoch_callback = function(epoch, n_epochs, coords, fixed_coords) {
        last_epoch <<- epoch
      }
    )
    expect_ok_matrix(embedding, nr = n)
    expect_equal(last_epoch, expected_epochs)
  }

  project(3, 100)
  project(10001, 30)
})

test_that("legacy single-block PCA models can select columns of stored scores", {
  m <- fit_small_model(pca = 2, pca_method = "svd")
  expected <- query_neighbors(m)
  m$metric <- list(euclidean = 1:2)
  expect_equal(query_neighbors(m), expected)
})

test_that("uniform density scaling stays finite in fitting and projection", {
  x <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
  fit <- function(ds, ret_model) {
    umap(
      x,
      nn_method = "annoy",
      n_neighbors = 2,
      init = "random",
      seed = 1,
      n_epochs = 5,
      n_threads = 0,
      dens_scale = ds,
      ret_model = ret_model
    )
  }
  for (ds in c(0, 1)) {
    expect_true(all(is.finite(fit(ds, FALSE))))
    m <- fit(ds, TRUE)
    expect_true(all(is.finite(m$embedding)))
    expect_equal(m$ai, rep(sqrt(m$a), nrow(x)))
    expect_equal(unname(m$rad_coeff), c(log(unname(m$a)), 0))
    expect_true(all(is.finite(umap_transform(
      x,
      m,
      n_epochs = 3,
      n_threads = 0
    ))))
  }
  expect_equal(scale_radii(1:3, 0, 4), rep(2, 3))
})

test_that("weighted projection handles tiny and all-zero query weights", {
  m <- fit_small_model(iris10[1:4, ], n_components = 1)
  m$embedding <- matrix(c(0, 10, 20, 30))
  for (distance in c(110, 1000)) {
    s <- Matrix::sparseMatrix(i = 2L, j = 1L, x = distance, dims = c(4L, 1L))
    q <- umap_transform(
      NULL,
      m,
      nn_method = s,
      n_epochs = 0,
      n_threads = 0,
      ret_extra = "fgraph"
    )
    expect_equal(as.numeric(q$embedding), 10)
    expect_identical(as.numeric(q$fgraph[1, 2]) > 0, distance == 110)
  }
  for (batch in c(FALSE, TRUE)) {
    final <- umap_transform(
      NULL,
      m,
      nn_method = s,
      n_epochs = 5,
      batch = batch,
      n_threads = 0
    )
    expect_equal(as.numeric(final), 10)
  }
  q <- init_transform_parallel(
    t(m$embedding),
    c(2L, 3L, 1L, 4L),
    c(0L, 2L, 4L),
    2,
    c(1e-320, 2e-320, 0, 0)
  )
  expect_equal(as.numeric(q), c(50 / 3, 15), tolerance = 1e-6)
})

test_that("precomputed queries reject invalid training references and empty queries", {
  model <- fit_small_model()
  project <- function(graph, ...) {
    umap_transform(
      model = model,
      nn_method = graph,
      n_epochs = 0,
      n_threads = 0,
      ...
    )
  }
  graph <- query_neighbors(model, iris10[1:3, ])
  expect_error(project(graph, init = matrix(0, 2, 2)), "wrong dimensions")
  expect_error(project(lapply(graph, function(x) x[, -1])), "number of columns")
  graph$idx[1, 1] <- nrow(iris10) + 1L
  expect_error(project(graph), "exceeds the number of training observations")
  graph$idx[1, 1] <- 0L
  expect_error(project(graph), "positive integers")
  graph$idx[1, 1] <- NA_integer_
  expect_error(project(graph), "missing data")

  sparse <- Matrix::sparseMatrix(
    i = 1:2,
    j = c(1, 1),
    x = c(1, 2),
    dims = c(nrow(iris10), 2L)
  )
  expect_error(project(sparse), "at least one neighbor")
  expect_error(project(sparse[, integer(), drop = FALSE]), "no observations")
})

test_that("older model settings reproduce the corresponding explicit projection", {
  model <- fit_small_model()
  old <- model
  old[c("rng_type", "batch", "opt_args", "binary_edge_weights")] <- NULL
  model$batch <- FALSE
  model$opt_args <- list()
  model$binary_edge_weights <- FALSE
  project <- function(m) {
    umap_transform(iris10[1:3, ], m, seed = 12, n_epochs = 4, n_threads = 0)
  }
  for (pcg in list(NULL, FALSE)) {
    old$pcg_rand <- pcg
    model$rng_type <- if (identical(pcg, FALSE)) "tausworthe" else "pcg"
    expect_equal(project(old), project(model))
  }
})

test_that("projection uses the numeric features selected before fitted scaling", {
  x <- data.frame(code = rep(1:2, 5), ui10, constant = 1)
  numeric_ids <- seq_len(ncol(ui10))
  for (columns in list(numeric_ids, colnames(ui10))) {
    for (scaling in list(FALSE, "Z")) {
      model <- fit_small_model(
        x,
        metric = list(euclidean = columns, categorical = "code"),
        scale = scaling
      )
      factor_x <- transform(x, code = factor(code))
      expected <- query_neighbors(model, factor_x)
      expect_equal(unname(expected$idx[, 1]), seq_len(nrow(x)))
      expect_equal(query_neighbors(model, x), expected)
      expect_equal(query_neighbors(model, as.matrix(x)), expected)
    }
  }
  model$numeric_columns <- NULL
  expect_equal(query_neighbors(model, factor_x), expected)
})

test_that("data frame projection matches the fitted numeric features by name", {
  x <- data.frame(ui10, group = factor(rep(1:2, 5)))
  model <- fit_small_model(
    x,
    metric = list(euclidean = 1:4, categorical = "group")
  )
  neighbor_ids <- function(query) query_neighbors(model, query)$idx
  expected <- neighbor_ids(x)
  expect_identical(neighbor_ids(x[c(4:1, 5)]), expected)
  expect_identical(neighbor_ids(x[c(5, 1:4)]), expected)

  renamed <- x
  names(renamed)[1] <- "unknown"
  expect_error(
    neighbor_ids(renamed),
    "Missing fitted numeric column.*Sepal.Length"
  )

  ambiguous <- x[c(5, 1:4)]
  names(ambiguous)[1] <- "Sepal.Length"
  expect_error(neighbor_ids(ambiguous), "unique column names")
})

test_that("Annoy projection inherits the fitted search_k unless overridden", {
  model <- fit_small_model(n_trees = 2, nn_args = list(search_k = 200))
  expect_equal(model$search_k, 200)
  search <- annoy_search
  seen <- numeric()
  local_mocked_bindings(annoy_search = function(..., search_k) {
    seen <<- c(seen, search_k)
    search(..., search_k = search_k)
  })
  query_neighbors(model)
  model$search_k <- 16 # Value stored by older fits before the nn_args override.
  query_neighbors(model)
  umap_transform(iris10, model, search_k = 50, n_epochs = 0, n_threads = 0)
  expect_equal(seen, c(200, 200, 50))
})

test_that("constant density preserves the specified fit and query repulsion", {
  fit <- function(ds) {
    fit_small_model(
      a = 1,
      b = 1,
      repulsion_strength = 0.25,
      dens_scale = ds,
      n_epochs = 8,
      seed = 42
    )
  }
  ordinary <- fit(NULL)
  density <- fit(0)
  expect_equal(density$embedding, ordinary$embedding)
  project <- function(m) {
    umap_transform(ui10, m, n_epochs = 8, n_threads = 0)
  }
  expect_equal(project(density), project(ordinary))
})

test_that("multiple metric blocks ignore density and diagnostics are independent", {
  blocks <- list(euclidean = 1:2, euclidean = 3:4)
  fit <- function(ds) {
    fit_small_model(metric = blocks, dens_scale = ds, n_epochs = 5, seed = 42)
  }
  ordinary <- fit(NULL)
  model <- fit(1)
  expect_equal(model$method, ordinary$method)
  expect_equal(model$embedding, ordinary$embedding)
  project <- function(m, extra = NULL) {
    umap_transform(ui10, m, n_epochs = 5, n_threads = 0, ret_extra = extra)
  }
  expect_equal(project(model), project(ordinary))
  radii <- project(model, "localr")
  bandwidth <- project(model, "sigma")
  expect_equal(radii$localr, bandwidth$sigma + bandwidth$rho)
  expect_length(radii$localr, nrow(ui10))
  expect_equal(radii$embedding, project(model))

  # Older fits used the first block's radii to fit the density curve.
  model$method <- "leopold"
  model$dens_scale <- 1
  model$ai <- rep(1, nrow(iris10))
  model$rad_coeff <- c(0, 1)
  expect_equal(project(model), project(model, "sigma")$embedding)
})

test_that("precomputed blocks preserve umap2's ordinary optimizer", {
  model <- umap2(
    X = NULL,
    nn_method = list(nn, nn),
    init = "random",
    a = 1,
    b = 1,
    dens_scale = 1,
    ret_model = TRUE,
    n_neighbors = 4,
    n_epochs = 4,
    n_threads = 0,
    seed = 42
  )
  expect_equal(model$method, "tumap")
})
