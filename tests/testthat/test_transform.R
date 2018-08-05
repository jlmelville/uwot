library(uwot)
context("Transform")

graph <- V_asymm + diag(1, nrow(V_asymm), ncol(V_asymm))
dV <- as.matrix(graph)
vdV <- as.vector(t(dV))
dgraph <- matrix(vdV[vdV > 0], byrow = TRUE, nrow = 10)
dgraph <- t(apply(dgraph, 1, function(x) { sort(x, decreasing = TRUE)}))

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
embedding <- init_new_embedding(train_embedding, nn, graph = NULL,
                                weighted = FALSE,
                                n_threads = 0,
                                verbose = FALSE)
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
embedding <- init_new_embedding(train_embedding, nn, graph = dgraph,
                                weighted = TRUE,
                                n_threads = 0,
                                verbose = FALSE)
expect_equal(embedding, wav, check.attributes = FALSE, tol = 1e-5)


# Check threaded code
RcppParallel::setThreadOptions(numThreads = 1)
embedding <- init_new_embedding(train_embedding, nn, graph = NULL,
                                weighted = FALSE,
                                n_threads = 1,
                                verbose = FALSE)
expect_equal(embedding, av, check.attributes = FALSE)

embedding <- init_new_embedding(train_embedding, nn, graph = dgraph,
                                weighted = TRUE,
                                n_threads = 1,
                                verbose = FALSE)
expect_equal(embedding, wav, check.attributes = FALSE, tol = 1e-5)


# av <- matrix(nrow = nrow(train_embedding), ncol = ncol(train_embedding))
# for (i in 1:nrow(train_embedding)) {
#   av[i, ] <- apply(train_embedding[nn$idx[i, ], ], 2, function(x) { weighted.mean(x, w = graph[nn$idx[i, ], i]) } )
# }
