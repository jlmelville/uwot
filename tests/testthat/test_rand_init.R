library(uwot)
context("Random initialization")

n_vertices <- 10
res_rand <- rand_init(n_vertices, ndim = 2, verbose = FALSE)

expect_ok_matrix(rand_init(n_vertices, ndim = 2, verbose = FALSE))
expect_ok_matrix(rand_init_lv(n_vertices, ndim = 2, verbose = FALSE))

expect_ok_matrix(rand_init(n_vertices, ndim = 1, verbose = FALSE))
expect_ok_matrix(rand_init_lv(n_vertices, ndim = 1, verbose = FALSE))
