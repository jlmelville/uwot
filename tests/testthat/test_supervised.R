context("Supervised")

# categorical y
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1, y = ycat
)
expect_ok_matrix(res)

# numeric y
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1, y = ynum
)
expect_ok_matrix(res)

# mixed categorical and numeric
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  y = data.frame(ycat, ynum)
)
expect_ok_matrix(res)

# multiple categorical y
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  y = data.frame(ycat, ycat2)
)
expect_ok_matrix(res)

# multiple numeric y
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  y = data.frame(ynum, ynum2)
)
expect_ok_matrix(res)

# multiple numeric and categorical
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  y = data.frame(ynum, ynum2, ycat, ycat2)
)
expect_ok_matrix(res)

# multiple numeric with different metrics and categorical
set.seed(1337)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1,
  target_metric = list("euclidean" = 1, "cosine" = 2),
  target_weight = 0.5,
  y = data.frame(ynum, ynum2, ycat, ycat2)
)
expect_ok_matrix(res)

sm <- Matrix::drop0(matrix(c(
  -0.9183907, -1.4071020, 0.70400164,
  0.4990913, -0.1631884, -0.03232201,
  0.2156861, 0.4341653, 0.92592670
), byrow = TRUE, nrow = 3))
# make matrix positive and symmetric like typical UMAP fuzzy graph
sms <- (Matrix::t(sm) + sm)^2

expected <- matrix(c(
  1,
  1, 0.43551409, 1, 0.24170431, 0.23371835,
  0.43551409, 0.23371835, 1
), byrow = TRUE, nrow = 3)

# checked against python version
expect_equal(as.matrix(reset_local_connectivity(sms)), expected,
  tol = 1e-7,
  check.attributes = FALSE
)

# tested on a modified python version with the effect n_neighbors changed
expected_reset_local_metric <-
  matrix(c(
    1, 1, 0.5972302,
    1, 0.44010492, 0.43783589,
    0.5972302, 0.43783589, 1
  ), byrow = TRUE, nrow = 3)
expect_equal(
  as.matrix(
    reset_local_connectivity(sms, reset_local_metric = TRUE, num_local_metric_neighbors = 3)
  ),
  expected_reset_local_metric,
  tol = 1e-7,
  check.attributes = FALSE
)

expect_equal(
  as.matrix(
    reset_local_connectivity(sms, reset_local_metric = TRUE, num_local_metric_neighbors = 3, n_threads = 2)
  ),
  expected_reset_local_metric,
  tol = 1e-7,
  check.attributes = FALSE
)

sparr <- new("dgCMatrix",
  i = c(0L, 2L, 0L, 1L, 2L, 0L, 1L),
  p = c(0L, 2L, 5L, 7L),
  Dim = c(3L, 3L), Dimnames = list(NULL, NULL),
  x = c(
    0.918390745913514, 0.215686070576616, 1.40710203887692,
    0.163188411813119, 0.434165332563817, 0.704001636268765,
    0.0323220081795518
  ), factors = list()
)

sparr2 <- new("dgCMatrix",
  i = c(0L, 1L, 2L, 1L, 2L, 0L, 1L),
  p = c(0L, 3L, 5L, 7L),
  Dim = c(3L, 3L), Dimnames = list(NULL, NULL),
  x = c(
    1.68463092, 2.91620546, 0.26469792,
    1.08820257, 0.96444675, 1.46399222, 2.72643589
  ),
  factors = list()
)


# Numbers taken from Python implementation
int09 <- general_simplicial_set_intersection(sparr, sparr2, 0.9)
res09 <- matrix(
  c(
    1.66877087146, 0.137467853888, 1.40799953091,
    1.84399206494, 0.889673751622, 1.86201852389,
    0.223218799442, 0.879058365893, 0.000000
  ),
  nrow = 3, byrow = TRUE
)
expect_equal(as.matrix(int09), res09, check.attributes = FALSE, tol = 1e-6)


int01 <- general_simplicial_set_intersection(sparr, sparr2, 0.1)
res01 <- matrix(
  c(
    0.97318335824, 1.12392924757, 0.734457833761,
    0.0182018202924, 0.164728272878, 0.0361324854953,
    0.186072986202, 0.432422466467, 0.000000
  ),
  nrow = 3, byrow = TRUE
)
expect_equal(as.matrix(int01), res01, check.attributes = FALSE, tol = 1e-6)

sp34 <- Matrix::drop0(matrix(nrow = 3, byrow = TRUE, c(
  0, 0.7403984, 0, 0.6574427,
  0, 0, 0.9472488, 0,
  0, 0.3039677, 0.2868714, 0
)))

expect_equal(colMaxs(sp34), c(0, 0.7403984, 0.9472488, 0.6574427))
