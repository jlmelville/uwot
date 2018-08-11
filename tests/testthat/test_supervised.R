context("Supervised")

y <- as.factor(c(levels(iris$Species)[c(rep(1, 3), rep(2, 3), rep(3, 3))], NA))
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, alpha = 0.5, min_dist = 0.001,
  init = "rand", verbose = FALSE, n_threads = 1, y = y
)
expect_ok_matrix(res)


sm <- Matrix::drop0(matrix(c(
  -0.9183907, -1.4071020, 0.70400164,
  0.4990913, -0.1631884, -0.03232201,
  0.2156861, 0.4341653, 0.92592670
), byrow = TRUE, nrow = 3))

expected <- matrix(c(
  -4.310855, 1, 1,
  1, -0.7608521, 0.4345031,
  1, 0.4345031, 1
), byrow = TRUE, nrow = 3)

expect_equal(as.matrix(reset_local_connectivity(sm)), expected,
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
res09 <- matrix(c(
  1.66877087146, 0.137467853888, 1.40799953091,
  1.84399206494, 0.889673751622, 1.86201852389,
  0.223218799442, 0.879058365893, 0.000000
),
nrow = 3, byrow = TRUE
)
expect_equal(as.matrix(int09), res09, check.attributes = FALSE, tol = 1e-6)


int01 <- general_simplicial_set_intersection(sparr, sparr2, 0.1)
res01 <- matrix(c(
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
