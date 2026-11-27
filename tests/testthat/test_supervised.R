test_that("categorical supervision removes edges between different labels", {
  z <- c(0, 1, 3, 6, 10, 15)
  x <- cbind(z, -z)
  labels <- factor(c("a", "a", "a", "b", "b", "b"))
  graph <- function(y = NULL) {
    similarity_graph(
      x,
      y = y,
      n_neighbors = 6,
      target_weight = 1,
      n_threads = 0
    )
  }
  between <- outer(labels, labels, "!=")
  expect_true(any(graph()[between] > 0))

  supervised <- graph(labels)
  expect_true(all(supervised[between] == 0))
  expect_true(any(supervised[!between] > 0))
})

test_that("both numeric and categorical targets influence a mixed target graph", {
  targets <- data.frame(
    score = ynum,
    other = ynum2,
    group = ycat,
    class = ycat2
  )
  graph <- function(y) {
    similarity_graph(
      iris10,
      y = y,
      target_metric = list(euclidean = "score", manhattan = "other"),
      nn_method = "annoy",
      n_neighbors = 4,
      n_threads = 0
    )
  }
  expected <- graph(targets)
  for (column in names(targets)) {
    changed <- targets
    changed[[column]] <- changed[[column]][c(2:10, 1)]
    expect_false(isTRUE(all.equal(graph(changed), expected)), label = column)
  }
})

test_that("local connectivity can be reset on a supervised graph", {
  sm <- Matrix::drop0(matrix(
    c(
      -0.9183907,
      -1.4071020,
      0.70400164,
      0.4990913,
      -0.1631884,
      -0.03232201,
      0.2156861,
      0.4341653,
      0.92592670
    ),
    byrow = TRUE,
    nrow = 3
  ))
  # make matrix positive and symmetric like typical UMAP fuzzy graph
  sms <- (Matrix::t(sm) + sm)^2

  expected <- matrix(
    c(
      1,
      1,
      0.43551409,
      1,
      0.24170431,
      0.23371835,
      0.43551409,
      0.23371835,
      1
    ),
    byrow = TRUE,
    nrow = 3
  )

  # checked against python version
  expect_equal(
    as.matrix(reset_local_connectivity(sms)),
    expected,
    tolerance = 1e-7
  )

  # tested on a modified python version with the effect n_neighbors changed
  expected_reset_local_metric <-
    matrix(
      c(
        1,
        1,
        0.5972302,
        1,
        0.44010492,
        0.43783589,
        0.5972302,
        0.43783589,
        1
      ),
      byrow = TRUE,
      nrow = 3
    )
  expect_equal(
    as.matrix(
      reset_local_connectivity(
        sms,
        reset_local_metric = TRUE,
        num_local_metric_neighbors = 3
      )
    ),
    expected_reset_local_metric,
    tolerance = 1e-7
  )

  expect_equal(
    as.matrix(
      reset_local_connectivity(
        sms,
        reset_local_metric = TRUE,
        num_local_metric_neighbors = 3,
        n_threads = 2
      )
    ),
    expected_reset_local_metric,
    tolerance = 1e-7
  )
})

test_that("weighted intersections agree with the reference graphs", {
  sparr <- new(
    "dgCMatrix",
    i = c(0L, 2L, 0L, 1L, 2L, 0L, 1L),
    p = c(0L, 2L, 5L, 7L),
    Dim = c(3L, 3L),
    Dimnames = list(NULL, NULL),
    x = c(
      0.918390745913514,
      0.215686070576616,
      1.40710203887692,
      0.163188411813119,
      0.434165332563817,
      0.704001636268765,
      0.0323220081795518
    ),
    factors = list()
  )

  sparr2 <- new(
    "dgCMatrix",
    i = c(0L, 1L, 2L, 1L, 2L, 0L, 1L),
    p = c(0L, 3L, 5L, 7L),
    Dim = c(3L, 3L),
    Dimnames = list(NULL, NULL),
    x = c(
      1.68463092,
      2.91620546,
      0.26469792,
      1.08820257,
      0.96444675,
      1.46399222,
      2.72643589
    ),
    factors = list()
  )

  # Numbers taken from Python implementation
  int09 <- general_simplicial_set_intersection(sparr, sparr2, 0.9)
  res09 <- matrix(
    c(
      1.66877087146,
      0.137467853888,
      1.40799953091,
      1.84399206494,
      0.889673751622,
      1.86201852389,
      0.223218799442,
      0.879058365893,
      0.000000
    ),
    nrow = 3,
    byrow = TRUE
  )
  expect_equal(
    as.matrix(int09),
    res09,
    tolerance = 1e-6
  )

  int01 <- general_simplicial_set_intersection(sparr, sparr2, 0.1)
  res01 <- matrix(
    c(
      0.97318335824,
      1.12392924757,
      0.734457833761,
      0.0182018202924,
      0.164728272878,
      0.0361324854953,
      0.186072986202,
      0.432422466467,
      0.000000
    ),
    nrow = 3,
    byrow = TRUE
  )
  expect_equal(
    as.matrix(int01),
    res01,
    tolerance = 1e-6
  )
})

test_that("sparse column maxima include empty columns", {
  sp34 <- Matrix::drop0(matrix(
    nrow = 3,
    byrow = TRUE,
    c(
      0,
      0.7403984,
      0,
      0.6574427,
      0,
      0,
      0.9472488,
      0,
      0,
      0.3039677,
      0.2868714,
      0
    )
  ))

  expect_equal(colMaxs(sp34), c(0, 0.7403984, 0.9472488, 0.6574427))
})

test_that("numeric targets calculate their own neighbors with precomputed input", {
  y <- 2^seq_len(nrow(iris10))
  graph <- function(target) {
    similarity_graph(
      NULL,
      nn_method = nn,
      y = target,
      target_n_neighbors = 3,
      n_threads = 0,
      verbose = FALSE
    )
  }
  ynn <- dist_nn(dist(y), 3)
  expect_equal(graph(y), graph(ynn))
  expect_gt(max(abs(graph(y) - graph(y[c(2:10, 1)]))), 0)
  m <- fit_small_model(y = ynn)
  expect_identical(m$nn_method, "annoy")
  expect_equal(m$n_neighbors, 4)
  expect_equal(ncol(query_neighbors(m)$idx), 4L)
  expect_silent(unload_uwot(m, cleanup = FALSE))
})

test_that("mixed supervision retains the name of a single numeric target", {
  y <- data.frame(group = factor(rep(1:2, 5)), score = seq_len(10))
  for (method in c("annoy", "hnsw")) {
    if (method == "hnsw" && !is_installed("RcppHNSW")) {
      next
    }
    graph <- function(cols) {
      similarity_graph(
        iris10,
        y = y,
        nn_method = method,
        n_neighbors = 4,
        target_metric = list(euclidean = cols, categorical = "group"),
        n_threads = 0
      )
    }
    expect_equal(graph("score"), graph(1L))
  }
})
