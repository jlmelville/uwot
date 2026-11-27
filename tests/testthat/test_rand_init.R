library(uwot)

test_that("uniform random initialization respects its bounds in one and two dimensions", {
  set.seed(42)
  for (ndim in 1:2) {
    uniform <- rand_init(10, ndim)
    expect_ok_matrix(uniform, nc = ndim)
    expect_true(all(abs(uniform) <= 10))
    expect_gt(diff(range(uniform)), 0)
    expect_ok_matrix(rand_init_lv(10, ndim), nc = ndim)
  }
})
