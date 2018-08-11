library(uwot)
context("Curve Parameters")

expect_equal(as.vector(find_ab_params(spread = 1, min_dist = 0.001)),
  c(1.929, 0.792),
  tol = 1e-3
)

expect_equal(as.vector(find_ab_params(spread = 1, min_dist = 0.1)),
  c(1.577, 0.895),
  tol = 1e-3
)
