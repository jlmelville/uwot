context("categorical")

sm <- Matrix::drop0(matrix(c(
  -0.9183907, -1.4071020,  0.70400164,
   0.4990913, -0.1631884, -0.03232201,
   0.2156861,  0.4341653,  0.92592670), byrow = TRUE, nrow = 3))

expected <- matrix(c(
  -4.310855, 1,         1,
  1,        -0.7608521, 0.4345031,
  1,         0.4345031, 1), byrow = TRUE, nrow = 3)

expect_equal(as.matrix(reset_local_connectivity(sm)), expected, tol = 1e-7,
             check.attributes = FALSE)

