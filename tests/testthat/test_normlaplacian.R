library(uwot)
library(RSpectra)
context("normalized laplacian")

# this exists as a separate file only because it's easier to comment out as
# part of temporarily removing any function calls to rspectra when using rhub
# with sanitizers and valgrind (probably the extended compilation time with
# eigen causes a preperror)

# These numbers come from running UMAP Python code:
# spectral_layout(pairwise_distances(iris.data[0:10, :]))
# NB:
# 1. iris data in scikit-learn is currently from UCI repo, which has errors
#   (although this doesn't affect the first ten entries)
# 2. eigenvector calculation is not that converged and specifies a starting
#   vector that we can't supply with either RSpectra or eigen.
# 3. The eigenvectors are only identical up to a sign, so we take the absolute
#   values.
abs_expected_norm_lap <-
  abs(
    c2y(
      0.7477, -0.1292, -0.03001, 0.02127, -0.563, -0.01149, 0.1402,
      -0.2725, -0.01241, 0.1084, -0.106, -0.5723, 0.2024, -0.3082,
      0.1642, -5.549e-05, -0.04843, -0.1747, 0.1684, 0.6611
    )
  )
sparse_m <- Matrix::drop0(x2d(iris[1:10, ]))

test_that("normalized laplacian", {
  res <- normalized_laplacian_init(sparse_m)
  expect_equal(abs(res), abs_expected_norm_lap, tolerance = 0.2)
})

test_that("irlba tsvd normalized", {
  res <- irlba_tsvd_normalized_laplacian_init(sparse_m)
  expect_equal(abs(res), abs_expected_norm_lap, tolerance = 0.2)
})

test_that("irlba normalized", {
  res <- irlba_normalized_laplacian_init(sparse_m)
  expect_equal(abs(res), abs_expected_norm_lap, tolerance = 0.2)
})

test_that("laplacian eigenmap", {
  # tested via sklearn
  # things to note:
  # 1. output eigenvectors are not scaled to 1, due to the D^-1/2 transformation
  #    from Lsym's eigenvectors back to Lrw
  # 2. Lsym is formed by calling
  #    scipy.sparse.csgraph.laplacian(normed=True) on the affinity matrix,
  #    which assumes the diagonal is zero.
  #
  #    symmetrized normalized graph laplacian
  # from sklearn.preprocessing import normalize
  # from sklearn.datasets import load_digits
  # from sklearn.manifold import SpectralEmbedding
  # X, _ = load_digits(return_X_y=True)
  # embedding = SpectralEmbedding(n_components=3, n_neighbors=4, affinity="rbf", gamma = 1e-3)
  # X_transformed = embedding.fit_transform(X[:10])
  # normalize(X_transformed, axis = 0)
  expected_lap_eig <- matrix(
    c(
      0.21050269, -0.07732118, 0.63486516,
      -0.33501476, 0.11755963, -0.40229306,
      -0.36728785, 0.38404235, 0.020391,
      0.20458482, -0.04123934, -0.44198941,
      -0.3841261, -0.47833969, 0.17196966,
      0.3883986, -0.03743132, -0.22790212,
      -0.36483447, -0.32492041, 0.01860336,
      -0.27419176, 0.68954246, 0.34392682,
      0.04537549, 0.14056785, 0.12175907,
      0.39742651, -0.00077821, 0.15609656
    ),
    ncol = 3, byrow = TRUE
  )


  A <- matrix(
    c(
      0, 0.0288109, 0.053397, 0.104038, 0.079341, 0.145439, 0.0946093,
      0.0434563, 0.139317, 0.189191, 0.0288109, 0, 0.176753, 0.126438,
      0.100761, 0.108501, 0.197306, 0.0744967, 0.0940433, 0.0614212,
      0.053397, 0.176753, 0, 0.0544213, 0.0662712, 0.0462358, 0.124431,
      0.0876854, 0.162838, 0.0494398, 0.104038, 0.126438, 0.0544213, 0,
      0.0725848, 0.322066, 0.107206, 0.0395971, 0.164969, 0.130029, 0.079341,
      0.100761, 0.0662712, 0.0725848, 0, 0.0532904, 0.255125, 0.0290714,
      0.0634819, 0.0482673, 0.145439, 0.108501, 0.0462358, 0.322066,
      0.0532904, 0, 0.0748701, 0.0202419, 0.187683, 0.380222, 0.0946093,
      0.197306, 0.124431, 0.107206, 0.255125, 0.0748701, 0, 0.0273237,
      0.142702, 0.0647643, 0.0434563, 0.0744967, 0.0876854, 0.0395971,
      0.0290714, 0.0202419, 0.0273237, 0, 0.0584841, 0.0431531, 0.139317,
      0.0940433, 0.162838, 0.164969, 0.0634819, 0.187683, 0.142702,
      0.0584841, 0, 0.158817, 0.189191, 0.0614212, 0.0494398, 0.130029,
      0.0482673, 0.380222, 0.0647643, 0.0431531, 0.158817, 0
    ),
    nrow = 10
  )
  res <- laplacian_eigenmap(A, ndim = 3)
  expect_equal(abs(res), abs(expected_lap_eig), tolerance = 1e-4)

  expect_equal(abs(irlba_laplacian_eigenmap(A, ndim = 3)), abs(expected_lap_eig), tolerance = 1e-4)
})
