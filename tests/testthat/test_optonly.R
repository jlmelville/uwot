library(uwot)
context("optimization only")

# this mainly exists to hedge against rhub timing out during tests and to still
# exercise the main umap code even with RSpectra and other dependencies
# temporarily removed due to excessive compilation times. Filter on
# this in testthat.R:
#
# test_check("uwot", filter = "optonly")
#
init_coords <- matrix(c(
  -0.286003508982688, 0.205935933716443, 0.212672369696097,
  0.318264664390379, -0.290855854751177, -0.84524521577413, 0.10829983500751,
  -0.163970086771776, 0.611654891362094, 0.12924697210725, 0.0469197151280184,
  0.226976224751362, -0.0547582725501509, -0.0373386048885834,
  -0.0879982948022033, -0.0168061906596455, -0.269189585603006,
  0.0524409053485183, -0.0619916567707883, 0.201745760046477
), ncol = 2)

nn <- list(
  euclidean =
    list(
      idx = matrix(
        c(
          1L, 2L, 3L, 4L, 5L, 6L,
          7L, 8L, 9L, 10L, 5L, 10L, 4L, 3L, 1L, 5L, 3L, 1L, 4L, 2L, 8L,
          3L, 7L, 9L, 8L, 1L, 4L, 5L, 3L, 3L, 10L, 4L, 2L, 10L, 7L, 8L,
          8L, 10L, 2L, 4L
        ),
        ncol = 4
      ),
      dist = matrix(
        c(
          0,
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0.141421356237309, 0.173205080756888,
          0.244948974278318, 0.244948974278318, 0.141421356237309, 0.616441400296898,
          0.264575131106459, 0.173205080756888, 0.3, 0.173205080756888,
          0.173205080756888, 0.3, 0.264575131106459, 0.3, 0.223606797749979,
          0.616441400296898, 0.33166247903554, 0.223606797749979, 0.435889894354067,
          0.316227766016838, 0.469041575982343, 0.331662479035541, 0.3,
          0.316227766016839, 0.458257569495584, 0.7, 0.424264068711929,
          0.33166247903554, 0.509901951359279, 0.316227766016839
        ),
        ncol = 4
      )
    )
)

res <- umap(iris10,
  n_neighbors = 4, n_epochs = 2, learning_rate = 0.5, min_dist = 0.001,
  init = init_coords, verbose = FALSE, n_threads = 0, nn_method = nn
)
expect_ok_matrix(res)

first_coords <- c()
test_callback <- function(epochs, n_epochs, coords) {
  first_coords <<- c(first_coords, coords[1, 1])
}
set.seed(42)
res <- umap(iris10,
  n_neighbors = 4, n_epochs = 10, learning_rate = 0.5, min_dist = 0.001,
  init = init_coords, verbose = FALSE, nn_method = nn,
  batch = TRUE, epoch_callback = test_callback, approx_pow = TRUE,
  n_threads = 2
)
expect_ok_matrix(res)
expect_equal(length(first_coords), 10)
