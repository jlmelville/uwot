library(uwot)
context("Smooth kNN distances")

flatmat <- function(x, nr) {
  as.vector(t(matrix(x, nrow = nr)))
}

### C++ tests
nn_8 <- find_nn(iris10, k = 8)
nbrs8 <- ncol(nn_8$dist)
target8 <- log2(nbrs8)
res <-
  smooth_knn_distances_parallel(
    as.vector(t(nn_8$dist)),
    nn_ptr = nbrs8,
    skip_first = TRUE,
    ret_sigma = TRUE,
    target = target8
  )
expect_equal(flatmat(res$matrix, nbrs8), c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 0.883551016667945, 0.563402698221087, 0.789087277089996,
  0.57607769298836, 0.750303047844192, 1, 0.656969510423194, 0.745156972400171,
  0.586089930619138, 0.481555035145846, 0.279104282300107, 0.488196761636568,
  0.514561802750023, 0.489639480054329, 0.330384991495572, 0.626571100714632,
  0.367875074735739, 0.39660773625396, 0.438113497020007, 0.481555035145844,
  0.238036308009313, 0.321078737378799, 0.423034729464725, 0.419490620521321,
  0.275816663602388, 0.120285582819758, 0.297337437716562, 0.247710312149843,
  0.377578194428495, 0.445038652274379, 0.229198240647999, 0.217928223724008,
  0.13265884755527, 0.419490620521318, 0.257869637066222, 0.110625826611838,
  0.260166429776407, 0.231017974667955, 0.364373813389398, 0.220580064268054,
  0.212929653970733, 0.217928223724007, 0.0998021268957899, 0.0776802085446247,
  0.195560609120723, 0.072176661510608, 0.215176296482231, 0.231017974667954,
  0.147146427255277, 0.209014049190051, 0.157184945393181, 0.191460580118967,
  0.0408496922133704, 0.0176222685661076, 0.190057981521641, 0.0703455098666948,
  0.202477876903057, 0.148483915739209, 0.086695317543654, 0.162252224543109
))
expect_equal(res$sigma, c(
  0.2567215, 0.22098923, 0.08285332, 0.09981823,
  0.28608322, 0.17873764, 0.15968704, 0.17134094,
  0.25434113, 0.19572449
))
expected_rho <- c(
  0.14142136, 0.17320508, 0.24494897, 0.24494897,
  0.14142136, 0.6164414, 0.26457513, 0.17320508,
  0.3, 0.17320508
)
expect_equal(res$rho, expected_rho)

nn_4 <- find_nn(iris10, k = 4)
nn4dist <- as.vector(t(nn_4$dist))
nbrs4 <- ncol(nn_4$dist)
target4 <- log2(nbrs4)
res <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    ret_sigma = TRUE,
    target = target4
  )
expect_equal(flatmat(res$matrix, nbrs4), c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 0.838084924053271, 0.538562488894191, 0.672032147261722,
  0.54465912156976, 0.719264468544344, 1, 0.646253111343435, 0.689408428587288,
  0.574825084073865, 0.499998079054375, 0.161908061792063, 0.461447148366392,
  0.327968484367062, 0.455344908480281, 0.280728069552432, 5.12868558931539e-10,
  0.35375192447642, 0.310590623941469, 0.425174782792857, 0.499998079054373
))
# nn4
expected_sigma4 <- c(
  0.17993927, 0.20488739, 0.0493803, 0.09060478,
  0.24940491, 0.00390625, 0.15367126, 0.13551712,
  0.24542618, 0.20633698
)
expect_equal(res$sigma, expected_sigma4)
expected_rho4 <- c(
  0.14142136, 0.17320508, 0.24494897, 0.24494897,
  0.14142136, 0.6164414, 0.26457513, 0.17320508, 0.3,
  0.17320508
)
expect_equal(res$rho, expected_rho4)

# explicitly provide pointers into distances
res <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = seq(from = 0, to = 40, by = 4),
    skip_first = TRUE,
    ret_sigma = TRUE,
    target = target4
  )
expect_equal(flatmat(res$matrix, nbrs4), c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 0.838084924053271, 0.538562488894191, 0.672032147261722,
  0.54465912156976, 0.719264468544344, 1, 0.646253111343435, 0.689408428587288,
  0.574825084073865, 0.499998079054375, 0.161908061792063, 0.461447148366392,
  0.327968484367062, 0.455344908480281, 0.280728069552432, 5.12868558931539e-10,
  0.35375192447642, 0.310590623941469, 0.425174782792857, 0.499998079054373
))
expect_equal(res$sigma, expected_sigma4)
expect_equal(res$rho, expected_rho4)

### Various fuzzy set matrices are defined in helper_fuzzy_sets.R
# unsymmetrized fuzzy set
res <- nng_to_sparse(t(nn_4$idx), res$matrix, self_nbr = TRUE, by_row = FALSE)
expect_equal(res, V_asymm, tol = 1e-4)

# Fuzzy Set Union
expect_equal(fuzzy_set_union(res), V_union, tol = 1e-4)

# mix intersection with union
expect_equal(fuzzy_set_union(res, set_op_mix_ratio = 0.5), V_mix,
  tol = 1e-4
)

# intersection
expect_equal(fuzzy_set_union(res, set_op_mix_ratio = 0), V_intersect,
  tol = 1e-4
)

res_cpp_conn1 <- smooth_knn_distances_parallel(
  nn4dist,
  nn_ptr = nbrs4,
  skip_first = TRUE,
  target = target4,
  n_iter = 64,
  local_connectivity = 1.0,
  tol = 1e-5,
  min_k_dist_scale = 1e-3,
  n_threads = 0
)

expect_equal(
  nng_to_sparse(nn_4$idx, flatmat(res_cpp_conn1$matrix, nbrs4),
    self_nbr = TRUE
  ),
  V_asymm,
  tol = 1e-4
)

res_cpp_conn1.5 <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    target = target4,
    n_iter = 64,
    local_connectivity = 1.5,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = 0
  )
expect_equal(
  nng_to_sparse(t(nn_4$idx), res_cpp_conn1.5$matrix,
    self_nbr = TRUE, by_row = FALSE
  ),
  V_asymm_local,
  tol = 1e-4
)


res_cpp_conn1 <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    target = target4,
    n_iter = 64,
    local_connectivity = 1.0,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = 1,
    grain_size = 1
  )
expect_equal(
  nng_to_sparse(t(nn_4$idx), res_cpp_conn1$matrix,
    self_nbr = TRUE, by_row = FALSE
  ),
  V_asymm,
  tol = 1e-4
)

res_cpp_conn1.5 <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    target = target4,
    n_iter = 64,
    local_connectivity = 1.5,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = 1,
    grain_size = 1
  )
expect_equal(
  nng_to_sparse(t(nn_4$idx), res_cpp_conn1.5$matrix,
    self_nbr = TRUE, by_row = FALSE
  ),
  V_asymm_local,
  tol = 1e-4
)


# Test cross-distances
V_asymm_local_cross <- V_asymm_local
diag(V_asymm_local_cross) <- 1
V_asymm_local_cross <- cbind(
  V_asymm_local_cross,
  matrix(0, nrow = 10, ncol = 2)
)

res_cpp_conn1.5_cross <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    target = target4,
    n_iter = 64,
    local_connectivity = 1.5,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = 0
  )
expect_equal(
  nng_to_sparse(
    t(nn_4$idx),
    res_cpp_conn1.5_cross$matrix,
    by_row = FALSE,
    self_nbr = FALSE,
    max_nbr_id = 12
  ),
  V_asymm_local_cross,
  tol = 1e-4
)

res_cpp_conn1.5_cross <-
  smooth_knn_distances_parallel(
    nn4dist,
    nn_ptr = nbrs4,
    skip_first = TRUE,
    target = target4,
    n_iter = 64,
    local_connectivity = 1.5,
    tol = 1e-5,
    min_k_dist_scale = 1e-3,
    n_threads = 1
  )
expect_equal(
  nng_to_sparse(
    t(nn_4$idx),
    res_cpp_conn1.5_cross$matrix,
    by_row = FALSE,
    self_nbr = FALSE,
    max_nbr_id = 12
  ),
  V_asymm_local_cross,
  tol = 1e-4
)

# smooth_knn_matrix

expected_sknn4m <- Matrix::drop0(matrix(
  c(
    0, 0, 0, 0, 1.0000000, 0, 0, 8.380849e-01, 0, 0.1619081,
    0, 0, 0.5385625, 0.4614471, 0, 0, 0, 0, 0, 1.0000000,
    0, 0.3279685, 0, 1.0000000, 0, 0, 0.6720321, 0, 0, 0,
    0, 0, 1.0000000, 0, 0, 0, 0, 0, 0.5446591, 0.4553449,
    1, 0, 0, 0, 0, 0, 0.2807281, 7.192645e-01, 0, 0,
    1, 0, 0, 0, 1.0000000, 0, 0, 5.128686e-10, 0, 0,
    0, 0, 1.0000000, 0.6462531, 0, 0, 0, 3.537519e-01, 0, 0,
    1, 0, 0, 0, 0.6894084, 0, 0, 0, 0, 0.3105906,
    0, 0.4251748, 0.5748251, 1.0000000, 0, 0, 0, 0, 0, 0,
    0, 1.0000000, 0.4999981, 0.4999981, 0, 0, 0, 0, 0, 0
  ),
  nrow = 10, byrow = TRUE
))

sknn4m <- smooth_knn_matrix(nn_4)$matrix
expect_equal(sknn4m@x, expected_sknn4m@x, tol = 1e-7)
expect_equal(sknn4m@i, expected_sknn4m@i)

nn4sp <- Matrix::drop0(matrix(
  c(
    0, 0, 0, 0, 0.1414214, 0.6164414, 0, 0.1732051, 0, 0,
    0, 0, 0.3000000, 0, 0, 0, 0, 0, 0.5099020, 0.1732051,
    0, 0.3000000, 0, 0.2449490, 0, 0, 0.2645751, 0, 0.4358899, 0.3162278,
    0, 0.3316625, 0.2449490, 0, 0, 0, 0.3316625, 0, 0.3000000, 0.3162278,
    0.1414214, 0, 0, 0, 0, 0.6164414, 0, 0.2236068, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0.2645751, 0, 0.4582576, 0, 0, 0, 0, 0,
    0.1732051, 0, 0, 0, 0.2236068, 0.7000000, 0.4242641, 0, 0, 0,
    0, 0, 0, 0.3000000, 0, 0, 0, 0, 0, 0,
    0.4690416, 0.1732051, 0, 0.3162278, 0, 0, 0, 0.3316625, 0, 0
  ),
  nrow = 10, byrow = TRUE
))

sknn4msp <- smooth_knn_matrix(nn4sp)$matrix
expect_equal(sknn4msp@x, expected_sknn4m@x, tol = 1e-6)
expect_equal(sknn4msp@i, expected_sknn4m@i)

nn3sp <- Matrix::drop0(matrix(c(
  0, 0, 0, 0, 0.1414212, 0.6164416, 0, 0.1732050, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1732050,
  0, 0.3000002, 0, 0.2449490, 0, 0, 0.2645751, 0, 0.4358897, 0,
  0, 0, 0.2449490, 0, 0, 0, 0.3316627, 0, 0.2999998, 0.3162279,
  0.1414212, 0, 0, 0, 0, 0.6164416, 0, 0.2236066, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0.2645751, 0, 0, 0, 0, 0, 0, 0,
  0.1732050, 0, 0, 0, 0.2236066, 0, 0, 0, 0, 0,
  0, 0, 0, 0.2999998, 0, 0, 0, 0, 0, 0,
  0, 0.1732050, 0, 0, 0, 0, 0, 0, 0, 0
), nrow = 10, byrow = TRUE))

expected_sknn3m <- Matrix::drop0(matrix(c(
  0, 0, 0,         0,         1.0000000, 0, 0,         0.5849702, 0,         0,
  0, 0, 0.5849609, 0,         0,         0, 0,         0,         0,         1,
  0, 0, 0,         1.0000000, 0,         0, 0.5849651, 0,         0,         0,
  0, 0, 1.0000000, 0,         0,         0, 0,         0,         0.5849684, 0,
  1, 0, 0,         0,         0,         0, 0,         0.5849684, 0,         0,
  1, 0, 0,         0,         1.0000000, 0, 0,         0,         0,         0,
  0, 0, 1.0000000, 0.5849615, 0,         0, 0,         0,         0,         0,
  1, 0, 0,         0,         0.5849545, 0, 0,         0,         0,         0,
  0, 0, 0.5849692, 1.0000000, 0,         0, 0,         0,         0,         0,
  0, 1, 0,         0.5849544, 0,         0, 0,         0,         0,         0
), nrow = 10, byrow = TRUE))
sknn3msp <- smooth_knn_matrix(nn3sp)$matrix
expect_equal(sknn3msp@x, expected_sknn3m@x, tol = 1e-6)
expect_equal(sknn3msp@i, expected_sknn3m@i)

nn34sp <- nn4sp
nn34sp[, c(2, 4, 6)] <- nn3sp[, c(2, 4, 6)]
expected_sknn34m <-
  Matrix::drop0(matrix(c(
    0, 0,         0,         0,         1.0000000, 0, 0,         0.8380850, 0,        0.1619081,
    0, 0,         0.5849609, 0,         0,         0, 0,         0,         0,        1.0000000,
    0, 0.3279687, 0,         1.0000000, 0,         0, 0.6720329, 0,         0,        0,
    0, 0,         1.0000000, 0,         0,         0, 0,         0,         0.584968, 0,
    1, 0,         0,         0,         0,         0, 0.2807281, 0.7192646, 0,        0,
    1, 0,         0,         0,         1.0000000, 0, 0,         0,         0,        0,
    0, 0,         1.0000000, 0.6462529, 0,         0, 0,         0.3537518, 0,        0,
    1, 0,         0,         0,         0.6894085, 0, 0,         0,         0,        0.3105906,
    0, 0.4251747, 0.5748251, 1.0000000, 0,         0, 0,         0,         0,        0,
    0, 1.0000000, 0.4999980, 0.4999980, 0,         0, 0,         0,         0,        0
  ), nrow = 10, byrow = TRUE))
sknn34msp <- smooth_knn_matrix(nn34sp, ret_sigma = TRUE)
expect_equal(sknn34msp$matrix@x, expected_sknn34m@x, tol = 1e-6)
expect_equal(sknn34msp$matrix@i, expected_sknn34m@i)
expect_equal(sknn34msp$n_failures, 1)
expect_equal(sknn34msp$sigma,
  c(0.1799393, 0.2364655, 0.0493803, 0.1026688, 0.2494049, 1.0000000, 0.1536713, 0.1355171, 0.2454262, 0.2063370),
  tol = 1e-7
)
expect_equal(sknn34msp$rho, expected_rho, tol = 1e-6)
