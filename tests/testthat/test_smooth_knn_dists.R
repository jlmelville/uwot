library(uwot)
context("Smooth kNN distances")

### C++ tests
nn_8 <- find_nn(iris10, k = 8)
res <- smooth_knn_distances_parallel(nn_8$dist)$matrix
expect_equal(as.vector(res), c(
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
# expect_equal(res$sigma, c(
#   0.2567215, 0.22098923, 0.08285332, 0.09981823,
#   0.28608322, 0.17873764, 0.15968704, 0.17134094,
#   0.25434113, 0.19572449
# ))
# expect_equal(res$rho, c(
#   0.14142136, 0.17320508, 0.24494897, 0.24494897,
#   0.14142136, 0.6164414, 0.26457513, 0.17320508,
#   0.3, 0.17320508
# ))


nn_4 <- find_nn(iris10, k = 4)
res <- smooth_knn_distances_parallel(nn_4$dist)$matrix
expect_equal(as.vector(res), c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 0.838084924053271, 0.538562488894191, 0.672032147261722,
  0.54465912156976, 0.719264468544344, 1, 0.646253111343435, 0.689408428587288,
  0.574825084073865, 0.499998079054375, 0.161908061792063, 0.461447148366392,
  0.327968484367062, 0.455344908480281, 0.280728069552432, 5.12868558931539e-10,
  0.35375192447642, 0.310590623941469, 0.425174782792857, 0.499998079054373
))
# nn4
# expect_equal(res$sigma, c(
#   0.17993927, 0.20488739, 0.0493803, 0.09060478,
#   0.24940491, 0.00390625, 0.15367126, 0.13551712,
#   0.24542618, 0.20633698
# ))
# expect_equal(res$rho, c(
#   0.14142136, 0.17320508, 0.24494897, 0.24494897,
#   0.14142136, 0.6164414, 0.26457513, 0.17320508, 0.3,
#   0.17320508
# ))

### Various fuzzy set matrices are defined in helper_fuzzy_sets.R
# unsymmetrized fuzzy set
res <- nn_to_sparse(nn_4$idx, as.vector(res), self_nbr = TRUE)
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


res_cpp_conn1 <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.0,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  n_threads = 0
)$matrix
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1),
  self_nbr = TRUE
), V_asymm, tol = 1e-4)

res_cpp_conn1.5 <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  n_threads = 0
)$matrix
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1.5),
  self_nbr = TRUE
), V_asymm_local, tol = 1e-4)


res_cpp_conn1 <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.0,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  n_threads = 1, grain_size = 1
)$matrix
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1),
  self_nbr = TRUE
), V_asymm, tol = 1e-4)

res_cpp_conn1.5 <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  n_threads = 1, grain_size = 1
)$matrix
expect_equal(nn_to_sparse(nn_4$idx, as.vector(res_cpp_conn1.5),
  self_nbr = TRUE
), V_asymm_local, tol = 1e-4)


# Test cross-distances
V_asymm_local_cross <- V_asymm_local
diag(V_asymm_local_cross) <- 1
V_asymm_local_cross <- cbind(
  V_asymm_local_cross,
  matrix(0, nrow = 10, ncol = 2)
)

res_cpp_conn1.5_cross <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3,
  n_threads = 0
)$matrix
expect_equal(nn_to_sparse(
  nn_4$idx, as.vector(res_cpp_conn1.5_cross),
  self_nbr = FALSE, max_nbr_id = 12
),
V_asymm_local_cross,
tol = 1e-4
)

res_cpp_conn1.5_cross <- smooth_knn_distances_parallel(nn_4$dist,
  n_iter = 64, local_connectivity = 1.5,
  bandwidth = 1.0, tol = 1e-5, min_k_dist_scale = 1e-3, 
  n_threads = 1
)$matrix
expect_equal(nn_to_sparse(
  nn_4$idx, as.vector(res_cpp_conn1.5_cross),
  self_nbr = FALSE, max_nbr_id = 12
),
V_asymm_local_cross,
tol = 1e-4
)
