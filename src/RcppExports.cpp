// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// connected_components_undirected
List connected_components_undirected(std::size_t N, const std::vector<int>& indices1, const std::vector<int>& indptr1, const std::vector<int>& indices2, const std::vector<int>& indptr2);
RcppExport SEXP _uwot_connected_components_undirected(SEXP NSEXP, SEXP indices1SEXP, SEXP indptr1SEXP, SEXP indices2SEXP, SEXP indptr2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::size_t >::type N(NSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type indices1(indices1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type indptr1(indptr1SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type indices2(indices2SEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type indptr2(indptr2SEXP);
    rcpp_result_gen = Rcpp::wrap(connected_components_undirected(N, indices1, indptr1, indices2, indptr2));
    return rcpp_result_gen;
END_RCPP
}
// annoy_search_parallel_cpp
List annoy_search_parallel_cpp(const std::string& index_name, NumericMatrix mat, std::size_t n_neighbors, std::size_t search_k, const std::string& metric, std::size_t n_threads, std::size_t grain_size);
RcppExport SEXP _uwot_annoy_search_parallel_cpp(SEXP index_nameSEXP, SEXP matSEXP, SEXP n_neighborsSEXP, SEXP search_kSEXP, SEXP metricSEXP, SEXP n_threadsSEXP, SEXP grain_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type index_name(index_nameSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_neighbors(n_neighborsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type search_k(search_kSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type grain_size(grain_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(annoy_search_parallel_cpp(index_name, mat, n_neighbors, search_k, metric, n_threads, grain_size));
    return rcpp_result_gen;
END_RCPP
}
// calc_row_probabilities_parallel
List calc_row_probabilities_parallel(NumericMatrix nn_dist, double perplexity, std::size_t n_iter, double tol, bool ret_sigma, std::size_t n_threads, std::size_t grain_size);
RcppExport SEXP _uwot_calc_row_probabilities_parallel(SEXP nn_distSEXP, SEXP perplexitySEXP, SEXP n_iterSEXP, SEXP tolSEXP, SEXP ret_sigmaSEXP, SEXP n_threadsSEXP, SEXP grain_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type nn_dist(nn_distSEXP);
    Rcpp::traits::input_parameter< double >::type perplexity(perplexitySEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_sigma(ret_sigmaSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type grain_size(grain_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_row_probabilities_parallel(nn_dist, perplexity, n_iter, tol, ret_sigma, n_threads, grain_size));
    return rcpp_result_gen;
END_RCPP
}
// optimize_layout_r
NumericMatrix optimize_layout_r(NumericMatrix head_embedding, Nullable<NumericMatrix> tail_embedding, const std::vector<unsigned int> positive_head, const std::vector<unsigned int> positive_tail, const std::vector<unsigned int> positive_ptr, unsigned int n_epochs, unsigned int n_head_vertices, unsigned int n_tail_vertices, const std::vector<float> epochs_per_sample, const std::string& method, List method_args, float initial_alpha, List opt_args, Nullable<Function> epoch_callback, float negative_sample_rate, bool pcg_rand, bool batch, std::size_t n_threads, std::size_t grain_size, bool move_other, bool verbose);
RcppExport SEXP _uwot_optimize_layout_r(SEXP head_embeddingSEXP, SEXP tail_embeddingSEXP, SEXP positive_headSEXP, SEXP positive_tailSEXP, SEXP positive_ptrSEXP, SEXP n_epochsSEXP, SEXP n_head_verticesSEXP, SEXP n_tail_verticesSEXP, SEXP epochs_per_sampleSEXP, SEXP methodSEXP, SEXP method_argsSEXP, SEXP initial_alphaSEXP, SEXP opt_argsSEXP, SEXP epoch_callbackSEXP, SEXP negative_sample_rateSEXP, SEXP pcg_randSEXP, SEXP batchSEXP, SEXP n_threadsSEXP, SEXP grain_sizeSEXP, SEXP move_otherSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type head_embedding(head_embeddingSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type tail_embedding(tail_embeddingSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int> >::type positive_head(positive_headSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int> >::type positive_tail(positive_tailSEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int> >::type positive_ptr(positive_ptrSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_epochs(n_epochsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_head_vertices(n_head_verticesSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_tail_vertices(n_tail_verticesSEXP);
    Rcpp::traits::input_parameter< const std::vector<float> >::type epochs_per_sample(epochs_per_sampleSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< List >::type method_args(method_argsSEXP);
    Rcpp::traits::input_parameter< float >::type initial_alpha(initial_alphaSEXP);
    Rcpp::traits::input_parameter< List >::type opt_args(opt_argsSEXP);
    Rcpp::traits::input_parameter< Nullable<Function> >::type epoch_callback(epoch_callbackSEXP);
    Rcpp::traits::input_parameter< float >::type negative_sample_rate(negative_sample_rateSEXP);
    Rcpp::traits::input_parameter< bool >::type pcg_rand(pcg_randSEXP);
    Rcpp::traits::input_parameter< bool >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type grain_size(grain_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type move_other(move_otherSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optimize_layout_r(head_embedding, tail_embedding, positive_head, positive_tail, positive_ptr, n_epochs, n_head_vertices, n_tail_vertices, epochs_per_sample, method, method_args, initial_alpha, opt_args, epoch_callback, negative_sample_rate, pcg_rand, batch, n_threads, grain_size, move_other, verbose));
    return rcpp_result_gen;
END_RCPP
}
// smooth_knn_distances_parallel
List smooth_knn_distances_parallel(NumericMatrix nn_dist, double target, std::size_t n_iter, double local_connectivity, double bandwidth, double tol, double min_k_dist_scale, bool ret_sigma, std::size_t n_threads, std::size_t grain_size);
RcppExport SEXP _uwot_smooth_knn_distances_parallel(SEXP nn_distSEXP, SEXP targetSEXP, SEXP n_iterSEXP, SEXP local_connectivitySEXP, SEXP bandwidthSEXP, SEXP tolSEXP, SEXP min_k_dist_scaleSEXP, SEXP ret_sigmaSEXP, SEXP n_threadsSEXP, SEXP grain_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type nn_dist(nn_distSEXP);
    Rcpp::traits::input_parameter< double >::type target(targetSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_iter(n_iterSEXP);
    Rcpp::traits::input_parameter< double >::type local_connectivity(local_connectivitySEXP);
    Rcpp::traits::input_parameter< double >::type bandwidth(bandwidthSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type min_k_dist_scale(min_k_dist_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type ret_sigma(ret_sigmaSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type grain_size(grain_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_knn_distances_parallel(nn_dist, target, n_iter, local_connectivity, bandwidth, tol, min_k_dist_scale, ret_sigma, n_threads, grain_size));
    return rcpp_result_gen;
END_RCPP
}
// fast_intersection_cpp
NumericVector fast_intersection_cpp(IntegerVector rows, IntegerVector cols, NumericVector values, IntegerVector target, double unknown_dist, double far_dist);
RcppExport SEXP _uwot_fast_intersection_cpp(SEXP rowsSEXP, SEXP colsSEXP, SEXP valuesSEXP, SEXP targetSEXP, SEXP unknown_distSEXP, SEXP far_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type rows(rowsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cols(colsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type target(targetSEXP);
    Rcpp::traits::input_parameter< double >::type unknown_dist(unknown_distSEXP);
    Rcpp::traits::input_parameter< double >::type far_dist(far_distSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_intersection_cpp(rows, cols, values, target, unknown_dist, far_dist));
    return rcpp_result_gen;
END_RCPP
}
// general_sset_intersection_cpp
NumericVector general_sset_intersection_cpp(IntegerVector indptr1, IntegerVector indices1, NumericVector data1, IntegerVector indptr2, IntegerVector indices2, NumericVector data2, IntegerVector result_row, IntegerVector result_col, NumericVector result_val, double mix_weight);
RcppExport SEXP _uwot_general_sset_intersection_cpp(SEXP indptr1SEXP, SEXP indices1SEXP, SEXP data1SEXP, SEXP indptr2SEXP, SEXP indices2SEXP, SEXP data2SEXP, SEXP result_rowSEXP, SEXP result_colSEXP, SEXP result_valSEXP, SEXP mix_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type indptr1(indptr1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices1(indices1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indptr2(indptr2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indices2(indices2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type result_row(result_rowSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type result_col(result_colSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type result_val(result_valSEXP);
    Rcpp::traits::input_parameter< double >::type mix_weight(mix_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(general_sset_intersection_cpp(indptr1, indices1, data1, indptr2, indices2, data2, result_row, result_col, result_val, mix_weight));
    return rcpp_result_gen;
END_RCPP
}
// hardware_concurrency
unsigned int hardware_concurrency();
RcppExport SEXP _uwot_hardware_concurrency() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(hardware_concurrency());
    return rcpp_result_gen;
END_RCPP
}
// init_transform_parallel
NumericMatrix init_transform_parallel(NumericMatrix train_embedding, IntegerMatrix nn_index, Nullable<NumericMatrix> nn_weights, std::size_t n_threads, std::size_t grain_size);
RcppExport SEXP _uwot_init_transform_parallel(SEXP train_embeddingSEXP, SEXP nn_indexSEXP, SEXP nn_weightsSEXP, SEXP n_threadsSEXP, SEXP grain_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type train_embedding(train_embeddingSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type nn_index(nn_indexSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type nn_weights(nn_weightsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< std::size_t >::type grain_size(grain_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(init_transform_parallel(train_embedding, nn_index, nn_weights, n_threads, grain_size));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_uwot_connected_components_undirected", (DL_FUNC) &_uwot_connected_components_undirected, 5},
    {"_uwot_annoy_search_parallel_cpp", (DL_FUNC) &_uwot_annoy_search_parallel_cpp, 7},
    {"_uwot_calc_row_probabilities_parallel", (DL_FUNC) &_uwot_calc_row_probabilities_parallel, 7},
    {"_uwot_optimize_layout_r", (DL_FUNC) &_uwot_optimize_layout_r, 21},
    {"_uwot_smooth_knn_distances_parallel", (DL_FUNC) &_uwot_smooth_knn_distances_parallel, 10},
    {"_uwot_fast_intersection_cpp", (DL_FUNC) &_uwot_fast_intersection_cpp, 6},
    {"_uwot_general_sset_intersection_cpp", (DL_FUNC) &_uwot_general_sset_intersection_cpp, 10},
    {"_uwot_hardware_concurrency", (DL_FUNC) &_uwot_hardware_concurrency, 0},
    {"_uwot_init_transform_parallel", (DL_FUNC) &_uwot_init_transform_parallel, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_uwot(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
