// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "MegaLMM_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rstdnorm_mat_f
MatrixXf rstdnorm_mat_f(int n, int p);
RcppExport SEXP _MegaLMM_rstdnorm_mat_f(SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(rstdnorm_mat_f(n, p));
    return rcpp_result_gen;
END_RCPP
}
// SingleSite_regression_sampler_parallel
Rcpp::List SingleSite_regression_sampler_parallel(MatrixXf Y, MatrixXf X1_base, Rcpp::List X1_list_, SEXP X2_, SEXP Vx_, Rcpp::IntegerVector h2s_index, Rcpp::List chol_V_list_, VectorXf Y_prec, VectorXf Y_prec_a0, VectorXf Y_prec_b0, MatrixXf prior_prec_alpha1, VectorXf prior_prec_alpha2, MatrixXf prior_mean_beta, MatrixXf prior_prec_beta, SEXP current_alpha1s_, SEXP current_alpha2s_, Rcpp::List BayesAlphabet_parms);
RcppExport SEXP _MegaLMM_SingleSite_regression_sampler_parallel(SEXP YSEXP, SEXP X1_baseSEXP, SEXP X1_list_SEXP, SEXP X2_SEXP, SEXP Vx_SEXP, SEXP h2s_indexSEXP, SEXP chol_V_list_SEXP, SEXP Y_precSEXP, SEXP Y_prec_a0SEXP, SEXP Y_prec_b0SEXP, SEXP prior_prec_alpha1SEXP, SEXP prior_prec_alpha2SEXP, SEXP prior_mean_betaSEXP, SEXP prior_prec_betaSEXP, SEXP current_alpha1s_SEXP, SEXP current_alpha2s_SEXP, SEXP BayesAlphabet_parmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXf >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MatrixXf >::type X1_base(X1_baseSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X1_list_(X1_list_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type X2_(X2_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Vx_(Vx_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type h2s_index(h2s_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_V_list_(chol_V_list_SEXP);
    Rcpp::traits::input_parameter< VectorXf >::type Y_prec(Y_precSEXP);
    Rcpp::traits::input_parameter< VectorXf >::type Y_prec_a0(Y_prec_a0SEXP);
    Rcpp::traits::input_parameter< VectorXf >::type Y_prec_b0(Y_prec_b0SEXP);
    Rcpp::traits::input_parameter< MatrixXf >::type prior_prec_alpha1(prior_prec_alpha1SEXP);
    Rcpp::traits::input_parameter< VectorXf >::type prior_prec_alpha2(prior_prec_alpha2SEXP);
    Rcpp::traits::input_parameter< MatrixXf >::type prior_mean_beta(prior_mean_betaSEXP);
    Rcpp::traits::input_parameter< MatrixXf >::type prior_prec_beta(prior_prec_betaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type current_alpha1s_(current_alpha1s_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type current_alpha2s_(current_alpha2s_SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type BayesAlphabet_parms(BayesAlphabet_parmsSEXP);
    rcpp_result_gen = Rcpp::wrap(SingleSite_regression_sampler_parallel(Y, X1_base, X1_list_, X2_, Vx_, h2s_index, chol_V_list_, Y_prec, Y_prec_a0, Y_prec_b0, prior_prec_alpha1, prior_prec_alpha2, prior_mean_beta, prior_prec_beta, current_alpha1s_, current_alpha2s_, BayesAlphabet_parms));
    return rcpp_result_gen;
END_RCPP
}
// LDLt
List LDLt(SEXP A_);
RcppExport SEXP _MegaLMM_LDLt(SEXP A_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A_(A_SEXP);
    rcpp_result_gen = Rcpp::wrap(LDLt(A_));
    return rcpp_result_gen;
END_RCPP
}
// make_chol_ZtZ_Kinv_list
Rcpp::List make_chol_ZtZ_Kinv_list(Rcpp::List chol_Ki_mats_, Map<MatrixXd> h2s_matrix, MSpMat ZtZ, double drop0_tol, bool verbose, SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar, int ncores);
RcppExport SEXP _MegaLMM_make_chol_ZtZ_Kinv_list(SEXP chol_Ki_mats_SEXP, SEXP h2s_matrixSEXP, SEXP ZtZSEXP, SEXP drop0_tolSEXP, SEXP verboseSEXP, SEXP pbSEXP, SEXP setTxtProgressBarSEXP, SEXP getTxtProgressBarSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_Ki_mats_(chol_Ki_mats_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type h2s_matrix(h2s_matrixSEXP);
    Rcpp::traits::input_parameter< MSpMat >::type ZtZ(ZtZSEXP);
    Rcpp::traits::input_parameter< double >::type drop0_tol(drop0_tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< Function >::type setTxtProgressBar(setTxtProgressBarSEXP);
    Rcpp::traits::input_parameter< Function >::type getTxtProgressBar(getTxtProgressBarSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(make_chol_ZtZ_Kinv_list(chol_Ki_mats_, h2s_matrix, ZtZ, drop0_tol, verbose, pb, setTxtProgressBar, getTxtProgressBar, ncores));
    return rcpp_result_gen;
END_RCPP
}
// make_chol_V_list
Rcpp::List make_chol_V_list(Rcpp::List ZKZts_, Map<MatrixXd> h2s_matrix, double drop0_tol, bool verbose, SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar, int ncores);
RcppExport SEXP _MegaLMM_make_chol_V_list(SEXP ZKZts_SEXP, SEXP h2s_matrixSEXP, SEXP drop0_tolSEXP, SEXP verboseSEXP, SEXP pbSEXP, SEXP setTxtProgressBarSEXP, SEXP getTxtProgressBarSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type ZKZts_(ZKZts_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type h2s_matrix(h2s_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type drop0_tol(drop0_tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< SEXP >::type pb(pbSEXP);
    Rcpp::traits::input_parameter< Function >::type setTxtProgressBar(setTxtProgressBarSEXP);
    Rcpp::traits::input_parameter< Function >::type getTxtProgressBar(getTxtProgressBarSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    rcpp_result_gen = Rcpp::wrap(make_chol_V_list(ZKZts_, h2s_matrix, drop0_tol, verbose, pb, setTxtProgressBar, getTxtProgressBar, ncores));
    return rcpp_result_gen;
END_RCPP
}
// record_sample_Posterior_array
void record_sample_Posterior_array(Map<MatrixXd> current_sample, Map<MatrixXd> Posterior_array_, int sp_num);
RcppExport SEXP _MegaLMM_record_sample_Posterior_array(SEXP current_sampleSEXP, SEXP Posterior_array_SEXP, SEXP sp_numSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type current_sample(current_sampleSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Posterior_array_(Posterior_array_SEXP);
    Rcpp::traits::input_parameter< int >::type sp_num(sp_numSEXP);
    record_sample_Posterior_array(current_sample, Posterior_array_, sp_num);
    return R_NilValue;
END_RCPP
}
// get_MegaLMM_nthreads
int get_MegaLMM_nthreads();
RcppExport SEXP _MegaLMM_get_MegaLMM_nthreads() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_MegaLMM_nthreads());
    return rcpp_result_gen;
END_RCPP
}
// set_MegaLMM_nthreads
void set_MegaLMM_nthreads(int n_threads);
RcppExport SEXP _MegaLMM_set_MegaLMM_nthreads(SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    set_MegaLMM_nthreads(n_threads);
    return R_NilValue;
END_RCPP
}
// get_omp_nthreads
void get_omp_nthreads();
RcppExport SEXP _MegaLMM_get_omp_nthreads() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    get_omp_nthreads();
    return R_NilValue;
END_RCPP
}
// set_omp_nthreads
void set_omp_nthreads(int threads);
RcppExport SEXP _MegaLMM_set_omp_nthreads(SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    set_omp_nthreads(threads);
    return R_NilValue;
END_RCPP
}
// matrix_multiply_toDense
MatrixXd matrix_multiply_toDense(SEXP X_, SEXP Y_);
RcppExport SEXP _MegaLMM_matrix_multiply_toDense(SEXP X_SEXP, SEXP Y_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Y_(Y_SEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_multiply_toDense(X_, Y_));
    return rcpp_result_gen;
END_RCPP
}
// rstdnorm_mat
MatrixXd rstdnorm_mat(int n, int p);
RcppExport SEXP _MegaLMM_rstdnorm_mat(SEXP nSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(rstdnorm_mat(n, p));
    return rcpp_result_gen;
END_RCPP
}
// find_candidate_states
VectorXd find_candidate_states(MatrixXd h2s_matrix, double step_size, int old_state);
RcppExport SEXP _MegaLMM_find_candidate_states(SEXP h2s_matrixSEXP, SEXP step_sizeSEXP, SEXP old_stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type h2s_matrix(h2s_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type old_state(old_stateSEXP);
    rcpp_result_gen = Rcpp::wrap(find_candidate_states(h2s_matrix, step_size, old_state));
    return rcpp_result_gen;
END_RCPP
}
// regression_sampler_parallel
Rcpp::List regression_sampler_parallel(Map<MatrixXd> Y, Map<MatrixXd> X1_base, Rcpp::List X1_list_, Map<MatrixXd> X2, SEXP Vx_, Rcpp::IntegerVector h2s_index, Rcpp::List chol_V_list_, VectorXd Y_prec, VectorXd Y_prec_a0, VectorXd Y_prec_b0, Map<MatrixXd> prior_prec_alpha1, VectorXd prior_prec_alpha2, Map<MatrixXd> prior_mean_beta, Map<MatrixXd> prior_prec_beta);
RcppExport SEXP _MegaLMM_regression_sampler_parallel(SEXP YSEXP, SEXP X1_baseSEXP, SEXP X1_list_SEXP, SEXP X2SEXP, SEXP Vx_SEXP, SEXP h2s_indexSEXP, SEXP chol_V_list_SEXP, SEXP Y_precSEXP, SEXP Y_prec_a0SEXP, SEXP Y_prec_b0SEXP, SEXP prior_prec_alpha1SEXP, SEXP prior_prec_alpha2SEXP, SEXP prior_mean_betaSEXP, SEXP prior_prec_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X1_base(X1_baseSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type X1_list_(X1_list_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Vx_(Vx_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type h2s_index(h2s_indexSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_V_list_(chol_V_list_SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec(Y_precSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_a0(Y_prec_a0SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_b0(Y_prec_b0SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_alpha1(prior_prec_alpha1SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type prior_prec_alpha2(prior_prec_alpha2SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_mean_beta(prior_mean_betaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_beta(prior_prec_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(regression_sampler_parallel(Y, X1_base, X1_list_, X2, Vx_, h2s_index, chol_V_list_, Y_prec, Y_prec_a0, Y_prec_b0, prior_prec_alpha1, prior_prec_alpha2, prior_mean_beta, prior_prec_beta));
    return rcpp_result_gen;
END_RCPP
}
// sample_MME_ZKZts_c
MatrixXd sample_MME_ZKZts_c(Map<MatrixXd> Y, SEXP Z_, Map<VectorXd> tot_Eta_prec, Rcpp::List chol_ZtZ_Kinv_list_, Map<MatrixXd> h2s, VectorXi h2s_index);
RcppExport SEXP _MegaLMM_sample_MME_ZKZts_c(SEXP YSEXP, SEXP Z_SEXP, SEXP tot_Eta_precSEXP, SEXP chol_ZtZ_Kinv_list_SEXP, SEXP h2sSEXP, SEXP h2s_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Z_(Z_SEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type tot_Eta_prec(tot_Eta_precSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_ZtZ_Kinv_list_(chol_ZtZ_Kinv_list_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type h2s(h2sSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type h2s_index(h2s_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_MME_ZKZts_c(Y, Z_, tot_Eta_prec, chol_ZtZ_Kinv_list_, h2s, h2s_index));
    return rcpp_result_gen;
END_RCPP
}
// log_p_h2s
MatrixXd log_p_h2s(Map<MatrixXd> Y, Map<VectorXd> tot_Eta_prec, Rcpp::List chol_V_list_, Map<VectorXd> discrete_priors);
RcppExport SEXP _MegaLMM_log_p_h2s(SEXP YSEXP, SEXP tot_Eta_precSEXP, SEXP chol_V_list_SEXP, SEXP discrete_priorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type tot_Eta_prec(tot_Eta_precSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_V_list_(chol_V_list_SEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type discrete_priors(discrete_priorsSEXP);
    rcpp_result_gen = Rcpp::wrap(log_p_h2s(Y, tot_Eta_prec, chol_V_list_, discrete_priors));
    return rcpp_result_gen;
END_RCPP
}
// sample_h2s
VectorXi sample_h2s(Map<ArrayXXd> log_ps);
RcppExport SEXP _MegaLMM_sample_h2s(SEXP log_psSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<ArrayXXd> >::type log_ps(log_psSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_h2s(log_ps));
    return rcpp_result_gen;
END_RCPP
}
// sample_h2s_discrete_MH_c
VectorXi sample_h2s_discrete_MH_c(Map<MatrixXd> Y, Map<VectorXd> tot_Eta_prec, Map<VectorXd> discrete_priors, VectorXi h2_index, Map<MatrixXd> h2s_matrix, Rcpp::List chol_V_list_, double step_size);
RcppExport SEXP _MegaLMM_sample_h2s_discrete_MH_c(SEXP YSEXP, SEXP tot_Eta_precSEXP, SEXP discrete_priorsSEXP, SEXP h2_indexSEXP, SEXP h2s_matrixSEXP, SEXP chol_V_list_SEXP, SEXP step_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type tot_Eta_prec(tot_Eta_precSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type discrete_priors(discrete_priorsSEXP);
    Rcpp::traits::input_parameter< VectorXi >::type h2_index(h2_indexSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type h2s_matrix(h2s_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type chol_V_list_(chol_V_list_SEXP);
    Rcpp::traits::input_parameter< double >::type step_size(step_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_h2s_discrete_MH_c(Y, tot_Eta_prec, discrete_priors, h2_index, h2s_matrix, chol_V_list_, step_size));
    return rcpp_result_gen;
END_RCPP
}
// sample_factors_scores_c
MatrixXd sample_factors_scores_c(Map<MatrixXd> Eta_tilde, Map<MatrixXd> prior_mean, Map<MatrixXd> Lambda, Map<VectorXd> resid_Eta_prec, Map<VectorXd> F_e_prec);
RcppExport SEXP _MegaLMM_sample_factors_scores_c(SEXP Eta_tildeSEXP, SEXP prior_meanSEXP, SEXP LambdaSEXP, SEXP resid_Eta_precSEXP, SEXP F_e_precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Eta_tilde(Eta_tildeSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_mean(prior_meanSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type resid_Eta_prec(resid_Eta_precSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type F_e_prec(F_e_precSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_factors_scores_c(Eta_tilde, prior_mean, Lambda, resid_Eta_prec, F_e_prec));
    return rcpp_result_gen;
END_RCPP
}
// sample_tau2_delta_c_Eigen_v2
Rcpp::List sample_tau2_delta_c_Eigen_v2(double tau2, double xi, VectorXd delta, Map<VectorXd> scores, double tau_0, double delta_shape, double delta_scale, int p, int times);
RcppExport SEXP _MegaLMM_sample_tau2_delta_c_Eigen_v2(SEXP tau2SEXP, SEXP xiSEXP, SEXP deltaSEXP, SEXP scoresSEXP, SEXP tau_0SEXP, SEXP delta_shapeSEXP, SEXP delta_scaleSEXP, SEXP pSEXP, SEXP timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type tau2(tau2SEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< double >::type tau_0(tau_0SEXP);
    Rcpp::traits::input_parameter< double >::type delta_shape(delta_shapeSEXP);
    Rcpp::traits::input_parameter< double >::type delta_scale(delta_scaleSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type times(timesSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_tau2_delta_c_Eigen_v2(tau2, xi, delta, scores, tau_0, delta_shape, delta_scale, p, times));
    return rcpp_result_gen;
END_RCPP
}
// sample_trunc_delta_c_Eigen
VectorXd sample_trunc_delta_c_Eigen(VectorXd delta, VectorXd tauh, Map<VectorXd> scores, Map<VectorXd> shapes, double delta_1_rate, double delta_2_rate, Map<MatrixXd> randu_draws, double trunc_point);
RcppExport SEXP _MegaLMM_sample_trunc_delta_c_Eigen(SEXP deltaSEXP, SEXP tauhSEXP, SEXP scoresSEXP, SEXP shapesSEXP, SEXP delta_1_rateSEXP, SEXP delta_2_rateSEXP, SEXP randu_drawsSEXP, SEXP trunc_pointSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type tauh(tauhSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type scores(scoresSEXP);
    Rcpp::traits::input_parameter< Map<VectorXd> >::type shapes(shapesSEXP);
    Rcpp::traits::input_parameter< double >::type delta_1_rate(delta_1_rateSEXP);
    Rcpp::traits::input_parameter< double >::type delta_2_rate(delta_2_rateSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type randu_draws(randu_drawsSEXP);
    Rcpp::traits::input_parameter< double >::type trunc_point(trunc_pointSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_trunc_delta_c_Eigen(delta, tauh, scores, shapes, delta_1_rate, delta_2_rate, randu_draws, trunc_point));
    return rcpp_result_gen;
END_RCPP
}
// sample_MME_single_diagK
VectorXd sample_MME_single_diagK(VectorXd y, MatrixXd X, VectorXd prior_mean, VectorXd prior_prec, MSpMat chol_R, double tot_Eta_prec, VectorXd randn_theta, VectorXd randn_e);
RcppExport SEXP _MegaLMM_sample_MME_single_diagK(SEXP ySEXP, SEXP XSEXP, SEXP prior_meanSEXP, SEXP prior_precSEXP, SEXP chol_RSEXP, SEXP tot_Eta_precSEXP, SEXP randn_thetaSEXP, SEXP randn_eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type prior_mean(prior_meanSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type prior_prec(prior_precSEXP);
    Rcpp::traits::input_parameter< MSpMat >::type chol_R(chol_RSEXP);
    Rcpp::traits::input_parameter< double >::type tot_Eta_prec(tot_Eta_precSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type randn_theta(randn_thetaSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type randn_e(randn_eSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_MME_single_diagK(y, X, prior_mean, prior_prec, chol_R, tot_Eta_prec, randn_theta, randn_e));
    return rcpp_result_gen;
END_RCPP
}
// sample_coefs_set_c
MatrixXd sample_coefs_set_c(Rcpp::List model_matrices, Map<MatrixXd> prior_mean, Map<MatrixXd> prior_prec);
RcppExport SEXP _MegaLMM_sample_coefs_set_c(SEXP model_matricesSEXP, SEXP prior_meanSEXP, SEXP prior_precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type model_matrices(model_matricesSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_mean(prior_meanSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec(prior_precSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_coefs_set_c(model_matrices, prior_mean, prior_prec));
    return rcpp_result_gen;
END_RCPP
}
// get_fitted_set_c
MatrixXd get_fitted_set_c(Rcpp::List model_matrices, Map<MatrixXd> coefs);
RcppExport SEXP _MegaLMM_get_fitted_set_c(SEXP model_matricesSEXP, SEXP coefsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type model_matrices(model_matricesSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type coefs(coefsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_fitted_set_c(model_matrices, coefs));
    return rcpp_result_gen;
END_RCPP
}
// parallel_block_regression_sampler
Rcpp::List parallel_block_regression_sampler(Map<MatrixXd> Y, Map<MatrixXd> X1, Map<MatrixXd> X2, SEXP V_, SEXP chol_V_, VectorXd Y_prec, VectorXd Y_prec_a0, VectorXd Y_prec_b0, Map<MatrixXd> prior_prec_alpha, Map<MatrixXd> prior_mean_beta, Map<MatrixXd> prior_prec_beta);
RcppExport SEXP _MegaLMM_parallel_block_regression_sampler(SEXP YSEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP V_SEXP, SEXP chol_V_SEXP, SEXP Y_precSEXP, SEXP Y_prec_a0SEXP, SEXP Y_prec_b0SEXP, SEXP prior_prec_alphaSEXP, SEXP prior_mean_betaSEXP, SEXP prior_prec_betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type V_(V_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type chol_V_(chol_V_SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec(Y_precSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_a0(Y_prec_a0SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_b0(Y_prec_b0SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_alpha(prior_prec_alphaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_mean_beta(prior_mean_betaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_beta(prior_prec_betaSEXP);
    rcpp_result_gen = Rcpp::wrap(parallel_block_regression_sampler(Y, X1, X2, V_, chol_V_, Y_prec, Y_prec_a0, Y_prec_b0, prior_prec_alpha, prior_mean_beta, prior_prec_beta));
    return rcpp_result_gen;
END_RCPP
}
// parallel_Single_regression_sampler
Rcpp::List parallel_Single_regression_sampler(Map<MatrixXd> Y_, Map<MatrixXd> X1_, Map<MatrixXd> X2_, SEXP V_, SEXP chol_V_, VectorXd Y_prec, VectorXd Y_prec_a0, VectorXd Y_prec_b0, Map<MatrixXd> prior_prec_alpha, Map<MatrixXd> prior_mean_beta, Map<MatrixXd> prior_prec_beta, Map<MatrixXd> current_alphas, Map<MatrixXd> betas_alpha, Map<MatrixXd> betas_beta, Map<MatrixXd> betas_pi, MatrixXi betas_delta, int run_sampler_times);
RcppExport SEXP _MegaLMM_parallel_Single_regression_sampler(SEXP Y_SEXP, SEXP X1_SEXP, SEXP X2_SEXP, SEXP V_SEXP, SEXP chol_V_SEXP, SEXP Y_precSEXP, SEXP Y_prec_a0SEXP, SEXP Y_prec_b0SEXP, SEXP prior_prec_alphaSEXP, SEXP prior_mean_betaSEXP, SEXP prior_prec_betaSEXP, SEXP current_alphasSEXP, SEXP betas_alphaSEXP, SEXP betas_betaSEXP, SEXP betas_piSEXP, SEXP betas_deltaSEXP, SEXP run_sampler_timesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X1_(X1_SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type X2_(X2_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type V_(V_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type chol_V_(chol_V_SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec(Y_precSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_a0(Y_prec_a0SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type Y_prec_b0(Y_prec_b0SEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_alpha(prior_prec_alphaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_mean_beta(prior_mean_betaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type prior_prec_beta(prior_prec_betaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type current_alphas(current_alphasSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type betas_alpha(betas_alphaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type betas_beta(betas_betaSEXP);
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type betas_pi(betas_piSEXP);
    Rcpp::traits::input_parameter< MatrixXi >::type betas_delta(betas_deltaSEXP);
    Rcpp::traits::input_parameter< int >::type run_sampler_times(run_sampler_timesSEXP);
    rcpp_result_gen = Rcpp::wrap(parallel_Single_regression_sampler(Y_, X1_, X2_, V_, chol_V_, Y_prec, Y_prec_a0, Y_prec_b0, prior_prec_alpha, prior_mean_beta, prior_prec_beta, current_alphas, betas_alpha, betas_beta, betas_pi, betas_delta, run_sampler_times));
    return rcpp_result_gen;
END_RCPP
}
// cbind_list_withName
void cbind_list_withName(Map<MatrixXd> combined, const Rcpp::List& list, std::string variable);
RcppExport SEXP _MegaLMM_cbind_list_withName(SEXP combinedSEXP, SEXP listSEXP, SEXP variableSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Map<MatrixXd> >::type combined(combinedSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type list(listSEXP);
    Rcpp::traits::input_parameter< std::string >::type variable(variableSEXP);
    cbind_list_withName(combined, list, variable);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MegaLMM_rstdnorm_mat_f", (DL_FUNC) &_MegaLMM_rstdnorm_mat_f, 2},
    {"_MegaLMM_SingleSite_regression_sampler_parallel", (DL_FUNC) &_MegaLMM_SingleSite_regression_sampler_parallel, 17},
    {"_MegaLMM_LDLt", (DL_FUNC) &_MegaLMM_LDLt, 1},
    {"_MegaLMM_make_chol_ZtZ_Kinv_list", (DL_FUNC) &_MegaLMM_make_chol_ZtZ_Kinv_list, 9},
    {"_MegaLMM_make_chol_V_list", (DL_FUNC) &_MegaLMM_make_chol_V_list, 8},
    {"_MegaLMM_record_sample_Posterior_array", (DL_FUNC) &_MegaLMM_record_sample_Posterior_array, 3},
    {"_MegaLMM_get_MegaLMM_nthreads", (DL_FUNC) &_MegaLMM_get_MegaLMM_nthreads, 0},
    {"_MegaLMM_set_MegaLMM_nthreads", (DL_FUNC) &_MegaLMM_set_MegaLMM_nthreads, 1},
    {"_MegaLMM_get_omp_nthreads", (DL_FUNC) &_MegaLMM_get_omp_nthreads, 0},
    {"_MegaLMM_set_omp_nthreads", (DL_FUNC) &_MegaLMM_set_omp_nthreads, 1},
    {"_MegaLMM_matrix_multiply_toDense", (DL_FUNC) &_MegaLMM_matrix_multiply_toDense, 2},
    {"_MegaLMM_rstdnorm_mat", (DL_FUNC) &_MegaLMM_rstdnorm_mat, 2},
    {"_MegaLMM_find_candidate_states", (DL_FUNC) &_MegaLMM_find_candidate_states, 3},
    {"_MegaLMM_regression_sampler_parallel", (DL_FUNC) &_MegaLMM_regression_sampler_parallel, 14},
    {"_MegaLMM_sample_MME_ZKZts_c", (DL_FUNC) &_MegaLMM_sample_MME_ZKZts_c, 6},
    {"_MegaLMM_log_p_h2s", (DL_FUNC) &_MegaLMM_log_p_h2s, 4},
    {"_MegaLMM_sample_h2s", (DL_FUNC) &_MegaLMM_sample_h2s, 1},
    {"_MegaLMM_sample_h2s_discrete_MH_c", (DL_FUNC) &_MegaLMM_sample_h2s_discrete_MH_c, 7},
    {"_MegaLMM_sample_factors_scores_c", (DL_FUNC) &_MegaLMM_sample_factors_scores_c, 5},
    {"_MegaLMM_sample_tau2_delta_c_Eigen_v2", (DL_FUNC) &_MegaLMM_sample_tau2_delta_c_Eigen_v2, 9},
    {"_MegaLMM_sample_trunc_delta_c_Eigen", (DL_FUNC) &_MegaLMM_sample_trunc_delta_c_Eigen, 8},
    {"_MegaLMM_sample_MME_single_diagK", (DL_FUNC) &_MegaLMM_sample_MME_single_diagK, 8},
    {"_MegaLMM_sample_coefs_set_c", (DL_FUNC) &_MegaLMM_sample_coefs_set_c, 3},
    {"_MegaLMM_get_fitted_set_c", (DL_FUNC) &_MegaLMM_get_fitted_set_c, 2},
    {"_MegaLMM_parallel_block_regression_sampler", (DL_FUNC) &_MegaLMM_parallel_block_regression_sampler, 11},
    {"_MegaLMM_parallel_Single_regression_sampler", (DL_FUNC) &_MegaLMM_parallel_Single_regression_sampler, 17},
    {"_MegaLMM_cbind_list_withName", (DL_FUNC) &_MegaLMM_cbind_list_withName, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MegaLMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
