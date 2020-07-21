#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _MegaLMM_find_candidate_states(SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_get_fitted_set_c(SEXP, SEXP);
extern SEXP _MegaLMM_LDLt(SEXP);
extern SEXP _MegaLMM_log_p_h2s(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_make_chol_V_list(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_make_chol_ZtZ_Kinv_list(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_matrix_multiply_toDense(SEXP, SEXP);
extern SEXP _MegaLMM_record_sample_Posterior_array(SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_regression_sampler_parallel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_rstdnorm_mat(SEXP, SEXP);
extern SEXP _MegaLMM_sample_coefs_set_c(SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_factors_scores_c(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_h2s(SEXP);
extern SEXP _MegaLMM_sample_h2s_discrete_MH_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_MME_single_diagK(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_MME_ZKZts_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_tau2_delta_c_Eigen_v2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_sample_trunc_delta_c_Eigen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MegaLMM_set_MegaLMM_nthreads(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_MegaLMM_find_candidate_states",         (DL_FUNC) &_MegaLMM_find_candidate_states,          3},
    {"_MegaLMM_get_fitted_set_c",              (DL_FUNC) &_MegaLMM_get_fitted_set_c,               2},
    {"_MegaLMM_LDLt",                          (DL_FUNC) &_MegaLMM_LDLt,                           1},
    {"_MegaLMM_log_p_h2s",                     (DL_FUNC) &_MegaLMM_log_p_h2s,                      4},
    {"_MegaLMM_make_chol_V_list",              (DL_FUNC) &_MegaLMM_make_chol_V_list,               7},
    {"_MegaLMM_make_chol_ZtZ_Kinv_list",       (DL_FUNC) &_MegaLMM_make_chol_ZtZ_Kinv_list,        8},
    {"_MegaLMM_matrix_multiply_toDense",       (DL_FUNC) &_MegaLMM_matrix_multiply_toDense,        2},
    {"_MegaLMM_record_sample_Posterior_array", (DL_FUNC) &_MegaLMM_record_sample_Posterior_array,  3},
    {"_MegaLMM_regression_sampler_parallel",   (DL_FUNC) &_MegaLMM_regression_sampler_parallel,   19},
    {"_MegaLMM_rstdnorm_mat",                  (DL_FUNC) &_MegaLMM_rstdnorm_mat,                   2},
    {"_MegaLMM_sample_coefs_set_c",            (DL_FUNC) &_MegaLMM_sample_coefs_set_c,             3},
    {"_MegaLMM_sample_factors_scores_c",       (DL_FUNC) &_MegaLMM_sample_factors_scores_c,        5},
    {"_MegaLMM_sample_h2s",                    (DL_FUNC) &_MegaLMM_sample_h2s,                     1},
    {"_MegaLMM_sample_h2s_discrete_MH_c",      (DL_FUNC) &_MegaLMM_sample_h2s_discrete_MH_c,       7},
    {"_MegaLMM_sample_MME_single_diagK",       (DL_FUNC) &_MegaLMM_sample_MME_single_diagK,        8},
    {"_MegaLMM_sample_MME_ZKZts_c",            (DL_FUNC) &_MegaLMM_sample_MME_ZKZts_c,             6},
    {"_MegaLMM_sample_tau2_delta_c_Eigen_v2",  (DL_FUNC) &_MegaLMM_sample_tau2_delta_c_Eigen_v2,   9},
    {"_MegaLMM_sample_trunc_delta_c_Eigen",    (DL_FUNC) &_MegaLMM_sample_trunc_delta_c_Eigen,     8},
    {"_MegaLMM_set_MegaLMM_nthreads",         (DL_FUNC) &_MegaLMM_set_MegaLMM_nthreads,          1},
    {NULL, NULL, 0}
};

void R_init_MegaLMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
