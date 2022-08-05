#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen; 

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]


// [[Rcpp::export()]]
void record_sample_Posterior_array(Map<MatrixXd> current_sample, Map<MatrixXd> Posterior_array_, int sp_num) {
  int nrow_sample = current_sample.rows();
  int ncol_sample = current_sample.cols();
  if(nrow_sample == 0) return;
  if(ncol_sample == 0) return;

  int length_array = Posterior_array_.size();
  int nr = length_array / ncol_sample;
  if(nr * ncol_sample != length_array) stop("Wrong dimensions of Posterior_array");
  Map<MatrixXd> Posterior_array(Posterior_array_.data(),nr,ncol_sample);

  int n_samples = nr / nrow_sample;
  if(sp_num < 1 || sp_num > n_samples) stop("Array out of bounds");

  for(int i = 0; i < nrow_sample; i++){
    Posterior_array.row(i*n_samples + sp_num-1) = current_sample.row(i);
  }
}

// [[Rcpp::export]]
void set_MegaLMM_nthreads(int threads) {
  if ( threads > 0 ) {
    REprintf("not currently implemented. Use omp_set_num_threads() from RhpcBLASctl package");
      // omp_set_num_threads( threads );
  }
  // REprintf("Number of threads=%i\n", omp_get_max_threads());
}
