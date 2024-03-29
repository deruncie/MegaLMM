#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen; 

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

// // [[Rcpp::export]]
// void set_MegaLMM_nthreads(int threads) {
//   if ( threads > 0 ) {
//     // REprintf("not currently implemented. Use omp_set_num_threads() from RhpcBLASctl package");
//     #ifdef _OPENMP
//       omp_set_num_threads( threads );
//       REprintf("Number of threads=%i\n", omp_get_max_threads());
//     #else
//       REprintf("no OMP found?");
//     #endif
//   }
//   // REprintf("Number of threads=%i\n", omp_get_max_threads());
// }

// [[Rcpp::export]]
int get_MegaLMM_nthreads() {
  return MegaLMM_nthreads;
}

// [[Rcpp::export]]
void set_MegaLMM_nthreads(int n_threads) {
  int max_threads = 1;
  #ifdef _OPENMP
    max_threads = omp_get_max_threads();
  #endif
  
  MegaLMM_nthreads = std::max(1,std::min(n_threads,max_threads));
}


// [[Rcpp::export]]
void get_omp_nthreads() {
#ifdef _OPENMP
  REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
}

// [[Rcpp::export]]
void set_omp_nthreads(int threads) {
#ifdef _OPENMP
  omp_set_num_threads( threads );
  set_MegaLMM_nthreads(get_MegaLMM_nthreads()); // ensure this is less than omp_get_max_threads()
  REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
}
