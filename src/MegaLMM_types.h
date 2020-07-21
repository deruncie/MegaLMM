// Copyright 2020 Daniel Runcie
// Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
// that can be found in the LICENSE file and available at
// https://polyformproject.org/licenses/noncommercial/1.0.0/

#include <RcppEigen.h>
#include <ZigguratR.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif


// #include <RcppParallel.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXf;                  // variable size matrix, double precision
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXf;                  // variable size vector, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int
using Eigen::ArrayXXf;                  // variable size array, double precision
using Eigen::ArrayXf;                  // variable size array, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::SparseMatrix<float> SpMat;

typedef Eigen::SparseMatrix<double> SpMatd;
typedef Eigen::Map<SpMatd> MSpMatd;
// typedef Eigen::SparseMatrix<float> SpMatf;
// typedef Eigen::Map<SpMatf> MSpMatf;
// typedef Eigen::SparseMatrix<bool> SpMatb;

using namespace Rcpp;
// using namespace RcppParallel;
using namespace Eigen;

static Ziggurat::R::ZigguratR ziggr;

VectorXd find_candidate_states(MatrixXf, double, int);
// MatrixXf uncorrelated_prec_mat(VectorXd,VectorXd,VectorXd);
MatrixXf rstdnorm_mat(int n,int p);
struct General_Matrix_f {
  MatrixXf dense;
  SpMat sparse;
  bool triangular;
  bool isDense;
  bool isNULL = false;
  General_Matrix_f(MatrixXf dense_, SpMat sparse_,bool triangular_, bool isDense_,bool isNULL_) : dense(dense_), sparse(sparse_), triangular(triangular_), isDense(isDense_), isNULL(isNULL_){}
  MatrixXf solve(MatrixXf) const;
  MatrixXf tsolve(MatrixXf) const;
  MatrixXf operator*(MatrixXf) const;
  MatrixXf crossprod(MatrixXf) const;
  float get_log_det();
  int rows() const;
  int cols() const;
};
void load_General_Matrix_f_list(const Rcpp::List, std::vector<General_Matrix_f>&, bool);
