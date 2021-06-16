#include <RcppEigen.h>
#include <ZigguratR.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif


// #include <RcppParallel.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int
using Eigen::ArrayXXd;                  // variable size array, double precision
using Eigen::ArrayXd;                  // variable size array, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

using namespace Rcpp;
// using namespace RcppParallel;
using namespace Eigen;

static Ziggurat::R::ZigguratR ziggr;

VectorXd find_candidate_states(MatrixXd, double, int);
// MatrixXd uncorrelated_prec_mat(VectorXd,VectorXd,VectorXd);
MatrixXd rstdnorm_mat(int n,int p);
struct R_matrix {
  Map<MatrixXd> dense;
  MSpMat sparse;
  bool isDense;
  R_matrix(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
};
void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix>& X_vector);
