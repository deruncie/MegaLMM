#include <RcppEigen.h>
#include <ZigguratR.h>
// [[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif


// #include <RcppParallel.h>

using Eigen::Upper;
using Eigen::Lower;

// double precision for original MegaLMM
using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int
using Eigen::ArrayXXd;                  // variable size array, double precision
using Eigen::ArrayXd;                  // variable size array, double precision
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;


MatrixXd rstdnorm_mat(int n,int p);
struct R_matrix {
  Map<MatrixXd> dense;
  MSpMat sparse;
  bool isDense;
  R_matrix(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
};
void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix>& X_vector);

// float for MegaBayesC
using Eigen::MatrixXf;                  // variable size matrix, double precision
using Eigen::VectorXf;                  // variable size vector, double precision
using Eigen::ArrayXXf;                  // variable size array, double precision
using Eigen::ArrayXf;                  // variable size array, double precision

typedef Eigen::SparseMatrix<float> SpMat_f;

MatrixXf rstdnorm_mat_f(int n,int p);
struct General_Matrix_f {
  MatrixXf dense;
  SpMat_f sparse;
  bool triangular;
  bool isDense;
  bool isNULL = false;
  General_Matrix_f(MatrixXf dense_, SpMat_f sparse_,bool triangular_, bool isDense_,bool isNULL_) : dense(dense_), sparse(sparse_), triangular(triangular_), isDense(isDense_), isNULL(isNULL_){}
  MatrixXf solve(MatrixXf) const;
  // General_Matrix_f solve(General_Matrix_f) const;
  MatrixXf tsolve(MatrixXf) const;
  MatrixXf operator*(MatrixXf) const;
  MatrixXf crossprod(MatrixXf) const;
  float get_log_det();
  int rows() const;
  int cols() const;
};
void load_General_Matrix_f_list(const Rcpp::List, std::vector<General_Matrix_f>&, bool);


using namespace Rcpp;
// using namespace RcppParallel;
using namespace Eigen;

static Ziggurat::R::ZigguratR ziggr;

VectorXd find_candidate_states(MatrixXd, double, int);
Rcpp::IntegerVector which(Rcpp::LogicalVector x);
// MatrixXd uncorrelated_prec_mat(VectorXd,VectorXd,VectorXd);
static int MegaLMM_nthreads = 1;

int get_MegaLMM_nthreads();
void set_MegaLMM_nthreads(int);



struct RMatrix {
  Map<MatrixXd> dense;
  MSpMat sparse;
  bool triangular;
  bool isDense;
  bool isNULL = false;
  RMatrix(Map<MatrixXd> dense_, MSpMat sparse_,bool triangular_, bool isDense_,bool isNULL_) : dense(dense_), sparse(sparse_), triangular(triangular_), isDense(isDense_), isNULL(isNULL_){}
  MatrixXd solve(MatrixXd) const;
  // RMatrix solve(RMatrix) const;
  MatrixXd tsolve(MatrixXd) const;
  MatrixXd operator*(MatrixXd) const;
  MatrixXd crossprod(MatrixXd) const;
  float get_log_det();
  int rows() const;
  int cols() const;
};
RMatrix load_RMatrix(const SEXP, bool);

