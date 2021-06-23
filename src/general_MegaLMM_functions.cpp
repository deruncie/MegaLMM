// Copyright 2020 Daniel Runcie
// Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
// that can be found in the LICENSE file and available at
// https://polyformproject.org/licenses/noncommercial/1.0.0/


#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

using namespace Eigen;


// -------------------------------------------- //
// ---------- helper functions --------- //
// -------------------------------------------- //
// functions to speed up sparse multiplication and conversion to dense matrices

//' Multiplies two matrices (sparse or dense by dense), returns the product as a dense matrix
//'
//' @param X_ First matrix (matrix or dgCMatrix)
//' @param Y_ Second matrix (matrix)
//' @return Product of X_ and Y_ as a dense matrix
// [[Rcpp::export()]] 
MatrixXf matrix_multiply_toDense(SEXP X_, SEXP Y_){
  if(Rf_isNull(X_)) return(as<MatrixXf >(Y_));
  if(Rf_isMatrix(X_)) {
    MatrixXf X = as<MatrixXf >(X_);
    if(Rf_isMatrix(Y_)) {
      MatrixXf Y = as<MatrixXf >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    } else{
      SpMat Y = as<SpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
  else {
    SpMat X = as<SpMat>(X_);
    if(Rf_isMatrix(Y_)) {
      MatrixXf Y = as<MatrixXf >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      setNbThreads(0);
      if(Eigen::nbThreads( ) > 1) {
        SparseMatrix<float,RowMajor> Xr = X;  // Convert to RowMajor so it can be parallelized
        return(Xr*Y);
      } else{
        return(X*Y);
      }
    } else{
      SpMat Y = as<SpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
}


//' Draws a matrix of standard normal variables
//' 
//' Uses the ziggr function of the RcppZiggurat package as the RNG 
//'
//' @param n number of rows of matrix
//' @param p number of columns of matrix
//' @return nxp matrix of independent std normal values
// [[Rcpp::export()]]
MatrixXf rstdnorm_mat(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat.cast<float>());
}

MatrixXd rstdnorm_mat_d(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat);
}


//' Finds the set of variance component proportions within a specified distance from a starting proportion
//'
//' @param h2s_matrix Mxl matrix of all valid variance component proportions for M random effects
//' @param step_size value measuring the maximum euclidean distance to a new set of variance component proportions
//' @param old_state index in h2s_matrix of the current value of the variance component proportions
//' @return vector of indices of \code{h2s_matrix} giving new candidate variance component proportions
// [[Rcpp::export()]]
VectorXf find_candidate_states(
    MatrixXf h2s_matrix,
    float step_size,
    int old_state
) {
  VectorXf dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
  VectorXf indices(dists.size());
  int count = 0;
  for(int i = 0; i < dists.size(); i++){
    if(dists[i] < step_size & dists[i] > 0) {
      indices[count] = i;
      count++;
    }
  }
  if(count == 0) {  // return all indices as candidates
    for(int i = 0; i < dists.size(); i++){
      indices[count] = i;
      count++;
    }
  }
  return indices.head(count);
}


// code to convert list of R matrices (sparse or dense) into a thread-safe object
// struct General_Matrix_f {
//   MatrixXf dense;
//   MSpMat sparse;
//   bool isDense;
//   General_Matrix_f(MatrixXf dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// };

// Loads a sparse or dense matrix passed from R into a \code{General_Matrix_f} object
General_Matrix_f load_General_Matrix_f(SEXP X_, bool triangular) {
  MatrixXf null_d = MatrixXf::Zero(0,0);
  if(Rf_isNull(X_)) {
    SpMat null_s = null_d.sparseView();
    General_Matrix_f Xm(null_d,null_s,triangular,false,true);
    return(Xm);
  } else if(Rf_isMatrix(X_)){
    MatrixXf X = as<MatrixXf >(X_);
    SpMat null_s = null_d.sparseView();
    General_Matrix_f Xm(X,null_s,triangular,true,false);
    return(Xm);
  } else{
    SpMat X = as<SpMat>(X_);
    General_Matrix_f Xm(null_d,X,triangular,false,false);
    return(Xm);
  }
}

MatrixXf General_Matrix_f::solve(MatrixXf Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.triangularView<Upper>().solve(Y));
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse.triangularView<Upper>().solve(Y));
    }
  } else{
    // Note: these Cholesky's could be stored so they could be re-used.
    if(isDense) {
      return(dense.ldlt().solve(Y));
    } else{
      Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
      return(ldlt.solve(Y));
    }
  }
}

// MatrixXf General_Matrix_f::solve(General_Matrix_f Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       if(Y.isDense) {
//         return(dense.triangularView<Upper>().solve(Y.dense));
//       } else {
//         MatrixXf Z(dense.rows,Y.cols);
//         for(int j = 0; j < Y.cols; j++) Z.col(j) = dense.triangularView<Upper>().solve(Y.sparse.col(j));
//       }
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       if(Y.isDense) {
//         return(sparse.triangularView<Upper>().solve(Y.dense));
//       } else {
//         return(sparse.triangularView<Upper>().solve(Y.sparse));
//       }
//     }
//   } else{
//     // Note: these Cholesky's could be stored so they could be re-used.
//     if(isDense) {
//       if(Y.isDense) {
//         return(dense.ldlt().solve(Y.dense));
//       } else {
//         return(dense.ldlt().solve(Y.sparse));
//       }
//     } else{
//       Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
//       if(Y.isDense) {
//         return(ldlt.solve(Y.dense));
//       } else {
//         return(ldlt.solve(Y.sparse));
//       }
//     }
//   }
// }
MatrixXf General_Matrix_f::tsolve(MatrixXf Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.transpose().triangularView<Lower>().solve(Y));
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse.transpose().triangularView<Lower>().solve(Y));
    }
  } else{
    // Note: these Cholesky's could be stored so they could be re-used.
    if(isDense) {
      return(dense.transpose().ldlt().solve(Y));
    } else{
      Eigen::SimplicialLDLT<SpMat> ldlt(sparse.transpose());
      return(ldlt.solve(Y));
    }
  }
}

MatrixXf General_Matrix_f::operator*(MatrixXf Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.triangularView<Upper>() * Y);
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse * Y);
    }
  } else {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense * Y);
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse * Y);
    }
  }
}
MatrixXf General_Matrix_f::crossprod(MatrixXf Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
      return(dense.transpose().triangularView<Lower>() * Y);
    } else{
      if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
      return(sparse.transpose().triangularView<Lower>() * Y);
    }
  } else{
    if(isDense) {
      if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
      return(dense.transpose() * Y);
    } else{
      if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
      return(sparse.transpose() * Y);
    }
  }
}
float General_Matrix_f::get_log_det() {
  if(isNULL) stop("no determinant of null matrix");
  if(!triangular) stop("not implemented for non-triangular matrix");
  if(isDense) {
    return(dense.diagonal().array().log().sum());
  } else {
    float log_det = 0;
    for(int j = 0; j < sparse.rows(); j++) {
      log_det += std::log(sparse.coeffRef(j,j));
    }
    return(log_det);
  }
}
int General_Matrix_f::rows() const {
  int r;
  if(isDense) {
    r = dense.rows(); 
  } else {
    r = sparse.rows();
  }
  return(r);
}
int General_Matrix_f::cols() const {
  int c;
  if(isDense) {
    c = dense.cols(); 
  } else {
    c = sparse.cols();
  }
  return(c);
}

// Code to convert list of R matrices (sparse or dense) into a thread-safe object
//
// @param X_list List of matrices (each can be dgCMatrix or matrix)
// @param X_vector \code{std::vector} of \code{General_Matrix_f} type which will be populated
void load_General_Matrix_f_list(const Rcpp::List X_list, std::vector<General_Matrix_f>& X_vector, bool triangular){
  int p = X_list.size();
  X_vector.reserve(p);
  for(int i = 0; i < p; i++){
    SEXP Xi_ = X_list[i];
    X_vector.push_back(load_General_Matrix_f(Xi_, triangular));
  }
}

// -------------------------------------------- //
// ---------- regression_sampler --------- //
// -------------------------------------------- //

// Replicates the \code{which} function from R
// 
// @param x Logical vector
// @return IntegerVector with indices of the TRUE values of \code{x}
Rcpp::IntegerVector which(Rcpp::LogicalVector x) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x];
}

// Replicates the \code{which} function from R
// 
// @param chol_R \code{General_Matrix_f} object with the upper-triangular Cholesky decomposition of a square nxn matrix R
// @param X nxp matrix
// @return solve(t(R),X) as a dense matrix
// MatrixXf get_RinvtX2(const General_Matrix_f& chol_V, MatrixXf X){
//   MatrixXf RinvtX2;
//   if(chol_V.isDense) {
//     RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   } else{
//     RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   }
//   return(RinvtX2);
// }

VectorXf regression_sampler_v1(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& X1,           // nxa
    const MatrixXf& RinvtX2,                // nxb
    const MatrixXf& C,                     // bxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const General_Matrix_f& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix. Also called R here with RtR = V
    float Y_prec,                    // float
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
    const float rgamma_1,
    const float Y_prec_b0
){
  int n = y.size();
  int a = X1.cols();
  int b = RinvtX2.cols();
  
  // Check inputs
  if(X1.rows() != n) stop("Wrong dimension of X1");
  if(RinvtX2.rows() != n) stop("Wrong dimension of X2");
  if(C.rows() != b || C.cols() != b) stop("Wrong dimension of C");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");
  
  // Calculate cholesky of A_beta
  // C = Xt(RtR)^-1X
  MatrixXf C_beta = C;
  C_beta.diagonal() += prior_prec_beta;
  LLT<MatrixXf> A_beta_llt;
  A_beta_llt.compute(C_beta);
  MatrixXf chol_A_beta = A_beta_llt.matrixU();
  // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // Calculate A_alpha = Y_prec*X1^T*V_beta^{-1}*X1 + D_alpha^{-1}
    // We don't need to actually calculate V_beta^{-1} directly.
    // Instead,
    // MatrixXf RinvtX1 = get_RinvtX2(chol_V,X1);  // n*n*a -> n x a
    MatrixXf RinvtX1 = chol_V.tsolve(X1);// n*n*a -> n x a
    MatrixXf X1tVinvX2 = RinvtX1.transpose() * RinvtX2; // a*n*b -> a*b
    MatrixXf cholAbetainvt_X2tVinvX1 = chol_A_beta.transpose().triangularView<Lower>().solve(X1tVinvX2.transpose()); // b*b*a -> b x a
    
    MatrixXf A_alpha = RinvtX1.transpose() * RinvtX1 - cholAbetainvt_X2tVinvX1.transpose() * cholAbetainvt_X2tVinvX1;
    A_alpha.diagonal() += prior_prec_alpha;
    
    VectorXf Rinvsqy = chol_V.tsolve(y); // n*n*q -> n x 1;
    VectorXf XtRinvy = RinvtX2.transpose() * Rinvsqy; // b*n*1 >- b x 1
    VectorXf invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
    
    VectorXf X1t_V_beta_inv_y = RinvtX1.transpose() * Rinvsqy - cholAbetainvt_X2tVinvX1.transpose() * invSqAbXtRinvy;
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    
    y_tilde = y - X1 * alpha;
  }
  
  // Step 2 - sample Y_prec
  // We don't need to actually calculate V_beta^{-1} directly.
  VectorXf RinvSqy = chol_V.tsolve(y_tilde);
  VectorXf XtRinvy = RinvtX2.transpose() * RinvSqy;
  VectorXf prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
  float score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  VectorXf XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
  VectorXf beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
  beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}


VectorXf regression_sampler_v2(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& X1,           // nxa
    const MatrixXf& X2,           // nxm or nxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const General_Matrix_f& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXf& V,
    float Y_prec,                    // float
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
    const Ref<const VectorXf>& randn_e,
    const float rgamma_1,
    const float Y_prec_b0
){
  int n = y.size();
  int a = X1.cols();
  int b = X2.cols();
  
  // Check inputs
  if(X1.rows() != n) stop("Wrong dimension of X1");
  if(X2.rows() != n) stop("Wrong dimension of X2");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");
  
  // Calculate inverse of V_beta
  MatrixXf Dbeta_X2t = prior_prec_beta.cwiseInverse().asDiagonal() * X2.transpose();
  MatrixXf V_beta = X2 * Dbeta_X2t + V;
  LDLT<MatrixXf> V_beta_ldlt;
  V_beta_ldlt.compute(V_beta);
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    MatrixXf V_beta_inv_X1 = V_beta_ldlt.solve(X1);
    MatrixXf A_alpha = V_beta_inv_X1.transpose() * X1;
    A_alpha.diagonal() += prior_prec_alpha;
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    VectorXf X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    y_tilde = y - X1 * alpha;
  }
  
  // Step 2 - sample Y_prec
  VectorXf e2 = y_tilde.transpose() * V_beta_ldlt.solve(y_tilde);
  float score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  // what about prior mean?
  VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXf v = sqrt(Y_prec) * X2 * u;
  v += chol_V.crossprod(randn_e);
  // if(chol_V.isDense) {
  //   v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
  // } else{
  //   v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
  // }
  VectorXf w = V_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
  VectorXf beta = u + Dbeta_X2t * w / sqrt(Y_prec);
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}


VectorXf regression_sampler_v3(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& X1,           // nxa
    const MatrixXf& Ux,           // nxm or nxb
    const MatrixXf& Vx,           // mxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const General_Matrix_f& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXf& Vinv,
    const MatrixXf& VinvUx,
    const MatrixXf& UtVinvU,
    float Y_prec,                    // float
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
    const Ref<const VectorXf>& randn_e,
    const float rgamma_1,
    const float Y_prec_b0
){
  int n = y.size();
  int a = X1.cols();
  if(Vx.rows() != Ux.cols()) stop("Wrong dimensions of Vx");
  MatrixXf X2 = Ux*Vx;
  int b = X2.cols();
  
  // Check inputs
  if(X1.rows() != n) stop("Wrong dimension of X1");
  if(X2.rows() != n) stop("Wrong dimension of X2");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");
  if(randn_e.size() != n) stop("Wrong length of randn_e");
  
  // Calculate inverse of V_beta
  // Using Ainv - Ainv * Ux * (I + BVAinvU)inv * BVAinv in case B = VxDbeta_Vxt is singular
  MatrixXf Dbeta_Vxt = prior_prec_beta.cwiseInverse().asDiagonal() * Vx.transpose();
  MatrixXf VxDbeta_Vxt = Vx * Dbeta_Vxt;
  if(VinvUx.rows() != n) stop("Wrong dimensions of VinvUx");
  MatrixXf inner = VxDbeta_Vxt * UtVinvU;
  inner.diagonal().array() += 1.0;
  LDLT<MatrixXf> inner_ldlt;
  inner_ldlt.compute(inner);
  // Sigma_beta_inv = Vinv - VinvUx * inner.ldlt().solve(VxDbeta_Vxt * VinvUx.transpose());  // Don't actually calculate this. Stay in mxm space
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // MatrixXf V_beta_inv_X1 = Sigma_beta_inv * X1;
    // MatrixXf A_alpha = Y_prec * V_beta_inv_X1.transpose() * X1;
    MatrixXf Vinv_X1 = Vinv * X1;
    MatrixXf Uxt_Vinv_X1 = Ux.transpose() * Vinv_X1;
    MatrixXf A_alpha = X1.transpose() * Vinv_X1 - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt) * Uxt_Vinv_X1;
    A_alpha.diagonal() += prior_prec_alpha;
    
    // VectorXf X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
    VectorXf Uxt_Vinv_y = VinvUx.transpose() * y;
    VectorXf X1t_V_beta_inv_y = Vinv_X1.transpose() * y - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    
    y_tilde = y - X1 * alpha;
  }
  
  // Step 2 - sample Y_prec
  // VectorXf e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
  VectorXf Vinv_y_tilde = Vinv * y_tilde;
  VectorXf Uxt_Vinv_y = VinvUx.transpose() * y_tilde;
  VectorXf e2 = y_tilde.transpose() * Vinv_y_tilde - Uxt_Vinv_y.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
  
  float score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  // what about prior mean?
  VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXf v = std::sqrt(Y_prec) * X2 * u;
  v += chol_V.crossprod(randn_e);
  // if(chol_V.isDense) {
  //   v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
  // } else{
  //   v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
  // }
  // VectorXf X1 = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
  VectorXf e = y_tilde * std::sqrt(Y_prec) - v;
  VectorXf Uxt_Vinv_e = VinvUx.transpose() * e;
  VectorXf w = Vinv * e - VinvUx * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_e);
  
  VectorXf beta = u + Dbeta_Vxt * (Ux.transpose() * w) / std::sqrt(Y_prec); //b*b*1 + b*n*1
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}


// Single-site updater
// BayesC implementation
VectorXf regression_sampler_v4(
    const Ref<const VectorXf>& y_,           // nx1
    const MatrixXf& X1_,           // nxa
    const MatrixXf& X2,
    const ArrayXf& diag_X2tX2,
    const General_Matrix_f& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    VectorXf a,
    VectorXf alpha,
    VectorXf beta,
    VectorXi delta,
    float invVarRes,
    const Ref<const VectorXf>& invVarEffects, // bx1
    const Ref<const VectorXf>& pi,
    const Ref<const VectorXf>& randn_a,
    const Ref<const VectorXf>& randn_beta,
    const Ref<const VectorXf>& rand_unif,
    const Ref<const VectorXf>& rgamma_1,
    float Y_prec_b0,
    int nIter
) {
  VectorXf y = chol_V.tsolve(y_);
  MatrixXf X1 = chol_V.tsolve(X1_);
  ArrayXf diag_X1tX1 = X1.cwiseProduct(X1).colwise().sum();
// VectorXf regression_sampler_v4(
//   const VectorXf& y,
//   const MatrixXf& X1,
//   const MatrixXf& X2,
//   const ArrayXf& diag_X1tX1,
//   const ArrayXf& diag_X2tX2,
//   VectorXf a,
//   VectorXf alpha,
//   VectorXf beta,
//   VectorXi delta,
//   float invVarRes,
//   const Ref<const VectorXf>& invVarEffects, // bx1
//   const Ref<const VectorXf>& pi,
//   const Ref<const VectorXf>& randn_a,
//   const Ref<const VectorXf>& randn_beta,
//   const Ref<const VectorXf>& rand_unif,
//   const Ref<const VectorXf>& rgamma_1,
//   float Y_prec_b0,
//   int nIter
// ) {
  // const Ref<const VectorXf>& invVarEffects, // bx1
  // const Ref<const VectorXf>& pi,
  // const Ref<const VectorXf>& randn_a,
  // const Ref<const VectorXf>& randn_beta,
  // const Ref<const VectorXf>& rand_unif,
  // const Ref<const VectorXf>& rgamma_1,
  
  int nMarkers      = alpha.size();

  ArrayXf logPi = pi.array().log();
  ArrayXf logPiComp = (1.0 - pi.array()).log();
  ArrayXf logDelta0 = logPi;
  ArrayXf logVarEffects = invVarEffects.array().inverse().log();
  // // if(invVarEffects.size() != nMarkers) stop("Wrong length of invVarEffects");
  // // if(randn_a.size() != X1.cols()* nIter) stop("Wrong length of a");
  // // if(pi.size() != nMarkers) stop("Wrong length of pi");
  // // if(randn_beta.size() != nMarkers * nIter) stop("Wrong length of randn_beta");
  // // if(rand_unif.size() != nMarkers * nIter) stop("Wrong length of rand_unif");
  // // if(rgamma_1.size() != nIter) stop("Wrong length of rgamma_1");

  VectorXf yCorr = y - X1*a;
  for(int j = 0; j < nMarkers; j++) {
    if(delta[j] != 0) yCorr -= X2.col(j)*alpha[j];
  }

  for(int i = 0; i < nIter; i++) {

    // Sample a
    for(int j = 0; j < a.size(); j++) {
      float rhs = (X1.col(j).dot(yCorr) + diag_X1tX1[j]*a[j])*invVarRes;
      float lhs = diag_X1tX1[j]*invVarRes;
      float invLhs = 1.0/lhs;
      float gHat = rhs * invLhs;
      float old_a = a[j];
      a[j] = gHat + randn_a(i*a.size() + j)*sqrt(invLhs);
      yCorr += X1.col(j) * (old_a - a[j]);
    }

    int nLoci         = 0;
    // Sample beta = alpha*delta
    for(int j = 0; j < nMarkers; j++) {
      float rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j])*invVarRes;
      float lhs = diag_X2tX2[j]*invVarRes + invVarEffects[j];
      float invLhs = 1.0/lhs;
      float gHat = rhs * invLhs;
      float logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp[j];
      float probDelta1 = 1.0 / (1.0 + exp(logDelta0[j] - logDelta1));
      float oldAlpha = alpha[j];

      float u = rand_unif(i*nMarkers + j);
      float r = randn_beta(i*nMarkers + j);
      if(u < probDelta1) {
        delta[j] = 1.0;
        beta[j] = gHat + r*sqrt(invLhs);
        alpha[j] = beta[j];
        yCorr += X2.col(j) * (oldAlpha - alpha[j]);
        nLoci++;
      } else {
        if(oldAlpha != 0) {
          yCorr += X2.col(j) * oldAlpha;
        }
        delta[j] = 0;
        beta[j] = r/sqrt(invVarEffects[j]);
        alpha[j] = 0;
      }
    }
    invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + Y_prec_b0);
  }
  VectorXf result(1+a.size() + 3*nMarkers);
  result << invVarRes,a,alpha,beta,delta.cast<float>();
  // VectorXf result = VectorXf::Zero(1+a.size() + 3*nMarkers);
  // VectorXf result = VectorXf::Zero(1);
  return(result);
}

//' Draws samples from all ``fixed" coefficients (fixed and random) of a set of parallel linear regression models, conditional on the variance components.
//'
//' The model is either: \itemize{
//' \item y_i = X1_base*alpha1 + X1_list_[i]*alpha2 + X2*beta + e, e ~ N(0,1/Y_prec[i]*V)
//' \item y_i = X1_base*alpha1 + X1_list_[i]*alpha2 + X2*V_*beta + e, e ~ N(0,1/Y_prec[i]*V)
//' }
//' Where \code{V = RtR}, priors on elements of alpha1, alpha2 and beta are independent.
//' Each column of Y is considered independent
//'
//' @param which_sampler int: \itemize{
//'   \item 1: block sampler: b < n
//'   \item 2: block sampler: n >= b and X doesn't factorize
//'   \item 3: block sampler: n >= b, but X factorizes into UxVx where Ux is n x m and Vx = m x b, and m << n <= b
//' }
//' @param Y n x p matrix of observations
//' @param X1_base n x a1 matrix of X1 covariates common to all p. Can be NULL
//' @param X1_list_ p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
//' @param X2 either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
//' @param V_ m x b matrix if X2 is Ux, otherwise NULL
//' @param beta1_list_ p-list of a2-vectors for X1 coefficients. Can be NULL
//' @param beta2 a b x p matrix of current values for beta
//' @param h2s_index p-vector of indices for to select appropriate V of each trait
//' @param chol_V_list_ list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//' @param Y_prec p-vector of Y current precisions
//' @param Y_prec_a0,Y_prec_b0 scalars giving the shape and rate of the Gamma distribution for the prior on Y_prec
//' @param prior_prec_alpha1 a1 x p matrix of prior precisions for alpha1
//' @param prior_prec_alpha2 p-vector of precision of alpha2s for each trait
//' @param prior_mean_beta b x p matrix of prior means of beta
//' @param prior_prec_beta b x p matrix of prior precisions of beta
//' @param beta2_alpha_ b x p matrix for BayesC priors for beta2. Can be NULL
//' @param beta2_delta_ b x p matrix for BayesC priors for beta2. Can be NULL,
//' @param beta2_p_i_ b x p matrix for BayesC priors for beta2. Can be NULL
//' @return List with elements: \itemize{
//'   \item alpha1 a1 x p matrix of alpha1
//'   \item alpha2 concatenated vector of alpha2 for all traits
//'   \item beta b x p matrix of beta
//'   \item Y_prec p x 1 vector of Y_prec
//'   \item beta2_alpha b x p matrix (optional)
//'   \item beta2_delta_ b x p matrix (optional)
//' }
// [[Rcpp::export]]
Rcpp::List regression_sampler_parallel(
    int which_sampler,        //
    MatrixXf Y,               //
    MatrixXf X1_base,          //
    Rcpp::List X1_list_,             // p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
    SEXP X2_,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
    SEXP Vx_,                       // m x b matrix if X2 is Ux
    Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
    Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
    VectorXf Y_prec,               // p-vector of Y current precisions
    VectorXf Y_prec_a0,
    VectorXf Y_prec_b0,
    MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
    VectorXf prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
    MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
    MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
    SEXP current_alpha1s_,
    SEXP current_alpha2s_,
    Rcpp::List BayesAlphabet_parms
) {
  
  // MatrixXf beta2,            // b x p matrix of previous iteration regression coefficients
  // SEXP beta2_alpha_, // b x p matrix alpha parameters for BayesC priors for beta2
  // SEXP beta2_delta_, // b x p matrix delta parameters for BayesC priors for beta2
  // SEXP beta2_p_i_ // b x p matrix pi parameters for BayesC priors for beta2
  
  if(which_sampler > 4) stop("sampler not implemented");
  if(which_sampler == 3 && Rf_isNull(Vx_)) stop("Vx_ needed for sampler v3");
  int run_sampler_times = 1;
  
  int n = Y.rows();
  int p = Y.cols();
  
  // X1_base
  if(X1_base.rows() != n) stop("Wrong dimension of X1_base");
  int a1 = X1_base.cols();
  
  // X1_list
  std::vector<General_Matrix_f> X1_list;
  load_General_Matrix_f_list(X1_list_, X1_list, true);
  if(X1_list.size() > 0) {
    if(X1_list.size() != p) stop("Wrong length of X1_list");
  }
  
  // X2
  int b = 0;
  MatrixXf Ux,X2;
  MatrixXf Vx; // = MatrixXf::Zero(0,0);
  std::vector<General_Matrix_f> X2_list, Vx_list;
  List X2__, Vx__;
  if(Rf_isList(X2_)) {
    X2__ = as<List>(X2_);
    load_General_Matrix_f_list(X2__,X2_list, false);
  }
  if(Rf_isNull(Vx_)) {
    if(Rf_isList(X2_)) {
      for(int i = 0; i < X2__.size(); i++) {
        b += X2_list[i].dense.cols();
      }
    } else{
      X2 = as<MatrixXf>(X2_);
      b = X2.cols();
    }
  } else{
    Ux = as<MatrixXf>(X2_);
    if(Rf_isList(Vx_)) {
      Vx__ = as<List>(Vx_);
      load_General_Matrix_f_list(Vx__,Vx_list, false);
      for(int i = 0; i < Vx__.size(); i++) {
        b += Vx_list[i].dense.cols();
      }
    } else{
      Vx = as<MatrixXf>(Vx_);
      b = Vx.cols();
    }
  }
  
  // for BayesAlphabet
  MatrixXf current_alpha1s, betas_alpha, betas_beta, betas_pi;
  MatrixXi betas_delta;
  std::vector<General_Matrix_f> current_alpha2s;
  if(which_sampler == 4) {
    current_alpha1s = as<MatrixXf>(current_alpha1s_);
    if(Rf_isList(current_alpha2s_)) {
      List temp = as<List>(current_alpha2s_);
      load_General_Matrix_f_list(temp,current_alpha2s, false);
    }
    betas_alpha = as<MatrixXf>(BayesAlphabet_parms["alpha"]);
    betas_beta = as<MatrixXf>(BayesAlphabet_parms["beta"]);
    betas_pi = as<MatrixXf>(BayesAlphabet_parms["pi"]); 
    betas_delta = as<MatrixXi>(BayesAlphabet_parms["delta"]); 
    run_sampler_times = as<int>(BayesAlphabet_parms["run_sampler_times"]); 
  }
  
  // chol_V_list
  std::vector<General_Matrix_f> chol_V_list;
  load_General_Matrix_f_list(chol_V_list_, chol_V_list, true);
  if(max(h2s_index) > chol_V_list.size()) {
    stop("max(h2s_index) > length(chol_V_list)");
  }
  
  // priors
  if(Y_prec.size() != p) {
    stop("Wrong length of Y_prec");
  }
  if(Y_prec_a0.size() != p) {
    stop("Wrong length of Y_prec_a0");
  }
  if(Y_prec_b0.size() != p) {
    stop("Wrong length of Y_prec_b0");
  }
  if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
  if(X1_list.size() > 0 && prior_prec_alpha2.size() != p) {
    stop("Wrong length of prior_prec_alpha2");
  }
  if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
  if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
  
  // generate random numbers
  MatrixXf randn_alpha1 = rstdnorm_mat(a1*run_sampler_times,p);
  std::vector<VectorXf> randn_alpha2;
  if(X1_list.size() > 0){
    for(int i = 0; i < p; i++){
      randn_alpha2.push_back(rstdnorm_mat(X1_list[i].dense.cols()*run_sampler_times,1));
    }
  }
  MatrixXf randn_beta = rstdnorm_mat(b*run_sampler_times,p);
  MatrixXf randn_e;
  if(which_sampler == 2 || which_sampler == 3) {
    randn_e = rstdnorm_mat(n*run_sampler_times,p);
  }
  MatrixXf rgamma_1(run_sampler_times,p);
  for(int i = 0; i < p; i++) {
    rgamma_1.col(i) = as<VectorXf>(rgamma(run_sampler_times,Y_prec_a0[i] + n/2.0,1.0));
  }
  MatrixXf rand_unif;
  if(which_sampler > 3) {
    rand_unif = MatrixXf::Zero(b*run_sampler_times,p);
    for(int i = 0; i < p; i++) {
      rand_unif.col(i) = as<VectorXf>(runif(b*run_sampler_times,0.0,1.0));
    }
  }
  
  // Results structures
  MatrixXf alpha1(a1,p);
  std::vector<VectorXf> alpha2;
  alpha2.reserve(X1_list.size());
  int alpha2_size = 0;
  if(X1_list.size() > 0){
    for(int i = 0; i < X1_list.size(); i++){
      int a2 = X1_list[i].dense.cols();
      alpha2.push_back(VectorXf::Zero(a2));
      alpha2_size += a2;
    }
  }
  MatrixXf beta(b,p);
  if(which_sampler == 4) {
    beta = MatrixXf::Zero(3*b,p);
  }
  
  // go through h2s indices and sample columns with same index as a set
  for(int i = min(h2s_index); i <= max(h2s_index); i++) {
    int h2_index = i;
    VectorXi trait_set = as<VectorXi>(which(h2s_index == h2_index));  // list of traits with same h2_index
    
    if(trait_set.size() > 0){
      // prepare matrices for sampler
      MatrixXf RinvtX2, C, V, Vinv, VinvUx, UtVinvU;//, RinvtX1;
      ArrayXf diag_X1tVinvX1, diag_X2tVinvX2;
      General_Matrix_f chol_V = chol_V_list[h2_index - 1];
      if(which_sampler == 1) {
        // use regression_sampler_v1
        RinvtX2 = chol_V.tsolve(X2);
        // if(chol_V.isDense) {
        //   RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X2);
        // } else{
        //   RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X2);
        // }
        C = RinvtX2.transpose() * RinvtX2;
      }
      else if(which_sampler == 2) {
        // use regression_sampler_v2
        V = chol_V.crossprod(chol_V.dense);
        // if(chol_V.isDense) {
        //   V = chol_V.dense.transpose().triangularView<Lower>() * chol_V.dense;
        // } else{
        //   V = chol_V.sparse.transpose().triangularView<Lower>() * chol_V.sparse;
        // }
      } else if(which_sampler == 3) {
        // use regression_sampler_v3
        Vinv = chol_V.solve(chol_V.tsolve(MatrixXf::Identity(n,n)));
        VinvUx = chol_V.solve(chol_V.tsolve(Ux));
        // if(chol_V.isDense) {
        //   Vinv = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
        //   VinvUx = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(Ux));
        // } else{
        //   Vinv = chol_V.sparse.triangularView<Upper>().solve(chol_V.sparse.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
        //   VinvUx = chol_V.sparse.triangularView<Upper>().solve(chol_V.sparse.transpose().triangularView<Lower>().solve(Ux));
        //   // VinvUx = Vinv * U;
        //   // Rcout << i << std::endl;
        //   // Rcout << Vinv.diagonal().transpose() << std::endl;
        // }
        UtVinvU = Ux.transpose() * VinvUx;
      } else if(which_sampler == 4) {
        // use regression_sampler_v4
        RinvtX2 = chol_V.tsolve(X2);
        // if(chol_V.isDense) {
        //   RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X2);
        // } else{
        //   RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X2);
        // }
        diag_X2tVinvX2 = RinvtX2.cwiseProduct(RinvtX2).colwise().sum();
      }
      
      // REprintf("Number of threads=%i\\n", omp_get_max_threads());
      // Rcout <<trait_set.size() << " " << omp_get_max_threads() << std::endl;
      #pragma omp parallel for
      for(int i = 0; i < trait_set.size(); i++){
        int j = trait_set[i];
        MatrixXf X1;
        int a;
        int a2 = 0;
        int b;
        VectorXf prior_prec_alpha;
        VectorXf randn_alpha;
        if(X1_list.size() == 0) {
          X1 = X1_base;
          a = a1;
          prior_prec_alpha = prior_prec_alpha1.col(j);
          randn_alpha = randn_alpha1.col(j);
        } else{
          MatrixXf X12 = X1_list[j].dense;
          a2 = X12.cols();
          a = a1+a2;
          X1 = MatrixXf(n,a);
          X1 << X1_base,X12;
          prior_prec_alpha = VectorXf(a);
          prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
          prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
          randn_alpha = VectorXf(a);
          randn_alpha.head(a1) = randn_alpha1.col(j);
          randn_alpha.tail(a2) = randn_alpha2[j];
        }
        
        VectorXf samples;// = VectorXf::Zero(1+a1+a2+beta.rows());
        if(which_sampler == 1) {
          b = RinvtX2.cols();
          samples = regression_sampler_v1(Y.col(j), X1, RinvtX2, C, prior_prec_alpha, prior_mean_beta.col(j),
                                           prior_prec_beta.col(j), chol_V, Y_prec[j], randn_alpha,
                                           randn_beta.col(j), rgamma_1(0,j),Y_prec_b0[j]);
        } else if(which_sampler == 2) {
          b = X2.cols();
          samples = regression_sampler_v2(Y.col(j), X1, X2, prior_prec_alpha, prior_mean_beta.col(j),
                                           prior_prec_beta.col(j), chol_V, V, Y_prec[j], randn_alpha,
                                           randn_beta.col(j), randn_e.col(j),rgamma_1(0,j),Y_prec_b0[j]);
        } else if(which_sampler == 3) {
          b = Vx.cols();
          samples = regression_sampler_v3(Y.col(j), X1, Ux, Vx, prior_prec_alpha, prior_mean_beta.col(j),
                                           prior_prec_beta.col(j), chol_V, Vinv, VinvUx, UtVinvU, Y_prec[j], randn_alpha,
                                           randn_beta.col(j), randn_e.col(j),rgamma_1(0,j),Y_prec_b0[j]);
        } else if(which_sampler == 4) {
          b = 3*RinvtX2.cols();
          // // RinvtX1 = chol_V.tsolve(X1);  // This causes malloc error. Why? Or sometimes it seems using this below.
          // VectorXf Rinvty = chol_V.tsolve(Y.col(j));
          // if(chol_V.isDense) {
          //   RinvtX1 = chol_V.dense.transpose().triangularView<Lower>().solve(X1);
          //   // Rinvty = chol_V.dense.transpose().triangularView<Lower>().solve(Y.col(j));
          // } else{
          //   RinvtX1 = chol_V.sparse.transpose().triangularView<Lower>().solve(X1);
          //   // Rinvty = chol_V.sparse.transpose().triangularView<Lower>().solve(Y.col(j));
          // }
          // diag_X1tVinvX1 = RinvtX1.cwiseProduct(RinvtX1).colwise().sum();
          // samples = regression_sampler_v4(Rinvty, RinvtX1, RinvtX2, diag_X1tVinvX1, diag_X2tVinvX2,
          //                                 current_alpha1s.col(j),
          //                                 betas_alpha.col(j), betas_beta.col(j), betas_delta.col(j),
          //                                 Y_prec[j], prior_prec_beta.col(j),
          //                                 betas_pi.col(j),
          //                                 randn_alpha,randn_beta.col(j), rand_unif.col(j),rgamma_1.col(j), Y_prec_b0[j],
          //                                 run_sampler_times);
          samples = regression_sampler_v4(Y.col(j), X1, RinvtX2, diag_X2tVinvX2,
                                          chol_V,
                                          current_alpha1s.col(j),
                                          betas_alpha.col(j), betas_beta.col(j), betas_delta.col(j),
                                          Y_prec[j], prior_prec_beta.col(j),
                                          betas_pi.col(j),
                                          randn_alpha,randn_beta.col(j), rand_unif.col(j),rgamma_1.col(j), Y_prec_b0[j],
                                          run_sampler_times);
        } else {
          stop("sampler not implemented");
        }
        if(samples.size() != 1 + a1 + a2 + b) stop("wrong length of samples");
        // samples = VectorXf::Ones(1+a1+a2+b);
        
        // float prec = samples[0];
        // extract samples
        Y_prec[j] = samples[0];
        if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
        if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
        if(b > 0) beta.col(j) = samples.segment(1+a1+a2,b);
      }
    }
  }
  
  // collect alpha2 into a vector
  VectorXf alpha2_vec(alpha2_size);
  if(X1_list.size() > 0){
    int index = 0;
    for(int i = 0; i < X1_list.size(); i++){
      alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
      index += alpha2[i].size();
    }
  }
  
  Rcpp::List result = Rcpp::List::create(
    Named("alpha1") = alpha1,
    Named("alpha2") = alpha2_vec,
    Named("beta") = beta,
    Named("Y_prec") = Y_prec
  );
  
  return(result);
}


// -------------------------------------------- //
// ------------ sample_MME_ZKZts -------------- //
// -------------------------------------------- //

// Samples from model:
// y = Zu + e
// u ~ N(0,K); solve(K) = t(chol_K_inv) %*% chol_K_inv
// e[i] ~ N(0,1/tot_Eta_prec)
// C = ZtRinvZ + diag(Kinv)
//// [[Rcpp::export()]]
// VectorXf sample_MME_single_diagR(
//     VectorXf y,           // nx1
//     General_Matrix_f Z,    // nxr dgCMatrix or dense
//     SpMat chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
//     float tot_Eta_prec,   // float
//     float pe,            // float
//     VectorXf randn_theta  // rx1
// ){
//   VectorXf b;
//   b = Z.crossprod(y)*pe;
//   // if(Rf_isMatrix(Z_)) {
//   //   MatrixXf Z = as<MatrixXf>(Z_);
//   //   b = Z.transpose() * y * pe;
//   // } else{
//   //   SpMat Z = as<SpMat>(Z_);
//   //   b = Z.transpose() * y * pe;
//   // }
//   // if(Z.isDense) {
//   //   b = Z.dense.transpose() * y * pe;
//   // } else{
//   //   b = Z.sparse.transpose() * y * pe;
//   // }
//   b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b / sqrt(tot_Eta_prec));
//   b += randn_theta;
//   b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
//   return(b);
// }
// 
// 
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// // E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// // For complete data, ie no missing obs.
// // [[Rcpp::export()]]
// MatrixXf sample_MME_ZKZts_c(
//     MatrixXf Y,                    // nxp
//     SEXP Z_,
//     VectorXf tot_Eta_prec,         // px1
//     Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
//     MatrixXf h2s,                  // n_RE x p
//     VectorXi h2s_index                 // px1
//     ) {
// 
//   General_Matrix_f Z = load_General_Matrix_f(Z_,false);
// 
//   int p = Y.cols();
//   int r = Z.cols();
//   // if(Z.isDense) {
//   //   r = Z.dense.cols();
//   // } else{
//   //   r = Z.sparse.cols();
//   // }
// 
//   MatrixXf randn_theta = rstdnorm_mat(r,p);
// 
//   std::vector<General_Matrix_f> chol_ZtZ_Kinv_list;
//   load_General_Matrix_f_list(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list, true);
// 
//   MatrixXf U(r,p);
//   ArrayXf h2_e = 1.0 - h2s.colwise().sum().array();
//   ArrayXf pes = tot_Eta_prec.array() / h2_e.array();
// 
//   #pragma omp parallel for
//   for(std::size_t j = 0; j < p; j++){
//     int h2_index = h2s_index[j] - 1;
//     // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
//     U.col(j) = sample_MME_single_diagR(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
//   }
// 
//   return(U);
// }
VectorXd sample_MME_single_diagR(
    VectorXd y,           // nx1
    MSpMatd Z,    // nxr dgCMatrix or dense
    MSpMatd chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double tot_Eta_prec,   // float
    double pe,            // float
    VectorXd randn_theta  // rx1
){
  VectorXd b;
  b = Z.transpose()* y *pe;
  // if(Rf_isMatrix(Z_)) {
  //   MatrixXf Z = as<MatrixXf>(Z_);
  //   b = Z.transpose() * y * pe;
  // } else{
  //   SpMat Z = as<SpMat>(Z_);
  //   b = Z.transpose() * y * pe;
  // }
  // if(Z.isDense) {
  //   b = Z.dense.transpose() * y * pe;
  // } else{
  //   b = Z.sparse.transpose() * y * pe;
  // }
  b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b / sqrt(tot_Eta_prec));
  b += randn_theta;
  b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
  return(b);
}

// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c(
    Map<MatrixXd> Y,                    // nxp
    MSpMatd Z,
    Map<VectorXd> tot_Eta_prec,         // px1
    Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
    Map<MatrixXd> h2s,                  // n_RE x p
    VectorXi h2s_index                 // px1
) {
  
  // General_Matrix_fq Z = load_General_Matrix_fq2(Z_,false);
  
  int p = Y.cols();
  int r = Z.cols();
  // if(Z.isDense) {
  //   r = Z.dense.cols();
  // } else{
  //   r = Z.sparse.cols();
  // }
  
  // MatrixXd randn_theta = rstdnorm_mat2(r,p).cast<double>();
  MatrixXd randn_theta = rstdnorm_mat_d(r,p);
  
  std::vector<MSpMatd> chol_ZtZ_Kinv_list;
  int m = chol_ZtZ_Kinv_list_.size();
  chol_ZtZ_Kinv_list.reserve(m);
  for(int i = 0; i < m; i++){
    SEXP Xi_ = chol_ZtZ_Kinv_list_[i];
    chol_ZtZ_Kinv_list.push_back(as<MSpMatd> (Xi_));
  }
  // load_General_Matrix_fq_list2(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list, true);
  
  MatrixXd U(r,p);
  ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
  ArrayXd pes = tot_Eta_prec.array() / h2_e.array();
  
#pragma omp parallel for
  for(std::size_t j = 0; j < p; j++){
    int h2_index = h2s_index[j] - 1;
    // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
    if(h2_e[h2_index] == 1.0) {
      U.col(j) = VectorXd::Zero(r);
    } else{
      U.col(j) = sample_MME_single_diagR(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index], tot_Eta_prec[j], pes[j],randn_theta.col(j));
    }
  }
  
  return(U);
}

// -------------------------------------------- //
// ---------------- sample h2s ---------------- //
// -------------------------------------------- //

// [[Rcpp::export()]]
MatrixXf log_p_h2s(
    MatrixXf Y,              // nxp
    VectorXf tot_Eta_prec,   // px1
    Rcpp::List chol_V_list_,       // List. Each element contains: R st RtR = V. may be sparse or dense
    VectorXf discrete_priors // n_h2 x 1
    )
{
  int b = discrete_priors.size();
  int p = Y.cols();
  int n = Y.rows();

  std::vector<General_Matrix_f> chol_V_list;
  load_General_Matrix_f_list(chol_V_list_, chol_V_list, true);

  MatrixXf log_ps(b,p);

  #pragma omp parallel for
  for(std::size_t i = 0; i < b; i++){
    General_Matrix_f chol_R = chol_V_list[i];
    MatrixXf y_std(n,b);
    float log_det_V;
    y_std = chol_R.tsolve(Y);
    log_det_V = 2*chol_R.get_log_det();
    // if(chol_R.isDense){
    //   y_std = chol_R.dense.transpose().triangularView<Lower>().solve(Y);
    //   log_det_V = 2*chol_R.dense.diagonal().array().log().sum();
    // } else{
    //   y_std = chol_R.sparse.transpose().triangularView<Lower>().solve(Y);
    //   log_det_V = 0;
    //   for(int j = 0; j < chol_R.sparse.rows(); j++) {
    //     log_det_V += 2*std::log(chol_R.sparse.coeffRef(j,j));
    //   }
    // }
    VectorXf scores2 = (y_std.transpose() * y_std).diagonal().array() * tot_Eta_prec.array();
    log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_V - n*tot_Eta_prec.array().log()) -
      0.5 * scores2.array() + log(discrete_priors[i]));
  }

    return(log_ps);
}

// [[Rcpp::export()]]
VectorXi sample_h2s(
    ArrayXXd log_ps
)
{

  int p = log_ps.cols();
  int b = log_ps.rows();
  VectorXd rs = as<VectorXd>(runif(p));

  VectorXi h2s_index(p);

  #pragma omp parallel for
  for(std::size_t j = 0; j < p; j++){
    // for(int j = 0; j < p; j++){
    double max_col = log_ps.col(j).maxCoeff();
    double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
    VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
    h2s_index[j] = 1;
    double cumsum = 0;
    for(int i = 0; i < b; i++){
      cumsum += ps_j[i];
      if(rs[j] > cumsum) {
        h2s_index[j] ++;
      }
    }
    if(h2s_index[j] > p) h2s_index[j] = p;
  }

  return h2s_index; // 1-based index
}


// -------------------------------------------- //
// --------------- sample h2s MH -------------- //
// -------------------------------------------- //

float log_prob_h2_c(
    const Ref<const VectorXf>& y,           // nx1
    General_Matrix_f chol_R,     // nxn upper-triangular. Dense or sparse
    int n,                // int
    float tot_Eta_prec,  // float
    float discrete_prior // float
){
  VectorXf y_std;
  float log_det_V;
  y_std = chol_R.tsolve(y);
  log_det_V = 2*chol_R.get_log_det();
  // if(chol_R.isDense){
  //   y_std = chol_R.dense.transpose().triangularView<Lower>().solve(y);
  //   log_det_V = 2*chol_R.dense.diagonal().array().log().sum();
  // } else{
  //   y_std = chol_R.sparse.transpose().triangularView<Lower>().solve(y);
  //   log_det_V = 0;
  //   for(int j = 0; j < chol_R.sparse.rows(); j++) {
  //     log_det_V += 2*std::log(chol_R.sparse.coeffRef(j,j));
  //   }
  // }
  float score2 = tot_Eta_prec * y_std.dot(y_std);

  float log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_V - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_c(
    MatrixXf Y,                // nxp
    VectorXf tot_Eta_prec,     // px1
    VectorXf discrete_priors,  // n_h2 x 1
    VectorXi h2_index,              // px1
    MatrixXf h2s_matrix,       // n_RE x n_h2
    Rcpp::List chol_V_list_,         // List of R st RtR = V, can be dense or sparse
    float step_size                // float
){

  int p = Y.cols();
  int n = Y.rows();

  VectorXf r_draws = as<VectorXf>(runif(p));
  VectorXf state_draws = as<VectorXf>(runif(p));

  VectorXi new_index(p);

  std::vector<General_Matrix_f> chol_V_list;
  load_General_Matrix_f_list(chol_V_list_, chol_V_list, true);


  #pragma omp parallel for
  for(std::size_t j = 0; j < p; j++){
    int old_state = h2_index[j] - 1;
    VectorXf candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
    int r = state_draws[j] * (candidate_new_states.size());
    int proposed_state = candidate_new_states[r];

    if(discrete_priors[proposed_state] == 0.0) {
      new_index[j] = old_state;  // don't bother with calculations if prior == 0.0
    } else{
      float old_log_p = log_prob_h2_c(Y.col(j),chol_V_list[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);
      float new_log_p = log_prob_h2_c(Y.col(j),chol_V_list[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

      VectorXf candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);

      float forward_prob = 1.0 / candidate_new_states.size();
      float back_prob = 1.0 / candidate_states_from_new_state.size();

      float log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);

      if(log(r_draws[j]) < log_MH_ratio) {
        new_index[j] = proposed_state;
      } else {
        new_index[j] = old_state;
      }
    }
    new_index[j] += 1;  // convert to 1-based index
  }

    return new_index;
}


// -------------------------------------------------- //
// -- Sample factor scores --- //
// -------------------------------------------------- //

// Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
// phenotype residuals
// Y - ZU = F * Lambda + E
// F[,k] ~ N(XFBF[,k] + ZUF[,k],1/F_e_prec[,j])
// E[,j] ~ N(0,1/resid_Eta_prec[,j])
// Sampling is done separately for each block of rows with the same pattern of missing observations
// [[Rcpp::export()]]
MatrixXf sample_factors_scores_c( // returns nxk matrix
    MatrixXf Eta_tilde,      // nxp
    MatrixXf prior_mean,     // nxk
    MatrixXf Lambda,         // kxp
    VectorXf resid_Eta_prec, // px1
    VectorXf F_e_prec        // kx1
) {
  int n = Eta_tilde.rows();
  int k = Lambda.rows();
  MatrixXf randn_draws = rstdnorm_mat(n,k);

  MatrixXf Lmsg = Lambda * resid_Eta_prec.asDiagonal();
  MatrixXf Sigma = Lmsg * Lambda.transpose();
  Sigma.diagonal() += F_e_prec;
  Eigen::LLT<MatrixXf> chol_Sigma;
  chol_Sigma.compute(Sigma);
  MatrixXf R = chol_Sigma.matrixU();

  MatrixXf Meta = R.transpose().triangularView<Lower>().solve((Eta_tilde * Lmsg.transpose() + prior_mean * F_e_prec.asDiagonal()).transpose());

  MatrixXf Ft = R.triangularView<Upper>().solve(Meta + randn_draws.transpose());

  return Ft.transpose();
}



// -------------------------------------------------- //
// -- Sample tau2 and delta scores --- //
// -------------------------------------------------- //

VectorXd cumprod(const VectorXd& x) {
  int n = x.size();
  VectorXd res(n);
  res[0] = x[0];
  if(n > 1) {
    for(int i = 1; i < n; i++){
      res[i] = res[i-1]*x[i];
    }
  }
  return(res);
}

// // [[Rcpp::export()]]
// Rcpp::List sample_tau2_delta_c_Eigen_v2(
//     float tau2,
//     float xi,
//     VectorXf delta,
//     VectorXf scores,
//     float tau_0,
//     float delta_shape,
//     float delta_rate,
//     int p,
//     int times
// ) {
//
//   int K = scores.size();
//   if(delta.size() != K) stop("Wrong size of delta");
//   float shape;
//   float scale;
//   VectorXf cumprod_delta = cumprod(delta);
//   for(int i = 0; i < times; i++){
//     // sample tau2
//     shape = (p*K + 1)/2.0;
//     scale = 1.0/xi + cumprod_delta.dot(scores);
//     tau2 = 1.0/R::rgamma(shape,1.0/scale);
//
//     // sample xi
//     shape = 1.0;
//     scale = 1.0/(tau_0*tau_0) + 1.0/tau2;
//     xi = 1.0/R::rgamma(shape,1.0/scale);
//
//     for(int h = 1; h < K; h++) {
//       // delta_h
//       shape = delta_shape + p*(K-h)/2.0;
//       scale = delta_rate + cumprod_delta.tail(K-h).dot(scores.tail(K-h)) / (tau2 * delta(h));
//       delta[h] = R::rgamma(shape,1.0/scale);
//       cumprod_delta = cumprod(delta);
//     }
//   }
//
//   return(Rcpp::List::create(Named("tau2") = tau2,
//                             Named("xi") = xi,
//                             Named("delta") = delta));
// }
// [[Rcpp::export()]]
Rcpp::List sample_tau2_delta_c_Eigen_v2(
    double tau2,
    double xi,
    VectorXd delta,
    VectorXd scores,
    double tau_0,
    double delta_shape,
    double delta_scale,  // shape and scale for inverse gamma distribution (shape and rate for Gamma)
    int p,
    int times
) {
  
  int K = scores.size();
  if(delta.size() != K) stop("Wrong size of delta");
  double shape;
  double scale;
  VectorXd cumprod_delta = cumprod(delta);
  for(int i = 0; i < times; i++){
    // sample tau2
    shape = (p*K + 1)/2.0;
    scale = 1.0/xi + cumprod_delta.cwiseInverse().dot(scores);
    tau2 = 1.0/R::rgamma(shape,1.0/scale);
    
    // sample xi
    shape = 1.0;
    scale = 1.0/(tau_0*tau_0) + 1.0/tau2;
    xi = 1.0/R::rgamma(shape,1.0/scale);
    
    for(int h = 1; h < K; h++) {
      // delta_h
      shape = delta_shape + p*(K-h)/2.0;
      scale = delta_scale + delta(h) * cumprod_delta.tail(K-h).cwiseInverse().dot(scores.tail(K-h)) / tau2;
      delta[h] = 1.0/R::rgamma(shape,1.0/scale);
      cumprod_delta = cumprod(delta);
    }
  }
  
  return(Rcpp::List::create(Named("tau2") = tau2,
                            Named("xi") = xi,
                            Named("delta") = delta));
}

// This didn't work with floats - must have been overflow issues. Maybe use logs instead?
// [[Rcpp::export()]]
VectorXd sample_trunc_delta_c_Eigen(
    VectorXd delta,
    VectorXd tauh,
    VectorXd scores,
    VectorXd shapes,
    double delta_1_rate,
    double delta_2_rate,
    MatrixXd randu_draws,
    double trunc_point
) {
  int times = randu_draws.rows();
  int k = tauh.size();
  double p,u;
  
  double rate,delta_old;
  for(int i = 0; i < times; i++){
    delta_old = delta(0);
    rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
    u = randu_draws(i,0); // don't truncate delta(0)
    delta(0) = R::qgamma(u,shapes(0),1.0/rate,1,0);
    // tauh = cumprod(delta);
    tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
    
    for(int h = 1; h < k; h++) {
      delta_old = delta(h);
      rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
      p = R::pgamma(trunc_point,shapes(h),1.0/rate,1,0);  // left-tuncate delta(h) at trunc_point
      if(p > 0.999) {
        // p = 0.999;  // prevent over-flow.
        delta(h) = trunc_point;
      } else {
        u = p + (1.0-p)*randu_draws(i,h);
        delta(h) = R::qgamma(u,shapes(h),1.0/rate,1,0);
      }
      // tauh = cumprod(delta);
      tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
      // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
    }
  }
  return(delta);
}



/// old functions to be replaced

// [[Rcpp::export()]]
VectorXf sample_MME_single_diagK(  // returns b x 1 vector
    VectorXf y,           // nx1
    MatrixXf X,           // nxb
    VectorXf prior_mean,  // bx1
    VectorXf prior_prec,  // bx1
    SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    float tot_Eta_prec, // float
    VectorXf randn_theta, // bx1
    VectorXf randn_e      // 0x1 or nx1. 0x1 if b<n
){
  if(randn_e.size() == 0){
    MatrixXf RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
    VectorXf XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    VectorXf XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
    MatrixXf C = RinvSqX.transpose() * RinvSqX;
    C.diagonal() += prior_prec;
    LLT<MatrixXf> C_llt;
    C_llt.compute(C);
    MatrixXf chol_C = C_llt.matrixU();

    VectorXf b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
    b += randn_theta;
    b = chol_C.triangularView<Upper>().solve(b);
    return(b);
  } else {
    // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
    MatrixXf Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
    VectorXf alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));

    VectorXf u = randn_theta.array() / prior_prec.cwiseSqrt().array();
    u += prior_mean;
    VectorXf v = Phi * u + randn_e;
    VectorXf alpha_v = alpha-v;

    MatrixXf D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
    MatrixXf cov = Phi * D_PhiT;
    cov.diagonal().array() += 1.0;

    VectorXf w = cov.ldlt().solve(alpha_v);

    VectorXf theta = u + D_PhiT * w;

    return(theta);
  }
}

// [[Rcpp::export()]]
MatrixXf sample_coefs_set_c(    // return pxn matrix
    Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x p), nonZero_cols_X.
    MatrixXf prior_mean,   // pxn
    MatrixXf prior_prec    // pxn
  ){

  int n = model_matrices.size();
  int p = prior_mean.rows();

  std::vector<MatrixXf> y_list;
  std::vector<MatrixXf> X_list;
  std::vector<MatrixXf> tot_Y_prec_list;
  std::vector<ArrayXi> nonZero_cols_X;
  std::vector<MatrixXf> randn_theta_list;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    y_list.push_back(Rcpp::as<MatrixXf>(model_matrix_i["y"]));                    // matrix of observations (n_i x t)
    X_list.push_back(Rcpp::as<MatrixXf>(model_matrix_i["X"]));                    // design matrix (n_i x b_X), only including columns that are non-zero
    tot_Y_prec_list.push_back(Rcpp::as<MatrixXf>(model_matrix_i["tot_Y_prec"]));                    // matrix of residual precisions (n_i x t)
    nonZero_cols_X.push_back(Rcpp::as<ArrayXi>(model_matrix_i["nonZero_cols_X"])); // list of which columns b_X correspond to in full X matrix

    int t = y_list[i].cols();
    MatrixXf randn_theta = rstdnorm_mat(p/t,t);
    randn_theta_list.push_back(randn_theta);
  }

  int n_traits = y_list[0].cols();

  MatrixXf coefs(p,n);


  #pragma omp parallel for
  for(std::size_t j = 0; j < n; j++){
    MatrixXf Y = y_list[j];
    MatrixXf X = X_list[j];
    MatrixXf tot_Y_prec = tot_Y_prec_list[j];
    int b = randn_theta_list[j].rows();
    int b_X = X.cols();
    MatrixXf randn_e = MatrixXf::Zero(0,n_traits);
    for(int t = 0; t < n_traits; t++) {
      // Create cholesky decomposition of the residual variance matrix using tot_Y_prec.
      SpMat chol_R = tot_Y_prec.col(t).cwiseSqrt().cwiseInverse().asDiagonal().toDenseMatrix().sparseView();
      // MSpMat chol_R(Rsqrt.rows(),Rsqrt.cols(), Rsqrt.nonZeros(),Rsqrt.outerIndexPtr(),Rsqrt.innerIndexPtr(),Rsqrt.valuePtr());

      // first assign the result vector to prior_mean + randn/sqrt(prec)
      // will then replace values with sampled values.
      VectorXf prior_mean_tj = prior_mean.block(t*b,j,b,1);
      VectorXf prior_prec_tj = prior_prec.block(t*b,j,b,1);
      VectorXf randn_theta_tj = randn_theta_list[j].col(t);
      coefs.block(t*b,j,b,1) = prior_mean_tj.array() + randn_theta_tj.array() / prior_prec_tj.array().sqrt();

      // now, pull out parameters for the coefficients corresponding to the columns of X
      VectorXf prior_mean_tj_X(b_X);
      VectorXf prior_prec_tj_X(b_X);
      VectorXf randn_theta_tj_X(b_X);
      for(int k = 0; k < b_X; k++){
        int element = nonZero_cols_X[j][k]-1;
        prior_mean_tj_X.coeffRef(k) = prior_mean_tj.coeffRef(element);
        prior_prec_tj_X.coeffRef(k) = prior_prec_tj.coeffRef(element);
        randn_theta_tj_X.coeffRef(k) = randn_theta_tj.coeffRef(element);
      }
      VectorXf coefs_X = sample_MME_single_diagK(Y.col(t), X,
                                                 prior_mean_tj_X, prior_prec_tj_X,
                                                 chol_R,1.0,
                                                 randn_theta_tj_X,randn_e.col(t));

      // now replace the values in coef with the corresponding ones in coefs_X
      for(int k = 0; k < b_X; k++){
        int element = nonZero_cols_X[j][k]-1;
        coefs.coeffRef(t*b + element,j) = coefs_X.coeffRef(k);
      }
    }
  }

  return(coefs);
}

// [[Rcpp::export()]]
MatrixXf get_fitted_set_c(  // returns n_tot x p matrix in same order as data
    Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x b), s (n_i x 1), position (n_i x 1)
    MatrixXf coefs  // b x n matrix
    ){

  std::vector<MatrixXf> X_list;
  std::vector<ArrayXi> position_list;
  std::vector<ArrayXi> nonZero_cols_X;
  int total_obs = 0;
  int n = model_matrices.size();
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    X_list.push_back(Rcpp::as<MatrixXf>(model_matrix_i["X"]));
    position_list.push_back(Rcpp::as<ArrayXi>(model_matrix_i["position"]));
    nonZero_cols_X.push_back(Rcpp::as<ArrayXi>(model_matrix_i["nonZero_cols_X"])); // list of which columns b_X correspond to in full X matrix

    int n_obs = X_list[i].rows();
    total_obs += n_obs;
  }

  Rcpp::List model_matrix_1 = Rcpp::as<Rcpp::List>(model_matrices[0]);
  MatrixXf Y = Rcpp::as<MatrixXf>(model_matrix_1["y"]);
  int n_traits = Y.cols();

  MatrixXf Y_fitted(total_obs,n_traits);

  #pragma omp parallel for
  for(std::size_t j = 0; j < n; j++){
    int b_X = X_list[j].cols();
    int b_tot = coefs.rows() / n_traits;
    VectorXf coefs_j(b_X * n_traits);
    for(int t = 0; t < n_traits; t++){
      for(int k = 0; k < b_X; k++){
        coefs_j[t*b_X+k] = coefs.coeffRef(t*b_tot+nonZero_cols_X[j][k]-1,j);
      }
    }
    Map<MatrixXf> Eta_i(coefs_j.data(),b_X,n_traits);
    MatrixXf Y_fitted_j = X_list[j] * Eta_i;
    for(int i = 0; i < position_list[j].size(); i++) {
      Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
    }
  }

  return(Y_fitted);
}


