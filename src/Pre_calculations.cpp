// Copyright 2020 Daniel Runcie
// Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
// that can be found in the LICENSE file and available at
// https://polyformproject.org/licenses/noncommercial/1.0.0/

#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export()]]
List LDLt(SEXP A_) {
  if(Rf_isMatrix(A_)){
    MatrixXd A = as<MatrixXd >(A_);
    Eigen::LDLT<MatrixXd> ldlt_A;
    ldlt_A.compute(A);
    MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
    MatrixXd P = ldlt_A.transpositionsP() * I;
    VectorXd d = ldlt_A.vectorD();
    MatrixXd L = ldlt_A.matrixL();
    SpMatd Lsp = L.sparseView();
    if(static_cast<float>(Lsp.nonZeros()) / I.size() > 0.25) {
      return(List::create(
          Named("P") = P.sparseView(),
          Named("L") = L,
          Named("d") = d
      ));
    } else{
      return(List::create(
          Named("P") = P.sparseView(),
          Named("L") = Lsp,
          Named("d") = d
      ));
    }
  } else{
    SpMatd A = as<SpMatd>(A_);
    Eigen::SimplicialLDLT<SpMatd> ldlt_A;
    ldlt_A.compute(A);
    if(ldlt_A.info() != Eigen::Success) {
      // LDLt failed? Try again as a dense matrix
      MatrixXd Ad = A.toDense();
      return(LDLt(wrap(Ad)));
    } else {
      MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
      MatrixXd P = ldlt_A.permutationP() * I;
      return(List::create(
          Named("P") = P.sparseView(),
          Named("L") = ldlt_A.matrixL(),
          Named("d") = ldlt_A.vectorD()
      ));
    }
  }
}


// SpMat make_chol_K_inv(const std::vector<General_Matrix_f>& chol_Ki_mats, VectorXf h2s,float tol){
//   int h = h2s.size();
//   VectorXf sizes(h);
//   int total_size = 0;
//   for(int i = 0; i < h; i++){
//     if(chol_Ki_mats[i].isDense) {
//       sizes(i) = chol_Ki_mats[i].dense.rows();
//     } else{
//       sizes(i) = chol_Ki_mats[i].sparse.rows();
//     }
//     total_size += sizes(i);
//   }
//   MatrixXf chol_K_inv_dense(total_size,total_size);
//   chol_K_inv_dense.setZero();
//   int curr_row = 0;
//   int curr_col = 0;
//   for(int i = 0; i < h; i++){
//     if(h2s[i] == 0) {
//       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal().setOnes();
//       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal() /= 0;
//     } else{
//       // SpMat chol_Ki = chol_Ki_mats[i].sparse;
//       if(chol_Ki_mats[i].isDense) {
//         chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki_mats[i].dense;
//       } else{
//         chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki_mats[i].sparse;
//       }
//       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) /= sqrt(h2s[i]);
//     }
//     curr_row += sizes(i);
//     curr_col += sizes(i);
//   }
//   SpMat chol_K_inv = chol_K_inv_dense.sparseView(0,tol);
//   return chol_K_inv;
// }
// 
// // [[Rcpp::export()]]
// Rcpp::List make_chol_ZtZ_Kinv_list(Rcpp::List chol_Ki_mats_,
//                                    MatrixXf h2s_matrix,
//                                    SpMat ZtZ,
//                                    float drop0_tol,
//                                    SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
//                                    int ncores) {
//   int s = h2s_matrix.cols();
// 
// 
//   std::vector<General_Matrix_f> chol_Ki_mats;
//   load_General_Matrix_f_list(chol_Ki_mats_,chol_Ki_mats,true);
// 
//   std::vector<SpMat> chol_ZtZ_Kinv_list;
//   chol_ZtZ_Kinv_list.reserve(s);
//   for(int i = 0; i < s; i++){
//     chol_ZtZ_Kinv_list.push_back(SpMat(0,0));
//   }
// 
//   int n_groups = s/ncores;
//   for(int i = 0; i <= n_groups; i++){
//     int start = ncores*i;
//     int end = std::min(static_cast<double>(s),ncores*(i+1.0));
// 
//     for(std::size_t j = start; j < end; j++){
//       VectorXf h2s = h2s_matrix.col(j);
//       SpMat chol_K_inv = make_chol_K_inv(chol_Ki_mats,h2s,drop0_tol);
//       MatrixXf ZtZ_Kinv = 1.0/(1.0 - h2s.sum()) * ZtZ + chol_K_inv.transpose() * chol_K_inv;
//       Eigen::LLT<MatrixXf> chol_ZtZ_Kinv(ZtZ_Kinv);
//       MatrixXf chol_ZtZ_Kinv_R = chol_ZtZ_Kinv.matrixU();
//       chol_ZtZ_Kinv_list[j] = chol_ZtZ_Kinv_R.sparseView(0,drop0_tol);
//     }
//     int pb_state = as<int>(getTxtProgressBar(pb));
//     setTxtProgressBar(pb,pb_state+(end-start));
//   }
// 
//   Rcpp::List chol_ZtZ_Kinv_list_out(s);
//   for(int i = 0; i < s; i++){
//     chol_ZtZ_Kinv_list_out[i] = chol_ZtZ_Kinv_list[i];
//   }
// 
//   return(chol_ZtZ_Kinv_list_out);
// }

// Switching out for dense matrices for these calculations. The sparse ones don't seem accurate enough for these calculations.
// Maybe could try LDLt instead of LLt?
struct R_matrix {
  Map<MatrixXd> dense;
  MSpMatd sparse;
  bool isDense;
  R_matrix(Map<MatrixXd> dense_, MSpMatd sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
};

R_matrix load_R_matrix(SEXP X_) {
  MatrixXd null_d = MatrixXd::Zero(0,0);
  if(Rf_isMatrix(X_)){
    Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
    SpMatd null_s = null_d.sparseView();
    MSpMatd M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
    R_matrix Xm(X,M_null_s,true);
    return(Xm);
  } else{
    MSpMatd X = as<MSpMatd>(X_);
    Map<MatrixXd> M_null_d(null_d.data(),0,0);
    R_matrix Xm(M_null_d,X,false);
    return(Xm);
  }
}

// Code to convert list of R matrices (sparse or dense) into a thread-safe object
//
// @param X_list List of matrices (each can be dgCMatrix or matrix)
// @param X_vector \code{std::vector} of \code{R_matrix} type which will be populated
void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix>& X_vector){
  // null_matrices
  MatrixXd null_d = MatrixXd::Zero(0,0);
  Map<MatrixXd> M_null_d(null_d.data(),0,0);
  SpMatd null_s = null_d.sparseView();
  MSpMatd M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
  
  int p = X_list.size();
  X_vector.reserve(p);
  for(int i = 0; i < p; i++){
    SEXP Xi_ = X_list[i];
    X_vector.push_back(load_R_matrix(Xi_));
    // if(Rf_isMatrix(Xi_)){
    //   Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
    //   R_matrix Xim(Xi,M_null_s,true);
    //   X_vector.push_back(Xim);
    // } else{
    //   MSpMatd Xi = as<MSpMatd>(Xi_);
    //   R_matrix Xim(M_null_d,Xi,false);
    //   X_vector.push_back(Xim);
    // }
  }
}


SpMatd make_chol_K_inv(const std::vector<R_matrix>& chol_Ki_mats, VectorXd h2s,double tol){
  int h = h2s.size();
  VectorXd sizes(h);
  int total_size = 0;
  for(int i = 0; i < h; i++){
    if(chol_Ki_mats[i].isDense) {
      sizes(i) = chol_Ki_mats[i].dense.rows();
    } else{
      sizes(i) = chol_Ki_mats[i].sparse.rows();
    }
    total_size += sizes(i);
  }
  MatrixXd chol_K_inv_dense(total_size,total_size);
  chol_K_inv_dense.setZero();
  int curr_row = 0;
  int curr_col = 0;
  for(int i = 0; i < h; i++){
    if(h2s[i] == 0) {
      chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal().setOnes();
      chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal() /= 0;
    } else{
      // MSpMatd chol_Ki = chol_Ki_mats[i].sparse;
      if(chol_Ki_mats[i].isDense) {
        chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki_mats[i].dense;
      } else{
        chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki_mats[i].sparse;
      }
      chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) /= sqrt(h2s[i]);
    }
    curr_row += sizes(i);
    curr_col += sizes(i);
  }
  SpMatd chol_K_inv = chol_K_inv_dense.sparseView(0,tol);
  return chol_K_inv;
}

// [[Rcpp::export()]]
Rcpp::List make_chol_ZtZ_Kinv_list(Rcpp::List chol_Ki_mats_,
                                   Map<MatrixXd> h2s_matrix,
                                   MSpMatd ZtZ,
                                   double drop0_tol,
                                   SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
                                   int ncores) {
  int s = h2s_matrix.cols();
  
  
  std::vector<R_matrix> chol_Ki_mats;
  load_R_matrices_list(chol_Ki_mats_,chol_Ki_mats);
  
  std::vector<SpMatd> chol_ZtZ_Kinv_list;
  chol_ZtZ_Kinv_list.reserve(s);
  for(int i = 0; i < s; i++){
    chol_ZtZ_Kinv_list.push_back(SpMatd(0,0));
  }
  
  int n_groups = s/ncores;
  for(int i = 0; i <= n_groups; i++){
    int start = ncores*i;
    int end = std::min(static_cast<double>(s),ncores*(i+1.0));
    
    for(std::size_t j = start; j < end; j++){
      VectorXd h2s = h2s_matrix.col(j);
      SpMatd chol_K_inv = make_chol_K_inv(chol_Ki_mats,h2s,drop0_tol);
      MatrixXd ZtZ_Kinv = 1.0/(1.0 - h2s.sum()) * ZtZ + chol_K_inv.transpose() * chol_K_inv;
      Eigen::LLT<MatrixXd> chol_ZtZ_Kinv(ZtZ_Kinv);
      MatrixXd chol_ZtZ_Kinv_R = chol_ZtZ_Kinv.matrixU();
      chol_ZtZ_Kinv_list[j] = chol_ZtZ_Kinv_R.sparseView(0,drop0_tol);
    }
    int pb_state = as<int>(getTxtProgressBar(pb));
    setTxtProgressBar(pb,pb_state+(end-start));
  }
  
  Rcpp::List chol_ZtZ_Kinv_list_out(s);
  for(int i = 0; i < s; i++){
    chol_ZtZ_Kinv_list_out[i] = chol_ZtZ_Kinv_list[i];
  }
  
  return(chol_ZtZ_Kinv_list_out);
}



SpMat make_chol_R(const std::vector<General_Matrix_f>& ZKZts, const VectorXf h2s, const float tol){  
  
  int n = ZKZts[0].rows();
  bool dense = ZKZts[0].isDense;
  int h = h2s.size();
  MatrixXf Rd(n,n);
  SpMat Rs(n,n);
  Rs.setZero();
  if(dense) {
    Rd.setZero();
  }
  for(int i = 0; i < h; i++){
    if(ZKZts[i].isDense) {
      if(!dense) {
        Rd = Rs.toDense();
        dense = true;
      }
      Rd += h2s[i] * ZKZts[i].dense;
    } else{
      if(dense) {
        Rd += h2s[i] * ZKZts[i].sparse;
      } else{
        Rs += h2s[i] * ZKZts[i].sparse;
      }
    }
  }
  if(dense) {
    Rd.diagonal().array() += (1.0-h2s.sum());
    Eigen::LLT<MatrixXf> chol_R(Rd);
    MatrixXf chol_R_U = chol_R.matrixU();
    return chol_R_U.sparseView(0,tol);
  } else{
    for(int i = 0; i < n; i++){
      Rs.coeffRef(i,i) += (1.0-h2s.sum());
    }
    // diagonal().array() += (1.0-h2s.sum());
    Eigen::SimplicialLLT<SpMat,Lower,NaturalOrdering<int> > chol_R(Rs);
    MatrixXf chol_R_U = chol_R.matrixU();
    return chol_R_U.sparseView(0,tol);
  }
}


// [[Rcpp::export()]]
Rcpp::List make_chol_V_list(Rcpp::List ZKZts_,
                            MatrixXf h2s_matrix,
                            float drop0_tol,
                            SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
                            int ncores) {
  int s = h2s_matrix.cols();

  std::vector<General_Matrix_f> ZKZts;
  load_General_Matrix_f_list(ZKZts_,ZKZts,false);

  std::vector<SpMat> chol_R_list;
  chol_R_list.reserve(s);
  for(int i = 0; i < s; i++){
    chol_R_list.push_back(SpMat(0,0));
  }

  int n_groups = s/ncores;
  for(int i = 0; i <= n_groups; i++){
    int start = ncores*i;
    int end = std::min(static_cast<double>(s),ncores*(i+1.0));
    for(std::size_t j = start; j < end; j++){
      chol_R_list[j] = make_chol_R(ZKZts, h2s_matrix.col(j), drop0_tol);
    }
    int pb_state = as<int>(getTxtProgressBar(pb));
    setTxtProgressBar(pb,pb_state+(end-start));
  }

  Rcpp::List chol_R_list_out(s);
  for(int i = 0; i < s; i++){
    chol_R_list_out[i] = chol_R_list[i];
  }

  return(chol_R_list_out);
}
