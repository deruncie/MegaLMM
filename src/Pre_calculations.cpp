#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export()]]
List LDLt(SEXP A_) {
  if(Rf_isMatrix(A_)){
    Map<MatrixXd> A = as<Map<MatrixXd> >(A_);
    Eigen::LDLT<MatrixXd> ldlt_A;
    ldlt_A.compute(A);
    MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
    MatrixXd P = ldlt_A.transpositionsP() * I;
    VectorXd d = ldlt_A.vectorD();
    MatrixXd L = ldlt_A.matrixL();
    SpMat Lsp = L.sparseView();
    if(static_cast<double>(Lsp.nonZeros()) / I.size() > 0.25) {
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
    MSpMat A = as<MSpMat>(A_);
    Eigen::SimplicialLDLT<SpMat> ldlt_A;
    ldlt_A.compute(A);
    MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
    MatrixXd P = ldlt_A.permutationP() * I;
    return(List::create(
        Named("P") = P.sparseView(),
        Named("L") = ldlt_A.matrixL(),
        Named("d") = ldlt_A.vectorD()
    ));
  }
}


SpMat make_chol_K_inv(const std::vector<R_matrix>& chol_Ki_mats, VectorXd h2s,double tol){
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
      // MSpMat chol_Ki = chol_Ki_mats[i].sparse;
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
  SpMat chol_K_inv = chol_K_inv_dense.sparseView(0,tol);
  return chol_K_inv;
}

// [[Rcpp::export()]]
Rcpp::List make_chol_ZtZ_Kinv_list(Rcpp::List chol_Ki_mats_,
                                   Map<MatrixXd> h2s_matrix,
                                   MSpMat ZtZ,
                                   double drop0_tol,
                                   bool verbose, SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
                                   int ncores) {
  int s = h2s_matrix.cols();


  std::vector<R_matrix> chol_Ki_mats;
  load_R_matrices_list(chol_Ki_mats_,chol_Ki_mats);

  std::vector<SpMat> chol_ZtZ_Kinv_list;
  chol_ZtZ_Kinv_list.reserve(s);
  for(int i = 0; i < s; i++){
    chol_ZtZ_Kinv_list.push_back(SpMat(0,0));
  }

  int n_groups = s/ncores;
  for(int i = 0; i <= n_groups; i++){
    int start = ncores*i;
    int end = std::min(static_cast<double>(s),ncores*(i+1.0));

    for(std::size_t j = start; j < end; j++){
      VectorXd h2s = h2s_matrix.col(j);
      SpMat chol_K_inv = make_chol_K_inv(chol_Ki_mats,h2s,drop0_tol);
      MatrixXd ZtZ_Kinv = 1.0/(1.0 - h2s.sum()) * ZtZ + chol_K_inv.transpose() * chol_K_inv;
      Eigen::LLT<MatrixXd> chol_ZtZ_Kinv(ZtZ_Kinv);
      MatrixXd chol_ZtZ_Kinv_R = chol_ZtZ_Kinv.matrixU();
      chol_ZtZ_Kinv_list[j] = chol_ZtZ_Kinv_R.sparseView(0,drop0_tol);
    }
    if(verbose) {
      int pb_state = as<int>(getTxtProgressBar(pb));
      setTxtProgressBar(pb,pb_state+(end-start));
    }
  }

  Rcpp::List chol_ZtZ_Kinv_list_out(s);
  for(int i = 0; i < s; i++){
    chol_ZtZ_Kinv_list_out[i] = chol_ZtZ_Kinv_list[i];
  }

  return(chol_ZtZ_Kinv_list_out);
}


SpMat make_chol_R(const std::vector<R_matrix>& ZKZts, const VectorXd h2s, const double tol){  //std::vector<Map<MatrixXd> > ZKZts
  // Map<MatrixXd> ZKZts_0 = as<Map<MatrixXd> >(ZKZts[0]);
  int n;
  bool dense = false;
  if(ZKZts[0].isDense) {
    n = ZKZts[0].dense.rows();
    dense = true;
  } else{
    n = ZKZts[0].sparse.rows();
  }
  int h = h2s.size();
  MatrixXd Rd(n,n);
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
    Eigen::LLT<MatrixXd> chol_R(Rd);
    MatrixXd chol_R_U = chol_R.matrixU();
    return chol_R_U.sparseView(0,tol);
  } else{
    for(int i = 0; i < n; i++){
      Rs.coeffRef(i,i) += (1.0-h2s.sum());
    }
    // diagonal().array() += (1.0-h2s.sum());
    Eigen::SimplicialLLT<SpMat,Lower,NaturalOrdering<int> > chol_R(Rs);
    MatrixXd chol_R_U = chol_R.matrixU();
    return chol_R_U.sparseView(0,tol);
  }
}


// [[Rcpp::export()]]
Rcpp::List make_chol_V_list(Rcpp::List ZKZts_,
                            Map<MatrixXd> h2s_matrix,
                            double drop0_tol,
                            bool verbose, SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
                            int ncores) {
  int s = h2s_matrix.cols();

  std::vector<R_matrix> ZKZts;
  load_R_matrices_list(ZKZts_,ZKZts);

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
    if(verbose) {
      int pb_state = as<int>(getTxtProgressBar(pb));
      setTxtProgressBar(pb,pb_state+(end-start));
    }
  }

  Rcpp::List chol_R_list_out(s);
  for(int i = 0; i < s; i++){
    chol_R_list_out[i] = chol_R_list[i];
  }

  return(chol_R_list_out);
}
