// Copyright 2020 Daniel Runcie
// Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
// that can be found in the LICENSE file and available at
// https://polyformproject.org/licenses/noncommercial/1.0.0/


// #include <math.h>
// #include <iostream>
// #include "MegaLMM_types.h"
// 
// using namespace Eigen;
// 
// // [[Rcpp::export()]]
// void load_sparse(MSpMat X, Map<MatrixXd> y){
//   MatrixXd z = X*y;
// }
// // [[Rcpp::export()]]
// void load_sparsef(SpMatf X, MatrixXf y){
//   MatrixXf z = X*y;
// }
// 
// // [[Rcpp::export()]]
// void load_sparsef2(MSpMat X_, MatrixXf y){
//   SpMatf X = X_.cast<float>();
//   MatrixXf z = X*y;
// }
// 
// // [[Rcpp::export()]]
// void load_dense(Map<MatrixXd> X, MatrixXd y){
//   // Spmatf X = X_.cast<float>();
//   MatrixXd z = X*y;
// }
// 
// // [[Rcpp::export()]]
// void load_dense2(MatrixXf X, MatrixXf y){
//   // Spmatf X = X_.cast<float>();
//   MatrixXf z = X*y;
// }


// // [[Rcpp::export()]]
// MatrixXf rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXf X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXf X_mat = Map<MatrixXf>(X_vec.data(),n,p);
//   return(X_mat);
// }
// 
// 
// struct R_matrix2 {
//   Map<MatrixXf> dense;
//   MSpMatf sparse;
//   bool isDense;
//   bool isNULL;
//   R_matrix2(Map<MatrixXf> dense_, MSpMatf sparse_,bool isDense_, bool isNULL_) : dense(dense_), sparse(sparse_), isDense(isDense_), isNULL(isNULL_){}
// };
// void load_R_matrices_list2(const Rcpp::List X_list, std::vector<R_matrix2>& X_vector);
// 
// // Loads a sparse or dense matrix passed from R into a \code{R_matrix2} object
// R_matrix2 load_R_matrix2(SEXP X_) {
//   MatrixXf null_d = MatrixXf::Zero(0,0);
//   if(X_ == R_NilValue) {
//     SpMatf null_s = null_d.sparseView();
//     MSpMatf M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//     Map<MatrixXf> M_null_d(null_d.data(),0,0);
//     R_matrix2 Xm(M_null_d,M_null_s,false,true);
//     return(Xm);
//   }
//   if(Rf_isMatrix(X_)){
//     Map<MatrixXf> X = as<Map<MatrixXf> >(X_);
//     SpMatf null_s = null_d.sparseView();
//     MSpMatf M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//     R_matrix2 Xm(X,M_null_s,true,false);
//     return(Xm);
//   } else{
//     MSpMatf X = as<MSpMatf>(X_);
//     Map<MatrixXf> M_null_d(null_d.data(),0,0);
//     R_matrix2 Xm(M_null_d,X,false,false);
//     return(Xm);
//   }
// }
// 
// // Code to convert list of R matrices (sparse or dense) into a thread-safe object
// //
// // @param X_list List of matrices (each can be dgCMatrix or matrix)
// // @param X_vector \code{std::vector} of \code{R_matrix2} type which will be populated
// void load_R_matrices_list2(const Rcpp::List X_list, std::vector<R_matrix2>& X_vector){
//   // null_matrices
//   MatrixXf null_d = MatrixXf::Zero(0,0);
//   Map<MatrixXf> M_null_d(null_d.data(),0,0);
//   SpMatf null_s = null_d.sparseView();
//   MSpMatf M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//   
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     X_vector.push_back(load_R_matrix2(Xi_));
//     // if(Rf_isMatrix(Xi_)){
//     //   Map<MatrixXf> Xi = as<Map<MatrixXf> >(Xi_);
//     //   R_matrix2 Xim(Xi,M_null_s,true);
//     //   X_vector.push_back(Xim);
//     // } else{
//     //   MSpMatf Xi = as<MSpMatf>(Xi_);
//     //   R_matrix2 Xim(M_null_d,Xi,false);
//     //   X_vector.push_back(Xim);
//     // }
//   }
// }
// 
// VectorXi which2(Rcpp::LogicalVector x) {
//   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
//   return as<VectorXi>(v[x]);
// }
// 
// 
// MatrixXf slice(MatrixXf& X,VectorXi rows, VectorXi cols) {
//   MatrixXf out(rows.size(),cols.size());
//   if(rows.maxCoeff() > X.rows() || rows.minCoeff() <= 0) stop("row index out of bounds");
//   if(cols.maxCoeff() > X.cols() || cols.minCoeff() <= 0) stop("col index out of bounds");
//   // #pragma omp parallel for
//   for(int i = 0; i < cols.size(); i++) {
//     VectorXf col_i = X.col(cols[i]-1);
//     for(int j = 0; j < rows.size(); j++) {
//       out.coeffRef(j,i) = col_i[rows[j]-1];
//     }
//   }
//   return(out);
// }
// 
// MatrixXf get_RinvtX2(const R_matrix2& chol_V, MatrixXf X){
//   MatrixXf RinvtX2;
//   if(chol_V.isDense) {
//     RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   } else{
//     RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   }
//   return(RinvtX2);
// }
// 
// VectorXf regression_sampler_v12(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
//     const Ref<const VectorXf>& y,           // nx1
//     const MatrixXf& X1,           // nxa
//     const MatrixXf& RinvtX2,                // nxb
//     const MatrixXf& C,                     // bxb
//     const VectorXf& prior_prec_alpha, // ax 1
//     const Ref<const VectorXf>& prior_mean_beta,  // bx1
//     const Ref<const VectorXf>& prior_prec_beta,  // bx1
//     const R_matrix2& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix. Also called R here with RtR = V
//     float Y_prec,                    // float
//     const VectorXf& randn_alpha,
//     const Ref<const VectorXf>& randn_beta,
//     const float rgamma_1,
//     const float Y_prec_b0
// ){
//   int n = y.size();
//   int a = X1.cols();
//   int b = RinvtX2.cols();
//   
//   // Check inputs
//   if(X1.rows() != n) stop("Wrong dimension of X1");
//   if(RinvtX2.rows() != n) stop("Wrong dimension of X2");
//   if(C.rows() != b || C.cols() != b) stop("Wrong dimension of C");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//   
//   // Calculate cholesky of A_beta
//   // C = Xt(RtR)^-1X
//   MatrixXf C_beta = C;
//   C_beta.diagonal() += prior_prec_beta;
//   LLT<MatrixXf> A_beta_llt;
//   A_beta_llt.compute(C_beta);
//   MatrixXf chol_A_beta = A_beta_llt.matrixU();
//   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
//   
//   // Step 1
//   VectorXf alpha(a);
//   VectorXf y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // Calculate A_alpha = Y_prec*X1^T*V_beta^{-1}*X1 + D_alpha^{-1}
//     // We don't need to actually calculate V_beta^{-1} directly.
//     // Instead,
//     MatrixXf RinvtX1 = get_RinvtX2(chol_V,X1);  // n*n*a -> n x a
//     MatrixXf X1tVinvX2 = RinvtX1.transpose() * RinvtX2; // a*n*b -> a*b
//     MatrixXf cholAbetainvt_X2tVinvX1 = chol_A_beta.transpose().triangularView<Lower>().solve(X1tVinvX2.transpose()); // b*b*a -> b x a
//     
//     MatrixXf A_alpha = RinvtX1.transpose() * RinvtX1 - cholAbetainvt_X2tVinvX1.transpose() * cholAbetainvt_X2tVinvX1;
//     A_alpha.diagonal() += prior_prec_alpha;
//     
//     VectorXf Rinvsqy = get_RinvtX2(chol_V,y); // n*n*q -> n x 1;
//     VectorXf XtRinvy = RinvtX2.transpose() * Rinvsqy; // b*n*1 >- b x 1
//     VectorXf invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
//     
//     VectorXf X1t_V_beta_inv_y = RinvtX1.transpose() * Rinvsqy - cholAbetainvt_X2tVinvX1.transpose() * invSqAbXtRinvy;
//     
//     LLT<MatrixXf> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
//     
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     
//     y_tilde = y - X1 * alpha;
//   }
//   
//   // Step 2 - sample Y_prec
//   // We don't need to actually calculate V_beta^{-1} directly.
//   VectorXf RinvSqy = get_RinvtX2(chol_V,y_tilde);
//   VectorXf XtRinvy = RinvtX2.transpose() * RinvSqy;
//   VectorXf prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
//   float score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
//   Y_prec = rgamma_1/score;
//   
//   // Step 3 - sample beta
//   VectorXf XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
//   VectorXf beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
//   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
//   
//   VectorXf result(1+a+b);
//   result << Y_prec,alpha,beta;
//   
//   return(result);
// }
// 
// 
// VectorXf regression_sampler_v22(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
//     const Ref<const VectorXf>& y,           // nx1
//     const MatrixXf& X1,           // nxa
//     const MatrixXf& X2,           // nxm or nxb
//     const VectorXf& prior_prec_alpha, // ax 1
//     const Ref<const VectorXf>& prior_mean_beta,  // bx1
//     const Ref<const VectorXf>& prior_prec_beta,  // bx1
//     const R_matrix2& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXf& V,
//     float Y_prec,                    // float
//     const VectorXf& randn_alpha,
//     const Ref<const VectorXf>& randn_beta,
//     const Ref<const VectorXf>& randn_e,
//     const float rgamma_1,
//     const float Y_prec_b0
// ){
//   int n = y.size();
//   int a = X1.cols();
//   int b = X2.cols();
//   
//   // Check inputs
//   if(X1.rows() != n) stop("Wrong dimension of X1");
//   if(X2.rows() != n) stop("Wrong dimension of X2");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//   
//   // Calculate inverse of V_beta
//   MatrixXf Dbeta_X2t = prior_prec_beta.cwiseInverse().asDiagonal() * X2.transpose();
//   MatrixXf V_beta = X2 * Dbeta_X2t + V;
//   LDLT<MatrixXf> V_beta_ldlt;
//   V_beta_ldlt.compute(V_beta);
//   
//   // Step 1
//   VectorXf alpha(a);
//   VectorXf y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     MatrixXf V_beta_inv_X1 = V_beta_ldlt.solve(X1);
//     MatrixXf A_alpha = V_beta_inv_X1.transpose() * X1;
//     A_alpha.diagonal() += prior_prec_alpha;
//     
//     LLT<MatrixXf> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
//     
//     VectorXf X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     y_tilde = y - X1 * alpha;
//   }
//   
//   // Step 2 - sample Y_prec
//   VectorXf e2 = y_tilde.transpose() * V_beta_ldlt.solve(y_tilde);
//   float score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
//   
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXf v = sqrt(Y_prec) * X2 * u;
//   if(chol_V.isDense) {
//     v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   VectorXf w = V_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
//   VectorXf beta = u + Dbeta_X2t * w / sqrt(Y_prec);
//   
//   VectorXf result(1+a+b);
//   result << Y_prec,alpha,beta;
//   
//   return(result);
// }
// 
// 
// VectorXf regression_sampler_v32(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
//     const Ref<const VectorXf>& y,           // nx1
//     const MatrixXf& X1,           // nxa
//     const MatrixXf& Ux,           // nxm or nxb
//     const MatrixXf& Vx,           // mxb
//     const VectorXf& prior_prec_alpha, // ax 1
//     const Ref<const VectorXf>& prior_mean_beta,  // bx1
//     const Ref<const VectorXf>& prior_prec_beta,  // bx1
//     const R_matrix2& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXf& Vinv,
//     const MatrixXf& VinvUx,
//     const MatrixXf& UtVinvU,
//     float Y_prec,                    // float
//     const VectorXf& randn_alpha,
//     const Ref<const VectorXf>& randn_beta,
//     const Ref<const VectorXf>& randn_e,
//     const float rgamma_1,
//     const float Y_prec_b0
// ){
//   int n = y.size();
//   int a = X1.cols();
//   if(Vx.rows() != Ux.cols()) stop("Wrong dimensions of Vx");
//   MatrixXf X2 = Ux*Vx;
//   int b = X2.cols();
//   
//   // Check inputs
//   if(X1.rows() != n) stop("Wrong dimension of X1");
//   if(X2.rows() != n) stop("Wrong dimension of X2");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//   if(randn_e.size() != n) stop("Wrong length of randn_e");
//   
//   // Calculate inverse of V_beta
//   // Using Ainv - Ainv * Ux * (I + BVAinvU)inv * BVAinv in case B = VxDbeta_Vxt is singular
//   MatrixXf Dbeta_Vxt = prior_prec_beta.cwiseInverse().asDiagonal() * Vx.transpose();
//   MatrixXf VxDbeta_Vxt = Vx * Dbeta_Vxt;
//   if(VinvUx.rows() != n) stop("Wrong dimensions of VinvUx");
//   MatrixXf inner = VxDbeta_Vxt * UtVinvU;
//   inner.diagonal().array() += 1.0;
//   LDLT<MatrixXf> inner_ldlt;
//   inner_ldlt.compute(inner);
//   // Sigma_beta_inv = Vinv - VinvUx * inner.ldlt().solve(VxDbeta_Vxt * VinvUx.transpose());  // Don't actually calculate this. Stay in mxm space
//   
//   // Step 1
//   VectorXf alpha(a);
//   VectorXf y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // MatrixXf V_beta_inv_X1 = Sigma_beta_inv * X1;
//     // MatrixXf A_alpha = Y_prec * V_beta_inv_X1.transpose() * X1;
//     MatrixXf Vinv_X1 = Vinv * X1;
//     MatrixXf Uxt_Vinv_X1 = Ux.transpose() * Vinv_X1;
//     MatrixXf A_alpha = X1.transpose() * Vinv_X1 - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt) * Uxt_Vinv_X1;
//     A_alpha.diagonal() += prior_prec_alpha;
//     
//     // VectorXf X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
//     VectorXf Uxt_Vinv_y = VinvUx.transpose() * y;
//     VectorXf X1t_V_beta_inv_y = Vinv_X1.transpose() * y - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
//     
//     LLT<MatrixXf> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
//     
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     
//     y_tilde = y - X1 * alpha;
//   }
//   
//   // Step 2 - sample Y_prec
//   // VectorXf e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
//   VectorXf Vinv_y_tilde = Vinv * y_tilde;
//   VectorXf Uxt_Vinv_y = VinvUx.transpose() * y_tilde;
//   VectorXf e2 = y_tilde.transpose() * Vinv_y_tilde - Uxt_Vinv_y.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
//   
//   float score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
//   
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXf v = std::sqrt(Y_prec) * X2 * u;
//   if(chol_V.isDense) {
//     v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   // VectorXf X1 = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
//   VectorXf e = y_tilde * std::sqrt(Y_prec) - v;
//   VectorXf Uxt_Vinv_e = VinvUx.transpose() * e;
//   VectorXf w = Vinv * e - VinvUx * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_e);
//   
//   VectorXf beta = u + Dbeta_Vxt * (Ux.transpose() * w) / std::sqrt(Y_prec); //b*b*1 + b*n*1
//   
//   VectorXf result(1+a+b);
//   result << Y_prec,alpha,beta;
//   
//   return(result);
// }
// 
// 
// struct regression_result {
//   MatrixXf alpha1;
//   VectorXf alpha2_vec;
//   MatrixXf beta;
//   VectorXf Y_prec;
//   regression_result(MatrixXf alpha1_, VectorXf alpha2_vec_, MatrixXf beta_, VectorXf Y_prec_)
//     alpha1(alpha1_), alpha2_vec(alpha2_vec_), beta(beta_), Y_prec(Y_prec_): {}
// };
// 
// //' Draws samples from all ``fixed" coefficients (fixed and random) of a set of parallel linear regression models, conditional on the variance components.
// //'
// //' The model is either: \itemize{
// //' \item y_i = X1_base*alpha1 + X1_list_[i]*alpha2 + X2*beta + e, e ~ N(0,1/Y_prec[i]*V)
// //' \item y_i = X1_base*alpha1 + X1_list_[i]*alpha2 + X2*V_*beta + e, e ~ N(0,1/Y_prec[i]*V)
// //' }
// //' Where \code{V = RtR}, priors on elements of alpha1, alpha2 and beta are independent.
// //' Each column of Y is considered independent
// //'
// //' @param Y n x p matrix of observations
// //' @param X1_base n x a1 matrix of X1 covariates common to all p. Can be NULL
// //' @param X1_list_ p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
// //' @param X2 either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
// //' @param V_ m x b matrix if X2 is Ux, otherwise NULL
// //' @param h2s_index p-vector of indices for to select appropriate V of each trait
// //' @param chol_V_list_ list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
// //' @param Y_prec p-vector of Y current precisions
// //' @param Y_prec_a0,Y_prec_b0 scalars giving the shape and rate of the Gamma distribution for the prior on Y_prec
// //' @param prior_prec_alpha1 a1 x p matrix of prior precisions for alpha1
// //' @param prior_prec_alpha2 p-vector of precision of alpha2s for each trait
// //' @param prior_mean_beta b x p matrix of prior means of beta
// //' @param prior_prec_beta b x p matrix of prior precisions of beta
// //' @return List with elements: \itemize{
// //'   \item alpha1 a1 x p matrix of alpha1
// //'   \item alpha2 concatenated vector of alpha2 for all traits
// //'   \item beta b x p matrix of beta
// //'   \item Y_prec p x 1 vector of Y_prec
// //' }
// regression_result regression_sampler_parallel2(
//     MatrixXf& Y,               //
//     MatrixXf& X1_base,          // std::vector<R_matrix2>& X1_list,             // p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
//     MatrixXf& X2,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
//     SEXP Vx_,                       // m x b matrix if X2 is Ux
//     VectorXi h2s_index, // p-vector of indices for appropriate V of each trait
//     std::vector<R_matrix2>& chol_V_list,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//     VectorXf Y_prec,               // p-vector of Y current precisions
//     VectorXf Y_prec_a0,
//     VectorXf Y_prec_b0,
//     MatrixXf& prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     VectorXf& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     MatrixXf& prior_mean_beta, // b x p matrix of prior means of beta
//     MatrixXf& prior_prec_beta // b x p matrix of prior precisions of beta
// ) {
//   
//   int n = Y.rows();
//   int p = Y.cols();
//   
//   // X1_base
//   if(X1_base.rows() != n) stop("Wrong dimension of X1_base");
//   int a1 = X1_base.cols();
//   
//   // X1_list
//   if(X1_list.size() > 0) {
//     if(X1_list.size() != p) stop("Wrong length of X1_list");
//   }
//   
//   // X2 or Ux and V
//   MatrixXf Ux = X2;
//   MatrixXf z = MatrixXf::Zero(0,0);
//   Map<MatrixXf> Vx(z.data(),0,0);
//   int b = X2.cols();
//   if(X2.rows() != n) stop("Wrong dimension of X2");
//   if(Rf_isMatrix(Vx_)) {
//     // Map<MatrixXf> V__ = as<Map<MatrixXf> >(Vx_);
//     // new (&v) Map<MatrixXf> (V__,V__.rows(),V__.cols());
//     new (&Vx) Map<MatrixXf> (as<Map<MatrixXf> >(Vx_));
//     if(Ux.cols() != Vx.rows()) stop("X2 and Vx_ have incompatible dimensions");
//     b = Vx.cols();
//   }
//   
//   // priors
//   if(Y_prec.size() != p) {
//     stop("Wrong length of Y_prec");
//   }
//   if(Y_prec_a0.size() != p) {
//     stop("Wrong length of Y_prec_a0");
//   }
//   if(Y_prec_b0.size() != p) {
//     stop("Wrong length of Y_prec_b0");
//   }
//   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
//   if(X1_list.size() > 0 && prior_prec_alpha2.size() != p) {
//     stop("Wrong length of prior_prec_alpha2");
//   }
//   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
//   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
//   
//   // generate random numbers
//   MatrixXf randn_alpha1 = rstdnorm_mat2(a1,p);
//   std::vector<VectorXf> randn_alpha2;
//   if(X1_list.size() > 0){
//     for(int i = 0; i < p; i++){
//       randn_alpha2.push_back(rstdnorm_mat2(X1_list[i].dense.cols(),1));
//     }
//   }
//   MatrixXf randn_beta = rstdnorm_mat2(b,p);
//   MatrixXf randn_e;
//   if(b > n) {
//     randn_e = rstdnorm_mat2(n,p);
//   }
//   VectorXf rgamma_1(p);
//   for(int i = 0; i < p; i++) {
//     rgamma_1(i) = rgamma(1,Y_prec_a0[i] + n/2.0,1.0)(0);
//   }
//   
//   // Results structures
//   MatrixXf alpha1(a1,p);
//   std::vector<VectorXf> alpha2;
//   alpha2.reserve(X1_list.size());
//   int alpha2_size = 0;
//   if(X1_list.size() > 0){
//     for(int i = 0; i < X1_list.size(); i++){
//       int a2 = X1_list[i].dense.cols();
//       alpha2.push_back(VectorXf::Zero(a2));
//       alpha2_size += a2;
//     }
//   }
//   MatrixXf beta(b,p);
//   
//   // go through h2s indices and sample columns with same index as a set
//   for(int i = h2s_index.minCoeff(); i <= h2s_index.maxCoeff(); i++) {
//     int h2_index = i;
//     VectorXi trait_set = which2(as<LogicalVector>(h2s_index.array() == h2_index));  // list of traits with same h2_index
//     
//     if(trait_set.size() > 0){
//       // prepare matrices for sampler
//       MatrixXf RinvtX2, C, V, Vinv, VinvUx, UtVinvU;
//       R_matrix2 chol_V = chol_V_list[h2_index - 1];
//       int which_sampler;
//       // Decide which sampler to use
//       if(b <= n) {
//         // use regression_sampler_v12
//         which_sampler = 1;
//         if(chol_V.isDense) {
//           RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X2);
//         } else{
//           RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X2);
//         }
//         C = RinvtX2.transpose() * RinvtX2;
//       }
//       else if(Vx.cols() == 0) {
//         // use regression_sampler_v22
//         which_sampler = 2;
//         if(chol_V.isDense) {
//           V = chol_V.dense.transpose().triangularView<Lower>() * chol_V.dense;
//         } else{
//           V = chol_V.sparse.transpose().triangularView<Lower>() * chol_V.sparse;
//         }
//       } else {
//         // use regression_sampler_v32
//         which_sampler = 3;
//         if(chol_V.isDense) {
//           Vinv = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
//           VinvUx = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(Ux));
//         } else{
//           Vinv = chol_V.sparse.triangularView<Upper>().solve(chol_V.sparse.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
//           VinvUx = chol_V.sparse.triangularView<Upper>().solve(chol_V.sparse.transpose().triangularView<Lower>().solve(Ux));
//           // VinvUx = Vinv * U;
//           // Rcout << i << std::endl;
//           // Rcout << Vinv.diagonal().transpose() << std::endl;
//         }
//         UtVinvU = Ux.transpose() * VinvUx;
//       }
//       // REprintf("Number of threads=%i\\n", omp_get_max_threads());
//       // Rcout <<trait_set.size() << " " << omp_get_max_threads() << std::endl;
//       #pragma omp parallel for
//       for(int i = 0; i < trait_set.size(); i++){
//         int j = trait_set[i];
//         MatrixXf X1;
//         int a;
//         int a2 = 0;
//         int b;
//         VectorXf prior_prec_alpha;
//         VectorXf randn_alpha;
//         if(X1_list.size() == 0) {
//           X1 = X1_base;
//           a = a1;
//           prior_prec_alpha = prior_prec_alpha1.col(j);
//           randn_alpha = randn_alpha1.col(j);
//         } else{
//           Map<MatrixXf> X12 = X1_list[j].dense;
//           a2 = X12.cols();
//           a = a1+a2;
//           X1 = MatrixXf(n,a);
//           X1 << X1_base,X12;
//           prior_prec_alpha = VectorXf(a);
//           prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
//           prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
//           randn_alpha = VectorXf(a);
//           randn_alpha.head(a1) = randn_alpha1.col(j);
//           randn_alpha.tail(a2) = randn_alpha2[j];
//         }
//         
//         VectorXf samples;
//         if(which_sampler == 1) {
//           b = RinvtX2.cols();
//           samples = regression_sampler_v12(Y.col(j), X1, RinvtX2, C, prior_prec_alpha, prior_mean_beta.col(j),
//                                           prior_prec_beta.col(j), chol_V, Y_prec[j], randn_alpha,
//                                           randn_beta.col(j), rgamma_1[j],Y_prec_b0[j]);
//         } else if(which_sampler == 2) {
//           b = X2.cols();
//           samples = regression_sampler_v22(Y.col(j), X1, X2, prior_prec_alpha, prior_mean_beta.col(j),
//                                           prior_prec_beta.col(j), chol_V, V, Y_prec[j], randn_alpha,
//                                           randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0[j]);
//         } else if(which_sampler == 3) {
//           b = Vx.cols();
//           samples = regression_sampler_v32(Y.col(j), X1, Ux, Vx, prior_prec_alpha, prior_mean_beta.col(j),
//                                           prior_prec_beta.col(j), chol_V, Vinv, VinvUx, UtVinvU, Y_prec[j], randn_alpha,
//                                           randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0[j]);
//           
//         } else {
//           stop("sampler not implemented");
//         }
//         
//         // extract samples
//         Y_prec[j] = samples[0];
//         if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
//         if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
//         if(b > 0) beta.col(j) = samples.tail(b);
//       }
//     }
//   }
//   
//   // collect alpha2 into a vector
//   VectorXf alpha2_vec(alpha2_size);
//   if(X1_list.size() > 0){
//     int index = 0;
//     for(int i = 0; i < X1_list.size(); i++){
//       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
//       index += alpha2[i].size();
//     }
//   }
//   regression_result result(alpha1,alpha2_vec,beta,Y_prec);
//   return(result);
//   // return(Rcpp::List::create(
//   //     Named("alpha1") = alpha1,
//   //     Named("alpha2") = alpha2_vec,
//   //     Named("beta") = beta,
//   //     Named("Y_prec") = Y_prec
//   // ));
// }
// 
// 
// // [[Rcpp::export]]
// List MegaLMM_sampler(
//   List current_state,
//   List data_matrices,
//   List prior_hyperparameters,
//   List run_parameters,
//   List run_variables,
//   int n_samples
// ) {
//   
//   // initialize everything!
//   
//   // data_matrices
//   MatrixXf X1 = as<MatrixXf>(data_matrices["X1"]);
//   MatrixXf X2_R = as<MatrixXf>(data_matrices["X2_R"]);
//   MatrixXf X2_F = as<MatrixXf>(data_matrices["X2_F"]);
//   MatrixXf U2_F = as<MatrixXf>(data_matrices["U2_F"]);
//   MatrixXf V2_F;
//   if(data_matrices["V2_F"] != R_NilValue) {
//     V2_F = as<MatrixXf>(data_matrices["V2_F"]);
//   } 
//   SpMatf ZL = as<SpMatf>(data_matrices["ZL"]);
//   SpMatf RE_L = as<SpMatf>(data_matrices["RE_L"]);
//   MatrixXf h2s_matrix = as<MatrixXf>(data_matrices["h2s_matrix"]);
//   MatrixXf Lambda_fixed;
//   if(data_matrices["Lambda_fixed"] != R_NilValue) {
//     Lambda_fixed = as<MatrixXf>(data_matrices["Lambda_fixed"]);
//   }
//   
//   // current_state
//   int nrun = as<int>(current_state["nrun"]);
//   float total_time = as<float>(current_state["total_time"]);
//   MatrixXf Eta, Lambda, B1, B2_R, B2_F, F, U_F, U_R, resid_h2, F_h2, tot_Eta_prec, tot_F_prec,delta;
//   VectorXi resid_h2_index, F_h2_index;
//   SpMatb Y_missing;
// 
//   Eta = as<MatrixXf>(current_state["Eta"]);
// 
//   Lambda = as<MatrixXf>(current_state["Lambda"]);
//   if(current_state["B1"] != R_NilValue) B1 = as<MatrixXf>(current_state["B1"]);
//   if(current_state["B2_R"] != R_NilValue) B2_R = as<MatrixXf>(current_state["B2_R"]);
//   if(current_state["B2_F"] != R_NilValue) B2_F = as<MatrixXf>(current_state["B2_F"]);
// 
//   F = as<MatrixXf>(current_state["F"]);
// 
//   U_F = as<MatrixXf>(current_state["U_F"]);
//   U_R = as<MatrixXf>(current_state["U_R"]);
// 
//   resid_h2 = as<MatrixXf>(current_state["resid_h2"]);
//   F_h2 = as<MatrixXf>(current_state["F_h2"]);
// 
//   tot_Eta_prec = as<MatrixXf>(current_state["tot_Eta_prec"]);
//   tot_F_prec = as<MatrixXf>(current_state["tot_F_prec"]);
// 
//   resid_h2_index = as<VectorXi>(current_state["resid_h2_index"]);
//   F_h2_index = as<VectorXi>(current_state["F_h2_index"]);
// 
//   delta = as<MatrixXf>(current_state["delta"]);
//   // 
//   // Y_missing = as<SpMatb>(current_state["Y_missing"]);
// 
//   // prior-specific values (not all will be initialized)
//   // if you need more variables for your priors, add them here
//   // also, be sure to add them to the constuctor and makeRList()
//     // Lambda priors:
//       // Horseshoe
//       MatrixXf Lambda_tau2, Lambda_xi, Lambda_phi2, Lambda_nu, Lambda_m_eff;
// 
//     // Beta priors:
//       // Horseshoe
//       MatrixXf B2_R_tau2, B2_R_xi, B2_R_phi2, B2_R_nu, B2_R_m_eff;
//       MatrixXf B2_F_tau2, B2_F_xi, B2_F_phi2, B2_F_nu, B2_F_m_eff;
// 
// 
//   // dimensions
//   int p = Eta.cols();
//   int n = Eta.rows();
//   int b1 = X1.cols();
//   int b2_R = X2_R.cols();
//   int b2_F = X2_F.cols();
//   int r = U_R.rows();
//   int K = F.cols();
//   // int Kr = as<int>(run_variables["Kr"]);
//   LogicalVector fixed_factors = as<LogicalVector>(run_variables["fixed_factors"]);
//   int Kr = sum(!fixed_factors);
// 
//   // control parameters
//   float h2_step_size = as<float>(run_parameters["h2_step_size"]);
//   int verbose = as<int>(run_parameters["verbose"]);
//   int burn = as<int>(run_parameters["burn"]);
//   int thin = as<int>(run_parameters["thin"]);
// 
//   // pre-calculations
//   List Missing_data_map = as<List>(run_variables["Missing_data_map"]);
//   List Missing_row_data_map = as<List>(run_variables["Missing_row_data_map"]);
// 
//   std::vector<R_matrix2> Qt_list,QtX1_list,QtX2_R_list;
//   std::vector<VectorXi> QtX1_keepColumns_list;
//   load_R_matrices_list2(as<List>(run_variables["Qt_list"]), Qt_list);
//   load_R_matrices_list2(as<List>(run_variables["QtX1_list"]), QtX1_list);
//   load_R_matrices_list2(as<List>(run_variables["QtX2_R_list"]), QtX2_R_list);
//   List QtX1_keepColumns_list_ = as<List>(run_variables["QtX1_keepColumns_list"]);
//   for(int i = 0; i < QtX1_keepColumns_list_.size(); i++) {
//     QtX1_keepColumns_list.push_back(as<VectorXi>(QtX1_keepColumns_list_[i]));
//   }
// 
//   MatrixXf Qt1_U2_F = as<MatrixXf>(run_variables["Qt1_U2_F"]);
//   MatrixXf Qt1_X2_F = as<MatrixXf>(run_variables["Qt1_X2_F"]);
//   std::vector<std::vector<R_matrix2>> chol_V_list_list, chol_ZtZ_Kinv_list_list;
//   List chol_V_list_list_ = as<List>(run_variables["chol_V_list_list"]);
//   for(int i = 0; i < Qt_list.size(); i++) {
//     std::vector<R_matrix2> temp_list;
//     load_R_matrices_list2(chol_V_list_list_[i], temp_list);
//     chol_V_list_list.push_back(temp_list);
//   }
//   List chol_ZtZ_Kinv_list_list_ = as<List>(run_variables["chol_ZtZ_Kinv_list_list"]);
//   for(int i = 0; i < Qt_list.size(); i++) {
//     std::vector<R_matrix2> temp_list;
//     load_R_matrices_list2(chol_ZtZ_Kinv_list_list_[i], temp_list);
//     chol_ZtZ_Kinv_list_list.push_back(temp_list);
//   }
//   
//   
//   // Gibbs loop
//   for(int I = 0; I < n_samples*thin; I++) {
//     // sample_latent_traits
//     int n_coefs = b2_R + Kr;
//     MatrixXf prior_mean = MatrixXf::Zero(n_coefs,p);
//     MatrixXf prior_prec(n_coefs,p);
//     if(b2_R > 0) {
//       prior_prec << B2_R_prec,Lambda_prec;
//     } else{
//       prior_prec = Lambda_prec;
//     }
//     for(int set = 0; set < Missing_data_map.size(); set++) {
//       List Missing_data_map_set = as<List>(Missing_data_map[set]);
//       VectorXi cols = as<VectorXi>(Missing_data_map_set["Y_cols"]);
//       VectorXi rows = as<VectorXi>(Missing_data_map_set["Y_obs"]);
// 
//       if(cols.size() == 0 || rows.size() == 0) continue;
// 
//       R_matrix2 Qt_set = Qt_list[set];
//       R_matrix2 QtX2_R_set = QtX2_R_list[set];
//       MatrixXf X_set(rows.size(),QtX2_R_set.dense.cols() + sum(!fixed_factors));
//       MatrixXf Eta_set(rows.size(),cols.size()),F_set(rows.size(),sum(!fixed_factors));
// 
//       Eta_set = slice(Eta,rows,cols);
//       F_set = slice(F,rows,which2(!fixed_factors));
//       MatrixXf Y_set;
// 
//       if(Qt_set.isNULL) {
//         Y_set = Eta_set;
//         X_set << QtX2_R_set.dense,F_set;
//       } else if(Qt_set.isDense) {
//         Y_set = Qt_set.dense * Eta_set;
//         X_set << QtX2_R_set.dense,F_set;
//       } else{
//         Y_set = Qt_set.sparse * Eta_set;
//         X_set << QtX2_R_set.dense,F_set;
//       }
//       
//       regression_result new_samples = regression_sampler_parallel2(
//         Y_set,
//         QtX1_list[set],
//         X_set,
//         NULL,
//         
//                  
//                  
//       )
//         
//     }
//     
//     
//     
//     
//     
//     
//     
//   }
//   
//   
//   
//   
//   
//   
//   
//   
//   // update current_state
//   current_state["nrun"] = nrun;
//   current_state["total_time"] = total_time;
//   current_state["Eta"] = Eta;
//   current_state["Lambda"] = Lambda;
//   current_state["B1"] = B1;
//   current_state["B2_R"] = B2_R;
//   current_state["F"] = F;
//   current_state["U_F"] = U_F;
//   current_state["U_R"] = U_R;
//   current_state["resid_h2"] = resid_h2;
//   current_state["F_h2"] = F_h2;
//   current_state["tot_Eta_prec"] = tot_Eta_prec;
//   current_state["tot_F_prec"] = tot_F_prec;
//   current_state["resid_h2_index"] = resid_h2_index;
//   current_state["F_h2_index"] = F_h2_index;
//   current_state["delta"] = delta;
//   
//   return(current_state);
// }








// class current_state{
// public:
//   int nrun;
//   float total_time;
//   MatrixXf Eta;
//   MatrixXf Lambda, B1, B2_R, B2_F;
//   MatrixXf F;
//   MatrixXf U_F, U_R;
//   MatrixXf resid_h2, F_h2;
//   MatrixXf tot_Eta_prec, tot_F_prec;
//   VectorXi F_h2_index, resid_h2_index;
//   MatrixXf delta;
//   MatrixXi Y_missing;
//   
//   List prior_hyperparameters;
//   
//   // calculated values
//   MatrixXf XB, Eta_mean, Lambda_prec, B2_R_prec, B2_F_prec;
//   VectorXf theta;
//   
//   // prior-specific values (not all will be initialized)
//   // if you need more variables for your priors, add them here
//   // also, be sure to add them to the constuctor and makeRList()
//     // Lambda priors:
//       // Horseshoe
//       MatrixXf Lambda_tau2, Lambda_xi, Lambda_phi2, Lambda_nu, Lambda_m_eff;
//       
//     // Beta priors:
//       // Horseshoe
//       MatrixXf B2_R_tau2, B2_R_xi, B2_R_phi2, B2_R_nu, B2_R_m_eff;
//       MatrixXf B2_F_tau2, B2_F_xi, B2_F_phi2, B2_F_nu, B2_F_m_eff;
//   
//   
//   List makeRList();
//   
//   current_state();
//   current_state(List current_state_, List prior_hyperparameters_);
//   
// };
// 
// current_state::current_state(List current_state_, List prior_hyperparameters_): prior_hyperparameters(prior_hyperparameters_) {
//   nrun = as<int>(current_state_["nrun"]);
//   total_time = as<float>(current_state_["total_time"]);
//   
//   Eta = as<MatrixXf>(current_state_["Eta"]);
//   
//   Lambda = as<MatrixXf>(current_state_["Lambda"]);
//   B1 = as<MatrixXf>(current_state_["B1"]);
//   B2_R = as<MatrixXf>(current_state_["B2_R"]);
//   B2_F = as<MatrixXf>(current_state_["B2_F"]);
//   
//   F = as<MatrixXf>(current_state_["F"]);
//   
//   U_F = as<MatrixXf>(current_state_["U_F"]);
//   U_R = as<MatrixXf>(current_state_["U_R"]);
//   
//   resid_h2 = as<MatrixXf>(current_state_["resid_h2"]);
//   F_h2 = as<MatrixXf>(current_state_["F_h2"]);
//   
//   tot_Eta_prec = as<MatrixXf>(current_state_["tot_Eta_prec"]);
//   tot_F_prec = as<MatrixXf>(current_state_["tot_F_prec"]);
//   
//   resid_h2_index = as<VectorXi>(current_state_["resid_h2_index"]);
//   F_h2_index = as<VectorXi>(current_state_["F_h2_index"]);
//   
//   delta = as<MatrixXf>(current_state_["delta"]);
//   
//   Y_missing = as<MatrixXi>(current_state_["Y_missing"]);
//   
//   // XB = 
//   
// }

// List current_state::makeRList
