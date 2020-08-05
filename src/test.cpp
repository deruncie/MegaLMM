// #include <math.h>
// #include <iostream>
// #include "MegaLMM_types.h"
// 
// 
// 
// using namespace Eigen;
// 
// struct General_Matrix_f2 {
//   MatrixXf dense;
//   SpMat sparse;
//   bool triangular;
//   bool isDense;
//   bool isNULL = false;
//   General_Matrix_f2(MatrixXf dense_, SpMat sparse_,bool triangular_, bool isDense_,bool isNULL_) : dense(dense_), sparse(sparse_), triangular(triangular_), isDense(isDense_), isNULL(isNULL_){}
//   MatrixXf solve(MatrixXf) const;
//   MatrixXf tsolve(MatrixXf) const;
//   MatrixXf operator*(MatrixXf) const;
//   MatrixXf crossprod(MatrixXf) const;
//   float get_log_det();
//   int rows() const;
//   int cols() const;
// };
// 
// // [[Rcpp::export()]]
// MatrixXf rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
//   return(X_mat.cast<float>());
// }
// // Loads a sparse or dense matrix passed from R into a \code{General_Matrix_f2} object
// General_Matrix_f2 load_General_Matrix_f2(SEXP X_, bool triangular) {
//   MatrixXf null_d = MatrixXf::Zero(0,0);
//   if(Rf_isNull(X_)) {
//     SpMat null_s = null_d.sparseView();
//     General_Matrix_f2 Xm(null_d,null_s,triangular,false,true);
//     return(Xm);
//   } else if(Rf_isMatrix(X_)){
//     MatrixXf X = as<MatrixXf >(X_);
//     SpMat null_s = null_d.sparseView();
//     General_Matrix_f2 Xm(X,null_s,triangular,true,false);
//     return(Xm);
//   } else{
//     SpMat X = as<SpMat>(X_);
//     General_Matrix_f2 Xm(null_d,X,triangular,false,false);
//     return(Xm);
//   }
// }
// 
// MatrixXf General_Matrix_f2::solve(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.triangularView<Upper>().solve(Y));
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse.triangularView<Upper>().solve(Y));
//     }
//   } else{
//     // Note: these Cholesky's could be stored so they could be re-used.
//     if(isDense) {
//       return(dense.ldlt().solve(Y));
//     } else{
//       Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
//       return(ldlt.solve(Y));
//     }
//   }
// }
// MatrixXf General_Matrix_f2::tsolve(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.transpose().triangularView<Lower>().solve(Y));
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse.transpose().triangularView<Lower>().solve(Y));
//     }
//   } else{
//     // Note: these Cholesky's could be stored so they could be re-used.
//     if(isDense) {
//       return(dense.transpose().ldlt().solve(Y));
//     } else{
//       Eigen::SimplicialLDLT<SpMat> ldlt(sparse.transpose());
//       return(ldlt.solve(Y));
//     }
//   }
// }
// 
// MatrixXf General_Matrix_f2::operator*(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.triangularView<Upper>() * Y);
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse * Y);
//     }
//   } else {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense * Y);
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse * Y);
//     }
//   }
// }
// MatrixXf General_Matrix_f2::crossprod(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
//       return(dense.transpose().triangularView<Lower>() * Y);
//     } else{
//       if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
//       return(sparse.transpose().triangularView<Lower>() * Y);
//     }
//   } else{
//     if(isDense) {
//       if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
//       return(dense.transpose() * Y);
//     } else{
//       if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
//       return(sparse.transpose() * Y);
//     }
//   }
// }
// float General_Matrix_f2::get_log_det() {
//   if(isNULL) stop("no determinant of null matrix");
//   if(!triangular) stop("not implemented for non-triangular matrix");
//   if(isDense) {
//     return(dense.diagonal().array().log().sum());
//   } else {
//     float log_det = 0;
//     for(int j = 0; j < sparse.rows(); j++) {
//       log_det += std::log(sparse.coeffRef(j,j));
//     }
//     return(log_det);
//   }
// }
// int General_Matrix_f2::rows() const {
//   int r;
//   if(isDense) {
//     r = dense.rows(); 
//   } else {
//     r = sparse.rows();
//   }
//   return(r);
// }
// int General_Matrix_f2::cols() const {
//   int c;
//   if(isDense) {
//     c = dense.cols(); 
//   } else {
//     c = sparse.cols();
//   }
//   return(c);
// }
// 
// void load_General_Matrix_f_list2(const Rcpp::List X_list, std::vector<General_Matrix_f2>& X_vector, bool triangular){
//   // null_matrices
//   // MatrixXf null_d = MatrixXf::Zero(0,0);
//   // MatrixXf M_null_d(null_d.data(),0,0);
//   // SpMat null_s = null_d.sparseView();
//   // MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//   
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     X_vector.push_back(load_General_Matrix_f2(Xi_, triangular));
//     // if(Rf_isMatrix(Xi_)){
//     //   MatrixXf Xi = as<MatrixXf >(Xi_);
//     //   General_Matrix_f2 Xim(Xi,M_null_s,true);
//     //   X_vector.push_back(Xim);
//     // } else{
//     //   MSpMat Xi = as<MSpMat>(Xi_);
//     //   General_Matrix_f2 Xim(M_null_d,Xi,false);
//     //   X_vector.push_back(Xim);
//     // }
//   }
// }
// VectorXf sample_MME_single_diagR2(
//     VectorXf y,           // nx1
//     General_Matrix_f2 Z,    // nxr dgCMatrix or dense
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
// // [[Rcpp::export()]]
// MatrixXf sample_MME_ZKZts_c2(
//     MatrixXf Y,                    // nxp
//     SEXP Z_,
//     VectorXf tot_Eta_prec,         // px1
//     Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
//     MatrixXf h2s,                  // n_RE x p
//     VectorXi h2s_index                 // px1
// ) {
//   
//   General_Matrix_f2 Z = load_General_Matrix_f2(Z_,false);
//   
//   int p = Y.cols();
//   int r = Z.cols();
//   // Rcout << Y.rows();
//   // Rcout << Z.rows();
//   // int r = Y.cols();
//   // if(Z.isDense) {
//   //   r = Z.dense.cols();
//   // } else{
//   //   r = Z.sparse.cols();
//   // }
//   
//   MatrixXf randn_theta = rstdnorm_mat2(r,p);
//   
//   std::vector<General_Matrix_f2> chol_ZtZ_Kinv_list;
//   load_General_Matrix_f_list2(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list, true);
//   
//   MatrixXf U = MatrixXf::Zero(r,p);
//   ArrayXf h2_e = 1.0 - h2s.colwise().sum().array();
//   ArrayXf pes = tot_Eta_prec.array() / h2_e.array();
//   
// #pragma omp parallel for
//   for(std::size_t j = 0; j < p; j++){
//     int h2_index = h2s_index[j] - 1;
//     // Rcout << h2_index << std::endl;
//     // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
//     U.col(j) = sample_MME_single_diagR2(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
//     // Rcout << a.size() << std::endl;
//     // U.col(j)
//   }
//   
//   return(U);
// }


// 
// struct Cr {
//   MatrixXf X;
//   Cr(MatrixXf X_): X(X_){}
//   MatrixXf times(MatrixXf Y);
//   MatrixXf operator*(MatrixXf Y);
// };
// MatrixXf Cr::times(MatrixXf Y) {
//   return(X*Y);
// }
// MatrixXf Cr::operator*(MatrixXf Y) {
//   return(X*Y);
// }
// 
// // [[Rcpp::export]]
// MatrixXf test_times(MatrixXf X_, MatrixXf Y) {
//   Cr X(X_);
//   return(X.times(Y));
// }
// 
// // [[Rcpp::export]]
// MatrixXf test_times2(MatrixXf X_, MatrixXf Y) {
//   Cr X(X_);
//   return(X * Y);
// }
// 
// // [[Rcpp::export]]
// MatrixXf test_times3(MatrixXf X, MatrixXf Y) {
//   return(X * Y);
// }
// 
// // [[Rcpp::export]]
// VectorXf regression_sampler_v4b(
//     VectorXf y,
//     MatrixXf& X1,
//     MatrixXf& X2,
//     ArrayXf& diag_X1tX1,
//     ArrayXf& diag_X2tX2,
//     VectorXf a,
//     VectorXf alpha,
//     VectorXf beta,
//     VectorXi delta,
//     float invVarRes,
//     VectorXf invVarEffects, // bx1
//     VectorXf pi,
//     VectorXf randn_a,
//     VectorXf randn_beta,
//     VectorXf rand_unif,
//     VectorXf rgamma_1,
//     float Y_prec_b0,
//     int nIter
// ) {
// 
//   ArrayXf logPi = pi.array().log();
//   ArrayXf logPiComp = (1.0 - pi.array()).log();
//   ArrayXf logDelta0 = logPi;
//   ArrayXf logVarEffects = invVarEffects.array().inverse().log();
//   int nMarkers      = alpha.size();
// 
//   VectorXf yCorr = y - X1*a;
//   for(int j = 0; j < nMarkers; j++) {
//     if(delta[j] != 0) yCorr -= X2.col(j)*alpha[j];
//   }
// 
//   for(int i = 0; i < nIter; i++) {
// 
//     // Sample a
//     for(int j = 0; j < a.size(); j++) {
//       float rhs = (X1.col(j).dot(yCorr) + diag_X1tX1[j]*a[j])*invVarRes;
//       float lhs = diag_X1tX1[j]*invVarRes;
//       float invLhs = 1.0/lhs;
//       float gHat = rhs * invLhs;
//       float old_a = a[j];
//       a[j] = gHat + randn_a(i*a.size() + j)*sqrt(invLhs);
//       yCorr += X1.col(j) * (old_a - a[j]);
//     }
// 
//     int nLoci         = 0;
//     // Sample beta = alpha*delta
//     for(int j = 0; j < nMarkers; j++) {
//       float rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j])*invVarRes;
//       float lhs = diag_X2tX2[j]*invVarRes + invVarEffects[j];
//       float invLhs = 1.0/lhs;
//       float gHat = rhs * invLhs;
//       float logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp[j];
//       float probDelta1 = 1.0 / (1.0 + exp(logDelta0[j] - logDelta1));
//       float oldAlpha = alpha[j];
// 
//       float u = rand_unif(i*nMarkers + j);
//       float r = randn_beta(i*nMarkers + j);
//       if(u < probDelta1) {
//         delta[j] = 1.0;
//         beta[j] = gHat + r*sqrt(invLhs);
//         alpha[j] = beta[j];
//         yCorr += X2.col(j) * (oldAlpha - alpha[j]);
//         nLoci++;
//       } else {
//         if(oldAlpha != 0) {
//           yCorr += X2.col(j) * oldAlpha;
//         }
//         delta[j] = 0;
//         beta[j] = r/sqrt(invVarEffects[j]);
//         alpha[j] = 0;
//       }
//     }
//     invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + Y_prec_b0);
//   }
//   VectorXf result(1+a.size() + 3*nMarkers);
//   result << invVarRes,a,alpha,beta,delta.cast<float>();
//   return(result);
// }
// void test_regression_sampler_parallel(
//   int which_sampler,        //
//   MatrixXf Y,               //
//   MatrixXf X1_base,          //
//   Rcpp::List X1_list_,             // p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
//   SEXP X2_,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
//   SEXP Vx_,                       // m x b matrix if X2 is Ux
//   Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
//   Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//   VectorXf Y_prec,               // p-vector of Y current precisions
//   VectorXf Y_prec_a0,
//   VectorXf Y_prec_b0,
//   MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//   VectorXf prior_prec_alpha2, // a1 x p matrix of prior precisions for alpha1
//   MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
//   MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
//   SEXP current_alpha1s_,
//   SEXP current_alpha2s_,
//   SEXP current_betas_,    // p-list of a2 x 1 vectors
//   Rcpp::List BayesAlphabet_parms
// )
// {
//   Rcout << "asdf" << std::endl;
//   // 
//   
// }

// 
// 
// // [[Rcpp::export]]
// VectorXd add1d(Map<VectorXd> x) {
//   x += VectorXd::Ones(x.size());
//   return(x);
// }
// // [[Rcpp::export]]
// List add1f(VectorXf x) {
//   x += VectorXf::Ones(x.size());
//   // VectorXd xd(x.data())
//   return(List::create(Named("x") = x));
//   // return(x);
// }
// // [[Rcpp::export]]
// List add1f2(List x_) {
//   VectorXf x = as<VectorXf>(x_["x"]);
//   x += VectorXf::Ones(x.size());
//   // VectorXd xd(x.data())
//   return(List::create(Named("x") = x));
//   // return(x);
// }
// 
// 
// // [[Rcpp::export]]
// double regression_sampler_v4(
//     Map<MatrixXd> X, // nxb
//     const ArrayXd& diag_XtX, // bx1
//     Map<VectorXd> yCorr, // nx1
//     Map<VectorXd> alpha, // bx1
//     Map<VectorXd> beta, // bx1
//     Map<VectorXd> delta, // bx1
//     double&  invVarRes,
//     const ArrayXd& varEffects, // bx1
//     double pi,
//     double df,
//     double scale
// ) {
//   double logPi         = log(pi);
//   double logPiComp     = log(1.0-pi);
//   double logDelta0     = logPi;
//   // double invVarRes     = 1.0 / vare;
//   ArrayXd invVarEffects = varEffects.inverse();
//   ArrayXd logVarEffects = varEffects.log();
//   int nLoci         = 0;
//   int nMarkers      = alpha.size();
// 
//   // Sample beta = alpha*delta
//   for(int j = 0; j < nMarkers; j++) {
//     double rhs = (X.col(j).dot(yCorr) + diag_XtX[j]*alpha[j])*invVarRes;
//     double lhs = diag_XtX[j]*invVarRes + invVarEffects[j];
//     double invLhs = 1.0/lhs;
//     double gHat = rhs * invLhs;
//     double logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp;
//     double probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1));
//     double oldAlpha = alpha[j];
// 
//     double u = R::runif(0,1);
//     double r = R::rnorm(0,1);
//     if(u < probDelta1) {
//       delta[j] = 1.0;
//       beta[j] = gHat + r*sqrt(invLhs);
//       alpha[j] = beta[j];
//       yCorr += X.col(j) * (oldAlpha - alpha[j]);
//       nLoci++;
//     } else {
//       if(oldAlpha != 0) {
//         yCorr += X.col(j) * oldAlpha;
//       }
//       delta[j] = 0;
//       beta[j] = r*sqrt(varEffects[j]);
//       alpha[j] = 0;
//     }
//   }
// 
//   // sample invVarRes
//   invVarRes = 1.0 / ((yCorr.dot(yCorr) + df*scale)/R::rchisq(yCorr.size()));
//   return(invVarRes);
// }
// 
// 
// // [[Rcpp::export]]
// float regression_sampler_v4b(
//     MatrixXf X, // nxb
//     const ArrayXf& diag_XtX, // bx1
//     VectorXf yCorr, // nx1
//     VectorXf alpha, // bx1
//     VectorXf beta, // bx1
//     VectorXf delta, // bx1
//     float&  invVarRes,
//     const ArrayXf& varEffects, // bx1
//     float pi,
//     float df,
//     float scale
// ) {
//   float logPi         = log(pi);
//   float logPiComp     = log(1.0-pi);
//   float logDelta0     = logPi;
//   // float invVarRes     = 1.0 / vare;
//   ArrayXf invVarEffects = varEffects.inverse();
//   ArrayXf logVarEffects = varEffects.log();
//   int nLoci         = 0;
//   int nMarkers      = alpha.size();
//   
//   // Sample beta = alpha*delta
//   for(int j = 0; j < nMarkers; j++) {
//     float rhs = (X.col(j).dot(yCorr) + diag_XtX[j]*alpha[j])*invVarRes;
//     float lhs = diag_XtX[j]*invVarRes + invVarEffects[j];
//     float invLhs = 1.0/lhs;
//     float gHat = rhs * invLhs;
//     float logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp;
//     float probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1));
//     float oldAlpha = alpha[j];
//     
//     float u = R::runif(0,1);
//     float r = R::rnorm(0,1);
//     if(u < probDelta1) {
//       delta[j] = 1.0;
//       beta[j] = gHat + r*sqrt(invLhs);
//       alpha[j] = beta[j];
//       yCorr += X.col(j) * (oldAlpha - alpha[j]);
//       nLoci++;
//     } else {
//       if(oldAlpha != 0) {
//         yCorr += X.col(j) * oldAlpha;
//       }
//       delta[j] = 0;
//       beta[j] = r*sqrt(varEffects[j]);
//       alpha[j] = 0;
//     }
//   }
//   
//   // sample invVarRes
//   invVarRes = 1.0 / ((yCorr.dot(yCorr) + df*scale)/R::rchisq(yCorr.size()));
//   return(invVarRes);
// }
// 
// 
// 
// 
// // [[Rcpp::export]]
// double regression_sampler_v4c(
//     Map<MatrixXd> X, // nxb
//     const ArrayXd& diag_XtX, // bx1
//     Map<VectorXd> y, // nx1
//     Map<VectorXd> alpha, // bx1
//     Map<VectorXd> beta, // bx1
//     Map<VectorXd> delta, // bx1
//     double&  invVarRes,
//     const ArrayXd& varEffects, // bx1
//     double pi,
//     double df,
//     double scale
// ) {
//   double logPi         = log(pi);
//   double logPiComp     = log(1.0-pi);
//   double logDelta0     = logPi;
//   // double invVarRes     = 1.0 / vare;
//   ArrayXd invVarEffects = varEffects.inverse();
//   ArrayXd logVarEffects = varEffects.log();
//   int nLoci         = 0;
//   int nMarkers      = alpha.size();
//   
//   VectorXd yCorr = y - X*alpha;
//   
//   // Sample beta = alpha*delta
//   for(int j = 0; j < nMarkers; j++) {
//     double rhs = (X.col(j).dot(yCorr) + diag_XtX[j]*alpha[j])*invVarRes;
//     double lhs = diag_XtX[j]*invVarRes + invVarEffects[j];
//     double invLhs = 1.0/lhs;
//     double gHat = rhs * invLhs;
//     double logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp;
//     double probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1));
//     double oldAlpha = alpha[j];
//     
//     double u = R::runif(0,1);
//     double r = R::rnorm(0,1);
//     if(u < probDelta1) {
//       delta[j] = 1.0;
//       beta[j] = gHat + r*sqrt(invLhs);
//       alpha[j] = beta[j];
//       yCorr += X.col(j) * (oldAlpha - alpha[j]);
//       nLoci++;
//     } else {
//       if(oldAlpha != 0) {
//         yCorr += X.col(j) * oldAlpha;
//       }
//       delta[j] = 0;
//       beta[j] = r*sqrt(varEffects[j]);
//       alpha[j] = 0;
//     }
//   }
//   
//   // sample invVarRes
//   invVarRes = 1.0 / ((yCorr.dot(yCorr) + df*scale)/R::rchisq(yCorr.size()));
//   return(invVarRes);
// }
// 
// 
// // [[Rcpp::export]]
// VectorXf regression_sampler_v5(
//     VectorXf& y,
//     MatrixXf& X1,
//     MatrixXf& X2,
//     ArrayXf& diag_X1tX1,
//     ArrayXf& diag_X2tX2,
//     VectorXf& a,
//     VectorXf& alpha,
//     VectorXf& beta,
//     VectorXi& delta,
//     float invVarRes,
//     const ArrayXf& invVarEffects, // bx1
//     const ArrayXf& pi,
//     MatrixXf& randn_a,
//     MatrixXf& randn_beta,
//     MatrixXf& rand_unif,
//     VectorXf& rgamma_1,
//     float Y_prec_b0,
//     int nIter
// ) {
//   
//   ArrayXf logPi = pi.log();
//   ArrayXf logPiComp = (1.0 - pi).log();
//   ArrayXf logDelta0 = logPi;
//   ArrayXf logVarEffects = invVarEffects.inverse().log();
//   int nMarkers      = alpha.size();
//   
//   VectorXf yCorr = y - X1*a;
//   for(int j = 0; j < nMarkers; j++) {
//     if(delta[j] != 0) yCorr -= X2.col(j)*alpha[j];
//   }
//   
//   for(int i = 0; i < nIter; i++) {
//     
//     // Sample a
//     for(int j = 0; j < a.size(); j++) {
//       float rhs = (X1.col(j).dot(yCorr) + diag_X1tX1[j]*a[j])*invVarRes;
//       float lhs = diag_X1tX1[j]*invVarRes;
//       float invLhs = 1.0/lhs;
//       float gHat = rhs * invLhs;
//       float old_a = a[j];
//       a[j] = gHat + randn_a(j,i)*sqrt(invLhs);
//       yCorr += X1.col(j) * (old_a - a[j]);
//     }
//     
//     int nLoci         = 0;
//     // Sample beta = alpha*delta
//     for(int j = 0; j < nMarkers; j++) {
//       float rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j])*invVarRes;
//       float lhs = diag_X2tX2[j]*invVarRes + invVarEffects[j];
//       float invLhs = 1.0/lhs;
//       float gHat = rhs * invLhs;
//       float logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp[j];
//       float probDelta1 = 1.0 / (1.0 + exp(logDelta0[j] - logDelta1));
//       float oldAlpha = alpha[j];
//       
//       float u = rand_unif(j,i);
//       float r = randn_beta(j,i);
//       if(u < probDelta1) {
//         delta[j] = 1.0;
//         beta[j] = gHat + r*sqrt(invLhs);
//         alpha[j] = beta[j];
//         yCorr += X2.col(j) * (oldAlpha - alpha[j]);
//         nLoci++;
//       } else {
//         if(oldAlpha != 0) {
//           yCorr += X2.col(j) * oldAlpha;
//         }
//         delta[j] = 0;
//         beta[j] = r/sqrt(invVarEffects[j]);
//         alpha[j] = 0;
//       }
//     }
//     invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + Y_prec_b0);
//   }
//   VectorXf result(1+a.size() + 3*nMarkers);
//   result << invVarRes,a,alpha,beta,delta.cast<float>();
//   return(result);
// }
// 
// // [[Rcpp::export]]
// VectorXd regression_sampler_v5b(
//     Map<VectorXd> y,
//     Map<MatrixXd> X1,
//     Map<MatrixXd> X2,
//     ArrayXd& diag_X1tX1,
//     ArrayXd& diag_X2tX2,
//     Map<VectorXd> a,
//     Map<VectorXd> alpha,
//     Map<VectorXd> beta,
//     VectorXi& delta,
//     double invVarRes,
//     const ArrayXd& invVarEffects, // bx1
//     const ArrayXd& pi,
//     Map<MatrixXd> randn_a,
//     Map<MatrixXd> randn_beta,
//     Map<MatrixXd> rand_unif,
//     Map<VectorXd> rgamma_1,
//     double Y_prec_b0,
//     int nIter
// ) {
//   
//   ArrayXd logPi = pi.log();
//   ArrayXd logPiComp = (1.0 - pi).log();
//   ArrayXd logDelta0 = logPi;
//   ArrayXd logVarEffects = invVarEffects.inverse().log();
//   int nMarkers      = alpha.size();
//   
//   VectorXd yCorr = y - X1*a - X2*alpha;
//   
//   for(int i = 0; i < nIter; i++) {
//     
//     // Sample a
//     for(int j = 0; j < a.size(); j++) {
//       double rhs = (X1.col(j).dot(yCorr) + diag_X1tX1[j]*a[j])*invVarRes;
//       double lhs = diag_X1tX1[j]*invVarRes;
//       double invLhs = 1.0/lhs;
//       double gHat = rhs * invLhs;
//       double old_a = a[j];
//       a[j] = gHat + randn_a(j,i)*sqrt(invLhs);
//       yCorr += X1.col(j) * (old_a - a[j]);
//     }
//     
//     int nLoci         = 0;
//     // Sample beta = alpha*delta
//     for(int j = 0; j < nMarkers; j++) {
//       double rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j])*invVarRes;
//       double lhs = diag_X2tX2[j]*invVarRes + invVarEffects[j];
//       double invLhs = 1.0/lhs;
//       double gHat = rhs * invLhs;
//       double logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp[j];
//       double probDelta1 = 1.0 / (1.0 + exp(logDelta0[j] - logDelta1));
//       double oldAlpha = alpha[j];
//       
//       double u = rand_unif(j,i);
//       double r = randn_beta(j,i);
//       if(u < probDelta1) {
//         delta[j] = 1.0;
//         beta[j] = gHat + r*sqrt(invLhs);
//         alpha[j] = beta[j];
//         yCorr += X2.col(j) * (oldAlpha - alpha[j]);
//         nLoci++;
//       } else {
//         if(oldAlpha != 0) {
//           yCorr += X2.col(j) * oldAlpha;
//         }
//         delta[j] = 0;
//         beta[j] = r/sqrt(invVarEffects[j]);
//         alpha[j] = 0;
//       }
//     }
//     invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + Y_prec_b0);
//   }
//   VectorXd result(1+a.size() + 3*nMarkers);
//   result << invVarRes,a,alpha,beta,delta.cast<double>();
//   return(result);
// }
// 
// // [[Rcpp::export]]
// VectorXd regression_sampler_v5c(
//     Map<VectorXd> y,
//     Map<MatrixXd> X1,
//     Map<MatrixXd> X2,
//     ArrayXd& diag_X1tX1,
//     ArrayXd& diag_X2tX2,
//     Map<VectorXd> a,
//     Map<VectorXd> alpha,
//     Map<VectorXd> beta,
//     VectorXi& delta,
//     double invVarRes,
//     const ArrayXd& invVarEffects, // bx1
//     const ArrayXd& pi,
//     Map<VectorXd> rgamma_1,
//     double Y_prec_b0,
//     int nIter
// ) {
//   
//   ArrayXd logPi = pi.log();
//   ArrayXd logPiComp = (1.0 - pi).log();
//   ArrayXd logDelta0 = logPi;
//   ArrayXd logVarEffects = invVarEffects.inverse().log();
//   int nMarkers      = alpha.size();
//   
//   VectorXd yCorr = y - X1*a - X2*alpha;
//   
//   for(int i = 0; i < nIter; i++) {
//     
//     // Sample a
//     for(int j = 0; j < a.size(); j++) {
//       double rhs = (X1.col(j).dot(yCorr) + diag_X1tX1[j]*a[j])*invVarRes;
//       double lhs = diag_X1tX1[j]*invVarRes;
//       double invLhs = 1.0/lhs;
//       double gHat = rhs * invLhs;
//       double old_a = a[j];
//       a[j] = gHat + R::rnorm(0,1)*sqrt(invLhs);
//       yCorr += X1.col(j) * (old_a - a[j]);
//     }
//     
//     int nLoci         = 0;
//     // Sample beta = alpha*delta
//     for(int j = 0; j < nMarkers; j++) {
//       double rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j])*invVarRes;
//       double lhs = diag_X2tX2[j]*invVarRes + invVarEffects[j];
//       double invLhs = 1.0/lhs;
//       double gHat = rhs * invLhs;
//       double logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp[j];
//       double probDelta1 = 1.0 / (1.0 + exp(logDelta0[j] - logDelta1));
//       double oldAlpha = alpha[j];
//       
//       double u =  R::runif(0,1);
//       double r = R::rnorm(0,1);
//       if(u < probDelta1) {
//         delta[j] = 1.0;
//         beta[j] = gHat + r*sqrt(invLhs);
//         alpha[j] = beta[j];
//         yCorr += X2.col(j) * (oldAlpha - alpha[j]);
//         nLoci++;
//       } else {
//         if(oldAlpha != 0) {
//           yCorr += X2.col(j) * oldAlpha;
//         }
//         delta[j] = 0;
//         beta[j] = r/sqrt(invVarEffects[j]);
//         alpha[j] = 0;
//       }
//     }
//     invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + Y_prec_b0);
//   }
//   VectorXd result(1+a.size() + 3*nMarkers);
//   result << invVarRes,a,alpha,beta,delta.cast<double>();
//   return(result);
// }
                            
  
  // // Sample beta = alpha*delta
  // for(int j = 0; j < nMarkers; j++) {
  //   float rhs = (X.col(j).dot(yCorr) + diag_XtX[j]*alpha[j])*invVarRes;
  //   float lhs = diag_XtX[j]*invVarRes + invVarEffects[j];
  //   float invLhs = 1.0/lhs;
  //   float gHat = rhs * invLhs;
  //   float logDelta1 = -0.5*(log(lhs) + logVarEffects[j] - gHat*rhs) + logPiComp;
  //   float probDelta1 = 1.0 / (1.0 + exp(logDelta0 - logDelta1));
  //   float oldAlpha = alpha[j];
    
    // float u = R::runif(0,1);
    // float r = R::rnorm(0,1);
  //   if(u < probDelta1) {
  //     delta[j] = 1.0;
  //     beta[j] = gHat + r*sqrt(invLhs);
  //     alpha[j] = beta[j];
  //     yCorr += X.col(j) * (oldAlpha - alpha[j]);
  //     nLoci++;
  //   } else {
  //     if(oldAlpha != 0) {
  //       yCorr += X.col(j) * oldAlpha;
  //     }
  //     delta[j] = 0;
  //     beta[j] = r*sqrt(varEffects[j]);
  //     alpha[j] = 0;
  //   }
  // }
  
  // sample invVarRes
//   invVarRes = 1.0 / ((yCorr.dot(yCorr) + df*scale)/R::rchisq(yCorr.size()));
//   return(invVarRes);
// }

//   
// 
// VectorXd regression_sampler_v12(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& X1,           // nxa
//     const MatrixXd& RinvtX2,                // nxb
//     const MatrixXd& C,                     // bxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix. Also called R here with RtR = V
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const double rgamma_1,
//     const double Y_prec_b0
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
//   MatrixXd C_beta = C;
//   C_beta.diagonal() += prior_prec_beta;
//   LLT<MatrixXd> A_beta_llt;
//   A_beta_llt.compute(C_beta);
//   MatrixXd chol_A_beta = A_beta_llt.matrixU();
//   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
//   
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // Calculate A_alpha = Y_prec*X1^T*V_beta^{-1}*X1 + D_alpha^{-1}
//     // We don't need to actually calculate V_beta^{-1} directly.
//     // Instead,
//     MatrixXd RinvtX1 = get_RinvtX22(chol_V,X1);  // n*n*a -> n x a
//     MatrixXd X1tVinvX2 = RinvtX1.transpose() * RinvtX2; // a*n*b -> a*b
//     MatrixXd cholAbetainvt_X2tVinvX1 = chol_A_beta.transpose().triangularView<Lower>().solve(X1tVinvX2.transpose()); // b*b*a -> b x a
//     
//     MatrixXd A_alpha = RinvtX1.transpose() * RinvtX1 - cholAbetainvt_X2tVinvX1.transpose() * cholAbetainvt_X2tVinvX1;
//     A_alpha.diagonal() += prior_prec_alpha;
//     
//     VectorXd Rinvsqy = get_RinvtX22(chol_V,y); // n*n*q -> n x 1;
//     VectorXd XtRinvy = RinvtX2.transpose() * Rinvsqy; // b*n*1 >- b x 1
//     VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
//     
//     VectorXd X1t_V_beta_inv_y = RinvtX1.transpose() * Rinvsqy - cholAbetainvt_X2tVinvX1.transpose() * invSqAbXtRinvy;
//     
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
//     
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     
//     y_tilde = y - X1 * alpha;
//   }
//   
//   // Step 2 - sample Y_prec
//   // We don't need to actually calculate V_beta^{-1} directly.
//   VectorXd RinvSqy = get_RinvtX22(chol_V,y_tilde);
//   VectorXd XtRinvy = RinvtX2.transpose() * RinvSqy;
//   VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
//   double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
//   Y_prec = rgamma_1/score;
//   
//   // Step 3 - sample beta
//   VectorXd XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
//   VectorXd beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
//   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
//   
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
//   
//   return(result);
// }



// 
// 
// 
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
//   return(X_mat);
// }
// 
// 
// 
// // code to convert list of R matrices (sparse or dense) into a thread-safe object
// // struct R_matrix {
// //   Map<MatrixXd> dense;
// //   MSpMat sparse;
// //   bool isDense;
// //   R_matrix(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// // };
// 
// // Loads a sparse or dense matrix passed from R into a \code{R_matrix} object
// R_matrix load_R_matrix2(SEXP X_) {
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   if(Rf_isMatrix(X_)){
//     Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
//     SpMat null_s = null_d.sparseView();
//     MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//     R_matrix Xm(X,M_null_s,true);
//     return(Xm);
//   } else{
//     MSpMat X = as<MSpMat>(X_);
//     Map<MatrixXd> M_null_d(null_d.data(),0,0);
//     R_matrix Xm(M_null_d,X,false);
//     return(Xm);
//   }
// }
// 
// // Code to convert list of R matrices (sparse or dense) into a thread-safe object
// //
// // @param X_list List of matrices (each can be dgCMatrix or matrix)
// // @param X_vector \code{std::vector} of \code{R_matrix} type which2 will be populated
// void load_R_matrices_list2(const Rcpp::List X_list, std::vector<R_matrix>& X_vector){
//   // null_matrices
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   Map<MatrixXd> M_null_d(null_d.data(),0,0);
//   SpMat null_s = null_d.sparseView();
//   MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
// 
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     X_vector.push_back(load_R_matrix2(Xi_));
//     // if(Rf_isMatrix(Xi_)){
//     //   Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
//     //   R_matrix Xim(Xi,M_null_s,true);
//     //   X_vector.push_back(Xim);
//     // } else{
//     //   MSpMat Xi = as<MSpMat>(Xi_);
//     //   R_matrix Xim(M_null_d,Xi,false);
//     //   X_vector.push_back(Xim);
//     // }
//   }
// }
// 
// // -------------------------------------------- //
// // ---------- regression_sampler --------- //
// // -------------------------------------------- //
// 
// // Replicates the \code{which2} function from R
// //
// // @param x Logical vector
// // @return IntegerVector with indices of the TRUE values of \code{x}
// Rcpp::IntegerVector which2(Rcpp::LogicalVector x) {
//   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
//   return v[x];
// }
// 
// // Replicates the \code{which2} function from R
// //
// // @param chol_V \code{R_matrix} object with the upper-triangular Cholesky decomposition of a square nxn matrix R
// // @param X nxp matrix
// // @return solve(t(R),X) as a dense matrix
// MatrixXd get_RinvtX22(const R_matrix& chol_V, MatrixXd X){
//   MatrixXd RinvtX2;
//   if(chol_V.isDense) {
//     RinvtX2 = chol_V.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   } else{
//     RinvtX2 = chol_V.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   }
//   return(RinvtX2);
// }
// 
// VectorXd regression_sampler_v12(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& X1,           // nxa
//     const MatrixXd& RinvtX2,                // nxb
//     const MatrixXd& C,                     // bxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix. Also called R here with RtR = V
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const double rgamma_1,
//     const double Y_prec_b0
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
//   MatrixXd C_beta = C;
//   C_beta.diagonal() += prior_prec_beta;
//   LLT<MatrixXd> A_beta_llt;
//   A_beta_llt.compute(C_beta);
//   MatrixXd chol_A_beta = A_beta_llt.matrixU();
//   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
// 
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // Calculate A_alpha = Y_prec*X1^T*V_beta^{-1}*X1 + D_alpha^{-1}
//     // We don't need to actually calculate V_beta^{-1} directly.
//     // Instead,
//     MatrixXd RinvtX1 = get_RinvtX22(chol_V,X1);  // n*n*a -> n x a
//     MatrixXd X1tVinvX2 = RinvtX1.transpose() * RinvtX2; // a*n*b -> a*b
//     MatrixXd cholAbetainvt_X2tVinvX1 = chol_A_beta.transpose().triangularView<Lower>().solve(X1tVinvX2.transpose()); // b*b*a -> b x a
// 
//     MatrixXd A_alpha = RinvtX1.transpose() * RinvtX1 - cholAbetainvt_X2tVinvX1.transpose() * cholAbetainvt_X2tVinvX1;
//     A_alpha.diagonal() += prior_prec_alpha;
// 
//     VectorXd Rinvsqy = get_RinvtX22(chol_V,y); // n*n*q -> n x 1;
//     VectorXd XtRinvy = RinvtX2.transpose() * Rinvsqy; // b*n*1 >- b x 1
//     VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
// 
//     VectorXd X1t_V_beta_inv_y = RinvtX1.transpose() * Rinvsqy - cholAbetainvt_X2tVinvX1.transpose() * invSqAbXtRinvy;
// 
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// 
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// 
//     y_tilde = y - X1 * alpha;
//   }
// 
//   // Step 2 - sample Y_prec
//   // We don't need to actually calculate V_beta^{-1} directly.
//   VectorXd RinvSqy = get_RinvtX22(chol_V,y_tilde);
//   VectorXd XtRinvy = RinvtX2.transpose() * RinvSqy;
//   VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
//   double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
//   Y_prec = rgamma_1/score;
// 
//   // Step 3 - sample beta
//   VectorXd XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
//   VectorXd beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
//   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
// 
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
// 
//   return(result);
// }
// 
// 
// VectorXd regression_sampler_v22(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& X1,           // nxa
//     const MatrixXd& X2,           // nxm or nxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& V,
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const Ref<const VectorXd>& randn_e,
//     const double rgamma_1,
//     const double Y_prec_b0
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
//   MatrixXd Dbeta_X2t = prior_prec_beta.cwiseInverse().asDiagonal() * X2.transpose();
//   MatrixXd V_beta = X2 * Dbeta_X2t + V;
//   LDLT<MatrixXd> V_beta_ldlt;
//   V_beta_ldlt.compute(V_beta);
// 
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     MatrixXd V_beta_inv_X1 = V_beta_ldlt.solve(X1);
//     MatrixXd A_alpha = V_beta_inv_X1.transpose() * X1;
//     A_alpha.diagonal() += prior_prec_alpha;
// 
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// 
//     VectorXd X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     y_tilde = y - X1 * alpha;
//   }
// 
//   // Step 2 - sample Y_prec
//   VectorXd e2 = y_tilde.transpose() * V_beta_ldlt.solve(y_tilde);
//   double score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
// 
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXd v = sqrt(Y_prec) * X2 * u;
//   if(chol_V.isDense) {
//     v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   VectorXd w = V_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
//   VectorXd beta = u + Dbeta_X2t * w / sqrt(Y_prec);
// 
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
// 
//   return(result);
// }
// 
// 
// VectorXd regression_sampler_v32(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& X1,           // nxa
//     const MatrixXd& Ux,           // nxm or nxb
//     const MatrixXd& Vx,           // mxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_V,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& Vinv,
//     const MatrixXd& VinvUx,
//     const MatrixXd& UtVinvU,
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const Ref<const VectorXd>& randn_e,
//     const double rgamma_1,
//     const double Y_prec_b0
// ){
//   int n = y.size();
//   int a = X1.cols();
//   if(Vx.rows() != Ux.cols()) stop("Wrong dimensions of Vx");
//   MatrixXd X2 = Ux*Vx;
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
//   MatrixXd Dbeta_Vxt = prior_prec_beta.cwiseInverse().asDiagonal() * Vx.transpose();
//   MatrixXd VxDbeta_Vxt = Vx * Dbeta_Vxt;
//   if(VinvUx.rows() != n) stop("Wrong dimensions of VinvUx");
//   MatrixXd inner = VxDbeta_Vxt * UtVinvU;
//   inner.diagonal().array() += 1.0;
//   LDLT<MatrixXd> inner_ldlt;
//   inner_ldlt.compute(inner);
//   // Sigma_beta_inv = Vinv - VinvUx * inner.ldlt().solve(VxDbeta_Vxt * VinvUx.transpose());  // Don't actually calculate this. Stay in mxm space
// 
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // MatrixXd V_beta_inv_X1 = Sigma_beta_inv * X1;
//     // MatrixXd A_alpha = Y_prec * V_beta_inv_X1.transpose() * X1;
//     MatrixXd Vinv_X1 = Vinv * X1;
//     MatrixXd Uxt_Vinv_X1 = Ux.transpose() * Vinv_X1;
//     MatrixXd A_alpha = X1.transpose() * Vinv_X1 - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt) * Uxt_Vinv_X1;
//     A_alpha.diagonal() += prior_prec_alpha;
// 
//     // VectorXd X1t_V_beta_inv_y = V_beta_inv_X1.transpose() * y;
//     VectorXd Uxt_Vinv_y = VinvUx.transpose() * y;
//     VectorXd X1t_V_beta_inv_y = Vinv_X1.transpose() * y - Uxt_Vinv_X1.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
// 
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// 
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec) * randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// 
//     y_tilde = y - X1 * alpha;
//   }
// 
//   // Step 2 - sample Y_prec
//   // VectorXd e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
//   VectorXd Vinv_y_tilde = Vinv * y_tilde;
//   VectorXd Uxt_Vinv_y = VinvUx.transpose() * y_tilde;
//   VectorXd e2 = y_tilde.transpose() * Vinv_y_tilde - Uxt_Vinv_y.transpose() * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_y);
// 
//   double score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
// 
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXd v = std::sqrt(Y_prec) * X2 * u;
//   if(chol_V.isDense) {
//     v += chol_V.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_V.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   // VectorXd X1 = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
//   VectorXd e = y_tilde * std::sqrt(Y_prec) - v;
//   VectorXd Uxt_Vinv_e = VinvUx.transpose() * e;
//   VectorXd w = Vinv * e - VinvUx * inner_ldlt.solve(VxDbeta_Vxt * Uxt_Vinv_e);
// 
//   VectorXd beta = u + Dbeta_Vxt * (Ux.transpose() * w) / std::sqrt(Y_prec); //b*b*1 + b*n*1
// 
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
// 
//   return(result);
// }
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
// // [[Rcpp::export]]
// Rcpp::List regression_sampler_parallel2(
//     Map<MatrixXd> Y,               //
//     Map<MatrixXd> X1_base,          //
//     Rcpp::List X1_list_,             // p-list of n x a2 matrices of X1 covariates unique to each p. Can be NULL
//     Map<MatrixXd> X2,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
//     SEXP Vx_,                       // m x b matrix if X2 is Ux
//     Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
//     Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//     VectorXd Y_prec,               // p-vector of Y current precisions
//     VectorXd Y_prec_a0,
//     VectorXd Y_prec_b0,
//     Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
//     Map<MatrixXd> prior_prec_beta // b x p matrix of prior precisions of beta
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
//   std::vector<R_matrix> X1_list;
//   load_R_matrices_list2(X1_list_, X1_list);
//   if(X1_list.size() > 0) {
//     if(X1_list.size() != p) stop("Wrong length of X1_list");
//   }
// 
//   // X2 or Ux and V
//   Map<MatrixXd> Ux = X2;
//   MatrixXd z = MatrixXd::Zero(0,0);
//   Map<MatrixXd> Vx(z.data(),0,0);
//   int b = X2.cols();
//   if(X2.rows() != n) stop("Wrong dimension of X2");
//   if(Rf_isMatrix(Vx_)) {
//     // Map<MatrixXd> V__ = as<Map<MatrixXd> >(Vx_);
//     // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
//     new (&Vx) Map<MatrixXd> (as<Map<MatrixXd> >(Vx_));
//     if(Ux.cols() != Vx.rows()) stop("X2 and Vx_ have incompatible dimensions");
//     b = Vx.cols();
//   }
//   
//   // chol_V_list
//   std::vector<R_matrix> chol_V_list;
//   load_R_matrices_list2(chol_V_list_, chol_V_list);
//   if(max(h2s_index) > chol_V_list.size()) {
//     stop("max(h2s_index) > length(chol_V_list)");
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
//   MatrixXd randn_alpha1 = rstdnorm_mat2(a1,p);
//   std::vector<VectorXd> randn_alpha2;
//   if(X1_list.size() > 0){
//     for(int i = 0; i < p; i++){
//       randn_alpha2.push_back(rstdnorm_mat2(X1_list[i].dense.cols(),1));
//     }
//   }
//   MatrixXd randn_beta = rstdnorm_mat2(b,p);
//   MatrixXd randn_e;
//   if(b > n) {
//     randn_e = rstdnorm_mat2(n,p);
//   }
//   VectorXd rgamma_1(p);
//   for(int i = 0; i < p; i++) {
//     rgamma_1(i) = rgamma(1,Y_prec_a0[i] + n/2.0,1.0)(0);
//   }
// 
//   // Results structures
//   MatrixXd alpha1(a1,p);
//   std::vector<VectorXd> alpha2;
//   alpha2.reserve(X1_list.size());
//   int alpha2_size = 0;
//   if(X1_list.size() > 0){
//     for(int i = 0; i < X1_list.size(); i++){
//       int a2 = X1_list[i].dense.cols();
//       alpha2.push_back(VectorXd::Zero(a2));
//       alpha2_size += a2;
//     }
//   }
//   MatrixXd beta(b,p);
//   
//   // go through h2s indices and sample columns with same index as a set
//   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//     int h2_index = i;
//     VectorXi trait_set = as<VectorXi>(which2(h2s_index == h2_index));  // list of traits with same h2_index
// 
//     if(trait_set.size() > 0){
//       // prepare matrices for sampler
//       MatrixXd RinvtX2, C, V, Vinv, VinvUx, UtVinvU;
//       R_matrix chol_V = chol_V_list[h2_index - 1];
//       int which_sampler;
//       // Decide which2 sampler to use
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
//           Vinv = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           VinvUx = chol_V.dense.triangularView<Upper>().solve(chol_V.dense.transpose().triangularView<Lower>().solve(Ux));
//         } else{
//           Vinv = chol_V.sparse.triangularView<Upper>().solve(chol_V.sparse.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
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
//         MatrixXd X1;
//         int a;
//         int a2 = 0;
//         int b;
//         VectorXd prior_prec_alpha;
//         VectorXd randn_alpha;
//         if(X1_list.size() == 0) {
//           X1 = X1_base;
//           a = a1;
//           prior_prec_alpha = prior_prec_alpha1.col(j);
//           randn_alpha = randn_alpha1.col(j);
//         } else{
//           Map<MatrixXd> X12 = X1_list[j].dense;
//           a2 = X12.cols();
//           a = a1+a2;
//           X1 = MatrixXd(n,a);
//           X1 << X1_base,X12;
//           prior_prec_alpha = VectorXd(a);
//           prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
//           prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
//           randn_alpha = VectorXd(a);
//           randn_alpha.head(a1) = randn_alpha1.col(j);
//           randn_alpha.tail(a2) = randn_alpha2[j];
//         }
// 
//         VectorXd samples;
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
//   VectorXd alpha2_vec(alpha2_size);
//   if(X1_list.size() > 0){
//     int index = 0;
//     for(int i = 0; i < X1_list.size(); i++){
//       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
//       index += alpha2[i].size();
//     }
//   }
// 
//   return(Rcpp::List::create(
//       Named("alpha1") = alpha1,
//       Named("alpha2") = alpha2_vec,
//       Named("beta") = beta,
//       Named("Y_prec") = Y_prec
//   ));
// }
// 
// 
// 
// 
// 
// 
// VectorXd regression_sampler_v33(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& W,           // nxa
//     const MatrixXd& U,           // nxm or nxb
//     const MatrixXd& V,           // mxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& Rinv,
//     const MatrixXd& RinvU,
//     const MatrixXd& UtRinvU,
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const Ref<const VectorXd>& randn_e,
//     const double rgamma_1,
//     const double Y_prec_b0
// ){
//   int n = y.size();
//   int a = W.cols();
//   if(V.rows() != U.cols()) stop("Wrong dimensions of V");
//   MatrixXd X = U*V;
//   int b = X.cols();
//   
//   // Check inputs
//   if(W.rows() != n) stop("Wrong dimension of W");
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//   if(randn_e.size() != n) stop("Wrong length of randn_e");
//   
//   // Calculate inverse of Sigma_beta
//   // MatrixXd Sigma_beta_inv;
//   // Using Ainv - Ainv * U * (I + BVAinvU)inv * BVAinv in case B = VDVt is singular
//   MatrixXd DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
//   MatrixXd VDVt = V * DVt;
//   if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
//   MatrixXd inner = VDVt * UtRinvU;
//   inner.diagonal().array() += 1.0;
//   LDLT<MatrixXd> inner_ldlt;
//   inner_ldlt.compute(inner);
//   // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(VDVt * RinvU.transpose());  // Don't actually calculate this. Stay in mxm space
//   
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // MatrixXd SbinvW = Sigma_beta_inv * W;
//     // MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
//     MatrixXd RinvW = Rinv * W;
//     MatrixXd UtRinvW = U.transpose() * RinvW;
//     MatrixXd A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvW));
//     A_alpha.diagonal() += prior_prec_alpha;
//     
//     // VectorXd WtSbinvy = SbinvW.transpose() * y;
//     VectorXd UtRinvy = RinvU.transpose() * y;
//     VectorXd WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
//     
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
//     
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     
//     y_tilde = y - W * alpha;
//   }
//   
//   // Step 2 - sample Y_prec
//   // VectorXd e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
//   VectorXd Rinv_y = Rinv * y_tilde;
//   VectorXd UtRinvy = RinvU.transpose() * y_tilde;
//   VectorXd e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
//   
//   double score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
//   
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXd v = std::sqrt(Y_prec) * X * u;
//   if(chol_R.isDense) {
//     v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   // VectorXd w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
//   VectorXd e = y_tilde * std::sqrt(Y_prec) - v;
//   VectorXd UtRinve = RinvU.transpose() * e;
//   VectorXd w = Rinv * e - RinvU * inner_ldlt.solve(VDVt * UtRinve);
//   
//   VectorXd beta = u + DVt * (U.transpose() * w) / std::sqrt(Y_prec); //b*b*1 + b*n*1
//   
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
//   
//   return(result);
// }
// 
// //' Draws samples from all ``fixed" coefficients (fixed and random) of a set of parallel linear regression models, conditional on the variance components.
// //' 
// //' The model is either: \itemize{
// //' \item y_i = W_base*alpha1 + W_list_[i]*alpha2 + X*beta + e, e ~ N(0,1/Y_prec[i]*V)
// //' \item y_i = W_base*alpha1 + W_list_[i]*alpha2 + X*V_*beta + e, e ~ N(0,1/Y_prec[i]*V)
// //' }
// //' where \code{V = RtR}, priors on elements of alpha1, alpha2 and beta are independent.
// //' Each column of Y is considered independent
// //' 
// //' @param Y n x p matrix of observations
// //' @param W_base n x a1 matrix of W covariates common to all p. Can be NULL
// //' @param W_list_ p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
// //' @param X either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
// //' @param V_ m x b matrix if X is U, otherwise NULL
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
// // [[Rcpp::export]]
// Rcpp::List regression_sampler_parallel3(
//     Map<MatrixXd> Y,               // 
//     Map<MatrixXd> W_base,          // 
//     Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
//     Map<MatrixXd> X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
//     SEXP V_,                       // m x b matrix if X is U
//     Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
//     Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//     VectorXd Y_prec,               // p-vector of Y current precisions
//     double Y_prec_a0,
//     double Y_prec_b0,
//     Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
//     Map<MatrixXd> prior_prec_beta // b x p matrix of prior precisions of beta
// ) {
//   
//   int n = Y.rows();
//   int p = Y.cols();
//   
//   // W_base
//   if(W_base.rows() != n) stop("Wrong dimension of W_base");
//   int a1 = W_base.cols();
//   
//   // W_list
//   std::vector<R_matrix> W_list;
//   load_R_matrices_list2(W_list_, W_list);
//   if(W_list.size() > 0) {
//     if(W_list.size() != p) stop("Wrong length of W_list");
//   }
//   
//   // X or U and V
//   Map<MatrixXd> U = X;
//   MatrixXd z = MatrixXd::Zero(0,0);
//   Map<MatrixXd> V(z.data(),0,0);
//   int b = X.cols();
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(Rf_isMatrix(V_)) {
//     // Map<MatrixXd> V__ = as<Map<MatrixXd> >(V_);
//     // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
//     new (&V) Map<MatrixXd> (as<Map<MatrixXd> >(V_));
//     if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
//     b = V.cols();
//   }
//   
//   // chol_V_list
//   std::vector<R_matrix> chol_V_list;
//   load_R_matrices_list2(chol_V_list_, chol_V_list);
//   if(max(h2s_index) > chol_V_list.size()) {
//     stop("max(h2s_index) > length(chol_V_list)");
//   }
//   
//   // priors
//   if(Y_prec.size() != p) {
//     stop("Wrong length of Y_prec");
//   }
//   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
//   if(W_list.size() > 0 && prior_prec_alpha2.size() != p) {
//     stop("Wrong length of prior_prec_alpha2");
//   }
//   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
//   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
//   
//   // generate random numbers
//   MatrixXd randn_alpha1 = rstdnorm_mat2(a1,p);
//   std::vector<VectorXd> randn_alpha2;
//   if(W_list.size() > 0){
//     for(int i = 0; i < p; i++){
//       randn_alpha2.push_back(rstdnorm_mat2(W_list[i].dense.cols(),1));
//     }
//   }
//   MatrixXd randn_beta = rstdnorm_mat2(b,p);
//   MatrixXd randn_e;
//   if(b > n) {
//     randn_e = rstdnorm_mat2(n,p);
//   }
//   VectorXd rgamma_1 = as<VectorXd>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
//   
//   // Results structures
//   MatrixXd alpha1(a1,p);
//   std::vector<VectorXd> alpha2;
//   alpha2.reserve(W_list.size());
//   int alpha2_size = 0;
//   if(W_list.size() > 0){
//     for(int i = 0; i < W_list.size(); i++){
//       int a2 = W_list[i].dense.cols();
//       alpha2.push_back(VectorXd::Zero(a2));
//       alpha2_size += a2;
//     }
//   }
//   MatrixXd beta(b,p);
//   
//   // go through h2s indices and sample columns with same index as a set
//   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//     int h2_index = i;
//     VectorXi trait_set = as<VectorXi>(which2(h2s_index == h2_index));  // list of traits with same h2_index
//     
//     if(trait_set.size() > 0){
//       // prepare matrices for sampler
//       MatrixXd RinvSqX, C, R, Rinv, RinvU, UtRinvU;
//       R_matrix chol_R = chol_V_list[h2_index - 1];
//       int which_sampler;
//       // Decide which2 sampler to use
//       if(b <= n) {
//         // use regression_sampler_v12
//         which_sampler = 1;
//         if(chol_R.isDense) {
//           RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
//         } else{
//           RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
//         }
//         C = RinvSqX.transpose() * RinvSqX;
//       }
//       else if(V.cols() == 0) {
//         // use regression_sampler_v22
//         which_sampler = 2;
//         if(chol_R.isDense) {
//           R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
//         } else{
//           R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
//         }
//       } else {
//         // use regression_sampler_v33
//         which_sampler = 3;
//         if(chol_R.isDense) {
//           Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
//         } else{
//           Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
//           // RinvU = Rinv * U;
//           // Rcout << i << std::endl;
//           // Rcout << Rinv.diagonal().transpose() << std::endl;
//         }
//         UtRinvU = U.transpose() * RinvU;
//       }
//       // REprintf("Number of threads=%i\\n", omp_get_max_threads());
//       // Rcout <<trait_set.size() << " " << omp_get_max_threads() << std::endl;
// #pragma omp parallel for
//       for(int i = 0; i < trait_set.size(); i++){
//         int j = trait_set[i];
//         MatrixXd W;
//         int a;
//         int a2 = 0;
//         int b;
//         VectorXd prior_prec_alpha;
//         VectorXd randn_alpha;
//         if(W_list.size() == 0) {
//           W = W_base;
//           a = a1;
//           prior_prec_alpha = prior_prec_alpha1.col(j);
//           randn_alpha = randn_alpha1.col(j);
//         } else{
//           Map<MatrixXd> W2 = W_list[j].dense;
//           a2 = W2.cols();
//           a = a1+a2;
//           W = MatrixXd(n,a);
//           W << W_base,W2;
//           prior_prec_alpha = VectorXd(a);
//           prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
//           prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
//           randn_alpha = VectorXd(a);
//           randn_alpha.head(a1) = randn_alpha1.col(j);
//           randn_alpha.tail(a2) = randn_alpha2[j];
//         }
//         
//         VectorXd samples;
//         if(which_sampler == 1) {
//           b = RinvSqX.cols();
//           // samples = regression_sampler_v12(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
//           //                                 prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
//           //                                 randn_beta.col(j), rgamma_1[j],Y_prec_b0);
//         } else if(which_sampler == 2) {
//           b = X.cols();
//           // samples = regression_sampler_v22(Y.col(j), W, X, prior_prec_alpha, prior_mean_beta.col(j),
//           //                                 prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
//           //                                 randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
//         } else if(which_sampler == 3) {
//           b = V.cols();
//           samples = regression_sampler_v33(Y.col(j), W, X, V, prior_prec_alpha, prior_mean_beta.col(j),
//                                           prior_prec_beta.col(j), chol_R, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
//                                           randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
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
//   VectorXd alpha2_vec(alpha2_size);
//   if(W_list.size() > 0){
//     int index = 0;
//     for(int i = 0; i < W_list.size(); i++){
//       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
//       index += alpha2[i].size();
//     }
//   }
//   
//   return(Rcpp::List::create(
//       Named("alpha1") = alpha1,
//       Named("alpha2") = alpha2_vec,
//       Named("beta") = beta,
//       Named("Y_prec") = Y_prec
//   ));
// }
// 
