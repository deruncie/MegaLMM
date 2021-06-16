// #include <math.h>
// #include <iostream>
// #include "MegaLMM_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
//
// // [[Rcpp::export()]]
// VectorXd sample_delta_c_Eigen(
//     VectorXd delta,
//     VectorXd tauh,
//     Map<VectorXd> scores,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randg_draws  // all done with rate = 1;
// ) {
//   int times = randg_draws.rows();
//   int k = tauh.size();
//
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
//     delta(0) = randg_draws(i,0) / rate;
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
//       delta(h) = randg_draws(i,h) / rate;
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//   }
//   return(delta);
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_trunc_delta_c_Eigen(
//     VectorXd delta,
//     VectorXd tauh,
//     Map<VectorXd> scores,
//     Map<VectorXd> shapes,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randu_draws,
//     double trunc_point
// ) {
//   int times = randu_draws.rows();
//   int k = tauh.size();
//   double p,u;
//
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
//     u = randu_draws(i,0); // don't truncate delta(0)
//     delta(0) = R::qgamma(u,shapes(0),1.0/rate,1,0);
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
//       p = R::pgamma(trunc_point,shapes(h),1.0/rate,1,0);  // left-tuncate delta(h) at trunc_point
//       if(p > 0.999) p = 0.999;  // prevent over-flow.
//       u = p + (1.0-p)*randu_draws(i,h);
//       delta(h) = R::qgamma(u,shapes(h),1.0/rate,1,0);
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//   }
//   return(delta);
// }
//
// // [[Rcpp::export()]]
// Rcpp::List sample_delta_omega_c_Eigen(
//     VectorXd delta,
//     VectorXd tauh,
//     double omega2,
//     double xi,
//     Map<VectorXd> scores,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randg_draws  // all done with rate = 1;
// ) {
//   int times = randg_draws.rows();
//   int k = tauh.size();
//
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//
//     rate = 1.0/xi + tauh.dot(scores);
//     omega2 = rate / randg_draws(i,0);
//
//     rate = 1.0 + 1.0 / omega2;
//     xi = rate / randg_draws(i,1);
//
//     VectorXd std_scores = scores / omega2;
//
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(std_scores);
//     delta(0) = randg_draws(i,2) / rate;
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(std_scores.tail(k-h));
//       delta(h) = randg_draws(i,2+h) / rate;
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//   }
//   return(Rcpp::List::create(omega2,xi,delta));
// }
//
//
// // [[Rcpp::export()]]
// Rcpp::List sample_trunc_delta_omega_c_Eigen(
//     VectorXd delta,
//     VectorXd tauh,
//     double omega2,
//     double xi,
//     Map<VectorXd> scores,
//     Map<VectorXd> shapes,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randu_draws,
//     double trunc_point
// ) {
//   int times = randu_draws.rows();
//   int k = tauh.size();
//
//   double rate,delta_old, u,p;
//   for(int i = 0; i < times; i++){
//
//     rate = 1.0/xi + tauh.dot(scores);
//     u = randu_draws(i,0); // don't truncate delta(0)
//     omega2 = 1.0 / R::qgamma(u,shapes(0),1.0/rate,1,0);
//
//     rate = 1.0 + 1.0 / omega2;
//     u = randu_draws(i,1); // don't truncate delta(0)
//     xi = 1.0 / R::qgamma(u,shapes(1),1.0/rate,1,0);
//
//     VectorXd std_scores = scores / omega2;
//
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(std_scores);
//     u = randu_draws(i,2); // don't truncate delta(0)
//     delta(0) = R::qgamma(u,shapes(2),1.0/rate,1,0);
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(std_scores.tail(k-h));
//       p = R::pgamma(trunc_point,shapes(2+h),1.0/rate,1,0);  // left-tuncate delta(h) at trunc_point
//       if(p > 0.999) p = 0.999;  // prevent over-flow.
//       u = p + (1.0-p)*randu_draws(i,2+h);
//       delta(h) = R::qgamma(u,shapes(2+h),1.0/rate,1,0);
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//   }
//   return(Rcpp::List::create(omega2,xi,delta));
// }
//
//
// // -------------------------- //
// // -------------------------- //
// // -------------------------- //
//
// // [[Rcpp::export()]]
// NumericVector rgig_multiple(
//   int n,
//   NumericVector lambda,
//   NumericVector chi,
//   NumericVector psi
// ){
//   NumericVector res(n);
//   SEXP (*fun)(SEXP,SEXP,SEXP,SEXP) = NULL;
//   if (!fun)
//     fun = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("GIGrvg", "rgig");
//   for(int i = 0; i < n; i++){
//     res[i] = as<double>(fun(wrap(1.0),wrap(lambda[i]),wrap(chi[i]),wrap(psi[i])));
//   }
//   return(res);
// }
//
//
// // -------------------------- //
// // ----Y = ZXB + E ---------- //
// // -------------------------- //
//
// VectorXd sample_MME_single_hierarchical_diagK(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MSpMat   Z,           // nxr   // use when r < n < b
//     MatrixXd X,           // rxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
//   // Phi = sqrt(tot_Eta_prec) * chol_R_invT * Z * X
//   // Phi = U * X = Vt * X, where U = Vt = sqrt(tot_Eta_prec) * chol_R_invT * Z
//   MatrixXd U = chol_R.transpose().triangularView<Lower>().solve(Z.toDense() * sqrt(tot_Eta_prec));
//   MatrixXd Phi = U * X;
//   MatrixXd V = U.transpose();
//   VectorXd alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//
//   VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
//   u += prior_mean;
//   VectorXd v = Phi * u + randn_e;
//   VectorXd alpha_v = alpha-v;
//
//   // Binomial inverse theorm: (A + UBV)^-1 = A^-1 - A^-1 * U (I + BVU)^-1*BVA^-1
//   MatrixXd B = X * prior_prec.cwiseInverse().asDiagonal() * X.transpose();
//   MatrixXd UB = U*B;
//   MatrixXd inner = B + UB.transpose()*UB;
//   VectorXd w = alpha_v - UB*inner.ldlt().solve(UB.transpose() * alpha_v);
//
//   VectorXd theta = u + prior_prec.cwiseInverse().asDiagonal() * (Phi.transpose() * w);
//   return(theta);
// }
//
// struct sample_MME_single_hierarchical_diagK_worker : public RcppParallel::Worker {
//   MatrixXd Y;
//   MSpMat Z;
//   MatrixXd X;
//   MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//   const std::vector<MSpMat> chol_R_list;
//   VectorXi h2s_index;
//   VectorXd tot_Eta_prec;
//   MatrixXd &coefs;
//
//   sample_MME_single_hierarchical_diagK_worker(
//     MatrixXd Y,           // nxp
//     MSpMat   Z,           // nxb
//     MatrixXd X,           // nxb
//     MatrixXd prior_mean,  // bxp
//     MatrixXd prior_prec,  // bxp
//     const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
//     VectorXi h2s_index,   // px1, 1-based index
//     VectorXd tot_Eta_prec,// px1
//     MatrixXd randn_theta, // bxp
//     MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
//     MatrixXd &coefs       // bxp
//   ):
//     Y(Y), Z(Z), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
//     chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       int h2_index = h2s_index[j] - 1;
//       MSpMat chol_R = chol_R_list[h2_index];
//       coefs.col(j) = sample_MME_single_hierarchical_diagK(Y.col(j), Z, X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
//     }
//   }
// };
//
// // Samples from B in model:
// // Y = ZXB + E
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_hierarchical_c(  // returns bxp matrix
//     Map<MatrixXd> Y,              // nxp
//     MSpMat        Z,              // nxr   // use when r < n < b
//     Map<MatrixXd> X,              // rxb
//     Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
//     VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
//     Map<VectorXd> tot_Eta_prec,   // px1
//     Map<MatrixXd> prior_mean,     // bxp
//     Map<MatrixXd> prior_prec,     // bxp
//     int grainSize) {
//
//   int b = X.cols();
//   int p = Y.cols();
//   int n = Y.rows();
//
//   MatrixXd randn_theta = rstdnorm_mat(b,p);
//   MatrixXd randn_e = rstdnorm_mat(n,p);
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < h2s_index.maxCoeff(); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
//
//   MatrixXd coefs(b,p);
//
//   sample_MME_single_hierarchical_diagK_worker sampler(Y,Z,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }
//
//
