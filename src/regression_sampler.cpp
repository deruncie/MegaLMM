#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List parallel_block_regression_sampler(
    Map<MatrixXd> Y,               //
    Map<MatrixXd> X1,
    Map<MatrixXd> X2,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
    SEXP V_,
    SEXP chol_V_,
    VectorXd Y_prec,               // p-vector of Y current precisions
    VectorXd Y_prec_a0,
    VectorXd Y_prec_b0,
    Map<MatrixXd> prior_prec_alpha, // a1 x p matrix of prior precisions for alpha1
    Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
    Map<MatrixXd> prior_prec_beta // b x p matrix of prior precisions of beta
) {

  
  int n = Y.rows();
  int p = Y.cols();
  
  // RMatrix X1 = load_RMatrix(X1_, false);
  RMatrix V = load_RMatrix(V_, false);
  RMatrix chol_V = load_RMatrix(chol_V_, true);
  
  int a = X1.cols();
  
  int b;
  if(V.isNULL) {
    b = X2.cols();
  } else{
    b = V.cols();
    if(b == 0) stop("Need V or X2 to have non-zero columns");  // is this true?
  }
  
  // generate random numbers
  MatrixXd randn_alpha = rstdnorm_mat(a,p);
  MatrixXd randn_beta = rstdnorm_mat(b,p);
  MatrixXd randn_e;
  if(b > n) {
    randn_e = rstdnorm_mat(n,p);
  }
  VectorXd rgamma_1(p);
  for(int i = 0; i < p; i++) {
    rgamma_1(i) = rgamma(1,Y_prec_a0[i] + n/2.0,1.0)(0);
  }
  
  // Results structures
  MatrixXd alpha(a,p);
  MatrixXd beta(b,p);
  
  // Remove prior_mean_beta from Y
  Y = Y - X2 * prior_mean_beta;
  
  // rotate by chol_V
  MatrixXd RinvtY = chol_V.tsolve(Y);
  MatrixXd RinvtX2 = chol_V.tsolve(X2);
  MatrixXd C = RinvtX2.transpose() * RinvtX2;
  MatrixXd XtRinvY = RinvtX2.transpose() * RinvtY;
  
  #pragma omp parallel for num_threads(get_MegaLMM_nthreads())  
  for(int j = 0; j < p; j++) {
    
    MatrixXd C_beta = C;
    C_beta.diagonal() += prior_prec_beta.col(j);
    LLT<MatrixXd> A_beta_llt;
    A_beta_llt.compute(C_beta);
    MatrixXd chol_A_beta = A_beta_llt.matrixU();
    
    // Step 1
    VectorXd Rinvty_tilde = RinvtY.col(j);
    if(a > 0) {
      // Sample alpha
      // Calculate A_alpha = Y_prec*X1^T*V_beta^{-1}*X1 + D_alpha^{-1}
      // We don't need to actually calculate V_beta^{-1} directly.
      // Instead,
      MatrixXd RinvtX1 = chol_V.tsolve(X1);  // n*n*a -> n x a
      MatrixXd X1tVinvX2 = RinvtX1.transpose() * RinvtX2; // a*n*b -> a*b
      MatrixXd cholAbetainvt_X2tVinvX1 = chol_A_beta.transpose().triangularView<Lower>().solve(X1tVinvX2.transpose()); // b*b*a -> b x a
      
      MatrixXd A_alpha = RinvtX1.transpose() * RinvtX1 - cholAbetainvt_X2tVinvX1.transpose() * cholAbetainvt_X2tVinvX1;
      A_alpha.diagonal() += prior_prec_alpha.col(j);
      
      // VectorXd Rinvsqy = get_RinvtX2(chol_V,y); // n*n*q -> n x 1;
      // VectorXd XtRinvy = RinvtX2.transpose() * Rinvsqy; // b*n*1 >- b x 1
      VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvY.col(j)); // b*b*1 -> b*1
      
      VectorXd X1t_V_beta_inv_y = RinvtX1.transpose() * RinvtY.col(j) - cholAbetainvt_X2tVinvX1.transpose() * invSqAbXtRinvy;
      
      LLT<MatrixXd> A_alpha_llt;
      A_alpha_llt.compute(A_alpha);
      MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
      
      alpha.col(j) = chol_A_alpha.transpose().triangularView<Lower>().solve(X1t_V_beta_inv_y) + 1.0/sqrt(Y_prec[j]) * randn_alpha.col(j);
      alpha.col(j) = chol_A_alpha.triangularView<Upper>().solve(alpha.col(j));
      
      Rinvty_tilde = Rinvty_tilde - RinvtX1 * alpha.col(j);
    }
    
    // Step 2 - sample Y_prec
    // We don't need to actually calculate V_beta^{-1} directly.
    // VectorXd RinvSqy = get_RinvtX2(chol_V,y_tilde);
    VectorXd XtRinvy = RinvtX2.transpose() * Rinvty_tilde;
    VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
    double score = Y_prec_b0[j] + (Rinvty_tilde.dot(Rinvty_tilde) - prod1.dot(prod1))/2;
    Y_prec[j] = rgamma_1[j]/score;
    
    // Step 3 - sample beta
    VectorXd XtRinvy_std_mu = XtRinvy*Y_prec[j]; // + prior_prec_beta.col(j).asDiagonal()*prior_mean_beta.col(j);
    beta.col(j) = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec[j]) + randn_beta.col(j);
    beta.col(j) = chol_A_beta.triangularView<Upper>().solve(beta.col(j)) / sqrt(Y_prec[j]);
    
    // Step 4 - add in prior_mean_beta
    beta.col(j) += prior_mean_beta.col(j);
  }
  
  return(Rcpp::List::create(
      Named("alpha") = alpha,
      Named("beta") = beta,
      Named("Y_prec") = Y_prec
  ));
}
   
  
  
  // [[Rcpp::export]]
Rcpp::List parallel_Single_regression_sampler(
      Map<MatrixXd> Y_,               //
      Map<MatrixXd> X1_,
      Map<MatrixXd> X2_,               // either X2, a n x b matrix, or Ux, a n x m matrix. If Ux, then V must be non-NULL
      SEXP V_,
      SEXP chol_V_,
      VectorXd Y_prec,               // p-vector of Y current precisions
      VectorXd Y_prec_a0,
      VectorXd Y_prec_b0,
      Map<MatrixXd> prior_prec_alpha, // a1 x p matrix of prior precisions for alpha1
      Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
      Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
      Map<MatrixXd> current_alphas,
      Map<MatrixXd> betas_alpha,
      Map<MatrixXd> betas_beta,
      Map<MatrixXd> betas_pi,
      MatrixXi betas_delta,
      int run_sampler_times
) {
    int n = Y_.rows();
    int p = Y_.cols();
    
    // RMatrix X1 = load_RMatrix(X1_, false);
    RMatrix V = load_RMatrix(V_, false);
    RMatrix chol_V = load_RMatrix(chol_V_, true);
    
    int a1 = X1_.cols();
    
    int b;
    if(V.isNULL) {
      b = X2_.cols();
    } else{
      b = V.cols();
      if(b == 0) stop("Need V or X2 to have non-zero columns");  // is this true?
    }
    
    // generate random numbers
    MatrixXd randn_alphas = rstdnorm_mat(a1*run_sampler_times,p);
    MatrixXd randn_betas = rstdnorm_mat(b*run_sampler_times,p);
    MatrixXd rgamma_1s(run_sampler_times,p);
    for(int i = 0; i < p; i++) {
      rgamma_1s.col(i) = as<VectorXd>(rgamma(run_sampler_times,Y_prec_a0[i] + n/2.0 + b/2.0,1.0));
    }
    MatrixXd rand_unifs;
    rand_unifs = MatrixXd::Zero(b*run_sampler_times,p);
    for(int i = 0; i < p; i++) {
      rand_unifs.col(i) = as<VectorXd>(runif(b*run_sampler_times,0.0,1.0));
    }
    
    // rotate by chol_V
    MatrixXd Y = chol_V.tsolve(Y_);
    MatrixXd X1 = chol_V.tsolve(X1_);
    MatrixXd X2 = chol_V.tsolve(X2_);

    ArrayXd diag_X2tX2 = X2.cwiseProduct(X2).colwise().sum();
    ArrayXd diag_X1tX1 = X1.cwiseProduct(X1).colwise().sum();
    
    #pragma omp parallel for num_threads(get_MegaLMM_nthreads())  
    for(int t = 0; t < p; t++) {
      VectorXd y = Y.col(t);
      
      VectorXd a = current_alphas.col(t);
      VectorXd alpha = betas_alpha.col(t);
      VectorXd beta = betas_beta.col(t);
      VectorXi delta = betas_delta.col(t);
      double invVarRes = Y_prec[t];
      VectorXd invVarEffects = prior_prec_beta.col(t);
      VectorXd pi = betas_pi.col(t);
      VectorXd randn_a = randn_alphas.col(t);
      VectorXd randn_beta = randn_betas.col(t);
      VectorXd rand_unif = rand_unifs.col(t);
      VectorXd rgamma_1 = rgamma_1s.col(t);
      // double Y_prec_b0 = Y_prec_b0[j];
      
      
      
      int nMarkers = b;
      
      ArrayXd logPi = pi.array().log();
      ArrayXd logPiComp = (1.0 - pi.array()).log();
      ArrayXd logDelta0 = logPi;
      
      VectorXd yCorr = y - X1*a;
      for(int j = 0; j < nMarkers; j++) {
        if(delta[j] != 0) yCorr -= X2.col(j)*alpha[j];
      }
      
      for(int i = 0; i < run_sampler_times; i++) {
        VectorXd invVarEffects = prior_prec_beta.col(t) * invVarRes;
        ArrayXd logVarEffects = invVarEffects.array().inverse().log();
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
        
        int nLoci = 0;
        // Sample beta = alpha*delta
        for(int j = 0; j < nMarkers; j++) {
          float rhs = (X2.col(j).dot(yCorr) + diag_X2tX2[j]*alpha[j] + prior_mean_beta(j,t)*invVarEffects[j])*invVarRes;
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
        invVarRes = rgamma_1[i]/(yCorr.dot(yCorr)/2.0 + beta.dot(beta)/2.0 + Y_prec_b0[t]);
      }
      current_alphas.col(t) = a;
      betas_alpha.col(t) = alpha;
      betas_beta.col(t) = beta;
      betas_delta.col(t) = delta;
      Y_prec[t] = invVarRes;
      
    }
    
    return(Rcpp::List::create(
        Named("alpha") = current_alphas,
        Named("betas_alpha") = betas_alpha,
        Named("betas_beta") = betas_beta,
        Named("betas_delta") = betas_delta,
        Named("Y_prec") = Y_prec
    ));
  }


