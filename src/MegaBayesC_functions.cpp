#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

using namespace Eigen;


// [[Rcpp::export()]]
MatrixXf rstdnorm_mat_f(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat.cast<float>());
}


General_Matrix_f load_General_Matrix_f(SEXP X_, bool triangular) {
  MatrixXf null_d = MatrixXf::Zero(0,0);
  if(Rf_isNull(X_)) {
    SpMat_f null_s = null_d.sparseView();
    General_Matrix_f Xm(null_d,null_s,triangular,false,true);
    return(Xm);
  } else if(Rf_isMatrix(X_)){
    MatrixXf X = as<MatrixXf >(X_);
    SpMat_f null_s = null_d.sparseView();
    General_Matrix_f Xm(X,null_s,triangular,true,false);
    return(Xm);
  } else{
    SpMat_f X = as<SpMat_f>(X_);
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
      Eigen::SimplicialLDLT<SpMat_f> ldlt(sparse);
      return(ldlt.solve(Y));
    }
  }
}

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
      Eigen::SimplicialLDLT<SpMat_f> ldlt(sparse.transpose());
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
Rcpp::List SingleSite_regression_sampler_parallel(
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
  
  int which_sampler = 4; // v4 sampler is the single site sampler
  
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
  MatrixXf randn_alpha1 = rstdnorm_mat_f(a1*run_sampler_times,p);
  std::vector<VectorXf> randn_alpha2;
  if(X1_list.size() > 0){
    for(int i = 0; i < p; i++){
      randn_alpha2.push_back(rstdnorm_mat_f(X1_list[i].dense.cols()*run_sampler_times,1));
    }
  }
  MatrixXf randn_beta = rstdnorm_mat_f(b*run_sampler_times,p);
  MatrixXf randn_e;
  if(which_sampler == 2 || which_sampler == 3) {
    randn_e = rstdnorm_mat_f(n*run_sampler_times,p);
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
          stop("Use 'regression_sampler_parallel' for which_sampler 1-3");
          // b = RinvtX2.cols();
          // samples = regression_sampler_v1(Y.col(j), X1, RinvtX2, C, prior_prec_alpha, prior_mean_beta.col(j),
          //                                 prior_prec_beta.col(j), chol_V, Y_prec[j], randn_alpha,
          //                                 randn_beta.col(j), rgamma_1(0,j),Y_prec_b0[j]);
        } else if(which_sampler == 2) {
          stop("Use 'regression_sampler_parallel' for which_sampler 1-3");
          // b = X2.cols();
          // samples = regression_sampler_v2(Y.col(j), X1, X2, prior_prec_alpha, prior_mean_beta.col(j),
          //                                 prior_prec_beta.col(j), chol_V, V, Y_prec[j], randn_alpha,
          //                                 randn_beta.col(j), randn_e.col(j),rgamma_1(0,j),Y_prec_b0[j]);
        } else if(which_sampler == 3) {
          stop("Use 'regression_sampler_parallel' for which_sampler 1-3");
          // b = Vx.cols();
          // samples = regression_sampler_v3(Y.col(j), X1, Ux, Vx, prior_prec_alpha, prior_mean_beta.col(j),
          //                                 prior_prec_beta.col(j), chol_V, Vinv, VinvUx, UtVinvU, Y_prec[j], randn_alpha,
          //                                 randn_beta.col(j), randn_e.col(j),rgamma_1(0,j),Y_prec_b0[j]);
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

