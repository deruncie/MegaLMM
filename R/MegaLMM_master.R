# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


#' Set MegaLMM run parameters
#'
#' Function to create run_parameters list for initializing MegaLMM model
#'
#' @param scale_Y Should the Y values be centered and scaled? Recommend, except for simulated data.
#' @param K number of factors
#' @param h2_divisions A scalar or vector of length equal to number of random effects. In MegaLMM, random
#'   effects are parameterized as proportions of the total variance of all random effects plus residuals. 
#'   The prior on the variance componets is discrete spanning the interval [0,1) over each varince component proportion
#'   with \code{h2_divisions} equally spaced values is constructed. If
#'   \code{h2_divisions} is a scalar, the prior for each variance component has this number of divisions. 
#'   If a vector, the length should equal the number of variance components, in the order of the random effects specified in the model
#' @param h2_step_size Either NULL, or a scaler in the range (0,1].
#'   If NULL, h2's will be sampled based on the marginal probability over all possible h2 vectors. 
#'   If a scalar, a Metropolis-Hastings update step will be used for each h2 vector.
#'   The trail value will be selected uniformly from all possible h2 vectors within this Euclidean distance from the current vector.
#' @param drop0_tol A scalar giving the a tolerance for the \code{drop0()} function that will be applied
#'     to various symmetric (possibly) sparse matrices to try to fix numerical errors and increase sparsity.
#' @param K_eigen_tol A scalar giving the minimum eigenvalue of a K matrix allowed. During pre-processing,
#'     eigenvalues of each K matrix will be calculated using \code{svd(K)}. Only eigenvectors of K with corresponding eigenvalues
#'     greater than this value will be kept. If smaller eigenvalues exist, the model will be transformed
#'     to reduce the rank of K, by multiplying Z by the remaining eigenvectors of K. This transformation
#'     is undone before posterior samples are recorded, so posterior samples of \code{U_F} and \code{U_R} are
#'     untransformed.
#' @param burn burnin length of the MCMC chain
#' @param thin thinning rate of the MCMC chain
#' @param max_NA_groups If 0, all NAs will be imputed during sampling. If Inf, all NAs will be marginalized over.
#'     If in (0,Inf), up to this many groups of columns will be separately sampled.
#'     The minimum number of NAs in each column not in one of these groups will be imputed.
#' @param svd_K If TRUE, the the diagonalization of ZKZt for the first random effect is accomplished using this algorithm:
#'     https://math.stackexchange.com/questions/67231/singular-value-decomposition-of-product-of-matrices which doesn't require forming ZKTt.
#'     If FALSE, the SVD of ZKZt for the first random effect is calculated directly. TRUE is generally faster if the same genomes are repeated several times.
#' @param verbose should progress during initiation and sampling be printed?
#' @param save_current_state should the current state of the sampler be saved every time the function \code{sample_MegaLMM} is called?
#' @seealso \code{\link{MegaLMM_init}}, \code{\link{sample_MegaLMM}}, \code{\link{print.MegaLMM_state}}
#'
MegaLMM_control = function(
                        which_sampler = list(Y = 1,F = 1),
                        run_sampler_times = 1,
                        scale_Y = c(T,F),
                        K = 20, h2_divisions = 100, h2_step_size = NULL,
                        drop0_tol = 1e-14, K_eigen_tol = 1e-10,
                        burn = 100,thin = 2,
                        max_NA_groups = Inf,
                        svd_K = TRUE,
                        verbose = TRUE,
                        save_current_state = TRUE,
                        ...
                        ) {
  formals_named = formals()
  formals_named = formals_named[names(formals_named) != '...']
  all_args = lapply(formals_named,function(x) {
    if(!is.list(eval(x))) {
      eval(x)[1]
    } else {
      eval(x)
    }})
  passed_args = lapply(as.list(match.call())[-1],eval)
  if(any(names(passed_args) %in% names(formals_named) == F)){
    unused_names = names(passed_args)[names(passed_args) %in% names(formals_named) == F]
    warning(sprintf('No argument(s) named %s',paste(unused_names,sep=', ')))
  }
  all_args[names(passed_args)] = passed_args
  return(all_args)
}



#' Set MegaLMM priors
#'
#' Function to create list of priors for MegaLMM model.
#'
#' Default values are provided, but any can be replaced. Note: \code{h2_priors_resids} and
#'     \code{h2_priors_factors} can be set after calling \code{MegaLMM_init} and before \code{sample_MegaLMM}
#'     if that is easier. See Vignette.
#'
#' @param tot_Y_var List of parameters of inverse gamma distribution for residual variances, specifically:
#'     \code{V} and \code{nu}, give shape = \code{nu-1} and scale = \code{1/(nu*V)}, so mean = \code{~V}
#' @param tot_F_var List of parameters of inverse gamma distribution for factor variances. See \code{tot_Y_var}.
#'     This parameter provides the parameter extension of Ghosh and Dunson (2009), but is removed
#'     from all Factor parameters before they are saved in Posterior
#' @param h2_priors_resids_fun function that returns prior probability for a given vector of variance component proportions for each random effect
#'     The function must take two arguments - a vector of variance component proportions: \code{h2},
#'     and \code{n} - the number of discrete levels of the prior.
#'     Alternatively, can be a scalar or vector of (relative) prior values for each value of the
#'     discrete prior.
#' @param h2_priors_factors_fun see \code{h2_priors_resids_fun}. Same, but for the h2s of the factors.
#' @param Lambda_prior A list with elements:
#'     1) \code{sampler}: a function that draws samples of the precision matrix for Lambda. Ex: \code{sample_Lambda_prec_reg_horseshoe}; 2)
#'     any other hyperparameters and control parameters for \code{sampler}
#' @param B_prior A list with elements:
#'     1) \code{sampler}: a function that draws samples of the precision matrix for B and B_F Ex: \code{sample_B2_prec_reg_horseshoe}; 2)
#'     any other hyperparameters and control parameters for \code{sampler}
#' @param cis_effects_prior Currently accepts a list with a single value giving the precision of each cis_effect variable
#'
#' @return a list with each of the prior components specified above.
#' @export
#'
MegaLMM_priors = function(
                        tot_Y_var = list(V = 0.5,   nu = 3),
                        tot_F_var = list(V = 18/20, nu = 20),
                        h2_priors_resids_fun = function(h2s, n) 1,
                        h2_priors_factors_fun = function(h2s, n) 1,
                        Lambda_prior = list(
                          sampler = sample_Lambda_prec_horseshoe,
                          prop_0 = 0.1,
                          delta_l = list(shape = 3, rate = 1),
                          delta_iterations_factor = 100
                        ),
                        B2_prior = list(
                          sampler = sample_B2_prec_horseshoe,
                          prop_0 = 0.1
                        ),
                        cis_effects_prior = list(
                          prec = 1
                        )

                    ) {
  passed_args = lapply(as.list(match.call())[-1],function(x) eval(x))
  default_args = formals()
  default_args = lapply(default_args[names(default_args) %in% names(passed_args) == F],function(x) eval(x))
  all_args = c(passed_args,default_args)
  return(all_args)
}


#' Set up a MegaLMM model
#'
#' Sets up the MegaLMM model, selects starting values, and pre-calculates matrices for the Gibbs
#' sampler.
#'
#' The first step in fitting a MegaLMM model. This function sets up the model matrices based on the
#' fixed and random effect formulas provided. This function must be followed by calls to
#' \link{set_priors_MegaLMM}, \link{initialize_variables_MegaLMM} and \link{initialize_MegaLMM}, before
#' the Gibbs sampler can be run with \link{sample_MegaLMM}.
#'
#' The model is specified as:
#'
#' y_i = g(eta_i)
#'
#' Eta = rbind(eta_1,...,eta_n) = X1*B1 + X2_R*B2_R + F*Lambda + Z*U_R + E_R
#'
#' F = X2_F * B2_F + Z*U_F + E_F
#'
#' For sampling, we reparameterize as:
#'
#' Qt*Eta = Qt*X*B + Qt*F*Lambda + Qt*ZL*U_R + Qt*E_R
#'
#' Qt*F = Qt*X_F * B_F + Qt*ZL*U_F + Qt*E_F
#'
#' where LTL = K and ZL = Z*L
#'
#' We sample the quantities Qt*Eta, Qt*F, B, Lambda, U_R, U_F. We then back-calculate Eta and F.
#' 
#' \strong{Note:} In the original \emph{MegaLMM} paper, we set \code{y_i = eta_i}, so replaced \code{Eta} with \code{Y} above.
#'
#' @param Y either a) a n x p matrix of data (n individuals x p traits), or b) a list describing
#'     the observation_model, data, and associated parameters. This list should contain:
#'     i) \code{observation_model}: a function modeled after \code{missing_observation_model}
#'         (see code by typing this into the console) that draws posterior samples of Eta conditional on
#'         the observations (Y) and the current state of the MegaLMM model (current_state).
#'         The function should have the same form as \link{missing_data},
#'         and must return Eta even if current_state is NULL.
#'     ii) \code{observations}: a data.frame containing the observaition-level data and associated covariates.
#'         This must include a column \code{ID} that is also present in \code{data}
#'     iii) any other parameters necessary for \code{observation_model}
#' @param formula RHS of a model. The syntax is similar to \link{lmer}.
#'     Random effects are specified by (1+factor | group), with the left side of the '|' a design
#'     matrix, and the right side a random effect factor (group). For each random effect factor, a
#'     covariance matrix (\code{K_mats}) or precision matrix (\code{K_inv_mats}) can be provided.
#'     Unlike in \code{lmer}, each variable or covariate in the design matrix is given an
#'     independent random effect (ie no covariance among random effects is modeled), so two bars '||'
#'     gives an identical model to one bar.
#'     Note: the speed of the model will decrease dramatically with the number of random effects (multiplied
#'     by h2_divisions).
#'     Fixed effects only apply to the model residuals (ie X1) below, and are not regularized (prior precision = 0)
#' @param extra_regressions Optional. A list including either:
#'     i) the matrix X (n x b) of regression coeffients, or
#'     ii) two matrices U (n x m) and V (m x b) such that X = U*V
#'     also, logical variables \code{resids} and \code{fixed} specify whether these coefficients apply to either or
#'     both of the model residuals or the factors
#' @param data data.frame with n rows containing columns corresponding to the fixed and random
#'   effects
#' @param relmat Optional. A list of covariance matrices for random effects. If none provided
#'     for any of the random effects, K is assumed to be the identity.
#' @param cis_genotypes Optional. A list of n x ci matrices of length p giving cis-effect coefficients for each trait
#' @param Lambda_fixed Optional. A matrix of the first k rows of Lambda that are fixed
#' @param run_parameters See \link{MegaLMM_control}
#' @param posteriorSample_params A character vector giving names of parameters to save all posterior samples
#' @param posteriorMean_params A character vector giving names of parameters to save only the posterior mean.
#' @param run_ID A unique identifier for this model. The code will create a folder with this name to hold all
#'     posterior samples and diagnostic information during the run.
#'
#' @return An object of class MegaLMM_state with components: \itemize{
#'     \item current_state: a list of parameters in the current iteration of the sampler. Initially empty
#'     \item Posterior: a list of arrays of posterior samples. Initially empty
#'     \item RNG: current state of R's Random number generator (for re-starting chaings)
#'     \item traitnames: vector of trait names (from colnames of Y)
#'     \item run_parameters, run_variables, data_matrices, priors: input data and parameters
#' }
#' @seealso \code{\link{MegaLMM_control}}, \code{\link{sample_MegaLMM}}, \code{\link{print.MegaLMM_state}}, \code{\link{plot.MegaLMM_state}}#'
#' @export
#'
setup_model_MegaLMM = function(Y,formula,extra_regressions=NULL,data,relmat=NULL, cis_genotypes = NULL, Lambda_fixed = NULL,
                            run_parameters = MegaLMM_control(),
                            posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec',
                                                       'resid_h2', 'B1', 'B2_F','B2_R','U_R','cis_effects','Lambda_m_eff',
                                                       'Lambda_pi','B2_R_pi','B2_F_pi'),
                            posteriorMean_params = c(),
                            run_ID = 'MegaLMM_run'){
  # creates model matrices, RE_setup, current_state
  # returns MegaLMM_state

  try(dir.create(run_ID),silent=T)

  # ----------------------------- #
  # -------- observation model ---------- #
  # ----------------------------- #

  if(is(Y,'list')){
    if(!'observation_model' %in% names(Y)) stop('observation_model not specified in Y')
    observation_model = Y$observation_model
    observation_model_parameters = Y[names(Y) != 'observation_model']
  } else{
    if(!is(Y,'matrix'))	Y = as.matrix(Y)
    if(nrow(Y) != nrow(data)) stop('Y and data have different numbers of rows')
    observation_model = missing_data_model
    observation_model_parameters = list(
      Y = Y,
      scale_Y = run_parameters$scale_Y
    )
  }

  # initialize observation_model
  observation_model_parameters$observation_setup = observation_model(observation_model_parameters,list(data_matrices = list(data = data)))
  n = nrow(data)
  p = observation_model_parameters$observation_setup$p

  traitnames = observation_model_parameters$observation_setup$traitnames
  if(is.null(traitnames)) traitnames = paste('trait',1:p,sep='_')
  if(is.null(observation_model_parameters$observation_setup$Y_missing)) {
    observation_model_parameters$observation_setup$Y_missing = matrix(F,n,p)
  }
  if(!is(observation_model_parameters$observation_setup$Y_missing,'lgTMatrix')){
    observation_model_parameters$observation_setup$Y_missing = as(observation_model_parameters$observation_setup$Y_missing,'lgTMatrix')
  }

  # ----------------------------- #
  # -------- model matrices ---------- #
  # ----------------------------- #

  model_setup = make_model_setup(formula,data,relmat)

  # X1 is the "fixed effects", un-shrunk covariates that only affect residuals
  X1 = model_setup$lmod$X
  b1 = ncol(X1)

  # -------- regressions ---------- #
  # X2_F and X2_R (and V_F and V_R) are the coefficient matrices for the regularized regression coefficients
  # these are specified as lists with either X (n x b) or {U(nxm),V(mxb)}
  X2_R = matrix(0,n,0)
  X2_F = matrix(0,n,0)
  U2_F = X2_F
  V2_F = NULL
  b2_R = 0
  b2_F = 0
  if(!is.null(extra_regressions)) {
    # project out intercept
    if('X' %in% names(extra_regressions)){
      M = diag(1,n) - matrix(1/n,n,n)
      extra_regressions$X = M %*% extra_regressions$X
      extra_regressions$X = extra_regressions$X[,colSums(abs(extra_regressions$X))>1e-10,drop=FALSE]
    } else if('V' %in% names(extra_regressions)){
      m = nrow(extra_regressions$V)
      M = diag(1,m) - matrix(1/m,m,m)
      extra_regressions$V = M %*% extra_regressions$V
      extra_regressions$V = extra_regressions$V[,colSums(abs(extra_regressions$V))>1e-10,drop=FALSE]
    }
    if(!is.null(extra_regressions$factors) && extra_regressions$factors == TRUE){
      if('X' %in% names(extra_regressions)){
        X2_F = extra_regressions$X
        U2_F = X2_F
        V2_F = NULL
        b2_F = ncol(X2_F)
      } else{
        if(all(c('U','V') %in% names(extra_regressions))){
          U2_F = extra_regressions$U
          V2_F = extra_regressions$V
          X2_F = U2_F %*% V2_F
          b2_F = ncol(V2_F)
        } else{
          stop("missing U or V in extra_regressions")
        }
      }
    }
    if(!is.null(extra_regressions$resids) && extra_regressions$resids == TRUE){
      if('X' %in% names(extra_regressions)){
        X2_R = extra_regressions$X
        b2_R = ncol(X2_R)
      } else{
        if(all(c('U','V') %in% names(extra_regressions))){
          X2_R = extra_regressions$U %*% extra_regressions$V
          b2_R = ncol(X2_R)
        } else{
          stop("missing U or V in extra_regressions")
        }
      }
    }
  }
  if(!nrow(X2_F) == n) stop("Wrong dimension of X2_F")
  if(!nrow(X2_R) == n) stop("Wrong dimension of X2_R")
  
  # -------- cis genotypes ---------- #
  if(is.null(cis_genotypes)){
    cis_genotypes = list()
    n_cis_effects = NULL
    cis_effects_index = NULL
  } else{
    n_cis_effects = sapply(cis_genotypes,ncol)
    cis_effects_index = rep(seq_len(p),n_cis_effects)
  }

  # -------- Random effects ---------- #

  # RE_setup is a list like from BGLR that describes the random effects (for both residuals and factors)
  RE_setup = model_setup$RE_setup

  # add names to RE_setup if needed
  n_RE = length(RE_setup)
  for(i in 1:n_RE){
    if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
      names(RE_setup)[i] = paste0('RE.',i)
    }
  }
  RE_names = names(RE_setup)

  # combine Z matrices
  Z = do.call(cbind,lapply(RE_setup,function(x) x$Z))
  Z = as(Z,'dgCMatrix')


  # find RE indices
  RE_lengths = sapply(RE_setup,function(x) ncol(x$Z))
  RE_starts = cumsum(c(0,RE_lengths)[1:n_RE])
  names(RE_starts) = RE_names
  RE_indices = lapply(RE_names,function(re) RE_starts[re] + 1:RE_lengths[re])
  names(RE_indices) = RE_names

  # function to ensure that covariance matrices are sparse and symmetric
  fix_K = function(x) forceSymmetric(drop0(x,tol = run_parameters$drop0_tol))

  # construct RE_L and ZL for each random effect
  for(i in 1:length(RE_setup)){
    re_name = names(RE_setup)[i]
    RE_setup[[i]] = within(RE_setup[[i]],{
      if(!'ZL' %in% ls()){
        if('K' %in% ls() && !is.null(K)){
          id_names = rownames(K)
          if(is(K,'Matrix') & isDiagonal(K)) {
            L = as(diag(1,nrow(K)),'dgCMatrix')
            K_inv = as(diag(1/diag(K)),'dgCMatrix')
          } else {
            ldl_k = LDLt(K)
            large_d = ldl_k$d > run_parameters$K_eigen_tol
            r_eff = sum(large_d)
            # if need to use reduced rank model, then use D of K in place of K and merge L into Z
            # otherwise, use original K, set L = Diagonal(1,r)
            if(r_eff < length(ldl_k$d)) {
              K = as(diag(ldl_k$d[large_d]),'dgCMatrix')
              K_inv = as(diag(1/ldl_k$d[large_d]),'dgCMatrix')
              L = t(ldl_k$P) %*% ldl_k$L[,large_d]
            } else{
              L = as(diag(1,nrow(K)),'dgCMatrix')
              K_inv = as(with(ldl_k,t(P) %*% crossprod(diag(1/sqrt(d)) %*% solve(L)) %*% P),'dgCMatrix')
            }
            rm(list=c('ldl_k','large_d','r_eff'))
          }
          if(is.null(rownames(K))) rownames(K) = 1:nrow(K)
          rownames(K_inv) = rownames(K)
        } else if ('K_inv' %in% ls() && !is.null(K_inv)){
          id_names = rownames(K_inv)
          if(is.null(rownames(K_inv))) rownames(K_inv) = 1:nrow(K_inv)
          K = solve(K_inv)
          rownames(K) = rownames(K_inv)
          L = as(diag(1,nrow(K)),'dgCMatrix')
        } else{
          K = as(diag(1,ncol(Z)),'dgCMatrix')
          rownames(K) = colnames(Z)
          id_names = rownames(K)
          K_inv = K
          L = as(diag(1,nrow(K)),'dgCMatrix')
        }
       # if(is.null(id_names)) id_names = 1:length(id_names)
        rownames(L) = paste(id_names,re_name,sep='::')
        K = fix_K(K)
        ZL = Z %*% L
      }
    })
  }
  ZL = do.call(cbind,lapply(RE_setup,function(re) re$ZL))
  if(nnzero(ZL)/length(ZL) < .25) {
    ZL = as(ZL,'dgCMatrix')
  } else{
    ZL = as.matrix(ZL)
  }


  if(length(RE_setup) > 1) {
    RE_L = do.call(bdiag,lapply(RE_setup,function(re) re$L))
    rownames(RE_L) = do.call(c,lapply(RE_setup,function(re) rownames(re$L)))
  } else{
    RE_L = RE_setup[[1]]$L
  }
  if(nnzero(RE_L)/length(RE_L) < 0.25) {
    RE_L = as(RE_L,'dgCMatrix')
  } else{
    RE_L = as.matrix(RE_L)
  }
  r_RE = sapply(RE_setup,function(re) ncol(re$ZL))

  h2_divisions = run_parameters$h2_divisions
  if(length(h2_divisions) < n_RE){
    if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
    h2_divisions = rep(h2_divisions,n_RE)
  }
  if(is.null(names(h2_divisions))) {
    names(h2_divisions) = RE_names
  }
  h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  colnames(h2s_matrix) = NULL

  # # -------------------------------------------#
  # # identify groups of traits with same pattern of missingness
  # # ideally, want to be able to restrict the number of sets. Should be possible to merge sets of columngs together.

  Missing_data_map = list(list(
    Y_obs = 1:n,
    Y_cols = 1:p
  ))
  Missing_row_data_map = list(list(
    Y_obs = 1:n,
    Y_cols = 1:p
  ))


  # # -------------------------------------------#
  # # Fixed factors
  if(is.null(Lambda_fixed)) {
    fixed_factors = rep(F,run_parameters$K)
  } else{
    Lambda_fixed = as.matrix(Lambda_fixed)
    if(ncol(Lambda_fixed) != p) stop("wrong dimensions of Lambda_fixed")
    if(nrow(Lambda_fixed) >= run_parameters$K) stop("nrow(Lambda_fixed) >= K")
    fixed_factors = rep(F,run_parameters$K)
    fixed_factors[1:nrow(Lambda_fixed)] = T
  }
  Kr = run_parameters$K - sum(fixed_factors)


  run_variables = list(
    p      = p,
    n      = n,
    r_RE   = r_RE,
    RE_names = RE_names,
    b1       = b1,
    b2_R = b2_R,
    b2_F = b2_F,
    n_cis_effects     = n_cis_effects,
    cis_effects_index = cis_effects_index,
    Missing_data_map      = Missing_data_map,
    Missing_row_data_map  = Missing_row_data_map,
    fixed_factors         = fixed_factors,
    Kr = Kr
  )

  data_matrices = list(
    X1           = X1,
    X2_F = X2_F,
    U2_F = U2_F,
    V2_F = V2_F,
    X2_R = X2_R,
    Z           = Z,
    ZL          = ZL,
    RE_setup    = RE_setup,
    RE_L        = RE_L,  # matrix necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
    RE_indices  = RE_indices,
    h2s_matrix  = h2s_matrix,
    cis_genotypes = cis_genotypes,
    Lambda_fixed = Lambda_fixed,
    data       = data
  )

  run_parameters$observation_model = observation_model
  run_parameters$observation_model_parameters = observation_model_parameters
  run_parameters$traitnames = traitnames


  # ----------------------------- #
  # -- create MegaLMM_state object - #
  # ----------------------------- #

  MegaLMM_state = list(
    current_state  = list(),
    run_ID         = run_ID,
    data_matrices  = data_matrices,
    priors         = list(),
    run_parameters = run_parameters,
    run_variables  = run_variables
  )
  class(MegaLMM_state) = append('MegaLMM_state',class(MegaLMM_state))


  MegaLMM_state$Posterior = list(
    posteriorSample_params = posteriorSample_params,
    posteriorMean_params = posteriorMean_params,
    total_samples = 0,
    folder = sprintf('%s/Posterior',run_ID),
    files = c()
  )

  MegaLMM_state
}

#' Set priors for MegaLMM model.
#'
#' See \link{MegaLMM_priors} for more information
#'
#' @param MegaLMM_state MegaLMM_state object as returned by \link{setup_model_MegaLMM}
#' @param priors List as returned by \link{MegaLMM_priors}
#'
#' @return MegaLMM_state object with prior information added.
#' @export
#'
set_priors_MegaLMM = function(MegaLMM_state,priors = MegaLMM_priors()) {
  # returns MegaLMM_state

  p = MegaLMM_state$run_variables$p
  K = MegaLMM_state$run_parameters$K
  h2s_matrix = MegaLMM_state$data_matrices$h2s_matrix

  # ----------------------------- #
  # ----- re-formulate priors --- #
  # ----------------------------- #
  # total precision
  if(length(priors$tot_Y_var$V) == 1) {
    priors$tot_Y_var$V = rep(priors$tot_Y_var$V,p)
    priors$tot_Y_var$nu = rep(priors$tot_Y_var$nu,p)
  }
  if(length(priors$tot_F_var$V) == 1) {
    priors$tot_F_var$V = rep(priors$tot_F_var$V,K)
    priors$tot_F_var$nu = rep(priors$tot_F_var$nu,K)
  }
  priors$tot_Eta_prec_rate   = with(priors$tot_Y_var,V * nu)
  priors$tot_Eta_prec_shape  = with(priors$tot_Y_var,nu - 1)
  priors$tot_F_prec_rate     = with(priors$tot_F_var,V * nu)
  priors$tot_F_prec_shape    = with(priors$tot_F_var,nu - 1)

  # h2_priors_resids
  if(exists('h2_priors_resids',priors)) {
    if(length(priors$h2_priors_resids) == 1) priors$h2_priors_resids = rep(priors$h2_priors_resids,ncol(h2s_matrix))
    if(!length(priors$h2_priors_resids) == ncol(h2s_matrix)) stop('wrong length of priors$h2_priors_resids')
  } else{
    if(!is(priors$h2_priors_resids_fun,'function')) stop('need to provide a priors$h2_priors_resids_fun() to specify discrete h2 prior for resids')
    priors$h2_priors_resids = apply(h2s_matrix,2,priors$h2_priors_resids_fun,n = ncol(h2s_matrix))
  }
  priors$h2_priors_resids = priors$h2_priors_resids/sum(priors$h2_priors_resids)
  # h2_priors_factors
  if(exists('h2_priors_factors',priors)) {
    if(length(priors$h2_priors_factors) == 1) priors$h2_priors_factors = rep(priors$h2_priors_factors,ncol(h2s_matrix))
    if(!length(priors$h2_priors_factors) == ncol(h2s_matrix)) stop('wrong length of priors$h2_priors_factors')
  } else{
    if(!is(priors$h2_priors_factors_fun,'function')) stop('need to provide a priors$h2_priors_factors_fun() to specify discrete h2 prior for factors')
    priors$h2_priors_factors = apply(h2s_matrix,2,priors$h2_priors_factors_fun,n = ncol(h2s_matrix))
  }
  priors$h2_priors_factors = priors$h2_priors_factors/sum(priors$h2_priors_factors)

  MegaLMM_state$priors = priors
  MegaLMM_state

}



#' Initialize MegaLMM variables
#'
#' Initializes all variables in MegaLMM model with random draws.
#' Most variables are drawn from N(0,1) distributions, or uniformly for h2 parameters
#'
#' @param MegaLMM_state MegaLMM_state object as returned by \link{setup_model_MegaLMM}
#' @param ... Currenlty not supported
#'
#' @return MegaLMM_state object with current_state initialized with all variables
#' @export
initialize_variables_MegaLMM = function(MegaLMM_state,...){
  run_parameters = MegaLMM_state$run_parameters
  run_variables = MegaLMM_state$run_variables
  data_matrices = MegaLMM_state$data_matrices
  priors = MegaLMM_state$priors

  MegaLMM_state$current_state = with(c(run_parameters,run_variables,data_matrices,priors),{

    # Factors loadings:

    Lambda_prec = matrix(1,K,p)

    # Lambda - factor loadings
    #   Prior: Normal distribution for each element.
    #       mu = 0
    #       sd = sqrt(1/Lambda_prec)
    Lambda = matrix(rnorm(K*p,0,sqrt(1/Lambda_prec)),nr = K,nc = p)
    colnames(Lambda) = traitnames
    Lambda[fixed_factors,] = Lambda_fixed

    # residuals
    # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
    #  Prior: Gamma distribution for each element
    #       shape = tot_Eta_prec_shape
    #       rate = tot_Eta_prec_rate
    tot_Eta_prec = matrix(rgamma(p,shape = tot_Eta_prec_shape,rate = tot_Eta_prec_rate),nrow = 1)
    colnames(tot_Eta_prec) = traitnames

    # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
    #  Prior: Gamma distribution for each element
    #       shape = tot_F_prec_shape
    #       rate = tot_F_prec_rate
    tot_F_prec = matrix(1,nrow=1,ncol=K)
    #with(priors,matrix(rgamma(K,shape = tot_F_prec_shape,rate = tot_F_prec_rate),nrow=1))

    # Factor scores:

    # Resid discrete variances
    # p-matrix of n_RE x p with
    resid_h2_index = sample(c(1:ncol(h2s_matrix))[h2_priors_resids>0],p,replace=T)
    resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]

    # Factor discrete variances
    # K-matrix of n_RE x K with
    F_h2_index = sample(c(1:ncol(h2s_matrix))[h2_priors_factors>0],K,replace=T)
    F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F = matrix(rnorm(sum(r_RE) * K, 0, sqrt(F_h2[1,] / tot_F_prec)),ncol = K, byrow = T)
    rownames(U_F) = colnames(ZL)

    U_R = matrix(rnorm(sum(r_RE) * p, 0, sqrt(resid_h2[1,] / tot_Eta_prec)),ncol = p, byrow = T)
    colnames(U_R) = traitnames
    rownames(U_R) = colnames(ZL)

    # Fixed effects
    B1 = matrix(rnorm(b1*p), ncol = p)
    colnames(B1) = traitnames

    B2_R = 0*matrix(rnorm(b2_R*p),b2_R,ncol = p)
    colnames(B2_R) = traitnames
    rownames(B2_R) = colnames(X2_R)

    # Factor fixed effects
    B2_F = 0*matrix(rnorm(b2_F * K),b2_F,K)
    rownames(B2_F) = colnames(X2_F)

    XB = X1 %**% B1 + X2_R %*% B2_R

    F = X2_F %*% B2_F + ZL %*% U_F + matrix(rnorm(n * K, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = K, byrow = T)
    F = as.matrix(F)

    cis_effects = matrix(rnorm(length(cis_effects_index)),1,length(cis_effects_index))

    # var_Eta
    if(!'var_Eta' %in% ls()) var_Eta = rep(1,p)


    # ----------------------- #
    # ---Save initial values- #
    # ----------------------- #
    current_state = list(
      K              = K,
      Lambda         = Lambda,
      tot_F_prec     = tot_F_prec,
      F_h2_index     = F_h2_index,
      F_h2           = F_h2,
      U_F            = U_F,
      F              = F,
      tot_Eta_prec   = tot_Eta_prec,
      resid_h2_index = resid_h2_index,
      resid_h2       = resid_h2,
      U_R            = U_R,
      B1              = B1,
      B2_R            = B2_R,
      B2_F            = B2_F,
      XB             = XB,
      cis_effects    = cis_effects,
      nrun           = 0,
      total_time     = 0
    )
    return(current_state)
  })

  # Initialize parameters for Lambda_prior, B_prior, and QTL_prior (may be model-specific)
  MegaLMM_state$current_state = MegaLMM_state$priors$Lambda_prior$sampler(MegaLMM_state)
  MegaLMM_state$current_state = MegaLMM_state$priors$B2_prior$sampler(MegaLMM_state)

  # Initialize Eta
  observation_model_state = run_parameters$observation_model(run_parameters$observation_model_parameters,MegaLMM_state)
  MegaLMM_state$current_state[names(observation_model_state$state)] = observation_model_state$state


  # save the initial RNG state
  MegaLMM_state$current_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )

  # ----------------------------- #
  # --- Initialize MegaLMM_state --- #
  # ----------------------------- #

  Posterior = reset_Posterior(MegaLMM_state$Posterior,MegaLMM_state)
  MegaLMM_state$Posterior = Posterior

  MegaLMM_state$Posterior$posteriorSample_params = unique(c(MegaLMM_state$Posterior$posteriorSample_params,observation_model_state$posteriorSample_params))
  MegaLMM_state$Posterior$posteriorMean_params = unique(c(MegaLMM_state$Posterior$posteriorMean_params,observation_model_state$posteriorMean_params))

  return(MegaLMM_state)
}



#' Initialized Gibbs sampler for MegaLMM model
#'
#' The pre-calculates a set of matrices that will be re-used through the Gibbs chains.
#' These calculations can be slow for large models, especially if n is large, the number of
#' random effects is > 1, or there are many groups of observations with different
#' missing data patterns.
#'
#' @param MegaLMM_state MegaLMM_state object as returned by \link{setup_model_MegaLMM}
#' @param ncores number of cores to use for parallel evaluations. Not really used as RcppParallel is used instead.
#'     Instead, we break up the computation into chunks of this size.
#' @param Qt_list Optionally, \code{Qt_list}, \code{chol_R_list} and \code{chol_ZKZt_list} can be provided
#'     from a previous MegaLMM_state object if the data and model is identical.
#' @param chol_R_list See \code{Qt_list}
#' @param chol_ZKZt_list See \code{Qt_list}
#'
#' @return MegaLMM_state object with \code{Qt_list}, \code{chol_R_list} and \code{chol_ZKZt_list} added to run_variables
#' @export
#'
initialize_MegaLMM = function(MegaLMM_state, ncores = my_detectCores(), Qt_list = NULL, chol_R_list = NULL, chol_ZKZt_list = NULL) {
  # calculates Qt_list, chol_R_list and chol_ZKZt_list
  # returns MegaLMM_state

  run_parameters = MegaLMM_state$run_parameters
  data_matrices = MegaLMM_state$data_matrices
  Y_missing = run_parameters$observation_model_parameters$observation_setup$Y_missing
  h2s_matrix = MegaLMM_state$data_matrices$h2s_matrix
  verbose = run_parameters$verbose

  RE_setup = data_matrices$RE_setup
  ZL = data_matrices$ZL

  Missing_data_map      = MegaLMM_state$run_variables$Missing_data_map
  Missing_row_data_map  = MegaLMM_state$run_variables$Missing_row_data_map


  n = MegaLMM_state$run_variables$n
  p = MegaLMM_state$run_variables$p
  n_RE = length(RE_setup)
  RE_names = names(RE_setup)


  # ------------------------------------ #
  # ----Precalculate ZKZts, chol_Ks ---- #
  # ------------------------------------ #


  # now, for each set of columns, pre-calculate a set of matrices, etc
  # do calculations in several chunks
  n_matrices = 2*ncol(h2s_matrix)

  if(verbose) {
    print(sprintf("Pre-calculating random effect inverse matrices for %d groups of traits and %d sets of random effect weights", length(Missing_data_map), ncol(h2s_matrix)))
    pb = txtProgressBar(min=0,max = length(Missing_data_map) * n_matrices,style=3)
  }

  X1   = data_matrices$X1
  X2_R = data_matrices$X2_R
  X2_F = data_matrices$X2_F
  U2_F = data_matrices$U2_F
  ZL   = data_matrices$ZL
  cis_genotypes = data_matrices$cis_genotypes


  # cholesky decompositions (RtR) of each K_inverse matrix
  # chol_Ki_mats = lapply(RE_setup,function(re) as(chol(as.matrix(re$K_inv)),'dgCMatrix'))
  chol_Ki_mats = lapply(RE_setup,function(re) {
    K_inv = re$K_inv
    if(is(K_inv,'Matrix')) {
      if(isDiagonal(K_inv)) {
        sparseMatrix(i=1:nrow(K_inv),j=1:nrow(K_inv),x=sqrt(diag(K_inv)))
      } else{
        as(chol(forceSymmetric(K_inv)),'dgCMatrix')
      }
    } else {
      as(chol(K_inv),'dgCMatrix')
    }
  })


  Qt_list = list()
  # QtZL_list = list()
  QtX1_list = list()
  QtX1_keepColumns_list = list()
  QtX2_R_list = list()
  Qt_cis_genotypes = list()

  chol_V_list_list = list()
  chol_ZtZ_Kinv_list_list = list()

  svd_K1 = NULL
  for(set in seq_along(Missing_data_map)){
    x = Missing_data_map[[set]]$Y_obs
    cols = Missing_data_map[[set]]$Y_cols
    if(length(x) == 0) next

    # find Qt = svd(ZLKLtZt)$u
    if(ncol(ZL) < nrow(ZL)) {
      # a faster way of taking the SVD of ZLKZLt, particularly if ncol(ZL) < nrow(ZL). Probably no benefit if ncol(K) > nrow(ZL)
      if(is.null(svd_K1)){
        svd_K1 = svd(RE_setup[[1]]$K)
      }
      qr_ZU = qr(RE_setup[[1]]$ZL[x,,drop=FALSE] %*% svd_K1$u)
      R_ZU = drop0(qr.R(qr_ZU,complete=F),tol=run_parameters$drop0_tol)
      Q_ZU = drop0(qr.Q(qr_ZU,complete=T),tol=run_parameters$drop0_tol)
      RKRt = R_ZU %*% diag(svd_K1$d) %*% t(R_ZU)
      svd_RKRt = svd(RKRt)
      RKRt_U = svd_RKRt$u
      if(ncol(Q_ZU) > ncol(RKRt_U)) RKRt_U = bdiag(RKRt_U,diag(1,ncol(Q_ZU)-ncol(RKRt_U)))
      Qt = t(Q_ZU %**% RKRt_U)
    } else{
      ZKZt = with(RE_setup[[1]],ZL[x,,drop=FALSE] %*% K %*% t(ZL[x,,drop=FALSE]))
      if(is(ZKZt,'Matrix') & isDiagonal(ZKZt)) {
        nn = nrow(ZKZt)
        result = list(d = diag(ZKZt),u = sparseMatrix(i=1:nn,j=1:nn,x=1))
      } else{
        result = svd(ZKZt)
      }
      Qt = t(result$u)
    }
    Qt = as(drop0(as(Qt,'dgCMatrix'),tol = run_parameters$drop0_tol),'dgCMatrix')
    if(nnzero(Qt)/length(Qt) > 0.5) Qt = as.matrix(Qt)  # only store as sparse if it is sparse

    QtZL_matrices_set = lapply(RE_setup,function(re) Qt %*% re$ZL[x,,drop=FALSE])
    # QtZL_set = do.call(cbind,QtZL_matrices_set[RE_names])
    # if(nnzero(QtZL_set)/length(QtZL_set) > 0.5)  QtZL_set = as(QtZL_set,'dgCMatrix')
    QtX1_set = Qt %**% X1[x,,drop=FALSE]
    QtX1_keepColumns = c(1:ncol(X1)) %in% caret::findLinearCombos(QtX1_set)$remove == F # assess which columns of X1 are identifiable in this data set
    QtX1_set = QtX1_set[,QtX1_keepColumns,drop=FALSE]  # drop extra columns

    QtX2_R_set = Qt %**% X2_R[x,,drop=FALSE]

    if(length(cis_genotypes) == p) {
      Qt_cis_genotypes[cols] = lapply(cis_genotypes[cols],function(X) Qt %**% X[x,,drop=FALSE])
    }

    Qt_list[[set]]   = Qt
    QtX1_list[[set]] = QtX1_set
    QtX1_keepColumns_list[[set]] = QtX1_keepColumns
    QtX2_R_list[[set]]  = QtX2_R_set

    ZKZts_set = list()
    for(i in 1:n_RE){
      ZKZts_set[[i]] = forceSymmetric(drop0(QtZL_matrices_set[[i]] %*% RE_setup[[i]]$K %*% t(QtZL_matrices_set[[i]]),tol = run_parameters$drop0_tol))
      ZKZts_set[[i]] = as(as(ZKZts_set[[i]],'CsparseMatrix'),'dgCMatrix')
      if(nnzero(ZKZts_set[[i]])/length(ZKZts_set[[i]]) > 0.5) {
        ZKZts_set[[i]] = as.matrix(ZKZts_set[[i]])
      }
    }

    ZtZ_set = as(forceSymmetric(drop0(crossprod(ZL[x,]),tol = run_parameters$drop0_tol)),'dgCMatrix')

    chol_V_list_list[[set]] = make_chol_V_list(ZKZts_set,h2s_matrix,run_parameters$drop0_tol,pb,setTxtProgressBar,getTxtProgressBar,ncores)
    # convert any to dense if possible
    for(i in 1:length(chol_V_list_list[[set]])){
      chol_V_list_list[[set]][[i]] = drop0(chol_V_list_list[[set]][[i]],tol = run_parameters$drop0_tol)
      if(nnzero(chol_V_list_list[[set]][[i]])/length(chol_V_list_list[[set]][[i]]) > 0.25){
        chol_V_list_list[[set]][[i]] = as.matrix(chol_V_list_list[[set]][[i]])
      }
    }
    if(length(RE_setup) == 1 & isDiagonal(ZtZ_set) & all(sapply(chol_Ki_mats,isDiagonal))) {
      # in the case of 1 random effect with diagonal covariance matrix, we can skip the expensive calculations
      chol_ZtZ_Kinv_list_list[[set]] = sapply(1:length(h2s_matrix),function(i) {
        nn = nrow(ZtZ_set)
        sparseMatrix(i=1:nn,j=1:nn,x = sqrt(1/(1-h2s_matrix[1,i])*diag(ZtZ_set) + 1/h2s_matrix[1,i]*diag(chol_Ki_mats[[1]])^2))
      })
    } else{
      chol_ZtZ_Kinv_list_list[[set]] = make_chol_ZtZ_Kinv_list(chol_Ki_mats,h2s_matrix,ZtZ_set,run_parameters$drop0_tol,pb,setTxtProgressBar,getTxtProgressBar,ncores)
    }
  }
  if(verbose) close(pb)

  # Qt matrices for factors are only used with row set 1
  x = Missing_data_map[[1]]$Y_obs
  Qt1_U2_F = Qt_list[[1]] %**% U2_F[x,,drop=FALSE]
  Qt1_X2_F = Qt_list[[1]] %**% X2_F[x,,drop=FALSE]


  MegaLMM_state$run_variables = c(MegaLMM_state$run_variables,
    list(
    Qt_list    = Qt_list,
    # QtZL_list   = QtZL_list,
    QtX1_list   = QtX1_list,
    QtX1_keepColumns_list = QtX1_keepColumns_list,
    QtX2_R_list = QtX2_R_list,
    Qt1_U2_F = Qt1_U2_F,
    Qt1_X2_F = Qt1_X2_F,
    Qt_cis_genotypes = Qt_cis_genotypes,
    chol_V_list_list          = chol_V_list_list,
    chol_ZtZ_Kinv_list_list = chol_ZtZ_Kinv_list_list
  ))

  return(MegaLMM_state)
}

#' Print more detailed statistics on current MegaLMM state
#'
#' Print more detailed statistics on current MegaLMM state
#' @seealso \code{\link{MegaLMM_control}}, \code{\link{sample_MegaLMM}}, \code{\link{MegaLMM_init}},
#'   \code{\link{print.MegaLMM_state}}, \code{\link{plot.MegaLMM_state}}
summary.MegaLMM_state = function(MegaLMM_state){
  with(MegaLMM_state,{
    cat(
      c(sprintf('\n MegaLMM_state object for data of size %d x %d \n',nrow(data_matrices$Y),ncol(data_matrices$Y))),
      c(sprintf('Model dimensions: factors = %d, fixed = %d, regression_R = %d, regression_F = %d, random = %d \n',
                current_state$K,
                ncol(data_matrices$X1),
                ncol(data_matrices$X2_R),ncol(data_matrices$X2_F),
                ncol(data_matrices$ZL))),
      c(sprintf('Sampler: %s \n',run_parameters$sampler)),
      c(sprintf('Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$total_samples)),
      c(sprintf('Total time: %s \n\n',format(current_state$total_time)))
    )
  })
}

#' Print statistics on current MegaLMM state
#'
#' Print statistics on current MegaLMM state
#' @seealso \code{\link{MegaLMM_control}}, \code{\link{sample_MegaLMM}}, \code{\link{MegaLMM_init}},
#'   \code{\link{summary.MegaLMM_state}}, \code{\link{plot.MegaLMM_state}}
#' @export
print.MegaLMM_state = function(MegaLMM_state){
  with(MegaLMM_state,{
    cat(
      c(sprintf('\n Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$total_samples)),
      c(sprintf('Total time: %s \n\n',format(current_state$total_time)))
    )
  })
}

#' Make plots of current MegaLMM state
#'
#' Make plots of current MegaLMM state
#' 
#' @param MegaLMM_state output of sample_MegaLMM
#' @param file Output file for pdf booklet
#' @param setup optional - a list of known values for Lambda (error_factor_lambda), h2, factor_h2s
#' @seealso \code{\link{MegaLMM_control}}, \code{\link{sample_MegaLMM}}, \code{\link{MegaLMM_init}},
#'   \code{\link{print.MegaLMM_state}}, \code{\link{summary.MegaLMM_state}}
plot.MegaLMM_state = function(MegaLMM_state,file = 'diagnostics_plots.pdf',setup = NULL){
  file = sprintf('%s/%s',MegaLMM_state$run_ID,file)
  tempfile = sprintf('%s.temp',file)
  pdf(tempfile)
  if(!is.null(setup)){
    MegaLMM_state$setup = setup
    plot_diagnostics_simulation(MegaLMM_state)
  } else{
    plot_diagnostics(MegaLMM_state)
  }
  dev.off()
  system(sprintf('mv %s %s',tempfile,file))
}
