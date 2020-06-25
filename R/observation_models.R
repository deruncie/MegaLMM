# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


#' Sample missing data
#'
#' Function to sample missing data given model parameters.
#'
#' This is also a template for \code{data_model} functions.
#'
#'    The function should draw a posterior sample for Eta (the (potentially) latent
#'    data on the linear scale) given the factor model and the observed data. It returns a list
#'    with Eta and any other variables associated with the data_model.
#'    These variables are added to current_state, but are not currently saved in Posterior
#'
#'    Initial values, hyperparameters, and any helper functions / parameters can be passed
#'    in \code{observation_model_parameters}.
#'
#'    When run with an empty \code{current_state} and \code{data_matrices == NULL}, it should return
#'    a matrix Eta of the correct dimensions, but the values are unimportant.
#'
#' @param Y data matrix n_Y x p_Y
#' @param observation_model_parameters List of parameters necessary for the data model.
#'      Here, a Matrix of coordinates of NAs in Y
#' @param MegaLMM_state a MegaLMM_state object. Generally, only current_state and data_matrices is used. If
#'    empty, will return a default set of parameters of the appropriate size for model initialization.
#' @return list of data_model variables including:
#' @return state a list of parameters associated with the data_model. Here, only the matrix Eta
#' @return posteriorSample_params a list of parameter names to record posterior samples of
#' @return posteriorMean_params a list of parameters to record the posterior mean, but not save individual
#'     posterior samples
missing_data_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices

  if(!'observation_setup' %in% names(observation_model_parameters)) {
    observation_setup = with(c(observation_model_parameters,data_matrices,current_state),{
      if(scale_Y){
        Mean_Y = colMeans(Y,na.rm=T)
        var_Eta = apply(Y,2,var,na.rm=T)
        Eta = sweep(Y,2,Mean_Y,'-')
        Eta = sweep(Eta,2,sqrt(var_Eta),'/')
      } else {
        Eta = Y
        p_Y = dim(Y)[2]
        Mean_Y = rep(0,p_Y)
        var_Eta = rep(1,p_Y)
      }
      Y_missing = as(is.na(Y),'lgTMatrix')# un-compressed logical sparse matrix
      return(list(
        Eta = Eta,
        n = nrow(Y),
        p = ncol(Y),
        traitnames = colnames(Y),
        Mean_Y = Mean_Y,
        var_Eta = var_Eta,
        Y_missing = Y_missing,
        n_missing = sum(Y_missing),
        missing_indices = which(Y_missing)
      ))
    })
    return(observation_setup)
  }

  observation_model_state = with(c(observation_model_parameters,observation_model_parameters$observation_setup,data_matrices,current_state),{
    Eta_mean = matrix(0,0,0)
    if(n_missing > 0){
      n = nrow(Y)
      p = ncol(Y)
      if(length(current_state) == 0) {
        Eta_mean = matrix(0,n,p)
        resids = rnorm(n_missing)
      } else{
        Eta_mean = XB + F %**% Lambda + ZL %**% U_R
        resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
        resids = rnorm(n_missing,0,sqrt(1/resid_Eta_prec[Y_missing@j+1]))  # sample resids from normal distribution with appropriate variance
      }
      Eta[missing_indices] = Eta_mean[missing_indices] + resids
    }
    return(list(Eta = Eta,Eta_mean = Eta_mean,Y_missing = Y_missing, var_Eta = var_Eta))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Eta'),
              posteriorMean_params = c('Eta_mean')
  ))
}



#' Sample Eta given a voom observation model
#'
#' \code{voom} or \code{voomWithQualityWeights} (limma) should be run on RNAseq data,
#'    resulting in logCPM (E) and inverseWeights (weights)
#'    for each observation. The observations (E) should be passed as Y, and the weights as
#'    observation_model_parameters$prec_Y.
#'
#' When running \code{voom}, a fully-specified fixed effect model for the data should be specified.
#'
#' @inheritParams missing_data_model
#' @references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014).
#'     voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
#'     Genome Biology, 15(2), R29. http://doi.org/10.1186/gb-2004-5-10-r80
voom_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices

  new_variables = c('Eta')
  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{
    Eta = Y
    if(length(current_state) > 0){
      n = nrow(Y)
      p = ncol(Y)
      Eta_mean = XB + F %*% Lambda + ZL %**% U_R
      resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
      prec = sweep(prec_Y,2,resid_Eta_prec,'+')
      # Y_std = Y * prec_Y. This is calculated once in the initialization.
      Eta_hat = (Y_std + sweep(Eta_mean,2,resid_Eta_prec,'*')) / prec
      resid = matrix(rnorm(n*p),n,p) / sqrt(prec)
      Eta = Eta_hat + resid
    }
    return(list(Eta = as.matrix(Eta)))
  })
  return(list(state = observation_model_state[new_variables],
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
             )
         )
}


#' B spline basis with option of centering and differencing coefficients
#'
#' parameters follow \link{bs}
#'
#' uses code from \link{bs}, but optionally transforms the basis into a difference
#' between spline parameters.
#'
#' Also, optionally centers the spline so that the spline average is zero.
#'
#'
#' @return list of data_model variables including: \itemize{
#'     \item state a list of parameters associated with the data_model. Here, only the matrix Eta
#'     \item posteriorSample_params a list of parameter names to record posterior samples of
#'     \item posteriorMean_params a list of parameters to record the posterior mean, but not save individual
#'     posterior samples
#' }
#' @export
#'
#' @references following code from https://github.com/SurajGupta/r-source/blob/master/src/library/splines/R/splines.R
#'
bs_diff = function(x, df = NULL, knots = NULL, degree = 3, intercept = TRUE,
                    Boundary.knots = range(x),
                    differences = 1,
                    periodic = FALSE,
                    center = TRUE
) {
  # following code from https://github.com/SurajGupta/r-source/blob/master/src/library/splines/R/splines.R
  if(periodic){
    if(is.null(knots)) {
      bs_X = pbs::pbs(x=x,df=df,degree=degree,intercept=intercept,Boundary.knots=Boundary.knots)
    } else {
      bs_X = pbs::pbs(x=x,knots=knots,degree=degree,intercept=intercept,Boundary.knots=Boundary.knots)
    }
  } else {
    bs_X = splines::bs(x,df,knots,degree,intercept,Boundary.knots)
  }
  X = bs_X
  diff = differences
  if(differences > 0) {
    # differences transformm the parameters into difference between consecutive parameters.
    # As we do this sequentially, it penalizes higher-order derivatives of the curve
    # we include the averages of the lower-order splines as predictors as well
    m = ncol(X)
    contr = MASS::contr.sdif(m-differences+1)
    if(differences > 1) {
      for(diff in (differences-1):1){
        contr = MASS::contr.sdif(m-diff+1) %*% cbind(1/m,contr)
      }
    }
    X = X %*% contr
  }
  bs_X_attributes = attributes(bs_X)
  bs_X_attributes = bs_X_attributes[names(bs_X_attributes) %in% c('dim','dimnames') == F]
  attributes(X) = c(attributes(X),bs_X_attributes)
  attr(X,'differences') = differences
  attr(X,'periodic') = periodic
  attr(X,'center') = center
  class(X) = c('bs_diff',class(X))
  X
}
makepredictcall.bs_diff <- function(var, call)
{
  if(as.character(call)[1L] != "bs_diff") return(call)
  at <- attributes(var)[c("degree", "knots", "Boundary.knots", "intercept","differences","periodic","center")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

#' Sample Eta given regression-splines individual-level model
#'
#' Eta is a matrix of individual-level parmaters for a regression-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique). During sampling, should sample regression coefficients \code{Eta} given the factor model state.
#'     Set up to allow for non iid errors, but currently not implemented.
#'
#' @param observation_model_parameters a list including:
#'     1) \code{observations}, a data.frame with observation-level data including columns \code{ID} and \code{Y}
#'     3) \code{individual_model} the model that should be applied to the data for each ID.
#' @param MegaLMM_state The current \code{MegaLMM_state}.
#'     For initialization, can be a list with \code{data_matrices = list(data=data)}
#'     Where \code{data} contains \code{ID} for ordering.
regression_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices

  if(!'observation_setup' %in% names(observation_model_parameters)) {
    observation_setup = with(c(observation_model_parameters,data_matrices,current_state),{
      if(!'ID' %in% colnames(data)) stop('ID column required in data')
      if(!length(unique(data$ID)) == nrow(data)) stop('duplicate IDs in data')

      # ensure ID is a character for matching
      observations$ID = as.character(observations$ID)
      data$ID = as.character(data$ID)
      observations = subset(observations,ID %in% data$ID)

      # extract model Terms and Y matrix
      mf = model.frame(individual_model,observations)
      mm = model.matrix(individual_model,observations)
      traits = all.vars(update(individual_model,'~.0'))
      Y = as.matrix(observations[,traits,drop=FALSE])
      Terms = delete.response(terms(mf))

      n_terms = ncol(mm)

      id_index = tapply(1:nrow(observations),observations$ID,function(x) x)
      model_matrices = lapply(data$ID,function(id) {
        x = id_index[[id]]
        if(length(x) > 0){
          X = mm[x,,drop=FALSE]
        } else{
          x = matrix(0,0,0)
          X = matrix(0,0,n_terms)
        }
        # keep only columns of X which are non-zero, but record which columns these were of original X.
        nonZero_cols_X = which(colSums(X!=0) > 0)
        X = X[,nonZero_cols_X,drop=FALSE]
        list(
          X = X,
          y = Y[x,,drop=FALSE],
          position = x,
          nonZero_cols_X = nonZero_cols_X
        )
      })
      names(model_matrices) = data$ID

      n = length(model_matrices)
      n_traits = ncol(model_matrices[[1]]$y) # number of traits
      p_trait = n_terms  # number of coefficients per trait
      p = p_trait * n_traits # total number of coefficients
      traits = colnames(model_matrices[[1]]$y)
      if(length(traits) > 1) {
        traitnames = paste(rep(traits,each = p_trait),rep(colnames(mm),length(traits)),sep='::')
      } else{
        traitnames = colnames(mm)
      }
      Eta_row_names = data$ID
      Y_missing = t(sapply(model_matrices,function(x) rep(!(seq_len(n_terms) %in% x$nonZero_cols_X),n_traits)))

      observation_setup = list(
        Y = Y,
        n = length(model_matrices),
        p = p,
        Terms = Terms,
        n_traits = n_traits,
        traitnames = traitnames,
        Eta_row_names = Eta_row_names,
        model_matrices = model_matrices,
        Y_missing = Y_missing
      )
      return(observation_setup)
    })
    return(observation_setup)
  }

  observation_model_state = with(c(observation_model_parameters,observation_model_parameters$observation_setup,data_matrices,current_state),{
    if(!'var_Eta' %in% ls()) var_Eta = rep(1,p)
    if(length(current_state) == 0){
      Eta_mean = matrix(0,n,p)
      resid_Eta_prec = matrix(1,1,p)
      resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate),nr=1) # only a single precision parameter for the data_model?
    } else{
      Eta_mean = XB + F %*% Lambda + ZL %**% U_R
      resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
      if(!'resid_Y_prec' %in% ls()) resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate),nr=1) # only a single precision parameter for the data_model?
    }

    # re-scale Eta_mean and resid_Eta_prec
    Eta_mean = sweep(Eta_mean,2,sqrt(var_Eta),'*')
    resid_Eta_prec[] = resid_Eta_prec / var_Eta
    for(i in 1:length(model_matrices)){
      model_matrices[[i]]$tot_Y_prec = with(model_matrices[[i]],matrix(resid_Y_prec,nr = nrow(y),nc = ncol(y),byrow=T))
    }
    coefs = sample_coefs_set_c(model_matrices,t(Eta_mean),matrix(resid_Eta_prec,length(resid_Eta_prec),n))

    Y_fitted = get_fitted_set_c(model_matrices,coefs)
    Eta = t(coefs)
    colnames(Eta) = traitnames
    rownames(Eta) = Eta_row_names

    # un-scale Eta
    Eta = sweep(Eta,2,sqrt(var_Eta),'/')

    Y_tilde = Y - Y_fitted

    resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape + 0.5*dim(Y_tilde)[1], rate = resid_Y_prec_rate + 0.5*colSums(Y_tilde^2)),nr=1)

    return(list(Eta = Eta, resid_Y_prec = resid_Y_prec, Y_fitted=Y_fitted, Y=Y, var_Eta = var_Eta))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Y_fitted','Eta','resid_Y_prec'),
              posteriorMean_params = c()
  )
  )
}


#' Sample cis_eQTL coefficients
#'
#' @param observation_model_parameters list with:
#'     \code{Y} gene expression data matrix n x p
#'     \code{cis_genotypes} a list of design matrices of length \code{p} (ie number of columns of \code{Eta})
#'     This is used to specify trait-specific fixed effects, such a cis-genotypes
#'
#' @param MegaLMM_state
#'
#' @return list including:
#'     \code{state} parameters to add to \code{current_state}
#'     \code{posteriorSample_params} character vector of parameter names to include in the list of parameters to record posterior samples
#'     \code{posteriorMean_params} character vector of parameter names to include in the list of parameters to record posterior mean
cis_eQTL_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices

  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{
    n = nrow(Y)
    p = ncol(Y)
    Ut_cis = as(diag(1,n),'dgCMatrix')
    s_cis = rep(1,n)
    if(!exists('Eta')) {
      Eta_mean = matrix(0,n,p)
    } else{
      Eta_mean = XB + F %*% Lambda + ZL %*% U_R
    }
    if(!exists('tot_Eta_prec')) {
      resid_Eta_prec = matrix(1,1,p)
    } else{
      resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
    }

    Y_tilde = Y - Eta_mean
    if(!exists('ncores')) ncores = 1
    cis_effects_list = mclapply(1:p,function(j) {
      cis_X_j = cis_genotypes[[j]]
      b_j = ncol(cis_X_j)
      if(b_j > 0) {
        prior_mean = matrix(0,b_j,1)
        prior_prec = matrix(1e-10,b_j,1)
        prior_prec[apply(cis_X_j,2,var)==0] = 1e10
        randn_theta = matrix(rnorm(b_j),b_j,1)
        randn_e = matrix(rnorm(n),n,1)
        coefs_j = sample_coefs_parallel_sparse_c_Eigen(Ut_cis,Y_tilde[,j],cis_X_j,
                                           0, resid_Eta_prec[,j],
                                           s_cis,prior_mean,prior_prec,
                                           randn_theta,randn_e,
                                           1)
        return(coefs_j)
      }
      return(NULL)
    },mc.cores = ncores)

    cis_fitted = do.call(cbind,mclapply(1:p,function(j) {
      cis_X_j = cis_genotypes[[j]]
      if(ncol(cis_X_j) > 0) {
        return(cis_X_j %*% cis_effects_list[[j]])
      }
      return(rep(0,n))
    },mc.cores = ncores))

    Eta = Y - cis_fitted
    cis_effects = matrix(do.call(c,cis_effects_list),nrow=1)
    return(list(Eta = Eta, cis_effects = cis_effects))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Eta','cis_effects'),
              posteriorMean_params = c()
  )
  )
}


#' Sample Eta given B-splines individual-level model with binomial observations
#'
#' Eta is a matrix of individual-level parmaters for a B-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique)
#'
#' @param Y a vector of observation. Pre-scaled and centered if desired.
#' @param observation_model_parameters a list including:
#'     \code{observations}, a data.frame with columns: \code{ID}, \code{N} and \code{covariate}, ordered by ID.
#'     and the variables df, knots, degree and intercept.
#'     Empty values will use the defaults for \code{bs()}.
bs_binomial_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices

  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{

  log_binom = function(beta,X,y,N,mu, sigma2){
    Xbeta = X %*% beta
    eXbeta = exp(Xbeta)
    p = eXbeta/(1+eXbeta)
    return( sum(y*log(p) + (N-y)*log(1-p))
            - sum((beta - mu)^2)/(2*sigma2)
    )
  }

  if(!exists('model_matrices')){
    if(!exists('df') || !is.numeric(df)) df = NULL
    if(!exists('knots') || !is.numeric(knots)) knots = NULL
    if(!exists('degree') || !is.numeric(degree)) degree = 3
    if(!exists('intercept') || !is.numeric(intercept)) intercept = FALSE
    covariate = observations$covariate
    coefficients = splines::bs(covariate,df = df,knots=knots,degree=degree,intercept = intercept)

    model_matrices = tapply(1:nrow(observations),observations$ID,function(x) {
      list(
        X = Matrix(coefficients[x,]),
        y = Y[x,],
        N = observations$N[x]
      )
    })
  }

  n = length(model_matrices)
  p = ncol(coefficients)

  if(length(current_state) == 0){
    Eta_mean = matrix(0,n,p)
    resid_Eta_prec = matrix(1e-10,1,p)
  } else{
    Eta_mean = X %*% B + F %*% Lambda + ZL %*% U_R
    resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
  }

  if(!exists(ncores)) ncores = 1
  Eta = do.call(rbind,mclapply(1:n,function(i) {
    X = model_matrices[[i]]$X
    y = model_matrices[[i]]$y
    N = model_matrices[[i]]$N
    eta_i = MfUSampler::MfU.Sample(Eta[i,], f=log_binom, uni.sampler="slice", X=X, y=y,N=N)
  },mc.cores = ncores))

    return(list(Eta = Eta, model_matrices = model_matrices, coefficients = coefficients))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
  )
  )
}

probe_gene_model = function(observation_model_parameters,MegaLMM_state = list()){
  current_state = MegaLMM_state$current_state
  data_matrices = MegaLMM_state$data_matrices
  new_variables = c('Eta','resid_Y_prec','mu_probe')

  if(!'observation_setup' %in% names(observation_model_parameters)) {
    observation_setup = with(c(observation_model_parameters,data_matrices,current_state),{
      list(n = nrow(Y),
           p = nrow(Z_Y))

    })
    return(observation_setup)
  }

  observation_model_state = with(c(observation_model_parameters,observation_model_parameters$observation_setup,data_matrices,current_state),{
    p = nrow(Z_Y)
    n = nrow(Y)

    if(!exists('Eta') || !all(dim(Eta) == c(n,p))) Eta = matrix(0,n,p)
    if(!exists('resid_Y_prec')) resid_Y_prec = rep(1,ncol(Y))

    Y_tilde = Y - Eta %*% Z_Y
    mu_probe = colMeans(Y_tilde) + rnorm(ncol(Y))/sqrt(n*resid_Y_prec)

    Y_tilde = sweep(Y,2,mu_probe,'-')
    Y_std = sweep(Y_tilde,2,resid_Y_prec,'*')

    if(!exists('tot_Eta_prec')){
      sum_Y_prec = Z_Y %*% resid_Y_prec
      prec = sum_Y_prec@x
      post_mean = sweep(Y_std %*% t(Z_Y),2,prec,'/')
      resid = sweep(matrix(rnorm(n*p),n,p),2,sqrt(prec),'/')
    } else{
      resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
      Eta_mean = XB + F %*% Lambda + ZL %**% U_R

      Eta_mean_std = sweep(Eta_mean,2,resid_Eta_prec,'*')
      sum_Y_prec = Z_Y %*% resid_Y_prec
      prec = sum_Y_prec@x + resid_Eta_prec

      post_mean = sweep(Y_std %*% t(Z_Y) + Eta_mean_std,2,prec,'/')
      resid = sweep(matrix(rnorm(n*p),n,p),2,sqrt(prec),'/')
    }
    Eta = as.matrix(post_mean) + resid

    Y_tilde = Y_tilde - Eta %*% Z_Y
    resid_Y_prec = rgamma(ncol(Y),shape = resid_Y_prec_shape + 0.5*nrow(Y),
                          rate = resid_Y_prec_rate + 0.5 * colSums(Y_tilde^2))
    return(list(Eta = Eta, mu_probe = mu_probe, resid_Y_prec = resid_Y_prec))
  })
  return(list(state = observation_model_state[new_variables],
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta','resid_Y_prec','mu_probe')
  )
  )
}



#' Sample Eta given a list of different observation_models
#'
#'
#' @param observation_model_parameters a list of observation_model lists
#' @param MegaLMM_state the current MegaLMM_state object. Can be just a list with \code{data}
#'
#' @return list of data_model variables including:
#' @return state a list of parameters associated with the data_model. Here, only the matrix Eta
#' @return posteriorSample_params a list of parameter names to record posterior samples of
#' @return posteriorMean_params a list of parameters to record the posterior mean, but not save individual
#'     posterior samples
#' @export
#'
combined_model = function(observation_model_parameters,MegaLMM_state = list()){

  ## NOTE: Somehow need to pass a subset of the factor model to each sub-model, or maybe calculate Eta_mean once and pass a subset to each model.

  sub_models = observation_model_parameters$sub_models
  n_models = length(sub_models)
  if(is.null(names(sub_models))) names(sub_models) = 1:n_models

  # run each of the separate observation models
  results = lapply(names(sub_models),function(i) {
    Y = sub_models[[i]]
    if(!'observation_model' %in% names(Y)) stop(sprintf('observation_model not specified in model %s',i))
    observation_model = Y$observation_model

    # adjust model-specific parameter names
    if(length(names(MegaLMM_state$current_state))>0){
      names(MegaLMM_state$current_state) = sub(sprintf('.%s',i),'',names(MegaLMM_state$current_state))
    }

    # run observation_model
    new_state = observation_model(Y[names(Y) != 'observation_model'],MegaLMM_state)

    # fix names of parameters
    names(new_state$state) = paste(names(new_state$state),i,sep='.')
    if(length(new_state$posteriorSample_params)>0) new_state$posteriorSample_params = paste(new_state$posteriorSample_params,i,sep='.')
    if(length(new_state$posteriorMean_params)>0) new_state$posteriorMean_params = paste(new_state$posteriorMean_params,i,sep='.')

    # return new_state
    new_state
  })
  names(results) = names(sub_models)

  # merge Etas
  Eta = do.call(cbind,lapply(names(results),function(x) results[[x]]$state[[paste('Eta',x,sep='.')]]))

  # merge Y_missing
  Y_missing = do.call(cbind,lapply(names(results),function(x) results[[x]]$state[[paste('Y_missing',x,sep='.')]]))

  # make combined new_state
  new_state = do.call('c',lapply(results,function(x) x$state))
  names(new_state) = do.call('c',lapply(results,function(x) names(x$state)))
  new_state$Eta = Eta
  new_state$Y_missing = Y_missing

  posteriorSample_params = do.call('c',lapply(results,function(x) x$posteriorSample_params))
  posteriorMean_params = do.call('c',lapply(results,function(x) x$posteriorMean_params))

  return(list(state = new_state,
              posteriorSample_params = posteriorSample_params,
              posteriorMean_params = posteriorMean_params
            )
        )
}




