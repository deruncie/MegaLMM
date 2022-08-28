load_simulation_data = function(file = NULL){
  if(is.null(file)){
    if(file.exists('../setup.RData')) {
      load('../setup.RData')
      for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    } else{
      if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("R.matlab needed to load setup.mat. Please install it.",
             call. = FALSE)
      }
      setup = readMat('../setup.mat')
      for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    }
  } else{
    load(file)
  }
  Y = setup$Y
  K = setup$A
  r = dim(K)[1]
  n = nrow(Y)
  if(dim(setup$X)[1] != n) setup$X = t(setup$X)
  data = data.frame(animal = gl(r,n/r))
  rownames(K) = data$animal
  if(is.null(colnames(Y))) colnames(Y) = paste('Trait',1:ncol(Y),sep='_')
  return(list(Y = Y, data = data, K_mats = list(animal = K),setup = setup))
}



#' Multiply matrices
#'
#' Multiplies two matrices. Like \code{\%*\%}, but always returns a base matrix, not a Matrix, even when
#' one or both of the two matrices are of class matrix.
#'
#' Uses RcppEigen when one or both matrices are of class dgCMatrix
#'
#' @param X1 matrix-like object
#' @param X2 matrix-like object
#'
#' @return matrix
#' @export
#'
`%**%` = function(X1,X2){
  if(is.null(X1)) return(X2)
  if(is.null(X2)) return(X1)
  result = matrix_multiply_toDense(X1,X2)
  rownames(result) = rownames(X1)
  colnames(result) = colnames(X2)
  return(result)
  # if(inherits(X1,'dgCMatrix') && inherits(X2,'matrix')) return(SxD(X1,X2))
  # if(inherits(X1,'dgCMatrix') && inherits(X2,'dgCMatrix')) return(SxS(X1,X2))
  # if(inherits(X1,'matrix') && inherits(X2,'matrix')) return(X1 %*% X2)
  # if(inherits(X1,'data.frame') || inherits(X2,'data.frame')) return(as.matrix(X1) %**% as.matrix(X2))
  # return(as.matrix(X1 %*% X2))
}


#' Modification to \code{image()}.
#'
#' Like `image,CHMfactor-method`, but with better defaults for correlation matrices especially.
#'
#' The default is to have the colors evenly spaced, instead of having the full range <0 and >0.
#'
#' @param X numeric matrix to plot.
#' @param dimnames should row/colnames of matrix be added to the plot
#' @param title Optional title
#' @param include_zero should the z-value zero be included in the scale?
#'
#' @return ggplot object
#' @export
#'
Image = function(X,dimnames=FALSE,title = NULL,include_zero = TRUE,...) {
  require(ggplot2)
  X = as.matrix(X)
  X[] = as.numeric(X)
  if(!dimnames) rownames(X) <- colnames(X) <- NULL
  if(length(unique(rownames(X))) < nrow(X)) rownames(X) <- NULL
  if(length(unique(colnames(X))) < ncol(X)) colnames(X) <- NULL

  X_tall = reshape2::melt(X)
  colnames(X_tall) = c('Var1','Var2','value')
  if(!is.null(colnames(X))) X_tall$Var2 = factor(X_tall$Var2,levels = colnames(X))
  if(!is.null(rownames(X))) X_tall$Var1 = factor(X_tall$Var1,levels = rev(rownames(X)))
  X_tall$value = as.numeric(X_tall$value)
  X_tall = X_tall[X_tall$value != 0,,drop=FALSE] # make it sparse again
  p <- ggplot(X_tall,aes(x=Var2,y=Var1,fill=value)) + geom_tile(alpha = 0.8,height=1,width=1) + xlab('') + ylab('') + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if(include_zero) p <- p + expand_limits(fill=0)
  if(!is.null(title)) p <- p + ggtitle(title)
  if(is.numeric(X_tall$Var2)) p <- p + xlim(0,ncol(X)+1)
  if(is.numeric(X_tall$Var1)) p <- p + ylim(nrow(X)+1,0)
  if(length(unique(X_tall$value))>2) p <- p + scale_fill_gradient2(na.value = "grey90",...)
  # print(p)
  p
}
Image2=function(x,zlim = NULL,breaks=20,colors = c('blue','white','red'),colorkey = TRUE,aspect=NULL,...){
  # if zlim not passed and the range of the data is outside of (-1,1), expands zlim range
  if(missing(zlim)){
    if(all(x[!is.na(x)]>0)) {
      zlim = c(0,max(x,na.rm=T))
    } else{
      zlim = c(-1,1)*max(abs(x),na.rm=T)
    }
  }
  zlim[is.na(zlim)] = 0
  if(missing(aspect)){
    if(ncol(x) > (2*nrow(x))) aspect = 1/1.5
    if(nrow(x) > (2*ncol(x))) aspect = 1.5
  }
  at = seq(zlim[1]-1e-10,zlim[2]+1e-10,length=breaks)
  colors = colorRampPalette(colors)(breaks)
  image(Matrix(x),at=at,col.regions=colors,colorkey=colorkey,aspect=aspect,...)
}

get_ki = function(MegaLMM_state){
  # calculation is from Piironen and Vehtari - 2.4
  # adjust for emperical variance of the column of F - since this should really be set to 1
  # could adjust for 1/tot_F_prec instead, since that's actually a parameter.
  # If sigma2 is in the prior for lambda, than I think it should cancel out here.
  get_k = function(n,sigma2,tau2,lambda2,s2=1){
    # 1/(1+n/sigma2 * tau2 * s2 * lambda2)
    1/(1+n * tau2 * s2 * lambda2)
  }
  current_state = MegaLMM_state$current_state
  current_state = within(current_state, {
    tauh = tauh / colMeans(F^2)
    p = ncol(Lambda)
    n = nrow(F)
  })

  ki = with(current_state,get_k(n,1/t(tot_Eta_prec[rep(1,k),]),
                                (Lambda_omega2[1]/tauh[rep(1,p),]),
                                Lambda_phi2)
  )
  ki
}
get_meff = function(MegaLMM_state){
  ki = get_ki(MegaLMM_state)
  meff = colSums(1-ki)
  meff
}

get_factor_variance = function(MegaLMM_state) {
  vars = with(MegaLMM_state$current_state, {
    apply(F,2,var)*rowMeans(Lambda^2)
  })
  vars
}

# change_K = function(MegaLMM_state,
#                     increment = 0.05*MegaLMM_state$run_variables$Kr,
#                     max_K = min(MegaLMM_state$run_variables$n/4,MegaLMM_state$run_variables$p),
#                     stat = NULL,
#                     threshold = 0.01)

#' Re-orders factors in decreasing order of magnitude
#' 
#' Re-orders factors in decreasing order of magnitude
#' 
#' @param MegaLMM_state output of \code{\link{sample_MegaLMM}}
#' @param factor_order optional: vector or new orders for the factors. If Null, will be calculated
#' @param drop_cor_threshold Factors correlated above this value (either F or U_F) will be dropped
#'
#' @seealso \code{\link{sample_MegaLMM}}, \code{\link{plot.MegaLMM_state}}
reorder_factors = function(MegaLMM_state,factor_order = NULL, drop_cor_threshold = 0.6){
  # first drop factors where F is correlated beyond a threshold
  MegaLMM_state = drop_correlated_factors(MegaLMM_state,cor_threshold = drop_cor_threshold)

  # Now, re-orders factors in decreasing size of F %*% lambda
  # based on current state
  # also re-orders Posterior

  current_state = MegaLMM_state$current_state

  # reorder factors based on var(Lambda) * var(F)
  Lambda = current_state$Lambda
  F = current_state$F
  fixed_factors = MegaLMM_state$run_variables$fixed_factors

  # size is sum lambda_ij^2 * var(F_i)
  if(!exists('factor_order') || is.null(factor_order)) {
    if('Lambda_m_eff' %in% names(current_state)) {
      factor_order = order(MegaLMM_state$current_state$Lambda_m_eff,decreasing = TRUE)
    } else{
      sizes = (rowSums(Lambda^2) * colMeans(F^2))[!fixed_factors]
      factor_order = order(sizes,decreasing=TRUE)
    }
    factor_order_withFixed = factor_order
    if(sum(fixed_factors)>0) {
      factor_order_withFixed = c(1:sum(fixed_factors),sum(fixed_factors)+factor_order)
    }
  }

  # reorder current state
  Lambda_params = c('Lambda','Lambda_prec','Lambda_pi','Lambda_beta','Lambda_delta')
  F_params = c('F','B2_F','U_F','F_h2','tot_F_prec','B2_F_prec','B2_F_pi','B2_F_delta','B2_F_beta')
  
  reorder_params = c('Lambda','Lambda_prec','Plam',
                     'delta',
                     'F','B2_F','U_F','F_h2','F_e_prec','tot_F_prec', 'B2_F_prec'
  )
  for(param in Lambda_params) {
    if(!param %in% names(current_state)) next
    if(ncol(current_state[[param]]) == nrow(Lambda)) {
      current_state[[param]] = current_state[[param]][,factor_order_withFixed,drop=FALSE]
    } else if(nrow(current_state[[param]]) == nrow(Lambda)){
      # for Lambda
      current_state[[param]] = current_state[[param]][factor_order_withFixed,,drop=FALSE]
    } else{
      # for Lambda parameters
      current_state[[param]] = current_state[[param]][factor_order,,drop=FALSE]
    }
  }
  for(param in F_params) {
    if(!param %in% names(current_state)) next
    current_state[[param]] = current_state[[param]][,factor_order_withFixed,drop=FALSE]
  }

  current_state$delta[1] = 1
  # current_state$delta = matrix(c(current_state$tauh[1],current_state$tauh[-1]/current_state$tauh[-length(current_state$tauh)]),nrow=1)
  MegaLMM_state$current_state = current_state

  # reset Lambda_prec and B_prec
  for(i in 1:10){
    # -----Sample Lambda_prec ------------- #
    MegaLMM_state$current_state = MegaLMM_state$priors$Lambda_prior$sampler(MegaLMM_state,1)

    # -----Sample B2_prec ------------- #
    MegaLMM_state$current_state = MegaLMM_state$priors$B2_prior$sampler(MegaLMM_state,1)
  }
  
  # not re-ordering Posterior. Should be cleared anyway.
  
  # # reorder Posterior
  # Posterior = MegaLMM_state$Posterior
  # 
  # for(param in reorder_params){
  #   if(! param %in% names(Posterior)) next
  #   if(dim(Posterior[[param]])[1] == 0) next
  #   Posterior[[param]] = Posterior[[param]][,,factor_order,drop=FALSE]
  # }
  # 
  # MegaLMM_state$Posterior = Posterior
  # just clear the Posterior
  MegaLMM_state = clear_Posterior(MegaLMM_state)

  return(MegaLMM_state)
}

drop_correlated_factors = function(MegaLMM_state,cor_threshold = 0.6){
  # drop factors correlated beyond abs(cor_threshold)
  # also set Lambda_m_eff to 0
  F = MegaLMM_state$current_state$F - MegaLMM_state$data_matrices$X2_F %*% MegaLMM_state$current_state$B2_F
  cor_F = abs(cor(F))
  cor_F[lower.tri(cor_F,diag = T)] = 0
  K = ncol(cor_F)
  for(i in 1:(K-1)) {
    drop_cols = which(cor_F[i,]>cor_threshold)
    if(length(drop_cols) > 0) {
      print(sprintf('dropping cols: %s',paste(drop_cols,collapse=',')))
      MegaLMM_state$current_state$F[,drop_cols] = 0
      if('Lambda_m_eff' %in% names(MegaLMM_state$current_state)) MegaLMM_state$current_state$Lambda_m_eff[drop_cols] = 0
      cor_F[,drop_cols] = 0
    }
  }

  # repeat for Z %*% U_F
  # drop factors correlated beyond abs(cor_threshold)
  # also set Lambda_m_eff to 0
  F = MegaLMM_state$data_matrices$ZL %**% MegaLMM_state$current_state$U_F
  cor_F = abs(cor(F))
  cor_F[lower.tri(cor_F,diag = T)] = 0
  K = ncol(cor_F)
  for(i in 1:(K-1)) {
    drop_cols = which(cor_F[i,]>cor_threshold)
    if(length(drop_cols) > 0) {
      print(sprintf('dropping cols: %s',paste(drop_cols,collapse=',')))
      MegaLMM_state$current_state$F[,drop_cols] = 0
      if('Lambda_m_eff' %in% names(MegaLMM_state$current_state)) MegaLMM_state$current_state$Lambda_m_eff[drop_cols] = 0
      cor_F[,drop_cols] = 0
    }
  }

  return(MegaLMM_state)
}


rescale_factors_F = function(MegaLMM_state){
  # rescale factors based on F
  MegaLMM_state$current_state = within(MegaLMM_state$current_state,{
    F_sizes = colMeans(F^2)
    F = sweep(F,2,sqrt(F_sizes),'/')
    B2_F = sweep(B2_F,2,sqrt(F_sizes),'/')
    U_F = sweep(U_F,2,sqrt(F_sizes),'/')
    B2_F_prec = sweep(B2_F_prec,2,F_sizes,'*')
    Lambda = sweep(Lambda,1,sqrt(F_sizes),'*')
    delta_factor = c(F_sizes[1],exp(diff(log(F_sizes))))
    delta[] = delta / delta_factor
    Lambda_prec[] = sweep(Lambda_prec,1,cumprod(delta),'*')
  })
  return(MegaLMM_state)
}

calc_functions = function(MegaLMM_state,functions) {
  # recover()
  MegaLMM_state$current_state = remove_nuisance_parameters(MegaLMM_state)
  terms = unique(unlist(lapply(functions,function(FUN) {
    # FUN = match.call()[[3]]
    if(is(FUN,'character')){
      FUN = parse(text=FUN)
    }
    terms = all.vars(FUN)
    terms
  })))
  extra_terms = terms[terms %in% with(MegaLMM_state,c(names(current_state),names(data_matrices),names(priors))) == F]
  extra_env = list()
  for(term in extra_terms){
    if(term %in% ls(parent.frame(2))) {
      extra_env[[term]] = parent.frame(2)[[term]]
    } else if(term %in% ls(parent.frame(3))) {
      extra_env[[term]] = parent.frame(3)[[term]]
    }
  }
  base_env = with(MegaLMM_state,c(data_matrices,priors,current_state))
  base_env = c(base_env,extra_env)
  env = c(MegaLMM_state$current_state,base_env)
  # if(!all(terms %in% names(env))) stop(sprintf('Terms %s not found',paste(terms[terms %in% names(env) == F],collapse=', ')))
  result = lapply(functions,function(FUN) {
    if(is(FUN,'character')){
      FUN = parse(text=FUN)
    }
    result = try(eval(FUN,envir = env))
    result = as.matrix(result)
  })
  names(result) = names(functions)
  result
}
  
remove_nuisance_parameters = function(MegaLMM_state) {
  Missing_data_map = MegaLMM_state$run_variables$Missing_data_map
  RE_L_list = MegaLMM_state$data_matrices$RE_L_list
  # S_list = MegaLMM_state$run_variables$S_list
  current_state = within(MegaLMM_state$current_state,{
    # re-transform random effects using RE_L (RE_L %*% diag(D) %*% t(RE_L) = bdiag(K_mats))
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      if(length(cols)) next
      U_R[,cols] = RE_L_list[[set]] %**% U_R[,cols,drop=FALSE]
    }
    U_F = RE_L_list[[1]] %**% U_F
    # transform variables so that the variance of each column of F is 1.
    F_var = 1/tot_F_prec
    U_F[] = sweep(U_F,2,sqrt(F_var),'/')
    B2_F[] = sweep(B2_F,2,sqrt(F_var),'/')
    F[] = sweep(F,2,sqrt(F_var),'/')
    Lambda[] = sweep(Lambda,1,sqrt(F_var),'*')
    
    # re-scale by var_Eta
    if(!'var_Eta' %in% ls()) var_Eta = rep(1,ncol(Lambda))
    U_R[] = sweep(U_R,2,sqrt(var_Eta),'*')
    B1[] = sweep(B1,2,sqrt(var_Eta),'*')
    B2_R[] = sweep(B2_R,2,sqrt(var_Eta),'*')
    Lambda[] = sweep(Lambda,2,sqrt(var_Eta),'*')
    Eta[] = sweep(Eta,2,sqrt(var_Eta),'*')
    tot_Eta_prec[] = tot_Eta_prec / var_Eta
    
    # add means to Eta/Eta_mean
    
    Eta[] = sweep(Eta,2,MegaLMM_state$run_parameters$observation_model_parameters$observation_setup$Mean_Y,'+')
    if(ncol(Eta_mean) == length(var_Eta)) {
      Eta_mean[] = sweep(Eta_mean,2,sqrt(var_Eta),'*')
      Eta_mean[] = sweep(Eta_mean,2,MegaLMM_state$run_parameters$observation_model_parameters$observation_setup$Mean_Y,'+')
    }

  })
  return(current_state)
}

#' Saves current state in Posterior
#'
#' Saves current state in Posterior
#' @seealso \code{\link{sample_MegaLMM}}, \code{\link{plot.MegaLMM_state}}
save_posterior_sample = function(MegaLMM_state) {
  # All parameters in current are matrices.
  # Posterior is a list of arrays.
  # All factor parameters are matrices with ncol == k
  # Posterior arrays are expanded / contracted as the number of factors changes (with update_k)
  # values are re-scaled by var_Eta so that they are on the scale of the original data

  current_state = remove_nuisance_parameters(MegaLMM_state)
  Posterior = MegaLMM_state$Posterior

  total_samples = Posterior$total_samples + 1
  sp_num = Posterior$sp_num + 1
  Posterior$total_samples = total_samples
  Posterior$sp_num = sp_num

  # current_state = within(current_state,{
  #   # re-transform random effects using RE_L (RE_L %*% diag(D) %*% t(RE_L) = bdiag(K_mats))
  #   U_R = MegaLMM_state$data_matrices$RE_L %**% U_R
  #   U_F = MegaLMM_state$data_matrices$RE_L %**% U_F
  #   # transform variables so that the variance of each column of F is 1.
  #   F_var = 1/tot_F_prec
  #   U_F[] = sweep(U_F,2,sqrt(F_var),'/')
  #   B2_F[] = sweep(B2_F,2,sqrt(F_var),'/')
  #   F[] = sweep(F,2,sqrt(F_var),'/')
  #   Lambda[] = sweep(Lambda,1,sqrt(F_var),'*')
  # 
  #   # re-scale by var_Eta
  #   if(!'var_Eta' %in% ls()) var_Eta = rep(1,ncol(Lambda))
  #   U_R[] = sweep(U_R,2,sqrt(var_Eta),'*')
  #   B1[] = sweep(B1,2,sqrt(var_Eta),'*')
  #   B2_R[] = sweep(B2_R,2,sqrt(var_Eta),'*')
  #   Lambda[] = sweep(Lambda,2,sqrt(var_Eta),'*')
  #   Eta[] = sweep(Eta,2,sqrt(var_Eta),'*')
  #   tot_Eta_prec[] = tot_Eta_prec / var_Eta
  # })

  sp = dim(Posterior$Lambda)[1]

  for(param in Posterior$posteriorSample_params){
    # parameters shouldn't change dimension here
    record_sample_Posterior_array(current_state[[param]],Posterior[[param]],sp_num)
  }

  for(param in Posterior$posteriorMean_params){
    Posterior[[param]] = (Posterior[[param]]*(total_samples - 1) + current_state[[param]])/total_samples
  }
  
  if(!is.null(Posterior$posteriorFunctions)) {
    # MegaLMM_state$current_state = current_state
    function_results = calc_functions(MegaLMM_state,Posterior$posteriorFunctions)
    for(param in names(function_results)) {
      # parameters shouldn't change dimension here
      record_sample_Posterior_array(function_results[[param]],Posterior[[param]],sp_num)
    }
  }

  return(Posterior)
}

reset_Posterior = function(Posterior,MegaLMM_state){
  Posterior = list(
    posteriorSample_params = Posterior$posteriorSample_params,
    posteriorMean_params = Posterior$posteriorMean_params,
    posteriorFunctions = Posterior$posteriorFunctions,
    total_samples = Posterior$total_samples,
    folder = sprintf('%s/Posterior',MegaLMM_state$run_ID),
    files = Posterior$files
  )
  
  current_state = remove_nuisance_parameters(MegaLMM_state)

  # re-transform random effects using RE_L
  # current_state$U_F = MegaLMM_state$data_matrices$RE_L %*% current_state$U_F
  # current_state$U_R = MegaLMM_state$data_matrices$RE_L %*% current_state$U_R

  for(param in Posterior$posteriorSample_params){
    if(param %in% names(current_state) && length(dim(current_state[[param]])) == 2) {
      Posterior[[param]] = array(0,dim = c(0,dim(current_state[[param]])))
      dimnames(Posterior[[param]])[2:3] = dimnames(current_state[[param]])
    } else{
      # drop param from Posterior$posteriorSample_params
      Posterior$posteriorSample_params = Posterior$posteriorSample_params[Posterior$posteriorSample_params != param]
      Posterior[[param]] = NULL
    }
  }
  for(param in Posterior$posteriorMean_params) {
    if(param %in% names(current_state) && length(dim(current_state[[param]])) == 2) {
      # if(Posterior$total_samples == 0) {  # I think this is wrong. We always want to set it back to 0.
        Posterior[[param]] = array(0,dim = dim(current_state[[param]]))
        dimnames(Posterior[[param]]) = dimnames(current_state[[param]])
      # }
    } else{
      # drop param from Posterior$posteriorMean_params
      Posterior$posteriorMean_params = Posterior$posteriorMean_params[Posterior$posteriorMean_params != param]
      Posterior[[param]] = NULL
    }
  }
  if(!is.null(Posterior$posteriorFunctions)) {
    function_results = calc_functions(MegaLMM_state,Posterior$posteriorFunctions)
    for(param in names(function_results)) {
      Posterior[[param]] = array(0,dim = c(0,dim(function_results[[param]])))
      dimnames(Posterior[[param]])[2:3] = dimnames(function_results[[param]])
    }
  }
  Posterior$sp_num = 0
  Posterior
}

expand_Posterior = function(Posterior,size){
  for(param in c(Posterior$posteriorSample_params,names(Posterior$posteriorFunctions))){
    Posterior[[param]] = abind(Posterior[[param]],array(NA,dim = c(size,dim(Posterior[[param]])[2:3])),along = 1)
  }
  Posterior
}

#' Resets Posterior samples
#'
#' Clears and resets the saved Posterior samples. Updates the burn parameter of
#'     \code{run_parameters} to reflect all previous samples in the chain now count as burnin
#' @seealso \code{\link{sample_MegaLMM}}, \code{\link{plot.MegaLMM_state}}
clear_Posterior = function(MegaLMM_state) {
  # resets Posterior samples if burnin was not sufficient
  Posterior = MegaLMM_state$Posterior
  run_parameters = MegaLMM_state$run_parameters

  run_parameters$burn = max(run_parameters$burn,MegaLMM_state$current_state$nrun)

  Posterior$total_samples = 0
  Posterior = reset_Posterior(Posterior,MegaLMM_state)

  if(length(list.files(path = Posterior$folder))>0) system(sprintf('rm %s/*',Posterior$folder))
  Posterior$files = c()

  MegaLMM_state$Posterior = Posterior
  MegaLMM_state$run_parameters = run_parameters
  return(MegaLMM_state)
}


#' Saves a chunk of posterior samples
#'
#' Saves a chunk of posterior samples, each paramter in its own RData file.
#' The complete chain for a particular parameter can be re-loaded with \link{load_posterior_param}
#'
#' @param Posterior a Posterior list from a MegaLMM_state object
#' @param folder the folder to save the RData files
#' @return Posterior
save_posterior_chunk = function(MegaLMM_state){
  Posterior = MegaLMM_state$Posterior
  folder = Posterior$folder
  if(!dir.exists(folder)) dir.create(folder)
  file_suffix = sprintf('%d.rds',Posterior$total_samples)
  res = sapply(c(Posterior$posteriorSample_params,Posterior$posteriorMean_params,names(Posterior$posteriorFunctions)),function(param) {
    file_name = sprintf('%s/%s_%s',folder,param,file_suffix)
    samples = Posterior[[param]]
    if(length(samples) > 0) {
      saveRDS(samples,file = file_name,compress = FALSE)
    }
  })
  # print(sprintf('%d files',length(grep(file_suffix,list.files(path = folder)))))
  if(length(grep(file_suffix,list.files(path = folder)))>0) {
    Posterior$files = unique(c(Posterior$files,file_suffix))
  }
  Posterior = reset_Posterior(Posterior,MegaLMM_state)
  MegaLMM_state$Posterior = Posterior
  saveRDS(Posterior,file = sprintf('%s/Posterior_base.rds',folder))
  return(MegaLMM_state)
}

#' load the posterior samples of a single parameter from all saved chunks
#'
#' @param folder folder to find Posterior chunk files
#' @param param Name of parameter to load
#' @param samples vector of sample indices to load. If NULL, all samples loaded
#' @return array with all chuncks of Posterior samples appended together
load_posterior_param = function(MegaLMM_state,param,chunks = NULL){
  if(length(MegaLMM_state$Posterior$files) == 0) return(c())
  if(grepl('RData',MegaLMM_state$Posterior$files[1])) {
    return(load_posterior_param_old(MegaLMM_state,param,chunks))
  }
  folder = MegaLMM_state$Posterior$folder
  param_files = paste0(folder,'/',param,'_',MegaLMM_state$Posterior$files)
  all_files = list.files(path=folder,full.names = T)
  param_files = param_files[param_files %in% all_files]
  n_files = length(param_files)
  if(is.null(chunks)) chunks = 1:n_files
  param_files = na.omit(param_files[chunks])
  if(length(param_files) == 0) return(c())

  samples = readRDS(param_files[1])
  samples_dim = dim(samples)
  current_row = 0
  if(param %in% MegaLMM_state$Posterior$posteriorMean_params) {
    all_samples = samples / length(param_files)
  } else{
    all_samples = array(0,dim = c(MegaLMM_state$Posterior$total_samples*length(chunks)/n_files,samples_dim[2:3]))
    dimnames(all_samples) = dimnames(samples)
    all_samples[1:samples_dim[1],,] = samples
    current_row = samples_dim[1]
    if(current_row == 0) return(c())
  }

  # load other files
  if(length(param_files) > 1) {
    for(i in 2:length(param_files)){
      samples = readRDS(param_files[i])
      if(length(samples_dim) == 2) {
        all_samples = all_samples + samples / length(param_files)
      } else{
        samples_dim = dim(samples)
        if(current_row+samples_dim[1] > dim(all_samples)[1]){
          all_samples = abind(all_samples,array(0,dim = c(current_row+samples_dim[1] - dim(all_samples)[1],dim(all_samples)[2:3])),along = 1)
        }
        if(samples_dim[2] > dim(all_samples)[2]){
          all_samples = abind(all_samples,array(0,dim = c(samples_dim[1],samples_dim[2]-dim(all_samples)[2],dim(all_samples)[3])),along = 2)
        }
        if(samples_dim[3] > dim(all_samples)[3]){
          all_samples = abind(all_samples,array(0,dim = c(dim(all_samples)[1:2],samples_dim[3]-dim(all_samples)[3])),along = 3)
        }
        all_samples[current_row+1:dim(samples)[1],1:samples_dim[2],1:samples_dim[3]] = samples
        current_row = current_row+dim(samples)[1]
      }
    }
  }

  return(all_samples)
}

#' load the posterior samples of a single parameter from all saved chunks
#'
#' @param folder folder to find Posterior chunk files
#' @param param Name of parameter to load
#' @return array with all chuncks of Posterior samples appended together
load_posterior_param_old = function(MegaLMM_state,param,chunks=NULL){
  if(length(MegaLMM_state$Posterior$files) == 0) return(c())
  folder = MegaLMM_state$Posterior$folder
  n_files = length(MegaLMM_state$Posterior$files)
  param_files = paste0(folder,'/',param,'_',MegaLMM_state$Posterior$files)
  all_files = list.files(path=folder,full.names = T)
  param_files = param_files[param_files %in% all_files]
  n_files = length(param_files)
  if(is.null(chunks)) chunks = 1:n_files
  param_files = na.omit(param_files[chunks])
  if(length(param_files) == 0) return(c())

  load(param_files[1])
  samples_dim = dim(samples)
  current_row = 0
  if(param %in% MegaLMM_state$Posterior$posteriorMean_params) {
    all_samples = samples / length(param_files)
  } else{
    all_samples = array(0,dim = c(MegaLMM_state$Posterior$total_samples*length(chunks)/n_files,samples_dim[2:3]))
    all_samples[1:samples_dim[1],,] = samples
    current_row = samples_dim[1]
    if(current_row == 0) return(c())
  }

  # load other files
  if(length(param_files) > 1) {
    for(i in 2:length(param_files)){
      load(param_files[i])
      if(length(samples_dim) == 2) {
        all_samples = all_samples + samples / length(param_files)
      } else{
        samples_dim = dim(samples)
        if(current_row+samples_dim[1] > dim(all_samples)[1]){
          all_samples = abind(all_samples,array(0,dim = c(current_row+samples_dim[1] - dim(all_samples)[1],dim(all_samples)[2:3])),along = 1)
        }
        if(samples_dim[2] > dim(all_samples)[2]){
          all_samples = abind(all_samples,array(0,dim = c(samples_dim[1],samples_dim[2]-dim(all_samples)[2],dim(all_samples)[3])),along = 2)
        }
        if(samples_dim[3] > dim(all_samples)[3]){
          all_samples = abind(all_samples,array(0,dim = c(dim(all_samples)[1:2],samples_dim[3]-dim(all_samples)[3])),along = 3)
        }
        all_samples[current_row+1:dim(samples)[1],1:samples_dim[2],1:samples_dim[3]] = samples
        current_row = current_row+dim(samples)[1]
      }
    }
  }

  return(all_samples)
}

#' Re-loads a full Posterior list with all parameters
#'
#' @param MegaLMM_state a MegaLMM_state object (with empty Posterior)
#' @param params list of parameters to load. If NULL, all parameters will be loaded
#' @return Posterior list, as part of a MegaLMM_state object
reload_Posterior = function(MegaLMM_state,params = NULL){
  folder = MegaLMM_state$Posterior$folder
  Posterior = readRDS(paste(folder,'Posterior_base.rds',sep='/'))
  MegaLMM_state$Posterior = Posterior
  if(is.null(params)) params = c(Posterior$posteriorSample_params,Posterior$posteriorMean_params)
  for(param in params){
    Posterior[[param]] = load_posterior_param(MegaLMM_state,param)
    try(dimnames(Posterior[[param]]) <- dimnames(MegaLMM_state$Posterior[[param]]),silent=T)
  }
  Posterior
}


# pull out specific parameters from a specific sample from Posterior. Only returns
# parameters in listed in \code{terms}, and only if they have posterior samples (not posterior means)
make_current_state = function(Posterior,sample,terms){
  sample_terms = terms[terms %in% Posterior$posteriorSample_params]
  current_state = lapply(sample_terms,function(x) array(Posterior[[x]][sample,,],dim = dim(Posterior[[x]])[-1],dimnames = dimnames(Posterior[[x]])[-1]))
  names(current_state) = sample_terms
  current_state
}

#' Calculates the posterior mean of a function of parameters
#'
#' This function will apply the supplied function to each posterior sample of the chain. Variables
#'    referenced in FUN will be selected from the following search locations (in this order):
#'    1) Posterior$posteriorSample_params
#'    2) data_matrices
#'    3) priors
#'    4) Posterior$posteriorMean_params
#'    5) current_state
#'    6) calling environment (ex. sapply)
#'    7) global environment
#'
#' If mc.cores > 1, the operation will be parallelized. To reduce memory requirements,
#' if mc.cores > 1, a PSOCK cluster is created. This cluster only copyies in data when it
#' is actually used. The code is written so that only necessary objects are used by this cluster,
#' not the whole user's environment. Inside this cluster, mclapply is run, which forks the process inside the PSOCK
#' cluster. Creating the cluster takes some time, so often mc.cores=1 is faster. If the operation
#' is large and Posterior has many samples, this can still create memory issues, so limiting mc.cores
#' might still be necessary
#'
#' @param MegaLMM_state A MegaLMM_state object including a re-loaded Posterior list
#' @param FUN Operations to be applied to each posterior sample. Write as if this were operating
#'     within current_state. Can use priors, data_matrices, and other elements of current_state
#' @param samples (optional) vector of sample indexes to use in the computation
#' @param mc.cores (optional) number of cores to use for computations. See note about memory requirements.
#'
#' @return array of n_samples x dim1 x dim2 where dim1 and dim2 are the dimensions of the calculated
#'     parameter per posterior sample
get_posterior_FUN = function(MegaLMM_state,FUN,samples = NULL,mc.cores = 1) {
  FUN = match.call()[[3]]
  if(is(FUN,'character')){
    FUN = parse(text=FUN)
  }
  terms = all.vars(FUN)
  extra_terms = terms[terms %in% with(MegaLMM_state,c(names(data_matrices),names(priors),names(Posterior))) == F]
  extra_env = list()
  for(term in extra_terms){
    if(term %in% ls(parent.frame(2))) {
      extra_env[[term]] = parent.frame(2)[[term]]
    } else if(term %in% ls(parent.frame(3))) {
      extra_env[[term]] = parent.frame(3)[[term]]
    }
  }
  if(is.null(samples)) {  # count # available samples for the first term in terms (assuming it is in Posterior)
    term1 = terms[terms %in% MegaLMM_state$Posterior$posteriorSample_params][1]
    if(!term1 %in% names(MegaLMM_state$Posterior)) {
      print(sprintf('%s not in Posterior',term1))
      return(NULL)
    }
    samples = 1:dim(MegaLMM_state$Posterior[[term1]])[1]
  }
  base_env = with(MegaLMM_state,c(data_matrices,priors,Posterior[MegaLMM_state$Posterior$posteriorMean_params]))
  base_env = c(base_env,extra_env)
  Posterior = MegaLMM_state$Posterior
  per_sample_fun = function(sample_index_i) {
    # get current sample of each of the terms in FUN
    current_sample = make_current_state(Posterior,sample_index_i,terms)
    # evaluate FUN in an environment constructed from current_sample, and MegaLMM_state, taking current_sample first
    env = c(current_sample,base_env)
    result = eval(FUN,envir = env)
    if(is(result,'Matrix')) result = as.matrix(result)
    result
  }
  sample_1_result <- tryCatch(per_sample_fun(1),
                              error = function(e) {
                                message(e)
                                return(NULL)
                              })
  if(is.null(sample_1_result)) return(NULL)
  dim_1 = dim(sample_1_result) # get the dimension of the returned value
  # calculate value for each sample

  if(mc.cores == 1) {
    res = do.call('c',lapply(samples,per_sample_fun))
  } else{
    cluster = parallel::makeCluster(1)
    doParallel::registerDoParallel(cluster)
    res = foreach::foreach(i=1,.export = c('make_current_state'),.combine = 'c') %dopar% {
      # This way, we make a new process, and only pass it the data needed for per_sample_fun
      do.call('c',parallel::mclapply(samples,per_sample_fun,mc.cores = mc.cores))
    }
    parallel::stopCluster(cluster)
    # res = do.call(c,mclapply(samples,per_sample_fun,mc.cores = mc.cores))
  }
  # re-formulate into an appropriate array with the first dimension as samples
  if(is.null(dim_1)) {
    dim_1 = length(sample_1_result)
    res = matrix(res,ncol = dim_1,byrow=T)
    colnames(res) = names(sample_1_result)
  } else {
    res = aperm(array(res,dim = c(dim_1,length(samples))),c(3,1,2))
    dimnames(res)[2:3] = dimnames(sample_1_result)
  }
  res
}


#' Calculates posterior mean of a function of parameters
#'
#' The synthetic parameter can be pre-calculated using \link{get_posterior_FUN}, or provided directly
#'     as an array with the first dimension the samples.
#'
#' @param X either an array of posterior samples (either a parameter from \code{Posterior}, or
#'     an object generated by \link{get_posterior_fun}), or the MegaLMM_state object with re-loaded Posterior
#' @param FUN (optional) if \code{X} is a MegaLMM_state object, the function to calculate the synthetic parameter
#' @param bychunk (optional) if \code{TRUE}, will re-load each chunk of the Posterior separately
#'     and calculate the posterior mean of \code{FUN} separately for each chunk, and then return
#'     the overall posterior mean. This saves memory because the whole posterior does not need to be
#'     loaded into memory at once. Only necessary terms from Posterior are loaded.
#'
#' @return posterior mean matrix
get_posterior_mean = function(X,FUN,bychunk = FALSE,mc.cores = 1,...){
  result = NULL
  if(!bychunk) {
    if(is(X,'MegaLMM_state')) {
      MegaLMM_state = X
      FUN = match.call()$FUN
      X = do.call(get_posterior_FUN,list(MegaLMM_state=MegaLMM_state,FUN=FUN,mc.cores=mc.cores))
    }
    if(length(dim(X)) == 3) {
      result = matrix(colMeans(matrix(X,nr = dim(X)[1])),nr = dim(X)[2])
      dimnames(result) = dimnames(X)[-1]
    }
    if(length(dim(X)) == 2) {
      result = colMeans(X)
      names(result) = colnames(X)
    }
  } else{
    if(!is(X,'MegaLMM_state')) stop('Provide a MegaLMM_state object as "X"')
    MegaLMM_state = X
    FUN = match.call()$FUN
    if(is(FUN,'character')){
      FUN = parse(text=FUN)
    }
    terms = all.vars(FUN)
    terms = terms[terms %in% MegaLMM_state$Posterior$posteriorSample_params]
    n_files = length(MegaLMM_state$Posterior$files)
    result = 0
    chunk=1
#
    pb = txtProgressBar(min=0,max = n_files,style=3)
    for(chunk in 1:n_files){
      for(term in terms){
        MegaLMM_state$Posterior[[term]] = load_posterior_param(MegaLMM_state,term,chunks = chunk)
      }
      if(length(MegaLMM_state$Posterior[[term]]) > 0) {
        samples = do.call(get_posterior_FUN,list(MegaLMM_state=MegaLMM_state,FUN=FUN,mc.cores=mc.cores))
        result = result + dim(samples)[1]*get_posterior_mean(samples,bychunk = FALSE)
        rm(samples)
        gc()
      }
      setTxtProgressBar(pb, chunk)
    }
    close(pb)
    result = result / MegaLMM_state$Posterior$total_samples
  }
  result
}



#' Calculates Highest Posteriod Density intervals of a function of parameters
#'
#' The synthetic parameter can be pre-calculated using \link{get_posterior_FUN}, or provided directly
#'     as an array with the first dimension the samples.#'
#' @param X either an array of posterior samples (either a parameter from \code{Posterior}, or
#'     an object generated by \link{get_posterior_FUN}), or the MegaLMM_state object with re-loaded Posterior
#' @param FUN (optional) if \code{X} is a MegaLMM_state object, the function to calculate the synthetic parameter
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the intervals.
#'     The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob
#' @param ... other parameters passed to \link{get_posterior_FUN}
#'
#' @return array with first dimension of length=2 giving the lower and upper limits of the parameters
#'     in the other two dimensions.
get_posterior_HPDinterval = function(X,FUN = NULL,prob = 0.95,...){
  if(is(X,'MegaLMM_state')) {
    MegaLMM_state = X
    FUN = match.call()[[3]]
    X = do.call(get_posterior_FUN,list(MegaLMM_state=MegaLMM_state,FUN=FUN))
  }
  dims = dim(X)
  if(length(dims) == 3) result = aperm(array(apply(X,3,function(x) HPDinterval(mcmc(x),prob=prob)),c(dims[2],2,dims[3])),c(2,1,3))
  if(length(dims) == 2) result = t(HPDinterval(mcmc(X),prob=prob))
  result
}


#' Converts dense Matrix to "matrix" faster!
#'
#' @param X a ddenseMatrix
#'
#' @return a matrix
#' @export
#'
toDense = function(X) {
  matrix(X@x,nrow(X))
}
