
partial_matrix_model_setup = function(Y, center = T,scale = T) {
  Y_mean = rep(0,ncol(Y))
  Y_scale = rep(1,ncol(Y))
  Eta = Y
  if(center) {
    Y_mean = colMeans(Y,na.rm=T)
    Eta = sweep(Y,2,Y_mean,'-')
  }
  if(scale) {
    var_Eta = apply(Y,2,var,na.rm=T)
    Eta = sweep(Eta,2,sqrt(var_Eta),'/')
  }
  Y_missing = as(as(as(is.na(Y), "lMatrix"), "generalMatrix"), "TsparseMatrix")
  n_missing = sum(Y_missing)
  return(list(Eta = Eta,parameters = list(Y_mean = Y_mean,var_Eta = var_Eta,Y_missing = Y_missing,n_missing=n_missing)))
}
partial_matrix_model_sampler = function(trait_sets_current_state,parameters) {
  new_state = with(c(trait_sets_current_state,parameters),{
    # recover()
    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
    resids = rnorm(n_missing,0,sqrt(1/resid_Eta_prec[Y_missing@j+1]))  # sample resids from normal distribution with appropriate variance
    Eta[which(Y_missing)] = Eta_mean[which(Y_missing)] + resids
    list(Eta = Eta)
  })
  trait_sets_current_state$Eta = new_state$Eta
  trait_sets_current_state
}



Lambda_prec_ARD_sampler = function(MegaLMM_state,...) {
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  data_matrices = MegaLMM_state$data_matrices
  run_parameters = MegaLMM_state$run_parameters
  current_state  = MegaLMM_state$current_state
  
  current_state$Lambda = assemble_matrix(current_state,matrix='Lambda')
  
  current_state = with(c(priors,run_variables,run_parameters,data_matrices),
                       with(Lambda_prior,{
                         if(which_sampler$Y == 4) stop("Y sampler must be 1-3 for the ARD Lambda prior")
                         
                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100
                         
                         delta_1_shape = delta_1$shape
                         delta_1_rate  = delta_1$rate
                         delta_2_shape = delta_2$shape
                         delta_2_rate  = delta_2$rate
                         
                         within(current_state,{
                           
                           # initialize variables if needed
                           if(!exists('delta')){
                             if(!exists('X_group',Lambda_prior)) X_group = rep(1,ncol(Lambda_X))
                             Lambda_beta = matrix(0,K,ncol(Lambda_X))
                             Lambda_beta_var = matrix(1,K,length(unique(X_group)))
                             delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),ncol=1))
                             delta[] = pmax(1,delta)
                             tauh  = matrix(cumprod(delta),ncol=1)
                             Lambda_phi = Lambda_prec = matrix(1,K,p)
                             trunc_point_delta = 1
                           }
                           
                           tauh = cumprod(delta)
                           
                           # sample Lambda_mean
                           if(ncol(Lambda_X) > 0 && (!exists('fit_X',Lambda_prior) || fit_X)) {
                             Lambda_beta = t(do.call(cbind,lapply(1:K,function(k) {
                               prec_e = tot_Eta_prec[1,]*Lambda_phi[k,]*tauh[k]
                               Lambda_X_std = sqrt(prec_e)*Lambda_X
                               lambda_std = sqrt(prec_e)*Lambda[k,]#*sqrt(var_Eta) # allow for variance scaling of Y
                               C = crossprod(Lambda_X_std)
                               diag(C) = diag(C) + 1/Lambda_beta_var[k,X_group]*tauh[k]
                               chol_C = chol(C)
                               backsolve(chol_C,forwardsolve(t(chol_C),crossprod(Lambda_X_std,lambda_std)) + rnorm(ncol(X)))
                             })))
                             Lambda_beta_var = do.call(cbind,
                                                       lapply(unique(X_group),function(i) {
                                                         beta_index = X_group==i
                                                         1/rgamma(K,
                                                                  shape = Lambda_beta_var_shape + 0.5*sum(beta_index),
                                                                  rate = Lambda_beta_var_rate + 0.5*rowSums(Lambda_beta[,beta_index,drop=FALSE]^2)*tauh) #
                                                       }))
                             Lambda_beta2_std = Lambda_beta^2 / Lambda_beta_var[,X_group,drop=FALSE]
                             Lambda_mean = Lambda_beta %*% t(X)
                             # Lambda_mean = sweep(Lambda_mean,2,sqrt(var_Eta),'/') # allow for variance scaling of Y
                           }
                           
                           Lambda2 = (Lambda - Lambda_mean)^2
                           Lambda2_std = sweep(Lambda2,2,tot_Eta_prec[1,],'*') #/ 2
                           
                           Lambda_phi[] = matrix(rgamma(K*p,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2_std,1,tauh,'*'))/2),nr = K,nc = p)
                           
                           # # -----Sample delta, update tauh------ #
                           scores = 0.5*rowSums(Lambda2_std*Lambda_phi)
                           shapes = c(delta_1_shape + 0.5*p*K,
                                      delta_2_shape + 0.5*p*((K-1):1))
                           if(ncol(Lambda_X) > 0 && (!exists('fit_X',Lambda_prior) || fit_X)) {
                             scores = scores + 0.5 * rowSums(Lambda_beta2_std)
                             shapes = shapes + c(0.5*ncol(Lambda_X)*K,
                                                 0.5*ncol(Lambda_X)*((K-1):1))
                           }
                           times = delta_iterations_factor
                           # randg_draws = matrix(rgamma(times*Kr,shape = shapes,rate = 1),nr=times,byrow=T)
                           # delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                           randu_draws = matrix(runif(times*Kr),nr=times)
                           delta[,1] = sample_trunc_delta_c_Eigen( delta,tauh,scores,shapes,delta_1_rate,delta_2_rate,randu_draws,trunc_point_delta)
                           if(max(delta[,1]) > 1e30) recover()
                           # delta[,1] = d2
                           tauh[]  = matrix(cumprod(delta),ncol=1)
                           # print(c(delta))
                           
                           Lambda_prec[] = sweep(Lambda_phi,1,tauh,'*')
                         })
                       }))
  for(set in seq_along(current_state$trait_sets_current_state)) {
    cols = current_state$trait_sets_current_state[[set]]$cols
    current_state$trait_sets_current_state[[set]]$Lambda_mean = current_state$Lambda_mean[,cols]
    current_state$trait_sets_current_state[[set]]$Lambda_prec = current_state$Lambda_prec[,cols]
  }
  return(current_state)
}



B2_prec_horseshoe_sampler = function(MegaLMM_state,...) {
  # sampling as described in Supplemental methods, except we multiply columns of Prec_lambda by delta
  # phi2 = \lambda^2 in methods
  # the delta sequence controls tau_k. We have tau_1~C+(0,tau_0), and tau_k = tau_1*prod_{l=1}^K(delta^{-1}_l)
  # delta_l controls the decrease in odds of inclusion of each element for the lth factor relative to the (l-1)th
  # Note:  Piironen and Vehtari do not include sigma^2 in prior, so I have removed
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  run_parameters = MegaLMM_state$run_parameters
  current_state  = MegaLMM_state$current_state
  current_state$B2_R = assemble_matrix(current_state,matrices = 'B2_R')
  
  current_state = with(c(priors,run_variables,run_parameters),
                       with(B2_prior,{
                         if(b2_R>0 & which_sampler$Y == 4) stop("Y sampler must be 1-3 for the horseshoe B2_R prior")
                         if(b2_F>0 & which_sampler$F == 4) stop("Y sampler must be 1-3 for the horseshoe B2_F prior")
                         
                         tau_0 = prop_0/(1-prop_0) * 1/sqrt(n)
                         
                         within(current_state,{
                           
                           p = ncol(B2_R) # note: not all traits might use B2_R.
                           
                           if(!any(c('B2_R_tau2','B2_F_tau2') %in% names(current_state))){
                             if(verbose) print('initializing B_prec horseshoe')
                             B2_R_xi = matrix(1/rgamma(p,shape=1/2,rate=1/tau_0^2),nr=1)
                             B2_R_tau2 = matrix(1/rgamma(p,shape = 1/2, rate = 1/B2_R_xi[1,]),nr=1)
                             B2_R_nu = matrix(1/rgamma(b2_R*p,shape = 1/2, rate = 1), nr = b2_R, nc = p)
                             B2_R_phi2 = matrix(1/rgamma(b2_R*p,shape = 1/2, rate = 1/B2_R_nu),nr=b2_R,nc = p)
                             B2_R_prec = 1 / sweep(B2_R_phi2,2,B2_R_tau2[1,],'*')
                             B2_R = rstdnorm_mat(b2_R,p)/sqrt(B2_R_prec)
                             
                             
                             B2_F_xi = matrix(1/rgamma(K,shape=1/2,rate=1/tau_0^2),nr=1)
                             B2_F_tau2 = matrix(1/rgamma(K,shape = 1/2, rate = 1/B2_F_xi[1,]),nr=1)
                             B2_F_nu = matrix(1/rgamma(b2_F*K,shape = 1/2, rate = 1), nr = b2_F, nc = K)
                             B2_F_phi2 = matrix(1/rgamma(b2_F*K,shape = 1/2, rate = 1/B2_F_nu),nr=b2_F,nc = K)
                             B2_F_prec = 1 / sweep(B2_F_phi2,2,B2_F_tau2[1,],'*')
                             B2_F = rstdnorm_mat(b2_F,K)/sqrt(B2_F_prec)
                           } else {
                             
                             B2_R_2 = B2_R^2
                             B2_R_2_std = sweep(B2_R_2,2,tot_Eta_prec[1,],'*')
                             B2_R_nu[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1 + 1/B2_R_phi2), nr = b2_R, nc = p)
                             B2_R_phi2[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1/B2_R_nu + sweep(B2_R_2_std,2,(2*B2_R_tau2),'/')),nr=b2_R,nc = p)
                             B2_R_xi[] = 1/rgamma(p,shape=1,rate=1/tau_0^2 + 1/B2_R_tau2[1,])
                             B2_R_tau2[] = 1/rgamma(p,shape = (b2_R + 1)/2, rate = 1/B2_R_xi[1,] + (colSums(B2_R_2_std/B2_R_phi2))/2)
                             
                             
                             B2_F_2 = B2_F^2
                             B2_F_2_std = sweep(B2_F_2,2,tot_F_prec[1,],'*')
                             B2_F_nu[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1 + 1/B2_F_phi2), nr = b2_F, nc = K)
                             B2_F_phi2[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1/B2_F_nu + sweep(B2_F_2_std,2,(2*B2_F_tau2),'/')),nr=b2_F,nc = K)
                             B2_F_xi[] = 1/rgamma(K,shape=1,rate=1/tau_0^2 + 1/B2_F_tau2[1,])
                             B2_F_tau2[] = 1/rgamma(K,shape = (b2_F + 1)/2, rate = 1/B2_F_xi[1,] + (colSums(B2_F_2_std/B2_F_phi2))/2)
                             
                             # -----Update Plam-------------------- #
                             B2_R_prec = 1 / sweep(B2_R_phi2,2,B2_R_tau2[1,],'*')
                             B2_F_prec = 1 / sweep(B2_F_phi2,2,B2_F_tau2[1,],'*')
                             
                             rm(list=c('B2_R_2','B2_R_2_std','B2_F_2','B2_F_2_std','p'))
                           }
                         })
                       }))
  for(set in seq_along(current_state$trait_sets_current_state)) {
    cols = match(colnames(current_state$B2_R),colnames(current_state$trait_sets_current_state[[set]]$Eta))
    if(length(cols) == 0) next
    current_state$trait_sets_current_state[[set]]$B2_R_prec = current_state$B2_R_prec[,cols]
    current_state$trait_sets_current_state[[set]]$Lambda_prec = current_state$Lambda_prec[,cols]
  }
  return(current_state)
}


B2_prec_BayesC_sampler = function(MegaLMM_state,...) {
  # sampling as described in Supplemental methods, except we multiply columns of Prec_lambda by delta
  # phi2 = \lambda^2 in methods
  # the delta sequence controls tau_k. We have tau_1~C+(0,tau_0), and tau_k = tau_1*prod_{l=1}^K(delta^{-1}_l)
  # delta_l controls the decrease in odds of inclusion of each element for the lth factor relative to the (l-1)th
  # Note:  Piironen and Vehtari do not include sigma^2 in prior, so I have removed
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  run_parameters = MegaLMM_state$run_parameters
  current_state  = MegaLMM_state$current_state
  
  current_state = with(c(priors,run_variables,run_parameters),
                       with(B2_prior,{
                         if(b2_R>0 & !which_sampler$Y == 4) stop("Y sampler must be 4 for the BayesC B2_R prior")
                         if(b2_F>0 & !which_sampler$F == 4) stop("Y sampler must be 4 for the BayesC B2_F prior")
                         if(!exists('fixed_pi')) fixed_pi = NULL
                         
                         within(current_state,{
                           
                           if(!any(c('B2_R_pi','B2_F_pi') %in% names(current_state))){
                             if(verbose) print('initializing B_prec BayesC')
                             B2_R_pi = matrix(1,1,p)
                             B2_R_delta = matrix(1,b2_R,p)
                             B2_R_beta = matrix(1,b2_R,p)
                             B2_R_prec = matrix(1,b2_R,p)
                             
                             B2_F_pi = matrix(1,1,K)
                             B2_F_delta = matrix(1,b2_F,K)
                             B2_F_beta = matrix(1,b2_F,K)
                             B2_F_prec = matrix(1,b2_F,K)
                           } else {    
                             if(b2_R > 0) {
                               nLoci = colSums(B2_R_delta)
                               if(is.null(fixed_pi)) {
                                 B2_R_pi = matrix(rbeta(p,b2_R-nLoci+1,nLoci+1),nrow = 1,ncol = p,byrow = TRUE)
                               } else{
                                 B2_R_pi = matrix(fixed_pi,nrow = 1,ncol = p)
                               }
                               B2_R_prec[] = 1/matrix((colSums(B2_R^2) + B2_R_df*B2_R_scale)/rchisq(p,nLoci+B2_R_df),nrow = b2_R,ncol = p,byrow=TRUE)
                             }
                             
                             if(b2_F > 0) {
                               nLoci = colSums(B2_F_delta)
                               if(is.null(fixed_pi)) {
                                 B2_F_pi = matrix(rbeta(K,b2_F-nLoci+1,nLoci+1),nrow = 1,ncol = K,byrow = TRUE)
                               } else{
                                 B2_F_pi = matrix(fixed_pi,nrow = 1,ncol = K)
                               }
                               B2_F_prec[] = 1/matrix((colSums(B2_F^2) + B2_F_df*B2_F_scale)/rchisq(K,nLoci+B2_F_df),nrow = b2_F,ncol = K,byrow=TRUE)
                             }
                             
                             rm(list=c('nLoci'))
                           }
                         })
                       }))
  return(current_state)
}


B2_prec_BayesB_sampler = function(MegaLMM_state,...) {
  # sampling as described in Supplemental methods, except we multiply columns of Prec_lambda by delta
  # phi2 = \lambda^2 in methods
  # the delta sequence controls tau_k. We have tau_1~C+(0,tau_0), and tau_k = tau_1*prod_{l=1}^K(delta^{-1}_l)
  # delta_l controls the decrease in odds of inclusion of each element for the lth factor relative to the (l-1)th
  # Note:  Piironen and Vehtari do not include sigma^2 in prior, so I have removed
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  run_parameters = MegaLMM_state$run_parameters
  current_state  = MegaLMM_state$current_state
  
  current_state = with(c(priors,run_variables,run_parameters),
                       with(B2_prior,{
                         if(b2_R>0 & !which_sampler$Y == 4) stop("Y sampler must be 4 for the BayesB B2_R prior")
                         if(b2_F>0 & !which_sampler$F == 4) stop("Y sampler must be 4 for the BayesB B2_F prior")
                         if(!exists('fixed_pi')) fixed_pi = NULL
                         
                         within(current_state,{
                           
                           if(!any(c('B2_R_pi','B2_F_pi') %in% names(current_state))){
                             if(verbose) print('initializing B_prec BayesB')
                             if(b2_R > 0) {
                               B2_R_pi = matrix(1,1,p)
                               B2_R_delta = matrix(1,b2_R,p)
                               B2_R_beta = matrix(1,b2_R,p)
                               B2_R_prec = matrix(1,b2_R,p)
                             }
                             if(b2_F > 0) {
                               B2_F_pi = matrix(1,1,K)
                               B2_F_delta = matrix(1,b2_F,K)
                               B2_F_beta = matrix(1,b2_F,K)
                               B2_F_prec = matrix(1,b2_F,K)
                             }
                           } else {    
                             if(b2_R > 0) {
                               nLoci = colSums(B2_R_delta)
                               if(is.null(fixed_pi)) {
                                 B2_R_pi = matrix(rbeta(p,b2_R-nLoci+1,nLoci+1),nrow = 1,ncol = p,byrow = TRUE)
                               } else{
                                 B2_R_pi = matrix(fixed_pi,nrow = 1,ncol = p)
                               }
                               B2_R_prec[] = 1/((B2_R_beta^2 + B2_R_df*B2_R_scale)/rchisq(b2_R*p,1+B2_R_df))
                             }
                             
                             if(b2_F > 0) {
                               nLoci = colSums(B2_F_delta)
                               if(is.null(fixed_pi)) {
                                 B2_F_pi = matrix(rbeta(K,b2_F-nLoci+1,nLoci+1),nrow = 1,ncol = K,byrow = TRUE)
                               } else{
                                 B2_F_pi = matrix(fixed_pi,nrow = 1,ncol = K)
                               }
                               B2_F_prec[] = 1/((B2_F_beta^2 + B2_F_df*B2_F_scale)/rchisq(b2_F*K,1+B2_F_df))
                             }
                             B2_F = delta*
                               rm(list=c('nLoci'))
                           }
                         })
                       }))
  return(current_state)
}
