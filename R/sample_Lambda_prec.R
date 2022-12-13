# note: Kr is number of non-fixed columns of Lambda. We specify the prior precision of the final Kr columns of Lambda = Lambda[,!fixed_factors]

sample_Lambda_prec_horseshoe = function(MegaLMM_state,...) {
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
                       with(Lambda_prior,{
                         if(which_sampler$Y == 4) stop("Y sampler must be 1-3 for the horseshoe Lambda prior")
                         
                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100

                         delta_shape = delta$shape
                         delta_scale  = delta$scale

                         tau_0 = prop_0/(1-prop_0) * 1/sqrt(n)

                         within(current_state,{
                           if(Kr == 0) {
                             Lambda_prec = matrix(0,nrow=Kr,ncol=p)
                             delta = matrix(1,1,1)
                             return()
                           }

                           # initialize variables if needed
                           if(!'Lambda_tau2' %in% names(current_state)){
                             if(verbose) print('initializing Lambda_prec horseshoe')
                             Lambda_mean = matrix(0,Kr,p)
                             Lambda_tau2 = matrix(1,1,1)
                             Lambda_xi = matrix(1,1,1)
                             Lambda_phi2 = matrix(1,Kr,p)
                             Lambda_nu = matrix(1,Kr,p)
                             delta = with(priors,matrix(c(1,1/rgamma(Kr-1,shape = delta_shape,rate = 1/delta_scale)),ncol=1))
                             Lambda_prec = matrix(1,Kr,p)
                             trunc_point_delta = 1
                             Lambda_m_eff = matrix(1,Kr,1)
                           }

                           Lambda2 = (Lambda - Lambda_mean)[!fixed_factors,,drop=FALSE]^2
                           Lambda2_std = sweep(Lambda2,2,tot_Eta_prec[1,]/2,'*')

                           Lambda_nu[] = matrix(1/rgamma(Kr*p,shape = 1, rate = 1 + 1/Lambda_phi2), nr = Kr, nc = p)
                           # Lambda2_std_delta = sweep(Lambda2_std,1, cumprod(delta),'*')  # with delta~Ga
                           Lambda2_std_delta = sweep(Lambda2_std,1, cumprod(delta),'/') # with delta~iG
                           Lambda_phi2[] = matrix(1/rgamma(Kr*p,shape = 1, rate = 1/Lambda_nu + Lambda2_std_delta / Lambda_tau2[1]),nr=Kr,nc = p)

                           scores = rowSums(Lambda2_std / Lambda_phi2)
                           # for(i in 1:delta_iterations_factor) {
                           #   cumprod_delta = cumprod(delta[,1])
                           #   Lambda_tau2[] = 1/rgamma(1,shape = (Kr*p+1)/2, rate = 1/Lambda_xi[1] + sum(cumprod_delta*scores))
                           #   Lambda_xi[] = 1/rgamma(1,shape = 1,rate = 1/tau_0^2 + 1/Lambda_tau2[1])
                           #   for(h in 2:Kr) {
                           #     delta[h] = rgamma(1,shape = delta_shape + p/2*(Kr-h+1),rate = delta_l_rate + sum(cumprod_delta[h:Kr]*scores[h:Kr])/(Lambda_tau2[1]*delta[h]))
                           #     cumprod_delta = cumprod(delta[1,])
                           #   }
                           # }
                           new_samples = sample_tau2_delta_c_Eigen_v2(Lambda_tau2[1],Lambda_xi[1],delta,scores,
                                                                      tau_0,delta_shape,delta_scale,
                                                                      p,delta_iterations_factor)

                           Lambda_tau2[] = new_samples$tau2
                           Lambda_xi[] = new_samples$xi
                           delta[,1] = new_samples$delta

                           # -----Update Plam-------------------- #
                           # Lambda_prec[] = 1/(Lambda_tau2[1] * sweep(Lambda_phi2,1,cumprod(delta),'/')) # with delta~Ga
                           Lambda_prec[] = 1/(Lambda_tau2[1] * sweep(Lambda_phi2,1,cumprod(delta),'*'))  # with delta~iG

                           # ----- Calculate m_eff -------------- #
                           kappa = 1/(1+n/Lambda_prec)
                           Lambda_m_eff[] = rowSums(1-kappa)

                           rm(list = c('Lambda2','Lambda2_std','Lambda2_std_delta','scores','new_samples','kappa'))
                         })
                       }))
  return(current_state)
}


sample_Lambda_prec_ARD = function(MegaLMM_state,...) {
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  run_parameters = MegaLMM_state$run_parameters
  current_state  = MegaLMM_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(Lambda_prior,{
                         if(which_sampler$Y == 4) stop("Y sampler must be 1-3 for the ARD Lambda prior")

                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100
                         
                         delta_1_shape = delta_1$shape
                         delta_1_rate  = delta_1$rate
                         delta_2_shape = delta_2$shape
                         delta_2_rate  = delta_2$rate

                         within(current_state,{
                           if(Kr == 0) {
                             Lambda_prec = matrix(Kr,p)
                             return()
                           }

                           # initialize variables if needed
                           if(!exists('delta')){
                             if(!exists('X',Lambda_prior)) X = matrix(0,0,0)
                             if(!exists('X_group',Lambda_prior)) X_group = rep(1,ncol(X))
                             Lambda_beta = matrix(0,Kr,ncol(X))
                             Lambda_beta_var = matrix(1,Kr,length(unique(X_group)))
                             delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(Kr-1,shape = delta_2_shape,rate = delta_2_rate)),ncol=1))
                             delta[] = pmax(1,delta)
                             tauh  = matrix(cumprod(delta),ncol=1)
                             Lambda_phi = Lambda_prec = matrix(1,Kr,p)
                             trunc_point_delta = 1
                           }
                           
                           tauh = cumprod(delta)
                           
                           # sample Lambda_mean
                           if(ncol(X) > 0 && (!exists('fit_X',Lambda_prior) || fit_X)) {
                             Lambda_beta = t(do.call(cbind,lapply(seq_along(which(!fixed_factors)),function(k) {
                               prec_e = tot_Eta_prec[1,]*Lambda_phi[k,]*tauh[k]
                               X_std = sqrt(prec_e)*X
                               lambda_std = sqrt(prec_e)*Lambda[k,]
                               C = crossprod(X_std)
                               diag(C) = diag(C) + 1/Lambda_beta_var[k,X_group]*tauh[k]
                               chol_C = chol(C)
                               backsolve(chol_C,forwardsolve(t(chol_C),crossprod(X_std,lambda_std)) + rnorm(ncol(X)))
                             })))
                             Lambda_beta_var = do.call(cbind,
                                                       lapply(unique(X_group),function(i) {
                                                         beta_index = X_group==i
                                                         1/rgamma(Kr,
                                                                  shape = Lambda_beta_var_shape + 0.5*sum(beta_index),
                                                                  rate = Lambda_beta_var_rate + 0.5*rowSums(Lambda_beta[,beta_index,drop=FALSE]^2)*tauh) #
                                                       }))
                             Lambda_beta2_std = Lambda_beta^2 / Lambda_beta_var[,X_group,drop=FALSE]
                             Lambda_mean = Lambda_beta %*% t(X)
                           }
                           
                           Lambda2 = (Lambda - Lambda_mean)[!fixed_factors,,drop=FALSE]^2
                           Lambda2_std = sweep(Lambda2,2,tot_Eta_prec[1,],'*') #/ 2

                           Lambda_phi[] = matrix(rgamma(Kr*p,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2_std,1,tauh,'*'))/2),nr = Kr,nc = p)

                           # # -----Sample delta, update tauh------ #
                           scores = 0.5*rowSums(Lambda2_std*Lambda_phi)
                           shapes = c(delta_1_shape + 0.5*p*Kr,
                                      delta_2_shape + 0.5*p*((Kr-1):1))
                           if(ncol(X) > 0 && (!exists('fit_X',Lambda_prior) || fit_X)) {
                             scores = scores + 0.5 * rowSums(Lambda_beta2_std)
                             shapes = shapes + c(0.5*ncol(X)*Kr,
                                                 0.5*ncol(X)*((Kr-1):1))
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
  return(current_state)
}

sample_Lambda_prec_BayesC = function(MegaLMM_state,...) {
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
                       with(Lambda_prior,{
                         if(!which_sampler$Y == 4) stop("Y sampler must be 4 for the BayesC Lambda prior")
                         
                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100
                         if(!exists('fixed_pi')) fixed_pi = NULL
                         
                         
                         delta_1_shape = delta_1$shape
                         delta_1_rate  = delta_1$rate
                         delta_2_shape = delta_2$shape
                         delta_2_rate  = delta_2$rate
                         
                         
                         within(current_state,{
                           if(Kr == 0) {
                             Lambda_prec = matrix(0,nrow=Kr,ncol=p)
                             delta = matrix(1,1,1)
                             return()
                           }
                           
                           # initialize variables if needed
                           if(!exists('delta')){
                             if(verbose) print('initializing Lambda_prec BayesC')
                             Lambda_mean = matrix(0,Kr,p)
                             Lambda_prec = matrix(1,Kr,p)
                             delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(Kr-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                             tauh  = matrix(cumprod(delta),nrow=1)
                             Lambda_pi = matrix(1,Kr,1)
                             Lambda_delta = matrix(1,Kr,p)
                             Lambda_beta = matrix(1,Kr,p)
                             trunc_point_delta = 0
                             varEffects = 1
                           } else{
                             # if(sum(Lambda != 0) == 0) recover()
                             nLoci = rowSums(Lambda_delta)
                             if(is.null(fixed_pi)) {
                               Lambda_pi = matrix(rbeta(Kr,p-nLoci+1,nLoci+1),nrow = Kr,ncol = 1)
                             } else{
                               Lambda_pi = matrix(fixed_pi,nrow = Kr, ncol = 1)
                             }
                             # print(cbind(Lambda_pi, 1-nLoci/p,rowSums(Lambda == 0)/p))
                             # print(rowSums(Lambda==0)/p)
                             
                             Lambda2 = (Lambda - Lambda_mean)[!fixed_factors,,drop=FALSE]^2
                             # varEffects = matrix((rowSums(sweep(Lambda2,1,tauh,'*')) + Lambda_df*Lambda_scale)/rchisq(Kr,nLoci + Lambda_df),nrow = Kr, ncol = p)
                             varEffects = (sum(sweep(Lambda2,1,tauh,'*')) + Lambda_df*Lambda_scale)/rchisq(1,sum(nLoci) + Lambda_df)
                             
                             # # -----Sample delta, update tauh------ #
                             scores = 0.5*rowSums(Lambda2 / varEffects)
                             # shapes = c(delta_1_shape + 0.5*p*Kr,
                             #            delta_2_shape + 0.5*p*((Kr-1):1))
                             shapes = c(delta_1_shape + 0.5*sum(Lambda_delta),
                                        delta_2_shape + 0.5*(sum(Lambda_delta)-cumsum(rowSums(Lambda_delta)))[-Kr])  # nLoci in all higher-order rows of Lambda
                             times = delta_iterations_factor
                             # randg_draws = matrix(rgamma(times*Kr,shape = shapes,rate = 1),nr=times,byrow=T)
                             # delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                             randu_draws = matrix(runif(times*Kr),nr=times)
                             delta[] = sample_trunc_delta_c_Eigen( delta,tauh,scores,shapes,delta_1_rate,delta_2_rate,randu_draws,trunc_point_delta)
                             tauh[]  = matrix(cumprod(delta),nrow=1)
                             rm(list = c('Lambda2','scores'))
                           }
                           
                           Lambda_prec[] = t(tauh)[,rep(1,p),drop=FALSE]/varEffects
                           # sweep(1/varEffects,1,tauh,'*')
                           
                           
                           # # ----- Calcualte m_eff -------------- #
                           # kappa = 1/(1+n/(Lambda_prec*Lambda_delta))
                           # Lambda_m_eff[] = rowSums(1-kappa)
                           
                         })
                       }))
  return(current_state)
}

sample_Lambda_prec_BayesC_onePi = function(MegaLMM_state,...) {
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
                       with(Lambda_prior,{
                         if(!which_sampler$Y == 4) stop("Y sampler must be 4 for the BayesC Lambda prior")
                         
                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100
                         if(!exists('fixed_pi')) fixed_pi = NULL
                         
                         
                         delta_1_shape = delta_1$shape
                         delta_1_rate  = delta_1$rate
                         delta_2_shape = delta_2$shape
                         delta_2_rate  = delta_2$rate
                         
                         
                         within(current_state,{
                           if(Kr == 0) {
                             Lambda_prec = matrix(0,nrow=Kr,ncol=p)
                             delta = matrix(1,1,1)
                             return()
                           }
                           
                           # initialize variables if needed
                           if(!exists('delta')){
                             if(verbose) print('initializing Lambda_prec BayesC')
                             Lambda_mean = matrix(0,Kr,p)
                             Lambda_prec = matrix(1,Kr,p)
                             delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(Kr-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                             tauh  = matrix(cumprod(delta),nrow=1)
                             Lambda_pi = matrix(1,Kr,1)
                             Lambda_delta = matrix(1,Kr,p)
                             Lambda_beta = matrix(1,Kr,p)
                             trunc_point_delta = 0
                             varEffects = 1
                           } else{
                             # if(sum(Lambda != 0) == 0) recover()
                             # nLoci = rowSums(Lambda_delta)
                             nLoci = sum(Lambda_delta)
                             if(is.null(fixed_pi)) {
                               # Lambda_pi = matrix(rbeta(Kr,p-nLoci+1,nLoci+1),nrow = Kr,ncol = 1)
                               Lambda_pi = matrix(rbeta(1,Kr*p-nLoci+1,nLoci+1),nrow = Kr,ncol = 1)
                             } else{
                               Lambda_pi = matrix(fixed_pi,nrow = Kr, ncol = 1)
                             }
                             # print(cbind(Lambda_pi, 1-nLoci/p,rowSums(Lambda == 0)/p))
                             # print(rowSums(Lambda==0)/p)
                             
                             Lambda2 = (Lambda - Lambda_mean)[!fixed_factors,,drop=FALSE]^2
                             # varEffects = matrix((rowSums(sweep(Lambda2,1,tauh,'*')) + Lambda_df*Lambda_scale)/rchisq(Kr,nLoci + Lambda_df),nrow = Kr, ncol = p)
                             varEffects = (sum(sweep(Lambda2,1,tauh,'*')) + Lambda_df*Lambda_scale)/rchisq(1,sum(nLoci) + Lambda_df)
                             
                             # # -----Sample delta, update tauh------ #
                             scores = 0.5*rowSums(Lambda2 / varEffects)
                             # shapes = c(delta_1_shape + 0.5*p*Kr,
                             #            delta_2_shape + 0.5*p*((Kr-1):1))
                             shapes = c(delta_1_shape + 0.5*sum(Lambda_delta),
                                        delta_2_shape + 0.5*(sum(Lambda_delta)-cumsum(rowSums(Lambda_delta)))[-Kr])  # nLoci in all higher-order rows of Lambda
                             times = delta_iterations_factor
                             # randg_draws = matrix(rgamma(times*Kr,shape = shapes,rate = 1),nr=times,byrow=T)
                             # delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                             randu_draws = matrix(runif(times*Kr),nr=times)
                             delta[] = sample_trunc_delta_c_Eigen( delta,tauh,scores,shapes,delta_1_rate,delta_2_rate,randu_draws,trunc_point_delta)
                             tauh[]  = matrix(cumprod(delta),nrow=1)
                             rm(list = c('Lambda2','scores'))
                           }
                           
                           Lambda_prec[] = t(tauh)[,rep(1,p),drop=FALSE]/varEffects
                           # sweep(1/varEffects,1,tauh,'*')
                           
                           
                           # # ----- Calcualte m_eff -------------- #
                           # kappa = 1/(1+n/(Lambda_prec*Lambda_delta))
                           # Lambda_m_eff[] = rowSums(1-kappa)
                           
                         })
                       }))
  return(current_state)
}

