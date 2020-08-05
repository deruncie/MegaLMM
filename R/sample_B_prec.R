# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


sample_B2_prec_horseshoe = function(MegaLMM_state,...) {
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
                         if(b2_R>0 & which_sampler$Y == 4) stop("Y sampler must be 1-3 for the horseshoe B2_R prior")
                         if(b2_F>0 & which_sampler$F == 4) stop("Y sampler must be 1-3 for the horseshoe B2_F prior")

                         tau_0 = prop_0/(1-prop_0) * 1/sqrt(n)

                         within(current_state,{

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

                             rm(list=c('B2_R_2','B2_R_2_std','B2_F_2','B2_F_2_std'))
                           }



                         })
                       }))
  return(current_state)
}


sample_B2_prec_BayesC = function(MegaLMM_state,...) {
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
                         
                         within(current_state,{
                           
                           if(!any(c('B2_R_pi','B2_F_pi') %in% names(current_state))){
                             if(verbose) print('initializing B_prec BayesC')
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
                               B2_R_pi = matrix(rbeta(p,b2_R-nLoci+1,nLoci+1),nrow = 1,ncol = p,byrow = TRUE)
                               B2_R_prec[] = 1/matrix((colSums(B2_R^2) + B2_R_df*B2_R_scale)/rchisq(p,nLoci+B2_R_df),nrow = b2_R,ncol = p,byrow=TRUE)
                             }
                             
                             if(b2_F > 0) {
                               nLoci = colSums(B2_F_delta)
                               B2_F_pi = matrix(rbeta(K,b2_F-nLoci+1,nLoci+1),nrow = 1,ncol = K,byrow = TRUE)
                               B2_F_prec[] = 1/matrix((colSums(B2_F^2) + B2_R_df*B2_R_scale)/rchisq(K,nLoci+B2_R_df),nrow = b2_F,ncol = K,byrow=TRUE)
                             }
                           
                             rm(list=c('nLoci'))
                           }
                           
                           
                           
                         })
                       }))
  return(current_state)
}


# old code currently non-functional!
# sample_B_prec_RE = function(MegaLMM_state,...){
#   # treats B as a random effect - no parameter-specific shrinkage
#   # error if ncol(B_F)>0
#   priors         = MegaLMM_state$priors
#   run_variables  = MegaLMM_state$run_variables
#   current_state  = MegaLMM_state$current_state
# 
#   current_state = with(c(priors,run_variables),
#                        with(B_prior,{
#                          #   list(
#                          #   # load priors
#                          #   B_df   = B_prior$B_df,
#                          #   B_F_df = B_prior$B_F_df,
#                          #   B_QTL_df   = B_prior$B_QTL_df,
#                          #   B_QTL_F_df = B_prior$B_QTL_F_df,
#                          #   separate_QTL_shrinkage = B_prior$separate_QTL_shrinkage
#                          # ),
#                          if(b > 0) {
#                            prec_shape = with(global, nu - 1)
#                            prec_rate = with(global, V * nu)
#                          }
#                          if(b_F > 0) {
#                            stop("sample_B_prec_RE doesn't work with B_F")
#                          }
#                          within(current_state,{
# 
#                            # initialize variables if needed
#                            if(!exists('B_prec')){
#                              if(b > 0) {
#                                B_prec = matrix(rgamma(b,shape = prec_shape,rate=prec_rate),nrow = b, ncol = p)
#                              } else{
#                                B_prec = matrix(0,nrow=0,ncol=p)
#                              }
#                            }
#                            B2 = B^2
# 
#                            B_prec[] = matrix(rgamma(b,shape = prec_shape + p/2,rate = prec_rate + rowSums(B2)/2),nrow = b, ncol = p)
# 
#                            if(length(resid_intercept) > 0){
#                              B_prec[1,resid_intercept] = 1e-10
#                            }
#                          })
#                        }))
#   return(current_state)
# }
# 
# 
# 
# sample_B_prec_ARD = function(MegaLMM_state,...){
#   priors         = MegaLMM_state$priors
#   run_variables  = MegaLMM_state$run_variables
#   current_state  = MegaLMM_state$current_state
# 
#   current_state = with(c(priors,run_variables),
#                        with(B_prior,{
#                        #   list(
#                        #   # load priors
#                        #   B_df   = B_prior$B_df,
#                        #   B_F_df = B_prior$B_F_df,
#                        #   B_QTL_df   = B_prior$B_QTL_df,
#                        #   B_QTL_F_df = B_prior$B_QTL_F_df,
#                        #   separate_QTL_shrinkage = B_prior$separate_QTL_shrinkage
#                        # ),
#                               if(b > 0) {
#                                 tau_shape = with(global, nu - 1)
#                                 tau_rate = with(global, V * nu)
#                               }
#                               if(b_F > 0) {
#                                 tau_F_shape = with(global_F, nu - 1)
#                                 tau_F_rate = with(global_F, V * nu)
#                               }
#                               within(current_state,{
# 
#                                  # initialize variables if needed
#                                  if(!exists('B_tau')){
#                                    if(b > 0) {
#                                      B_tau = matrix(c(1e-10,rgamma(b-1,shape = tau_shape, rate = tau_rate)),nrow=1)
#                                    } else{
#                                      B_tau = matrix(0,nrow=1,ncol=0)
#                                    }
#                                    if(b_F > 0) {
#                                      B_F_tau = matrix(rgamma(b_F,shape = tau_F_shape, rate = tau_F_rate),nrow=1)
#                                    } else{
#                                      B_F_tau = matrix(0,nrow=1,ncol=0)
#                                    }
# 
#                                    B_prec = matrix(B_tau,nrow = b, ncol = p)
#                                    B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
#                                  }
#                                  B2 = B^2
#                                  B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec
# 
#                                  if(ncol(fixed_effects_common)>0){
#                                    ib = fixed_effects_common[1,]
#                                    ib_F = fixed_effects_common[2,]
#                                    tau = rgamma(length(ib),
#                                                 shape = tau_shape + ncol(B2)/2 + ncol(B_F2)/2,
#                                                 rate = tau_rate +
#                                                   rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2 +
#                                                   rowSums(B_F2[ib_F,,drop=FALSE] * B_F_prec[ib_F,,drop=FALSE]/c(B_F_tau)[ib_F])/2)
#                                    B_tau[1,ib] = tau
#                                    B_F_tau[1,ib_F] = tau
#                                  }
#                                  if(length(fixed_effects_only_resid) > 0){
#                                    ib = fixed_effects_only_resid
#                                    tau = rgamma(length(ib),
#                                                 shape = tau_shape + ncol(B2)/2,
#                                                 rate = tau_rate +
#                                                   rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2)
#                                    B_tau[1,ib] = tau
# 
#                                  }
#                                  if(length(fixed_effects_only_factors) > 0){
#                                    ib = fixed_effects_only_factors
#                                    tau = rgamma(length(ib),
#                                                 shape = tau_shape + ncol(B_F2)/2,
#                                                 rate = tau_rate +
#                                                   rowSums(B_F2[ib,,drop=FALSE] * B_F_prec[ib,,drop=FALSE]/c(B_F_tau)[ib])/2)
#                                    B_F_tau[1,ib] = tau
# 
#                                  }
#                                  B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
#                                  B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*t(B_F_tau)[,rep(1,k)])/2),nr = b_F,nc = k)
#                                  B_prec[] = B_prec*c(B_tau)
#                                  if(resid_intercept){
#                                    B_tau[1,1] = 1e-10
#                                    B_prec[1,] = 1e-10
#                                  }
#                                  B_F_prec[] = B_F_prec*t(B_F_tau)[,rep(1,k)]
#                                  B_F_tau[1,X_F_zero_variance] = 1e10
#                                  B_F_prec[X_F_zero_variance,] = 1e10
#                                })
#                          }))
#   return(current_state)
# }
# 
