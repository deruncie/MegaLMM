# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


sample_QTL_prec_ARD = function(MegaLMM_state,...){
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  current_state  = MegaLMM_state$current_state

  current_state = with(c(priors,run_variables),
                       with(QTL_prior,{
                         if(b_QTL + b_QTL_F == 0) return(current_state)

                         tau_shape = with(global, nu - 1)
                         tau_rate = with(global, V * nu)

                         within(current_state,{

                           B_QTL2 = B_QTL^2 * tot_Eta_prec[rep(1,b_QTL),]
                           B_QTL_F2 = B_QTL_F^2 * tot_F_prec[rep(1,b_QTL_F),]

                           if(is.null(separate_QTL_shrinkage) || separate_QTL_shrinkage == FALSE){
                             print('not separate_QTL_shrinkage')
                             if(b_QTL + b_QTL_F > 0){
                               if(!exists('B_QTL_tau')){
                                 B_QTL_tau = matrix(1,1,b_QTL)
                                 B_QTL_F_tau = matrix(1,1,b_QTL_F)
                                 B_QTL_prec = matrix(1,b_QTL,p)
                                 B_QTL_F_prec = matrix(1,b_QTL_F,k)
                               }

                               tau = rgamma(1,
                                            shape = tau_shape + length(B_QTL2)/2 + length(B_QTL_F2)/2,
                                            rate = tau_rate +
                                              sum(B_QTL2 * B_QTL_prec/c(B_QTL_tau))/2 +
                                              sum(B_QTL_F2 * B_QTL_F_prec/c(B_QTL_F_tau))/2)
                               B_QTL_tau[] = tau
                               B_QTL_F_tau[] = tau

                               B_QTL_prec[] = matrix(rgamma(b_QTL*p,shape = (QTL_df + 1)/2,rate = (QTL_df + B_QTL2*c(B_QTL_tau))/2),nr = b_QTL,nc = p)
                               B_QTL_F_prec[] = matrix(rgamma(b_QTL_F*k,shape = (QTL_df + 1)/2,rate = (QTL_df + B_QTL_F2*t(B_QTL_F_tau)[,rep(1,k)])/2),nr = b_QTL_F,nc = k)
                               B_QTL_prec[] = B_QTL_prec * c(B_QTL_tau)
                               B_QTL_F_prec[] = B_QTL_F_prec * t(B_QTL_F_tau)[,rep(1,k)]
                             }
                           } else {
                             if(b_QTL > 0){
                               if(!exists('B_QTL_tau')){
                                 B_QTL_tau = matrix(1,1,p)
                                 B_QTL_prec = matrix(1,b_QTL,p)
                               }
                               tau = rgamma(p,
                                            shape = tau_shape + dim(B_QTL2)[1]/2,
                                            rate = tau_rate +
                                              colSums(B_QTL2 * B_QTL_prec/B_QTL_tau[rep(1,b_QTL),])/2)
                               B_QTL_tau[] = tau

                               B_QTL_prec[] = matrix(rgamma(b_QTL*p,shape = (QTL_df + 1)/2,rate = (QTL_df + B_QTL2*B_QTL_tau[rep(1,b_QTL),])/2),nr = b_QTL,nc = p)
                               B_QTL_prec[] = B_QTL_prec * B_QTL_tau[rep(1,b_QTL),]
                             }
                             if(b_QTL_F > 0){
                               if(!exists('B_QTL_F_tau')){
                                 B_QTL_F_tau = matrix(1,1,k)
                                 B_QTL_F_prec = matrix(1,b_QTL_F,k)
                               }
                               tau = rgamma(k,
                                            shape = tau_shape + dim(B_QTL_F2)[1]/2,
                                            rate = tau_rate +
                                              colSums(B_QTL_F2 * B_QTL_F_prec/B_QTL_F_tau[rep(1,b_QTL_F),])/2)
                               B_QTL_F_tau[] = tau

                               B_QTL_F_prec[] = matrix(rgamma(b_QTL_F*k,shape = (QTL_df + 1)/2,rate = (QTL_df + B_QTL_F2*B_QTL_F_tau[rep(1,b_QTL_F),])/2),nr = b_QTL_F,nc = k)
                               B_QTL_F_prec[] = B_QTL_F_prec * B_QTL_F_tau[rep(1,b_QTL_F),]
                             }
                           }
                       })
                     }))
  return(current_state)
}


sample_QTL_prec_horseshoe = function(MegaLMM_state,...){
  # using algorithm from Makalic and Schmidt (2015)
  priors         = MegaLMM_state$priors
  run_variables  = MegaLMM_state$run_variables
  current_state  = MegaLMM_state$current_state

  current_state = with(c(priors,run_variables),
                       with(QTL_prior,{
                         if(b_QTL + b_QTL_F == 0) return(current_state)

                         within(current_state,{

                           if(b_QTL > 0)   B_QTL2 = B_QTL^2# * tot_Eta_prec[rep(1,b_QTL),]
                           if(b_QTL_F > 0) B_QTL_F2 = B_QTL_F^2# * tot_F_prec[rep(1,b_QTL_F),]

                           if(!exists('cauchy_iteractions_factor')) cauchy_iteractions_factor = 1

                           for(rep in 1:cauchy_iteractions_factor) {

                             if(is.null(separate_QTL_shrinkage) || separate_QTL_shrinkage == FALSE){
                               print('not separate_QTL_shrinkage')
                               if(b_QTL + b_QTL_F > 0){
                                 if(!exists('B_QTL_tau2')){
                                   B_QTL_xi = matrix(1,1,1)
                                   B_QTL_tau2 = matrix(1,1,1)
                                   B_QTL_nu = matrix(1,b_QTL,p)
                                   B_QTL_prec = matrix(1,b_QTL,p)
                                   B_QTL_F_tau2 = matrix(1,1,1)
                                   B_QTL_F_nu = matrix(1,b_QTL_F,k)
                                   B_QTL_F_prec = matrix(1,b_QTL_F,k)
                                 }

                                 tau2 = 1/rgamma(1,
                                              shape = (length(B_QTL2) + length(B_QTL_F2) + 1) / 2,
                                              rate = 1/B_QTL_xi +
                                                sum(B_QTL2 * B_QTL_prec*c(B_QTL_tau2))/2 +
                                                sum(B_QTL_F2 * B_QTL_F_prec*c(B_QTL_F_tau2))/2
                                              )
                                 B_QTL_xi[]  = 1/rgamma(1,shape = 1,rate = 1+1/tau2)

                                 B_QTL_tau2[] = tau2  # check this!
                                 B_QTL_F_tau2[] = tau2

                                 B_QTL_prec[] = matrix(rgamma(b_QTL*p,shape = 1,rate = 1/B_QTL_nu + B_QTL2/(2*c(B_QTL_tau2))),nr = b_QTL,nc = p)
                                 B_QTL_F_prec[] = matrix(rgamma(b_QTL_F*k,shape = 1,rate = 1/B_QTL_F_nu + B_QTL_F2/(2*c(B_QTL_F_tau2))),nr = b_QTL_F,nc = k)

                                 B_QTL_nu[] = matrix(1/rgamma(b_QTL*p,shape = 1,rate = 1 + B_QTL_prec[]),nr = b_QTL,nc = p)
                                 B_QTL_F_nu[] = matrix(1/rgamma(b_QTL_F*k,shape = 1,rate = 1 + B_QTL_F_prec[]),nr = b_QTL_F,nc = k)

                                 B_QTL_prec[] = B_QTL_prec / c(B_QTL_tau2)
                                 B_QTL_F_prec[] = B_QTL_F_prec / c(B_QTL_F_tau2)
                               }
                             } else {
                               if(b_QTL > 0){
                                 if(!exists('B_QTL_tau2')){
                                   B_QTL_xi = matrix(1,1,p)
                                   B_QTL_tau2 = matrix(1,1,p)
                                   B_QTL_nu = matrix(1,b_QTL,p)
                                   B_QTL_prec = matrix(1,b_QTL,p)
                                 }
                                 if(any(B_QTL != 0)) {
                                   B_QTL_tau2[] = 1/rgamma(p,
                                                          shape = (b_QTL+1)/2,
                                                          rate = 1/B_QTL_xi +
                                                            colSums(B_QTL2 * B_QTL_prec*B_QTL_tau2[rep(1,b_QTL),])/2)

                                   B_QTL_xi[]  = 1/rgamma(p,shape = 1,rate = 1+1/B_QTL_tau2)

                                   B_QTL_prec[] = matrix(rgamma(b_QTL*p,shape = 1,rate = 1/B_QTL_nu + B_QTL2 / (2*B_QTL_tau2)[rep(1,b_QTL),]),nr = b_QTL,nc = p)
                                   B_QTL_nu[] = matrix(1/rgamma(b_QTL*p,shape = 1,rate = 1 + B_QTL_prec),nr = b_QTL,nc = p)

                                   B_QTL_prec[] = B_QTL_prec# / B_QTL_tau2[rep(1,b_QTL),]
                                 }
                               }
                               if(b_QTL_F > 0){
                                 if(!exists('B_QTL_F_tau2')){
                                   tau_0 = p0/(1-p0)/sqrt(nrow(F))
                                   B_QTL_F_xi = matrix(1,1,k)
                                   B_QTL_F_tau2 = matrix(1,1,k)
                                   B_QTL_F_nu = matrix(1,b_QTL_F,k)
                                   B_QTL_F_prec = matrix(1,b_QTL_F,k)
                                 }
                                 if(any(B_QTL_F != 0)) {
                                   B_QTL_F_tau2[] = 1/rgamma(k,
                                                          shape = (b_QTL_F+1)/2,
                                                          rate = 1/B_QTL_F_xi +
                                                            colSums(B_QTL_F2 * B_QTL_F_prec * B_QTL_F_tau2[rep(1,b_QTL_F),])/2)
                                   B_QTL_F_xi[]  = 1/rgamma(k,shape = 1,rate = 1+1/B_QTL_F_tau2)

                                   B_QTL_F_prec[] = matrix(rgamma(b_QTL_F*k,shape = 1,rate = 1/B_QTL_F_nu + B_QTL_F2 / (2*B_QTL_F_tau2 * tau_0^2)[rep(1,b_QTL_F),]),nr = b_QTL_F,nc = k)
                                   B_QTL_F_nu[] = matrix(1/rgamma(b_QTL_F*k,shape = 1,rate = 1 + B_QTL_F_prec),nr = b_QTL_F,nc = k)

                                   B_QTL_F_prec[] = B_QTL_F_prec / B_QTL_F_tau2[rep(1,b_QTL_F),] / tau_0^2
                                 }
                               }
                             }
                           }
                         })
                       }))
  return(current_state)
}

