# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/


sample_latent_traits = function(MegaLMM_state,...) {
  data_matrices  = MegaLMM_state$data_matrices
  priors         = MegaLMM_state$priors
  run_parameters = MegaLMM_state$run_parameters
  run_variables  = MegaLMM_state$run_variables
  current_state  = MegaLMM_state$current_state

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    # k = nrow(Lambda)
    # recover()

    # ------------------------------------------------------------#
    # ----- Sample Lambda, Beta, cis_effects, tot_Eta_prec -------#
    # ------------------------------------------------------------#

    # conditioning on F and prior precisions
    # do sampling in groups of columns with same patterns of missing data

    n_coefs = b2_R + Kr
    prior_mean = matrix(0,n_coefs,p)
    if(b2_R > 0) {
      prior_prec = rbind(B2_R_prec,Lambda_prec)
    } else{ # b == 0
      prior_prec = Lambda_prec
    }
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      Y_set = Qt_list[[set]] %**% Eta[rows,cols,drop=FALSE]
      X_set = cbind(QtX2_R_list[[set]],Qt_list[[set]] %**% F[rows,!fixed_factors,drop=FALSE])
      if(length(Qt_cis_genotypes) == p) {
        Qt_cis_genotypes_set = Qt_cis_genotypes[cols]
      } else{
        Qt_cis_genotypes_set = list()
      }

      if(!which_sampler$Y == 4) {
        new_samples = regression_sampler_parallel(
          which_sampler = which_sampler$Y,
          Y = Y_set,
          X1_base = QtX1_list[[set]],
          X1_list_ = Qt_cis_genotypes_set,
          X2_ = X_set,
          Vx_ = NULL,
          h2s_index = resid_h2_index[cols],
          chol_V_list_ = chol_V_list_list[[set]],
          Y_prec = tot_Eta_prec[cols],
          Y_prec_a0 = tot_Eta_prec_shape[cols], Y_prec_b0 = tot_Eta_prec_rate[cols],
          prior_prec_alpha1 = matrix(0,nr = sum(QtX1_keepColumns_list[[set]]),ncol = length(cols)),
          prior_prec_alpha2 = rep(priors$cis_effects_prior$prec,sum(cis_effects_index %in% cols)),
          prior_mean_beta = prior_mean[,cols,drop=FALSE],
          prior_prec_beta = prior_prec[,cols,drop=FALSE],
          NULL,NULL,NULL
        )
        # extract samples
        
        # alpha1 -> un-shunk rows of B
        B1[QtX1_keepColumns_list[[set]],cols] = new_samples$alpha1
        B1[!QtX1_keepColumns_list[[set]],cols] = 0
        
        # alpha2 -> cis_effects
        if(!is.null(cis_effects_index)) {
          cis_effects[1,cols %in% cis_effects_index] = new_samples$alpha2
        }
        # beta -> B2_R and Lambda
        if(b2_R > 0) {
          B2_R[,cols] = new_samples$beta[1:b2_R,]
        }
        Lambda[!fixed_factors,cols] = new_samples$beta[b2_R + 1:Kr,]
        
        # Y_prec -> tot_Eta_prec
        if(any(is.na(new_samples$Y_prec))) break
        tot_Eta_prec[cols] = new_samples$Y_prec
        
      } else{
        current_alpha1s = B1[QtX1_keepColumns_list[[set]],cols,drop=FALSE]
        BayesAlphabet_parms = list(
          alpha = Lambda[!fixed_factors,cols,drop=FALSE],
          beta = Lambda_beta[,cols,drop=FALSE],
          pi = matrix(Lambda_pi,nrow = Kr,ncol = length(cols)),
          delta = Lambda_delta[,cols,drop=FALSE],
          run_sampler_times = run_sampler_times
        )
        if(b2_R > 0) {
          BayesAlphabet_parms$alpha = rbind(B2_R[,cols,drop=FALSE],Lambda[!fixed_factors,cols,drop=FALSE])
          BayesAlphabet_parms$beta = rbind(B2_R_beta[,cols,drop=FALSE],Lambda_beta[,cols,drop=FALSE])
          BayesAlphabet_parms$pi = rbind(matrix(B2_R_pi[cols],b2_R,length(cols),byrow=T),matrix(Lambda_pi,nrow = Kr,ncol = length(cols)))
          BayesAlphabet_parms$delta = rbind(B2_R_delta[,cols,drop=FALSE],Lambda_delta[,cols,drop=FALSE])
        }
        # recover()
        # asdf=0
        # while(sum(is.na(new_samples$Y_prec)) == 0) {
        # while(T) {
        # asdf = asdf+1
        new_samples = regression_sampler_parallel(
          which_sampler = which_sampler$Y,
          Y = Y_set,
          X1_base = QtX1_list[[set]],
          X1_list_ = Qt_cis_genotypes_set,
          X2_ = X_set,
          Vx_ = NULL,
          h2s_index = resid_h2_index[cols],
          chol_V_list_ = chol_V_list_list[[set]],
          Y_prec = tot_Eta_prec[cols],
          Y_prec_a0 = tot_Eta_prec_shape[cols], Y_prec_b0 = tot_Eta_prec_rate[cols],
          prior_prec_alpha1 = matrix(0,nr = sum(QtX1_keepColumns_list[[set]]),ncol = length(cols)),
          prior_prec_alpha2 = rep(priors$cis_effects_prior$prec,sum(cis_effects_index %in% cols)),
          prior_mean_beta = prior_mean[,cols,drop=FALSE],
          prior_prec_beta = prior_prec[,cols,drop=FALSE],
          current_alpha1s_ = current_alpha1s,
          current_alpha2s_ = NULL,
          BayesAlphabet_parms = BayesAlphabet_parms
        )
        if(sum(is.na(new_samples$alpha1)) + sum(is.na(new_samples$beta)) + sum(is.na(new_samples$Y_prec)) > 0) break
        # }
        # }
        # if(sum(new_samples$beta[b2_R + 1:Kr,]) == 0) recover()
        if(sum(is.na(new_samples$alpha1)) + sum(is.na(new_samples$beta)) + sum(is.na(new_samples$Y_prec)) > 0) recover()
        
        # extract samples
        
        # alpha1 -> un-shunk rows of B
        B1[QtX1_keepColumns_list[[set]],cols] = new_samples$alpha1
        B1[!QtX1_keepColumns_list[[set]],cols] = 0
        
        # alpha2 -> cis_effects
        if(!is.null(cis_effects_index)) {
          cis_effects[1,cols %in% cis_effects_index] = new_samples$alpha2
        }
        # beta -> B2_R and Lambda
        if(b2_R > 0) {
          B2_R[,cols] = new_samples$beta[1:b2_R,]
          B2_R_beta[,cols] = new_samples$beta[b2_R+Kr+1:b2_R,]
          B2_R_delta[,cols] = new_samples$beta[2*(b2_R+Kr)+1:b2_R,]
        }
        Lambda[!fixed_factors,cols] = new_samples$beta[b2_R + 1:Kr,]
        Lambda_beta[,cols] = new_samples$beta[b2_R+Kr+b2_R + 1:Kr,]
        Lambda_delta[,cols] = new_samples$beta[2*(b2_R+Kr)+b2_R + 1:Kr,]
        
        # Y_prec -> tot_Eta_prec
        if(any(is.na(new_samples$Y_prec))) break
        tot_Eta_prec[cols] = new_samples$Y_prec

      }
      
    }

    XB = X1 %**% B1 + X2_R %**% B2_R
    # remove cis effects
    if(length(cis_genotypes) > 0){
      for(j in 1:p){
        if(n_cis_effects[j] > 0){
          cis_X_j = cis_genotypes[[j]]
          XB[,j] = XB[,j] + cis_X_j %*% cis_effects[cis_effects_index == j]
        }
      }
    }
    Eta_tilde = Eta - XB - F %**% Lambda


    # ------------------------------------------------------------#
    # ----- random effects: resid_h2 and U_R ---------------------#
    # ------------------------------------------------------------#
    
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      # sample resid_h2_index

      QtEta_tilde_set = Qt_list[[set]] %**% Eta_tilde[rows,cols,drop=FALSE]

      if(!length(h2_priors_resids) == ncol(h2s_matrix)) stop('wrong length of h2_priors_resids')
      if(is.null(h2_step_size)) {
        log_ps = log_p_h2s(QtEta_tilde_set,
                           tot_Eta_prec[cols],
                           chol_V_list_list[[set]],
                           h2_priors_resids)
        resid_h2_index[cols] = sample_h2s(log_ps)
      } else{
        resid_h2_index[cols] = sample_h2s_discrete_MH_c(QtEta_tilde_set,
                                                        tot_Eta_prec[cols],
                                                        h2_priors_resids,
                                                        resid_h2_index[cols],
                                                        h2s_matrix,
                                                        chol_V_list_list[[set]],
                                                        h2_step_size)
      }
      resid_h2[,cols] = h2s_matrix[,resid_h2_index[cols],drop=FALSE]


      U_R[,cols] = sample_MME_ZKZts_c(Eta_tilde[rows,cols,drop=FALSE],
                                      ZL[rows,,drop=FALSE],
                                      tot_Eta_prec[cols,drop=FALSE],
                                      chol_ZtZ_Kinv_list_list[[set]],
                                      resid_h2[,cols,drop=FALSE],
                                      resid_h2_index[cols,drop=FALSE])
    }


    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample F variables ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    # only use rows in Missing_data_map[[1]]. These are all the rows with non-missing data
    rows = Missing_data_map[[1]]$Y_obs

    prior_mean = matrix(0,b2_F,K)

    Qt_F = Qt_list[[1]] %**% F[rows,,drop=FALSE]
    if(!which_sampler$F == 4) {
      new_samples = regression_sampler_parallel(
        which_sampler = which_sampler$F,
        Y = Qt_F,
        X1_base = matrix(0,length(rows),0),
        X1_list_ = NULL,
        X2_ = Qt1_U2_F,
        Vx_ = V2_F,
        h2s_index = F_h2_index,
        chol_V_list_ = chol_V_list_list[[1]],
        Y_prec = tot_F_prec,
        Y_prec_a0 = tot_F_prec_shape, Y_prec_b0 = tot_F_prec_rate,
        prior_prec_alpha1 = matrix(0,nr = 0,ncol = K),
        prior_prec_alpha2 = rep(0,0),
        prior_mean_beta = prior_mean,
        prior_prec_beta = B2_F_prec,
        NULL,NULL,NULL
      )
      
      # extract samples
      
      # alpha1 -> NULL
      # alpha2 -> NULL
      # beta -> B2_F
      XFBF = 0
      F_tilde = F
      Qt_F_tilde = Qt_F
      if(b2_F > 0) {
        B2_F = new_samples$beta
        XFBF = X2_F %**% B2_F
        F_tilde = F - XFBF
        if( b2_F > length(rows)) {
          Qt_F_tilde = Qt_list[[1]] %**% F_tilde[rows,,drop=FALSE]
        } else{
          Qt_F_tilde = Qt_F - Qt1_X2_F %**% B2_F
        }
      }
      
      # Y_prec -> tot_F_prec
      tot_F_prec[] = new_samples$Y_prec
      
    } else{
      current_alpha1s = matrix(0,0,0)
      BayesAlphabet_parms = list(
        alpha = B2_F,
        beta = B2_F_beta,
        pi = B2_F_pi[rep(1,b2_F),,drop=FALSE],
        delta = B2_F_delta,
        run_sampler_times = run_sampler_times
      )
      new_samples = regression_sampler_parallel(
        which_sampler = which_sampler$F,
        Y = Qt_F,
        X1_base = matrix(0,length(rows),0),
        X1_list_ = NULL,
        X2_ = Qt1_U2_F,
        Vx_ = V2_F,
        h2s_index = F_h2_index,
        chol_V_list_ = chol_V_list_list[[1]],
        Y_prec = tot_F_prec,
        Y_prec_a0 = tot_F_prec_shape, Y_prec_b0 = tot_F_prec_rate,
        prior_prec_alpha1 = matrix(0,nr = 0,ncol = K),
        prior_prec_alpha2 = rep(0,0),
        prior_mean_beta = prior_mean,
        prior_prec_beta = B2_F_prec,
        current_alpha1s_ = current_alpha1s,
        current_alpha2s_ = NULL,
        BayesAlphabet_parms = BayesAlphabet_parms
      )
      
      # extract samples
      
      # alpha1 -> NULL
      # alpha2 -> NULL
      # beta -> B2_F
      XFBF = 0
      F_tilde = F
      Qt_F_tilde = Qt_F
      if(b2_F > 0) {
        B2_F = new_samples$beta[1:b2_F,,drop=FALSE]
        B2_F_beta = new_samples$beta[b2_F+1:b2_F,,drop=FALSE]
        B2_F_delta = new_samples$beta[2*b2_F+1:b2_F,,drop=FALSE]
        XFBF = X2_F %**% B2_F
        F_tilde = F - XFBF
        if( b2_F > length(rows)) {
          Qt_F_tilde = Qt_list[[1]] %**% F_tilde[rows,,drop=FALSE]
        } else{
          Qt_F_tilde = Qt_F - Qt1_X2_F %**% B2_F
        }
      }
      
      # Y_prec -> tot_F_prec
      tot_F_prec[] = new_samples$Y_prec
    }


    if(!length(h2_priors_factors) == ncol(h2s_matrix)) stop('wrong length of h2_priors_factors')
    if(is.null(h2_step_size)) {
      log_ps = log_p_h2s(Qt_F_tilde,
                         tot_F_prec,
                         chol_V_list_list[[1]],
                         h2_priors_factors)
      F_h2_index = sample_h2s(log_ps)
    } else{
      F_h2_index = sample_h2s_discrete_MH_c(Qt_F_tilde,
                                          tot_F_prec,
                                          h2_priors_factors,
                                          F_h2_index,
                                          h2s_matrix,
                                          chol_V_list_list[[1]],
                                          h2_step_size)
    }
    F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_MME_ZKZts_c(F_tilde[rows,,drop=FALSE], ZL[rows,,drop=FALSE], tot_F_prec, chol_ZtZ_Kinv_list_list[[1]], F_h2, F_h2_index)

    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,Lambda, F_h2
    Eta_tilde = Eta - XB - ZL %**% U_R
    F_e_prec = tot_F_prec / (1-colSums(F_h2))
    prior_mean = ZL %**% U_F + XFBF

    for(set in seq_along(Missing_row_data_map)){
      cols = Missing_row_data_map[[set]]$Y_cols
      rows = Missing_row_data_map[[set]]$Y_obs
      if(length(rows) == 0) next
      if(length(cols) > 0) {
        F[rows,] = sample_factors_scores_c(Eta_tilde[rows,cols,drop=FALSE],
                                           prior_mean[rows,,drop=FALSE],
                                           Lambda[,cols,drop=FALSE],
                                           resid_Eta_prec[cols],
                                           F_e_prec)
      } else{
        F[rows,] = prior_mean[rows,,drop=FALSE] + sweep(rstdnorm_mat(length(rows),K),2,sqrt(F_e_prec[1,]),'/')
      }
    }

  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

