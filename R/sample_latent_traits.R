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
    if(b2_R > 0) {
      prior_mean = rbind(matrix(0,b2_R,p),Lambda_mean)
      prior_prec = rbind(B2_R_prec,Lambda_prec)
    } else{ # b == 0
      prior_mean = Lambda_mean
      prior_prec = Lambda_prec
    }

    ## 5/31/2022
    ## I think the following lines are needed and were missing
    Eta_resid = Eta
    if(any(fixed_factors)) Eta_resid = Eta_resid - F[,fixed_factors,drop=FALSE] %*% Lambda[fixed_factors,,drop=FALSE]

    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      Y_set = Qt_list[[set]] %**% Eta_resid[rows,cols,drop=FALSE]
      X_set = cbind(QtX2_R_list[[set]],Qt_list[[set]] %**% F[rows,!fixed_factors,drop=FALSE])
      if(length(Qt_cis_genotypes) == p) {
        Qt_cis_genotypes_set = Qt_cis_genotypes[cols]
      } else{
        Qt_cis_genotypes_set = list()
      }

      if(which_sampler$Y == 1) {
        new_samples = regression_sampler_parallel(
          Y_set,
          QtX1_list[[set]]$X1,
          Qt_cis_genotypes_set,
          X_set,
          NULL,
          resid_h2_index[cols],
          chol_V_list_list[[set]],
          tot_Eta_prec[cols],
          tot_Eta_prec_shape[cols], tot_Eta_prec_rate[cols],
          matrix(0,nr = sum(QtX1_list[[set]]$keepColumns),ncol = length(cols)),
          rep(priors$cis_effects_prior$prec,sum(cis_effects_index %in% cols)),
          prior_mean[,cols,drop=FALSE],
          prior_prec[,cols,drop=FALSE]
        )
        # extract samples
        
        # alpha1 -> un-shunk rows of B
        B1[QtX1_list[[set]]$keepColumns,cols] = new_samples$alpha1
        B1[!QtX1_list[[set]]$keepColumns,cols] = 0
        
        # alpha2 -> cis_effects
        if(!is.null(cis_effects_index)) {
          cis_effects[1,cols %in% cis_effects_index] = new_samples$alpha2
        }
        # beta -> B2_R and Lambda
        if(b2_R > 0) {
          B2_R[,cols] = new_samples$beta[1:b2_R,]
        }
        Lambda[!fixed_factors,cols] = new_samples$beta[b2_R + 1:Kr,]# + Lambda_mean[,cols]
        
        # Y_prec -> tot_Eta_prec
        if(any(is.na(new_samples$Y_prec))) break
        tot_Eta_prec[cols] = new_samples$Y_prec 
      } else if(which_sampler$Y == 2) {
        unique_h2s = unique(resid_h2_index[cols])
        for(h2_i in unique_h2s) {
          cols_i = cols[resid_h2_index[cols] == h2_i]
          new_samples = parallel_block_regression_sampler(
            Y_set[,resid_h2_index[cols] == h2_i,drop=FALSE],
            QtX1_list[[set]]$X1,
            X_set,
            NULL,
            chol_V_list_list[[set]][[resid_h2_index[cols_i][1]]],
            tot_Eta_prec[cols_i],
            tot_Eta_prec_shape[cols_i], tot_Eta_prec_rate[cols_i],
            matrix(0,nr = sum(QtX1_list[[set]]$keepColumns),ncol = sum(cols_i)),
            prior_mean[,cols_i,drop=FALSE],
            prior_prec[,cols_i,drop=FALSE]
          )
          # alpha1 -> un-shunk rows of B
          B1[QtX1_list[[set]]$keepColumns,cols_i] = new_samples$alpha
          B1[!QtX1_list[[set]]$keepColumns,cols_i] = 0
          
          # beta -> B2_R and Lambda
          if(b2_R > 0) {
            B2_R[,cols_i] = new_samples$beta[1:b2_R,]
          }
          Lambda[!fixed_factors,cols_i] = new_samples$beta[b2_R + 1:Kr,]
          
          # Y_prec -> tot_Eta_prec
          if(any(is.na(new_samples$Y_prec))) break
          tot_Eta_prec[cols_i] = new_samples$Y_prec 
        }
      } else if(which_sampler$Y == 3) {
        unique_h2s = unique(resid_h2_index[cols])
        for(h2_i in unique_h2s) {
          cols_i = cols[resid_h2_index[cols] == h2_i]
          new_samples = parallel_Single_regression_sampler(
            Y_set[,resid_h2_index[cols] == h2_i,drop=FALSE],
            QtX1_list[[set]]$X1,
            X_set,
            NULL,
            chol_V_list_list[[set]][[resid_h2_index[cols_i][1]]],
            tot_Eta_prec[cols_i],
            tot_Eta_prec_shape[cols_i], tot_Eta_prec_rate[cols_i],
            matrix(0,nr = sum(QtX1_list[[set]]$keepColumns),ncol = length(cols_i)),
            prior_mean[,cols_i,drop=FALSE],
            prior_prec[,cols_i,drop=FALSE],
            B1[QtX1_list[[set]]$keepColumns,cols_i,drop=FALSE],
            Lambda[!fixed_factors,cols_i,drop=FALSE],
            Lambda[!fixed_factors,cols_i,drop=FALSE],
            matrix(0,sum(!fixed_factors),length(cols_i)),
            matrix(1,sum(!fixed_factors),length(cols_i)),
            1
          )
          if(sum(is.na(unlist(new_samples)))>0) recover()

          # alpha1 -> un-shunk rows of B
          B1[QtX1_list[[set]]$keepColumns,cols_i] = new_samples$alpha
          B1[!QtX1_list[[set]]$keepColumns,cols_i] = 0

          # beta -> B2_R and Lambda
          if(b2_R > 0) {
            B2_R[,cols_i] = new_samples$betas_beta[1:b2_R,]
          }
          Lambda[!fixed_factors,cols_i] = new_samples$betas_beta[b2_R + 1:Kr,]

          # Y_prec -> tot_Eta_prec
          if(any(is.na(new_samples$Y_prec))) break
          tot_Eta_prec[cols_i] = new_samples$Y_prec
        }
      } else if(which_sampler$Y == 4){
        current_alpha1s = B1[QtX1_list[[set]]$keepColumns,cols,drop=FALSE]
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
        new_samples = SingleSite_regression_sampler_parallel(
          Y = Y_set,
          X1_base = QtX1_list[[set]]$X1,
          X1_list_ = Qt_cis_genotypes_set,
          X2_ = X_set,
          Vx_ = NULL,
          h2s_index = resid_h2_index[cols],
          chol_V_list_ = chol_V_list_list[[set]],
          Y_prec = tot_Eta_prec[cols],
          Y_prec_a0 = tot_Eta_prec_shape[cols], Y_prec_b0 = tot_Eta_prec_rate[cols],
          prior_prec_alpha1 = matrix(0,nr = sum(QtX1_list[[set]]$keepColumns),ncol = length(cols)),
          prior_prec_alpha2 = rep(priors$cis_effects_prior$prec,sum(cis_effects_index %in% cols)),
          prior_mean_beta = prior_mean[,cols,drop=FALSE],
          prior_prec_beta = prior_prec[,cols,drop=FALSE],
          current_alpha1s_ = current_alpha1s,
          current_alpha2s_ = NULL,
          BayesAlphabet_parms = BayesAlphabet_parms
        )
        # extract samples
        
        # alpha1 -> un-shunk rows of B
        B1[QtX1_list[[set]]$keepColumns,cols] = new_samples$alpha1
        B1[!QtX1_list[[set]]$keepColumns,cols] = 0
        
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
        
      } else {
        stop("which_sampler$Y > 4 not implemented")
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
    
    F_Lambda = F %**% Lambda
    Eta_tilde = Eta - XB - F_Lambda


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
                                      ZL_list[[set]][rows,,drop=FALSE],
                                      tot_Eta_prec[cols,drop=FALSE],
                                      chol_ZtZ_Kinv_list_list[[set]],
                                      resid_h2[,cols,drop=FALSE],
                                      resid_h2_index[cols,drop=FALSE])
      
      # prepare Eta_tilde for sampling F below
      Eta_tilde[,cols] = Eta_tilde[,cols] + F_Lambda[,cols] - ZL_list[[set]] %*% U_R[,cols,drop=FALSE] 
    }


    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample F variables ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    # only use rows in Missing_data_map[[1]]. These are all the rows with non-missing data
    rows = Missing_data_map[[1]]$Y_obs

    prior_mean = matrix(0,b2_F,K)

    Qt_F = Qt_list[[1]] %**% F[rows,,drop=FALSE]
    if(which_sampler$F %in% 1:3) {
      new_samples = regression_sampler_parallel(
        Qt_F,
        matrix(0,length(rows),0),
        NULL,
        Qt1_U2_F,
        V2_F,
        F_h2_index,
        chol_V_list_list[[1]],
        tot_F_prec,
        tot_F_prec_shape, tot_F_prec_rate,
        matrix(0,nr = 0,ncol = K),
        rep(0,0),
        prior_mean,
        B2_F_prec
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
    } else if(which_sampler$F == 4) {
      current_alpha1s = matrix(0,0,0)
      BayesAlphabet_parms = list(
        alpha = B2_F,
        beta = B2_F_beta,
        pi = B2_F_pi[rep(1,b2_F),,drop=FALSE],
        delta = B2_F_delta,
        run_sampler_times = run_sampler_times
      )
      new_samples = SingleSite_regression_sampler_parallel(
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
      
    } else {
      stop("which_sampler$F > 4 not implemented")
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
    
    U_F[] = sample_MME_ZKZts_c(F_tilde[rows,,drop=FALSE], ZL_list[[1]][rows,,drop=FALSE], tot_F_prec, chol_ZtZ_Kinv_list_list[[1]], F_h2, F_h2_index)

    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,Lambda, F_h2
    # Eta_tilde = Eta - XB - ZL_list[[1]] %**% U_R  # This was done above
    F_e_prec = tot_F_prec / (1-colSums(F_h2))
    prior_mean = ZL_list[[1]] %**% U_F + XFBF

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
    
    # optional - rescale F to unit variance. Not "right" because this is not a parameter. But might help?
    # if(exists('forceF_var1') & forceF_var1) {
    #   var_F = apply(F,2,var)
    #   # print(var_F)
    #   if(min(var_F) < 1e-10) recover()
    #   F = sweep(F,2,sqrt(var_F),'/')
    #   B2_F = sweep(B2_F,2,sqrt(var_F),'/')
    #   U_F = sweep(U_F,2,sqrt(var_F),'/')
    #   F_e_prec = F_e_prec / var_F
    #   # Lambda = sweep(Lambda,1,sqrt(var_F),'*')
    #   # Lambda_beta = sweep(Lambda_beta,1,sqrt(var_F),'*')
    #   # Lambda_beta_var = sweep(Lambda_beta_var,1,var_F,'*')
    #   # Lambda_mean = sweep(Lambda_mean,1,sqrt(var_F),'*')
    #   # delta = delta/var_F
    #   # tauh[] = cumprod(delta)
    # }

  }))
  current_state = current_state[current_state_names]
  # if(any(lapply(current_state,function(x) sum(is.na(x))>0))) recover()

  return(current_state)
}

