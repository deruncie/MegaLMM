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
      Lambda[!fixed_factors,cols] = new_samples$beta[b2_R + 1:Kr,]

      # Y_prec -> tot_Eta_prec
      if(any(is.na(new_samples$Y_prec))) break
      tot_Eta_prec[cols] = new_samples$Y_prec
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
                                      tot_Eta_prec[cols],
                                      chol_ZtZ_Kinv_list_list[[set]],
                                      resid_h2[,cols,drop=FALSE],
                                      resid_h2_index[cols])
    }


    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample F variables ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    # only use rows in Missing_data_map[[1]]. These are all the rows with non-missing data
    rows = Missing_data_map[[1]]$Y_obs

    prior_mean = matrix(0,b2_F,K)

    Qt_F = Qt_list[[1]] %**% F[rows,,drop=FALSE]
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

