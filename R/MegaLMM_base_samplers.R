#' Run MegaLMM Gibbs sampler
#'
#' Run MCMC chain for a specified number of iterations
#'
#' @param MegaLMM_state MegaLMM_state object of current chain
#' @param n_iter Number of iterations to add to the chain (not number of posterior samples to draw.
#'     This is determined by n_iter / thin)
sample_MegaLMM = function(MegaLMM_state,n_iter,grainSize = 1,verbose=TRUE,...) {
  data_matrices  = MegaLMM_state$data_matrices
  priors         = MegaLMM_state$priors
  run_parameters = MegaLMM_state$run_parameters
  run_variables  = MegaLMM_state$run_variables

  # ----------------------------------------------- #
  # -----------Reset Global Random Number Stream--- #
  # ----------------------------------------------- #
  RNG = MegaLMM_state$current_state$RNG
  do.call("RNGkind",as.list(RNG$RNGkind))  ## must be first!
  assign(".Random.seed", RNG$Random.seed, .GlobalEnv)

  # ----------------------------------------------- #
  # ----------------Set up run--------------------- #
  # ----------------------------------------------- #
  save_freq    = run_parameters$save_freq
  burn         = run_parameters$burn
  thin         = run_parameters$thin
  start_i      = MegaLMM_state$current_state$nrun

  # ----------------------------------------------- #
  # ---Extend posterior matrices for new samples--- #
  # ----------------------------------------------- #

  sp = (start_i + n_iter - burn)/thin - MegaLMM_state$Posterior$total_samples
  MegaLMM_state$Posterior = expand_Posterior(MegaLMM_state$Posterior,max(0,sp))

  # ----------------------------------------------- #
  # --------------start gibbs sampling------------- #
  # ----------------------------------------------- #

  if(verbose) pb = txtProgressBar(min=start_i,max = start_i+n_iter,style=3)
  start_time = Sys.time()
  for(i in start_i+(1:n_iter)){
    MegaLMM_state$current_state$nrun = i
    MegaLMM_state$current_state = MegaLMM_state$current_state[!sapply(MegaLMM_state$current_state,is.null)]

    # ----- Sample model parameters  except precisions ---------------- #
    MegaLMM_state$current_state = sample_latent_traits(MegaLMM_state,...)

    # -----Sample Lambda_prec ------------- #
    MegaLMM_state$current_state = MegaLMM_state$priors$Lambda_prior$sampler(MegaLMM_state,...)

    # -----Sample B2_prec ------------- #
    MegaLMM_state$current_state = MegaLMM_state$priors$B2_prior$sampler(MegaLMM_state,...)

    # ----- sample Eta ----- #
    observation_model_state = run_parameters$observation_model(run_parameters$observation_model_parameters,MegaLMM_state)$state
    MegaLMM_state$current_state[names(observation_model_state)] = observation_model_state

    # -- adapt number of factors to samples ---#
    # if(i > 200 && i < burn && runif(1) < with(MegaLMM_state$run_parameters,1/exp(b0 + b1*i))){  # adapt with decreasing probability per iteration
    #   MegaLMM_state$current_state = update_k(MegaLMM_state)
    # }

    # -- save sampled values (after thinning) -- #
    if( (i-burn) %% thin == 0 && i > burn) {
      MegaLMM_state$Posterior = save_posterior_sample(MegaLMM_state)
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  end_time = Sys.time()
  if(verbose) close(pb)
  print(end_time - start_time)
  MegaLMM_state$current_state$total_time = MegaLMM_state$current_state$total_time + end_time - start_time

  # ----------------------------------------------- #
  # ------------Save state for restart------------- #
  # ----------------------------------------------- #
  MegaLMM_state$current_state$RNG = list(
      Random.seed = .Random.seed,
      RNGkind = RNGkind()
    )

  current_state = MegaLMM_state$current_state
  if(run_parameters$save_current_state) saveRDS(current_state,file=sprintf('%s/current_state.rds',MegaLMM_state$run_ID))

  return(MegaLMM_state)
}
