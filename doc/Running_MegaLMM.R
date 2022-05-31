## ------------------------------------------------------------------------
library(MegaLMM)

## ------------------------------------------------------------------------
seed = 1  # for reproducibility
nSire = 50 # simulation design is a half-sib design: each of nSire fathers has nRep children, each with a different female
nRep = 10
nTraits = 100
nFixedEffects = 2 # fixed effects are effects on the factors
nFactors = 10
factor_h2s = c(rep(0,nFactors/2),rep(0.9,nFactors/2))  # factor_h2 is the % of variance in each factor trait that is explained by additive genetic variation.
# The is after accounting for the fixed effects on the factors
Va = 2 # residual genetic variance in each of the observed traits after accounting for the factors
Ve = 2 # residual microenvironmental variance in each of the observed traits after accounting for the factors
Vb = 0 # magnitude of the fixed effects (just factors)
new_halfSib_simulation('Sim_FE_1', nSire=nSire,nRep=nRep,p=nTraits, b=nFixedEffects, factor_h2s= factor_h2s,Va = Va, Ve = Ve,Vb = Vb)

## ------------------------------------------------------------------------
load('setup.RData')
Y = setup$Y
data = setup$data
K = setup$K

## ------------------------------------------------------------------------
run_parameters = MegaLMM_control(
  max_NA_groups = 3,
  scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 00,  # number of burn in samples before saving posterior samples
  K = 15 # number of factors
)

## ------------------------------------------------------------------------
priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      # Prior variance of trait residuals after accounting for fixed effects and factors
  tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
  Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe, # function that implements the horseshoe-based Lambda prior described in Runcie et al 2020. See code to see requirements for this function.
    prop_0 = 0.1,    # prior guess at the number of non-zero loadings in the first and most important factor
    delta = list(shape = 3, scale = 1),    # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
    delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  ),
  h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
)

## ------------------------------------------------------------------------
MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                              ~Fixed1 + (1|animal),  # RHS of base model for factors and residuals. Fixed effects defined here only apply to the factor residuals.
                              data = data,         # the data.frame with information for constructing the model matrices
                              relmat = list(animal = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                              run_parameters=run_parameters,
                              run_ID = 'MegaLMM_example'
                                )
maps = make_Missing_data_map(MegaLMM_state)
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)

MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)  # apply the priors
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state) # initialize the model
MegaLMM_state = initialize_MegaLMM(MegaLMM_state) # run the initial calculations
MegaLMM_state = clear_Posterior(MegaLMM_state) # prepare the output directories


## ----results="hide"------------------------------------------------------
# The following code is optional, but tries to guess for your system how many CPUs to use for fastest processing
(n_threads = optimize_n_threads(MegaLMM_state,seq(1,RcppParallel::defaultNumThreads(),by=1),times=2))
set_MegaLMM_nthreads(n_threads$optim)
# now do sampling is smallish chunks
n_iter = 100;  # how many samples to collect at once?
for(i  in 1:5) {
  print(sprintf('Run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_samples iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
  
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
  print(MegaLMM_state) # print status of current chain
  plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf

  # set of commands to run during burn-in period to help chain converge
  if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i < 3) {
    MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    print(MegaLMM_state$run_parameters$burn)
  }
}

## ------------------------------------------------------------------------
# Because this was a simulation, we can make some special diagnostic plots by passing in the true values
plot(MegaLMM_state,setup = setup) 

## ------------------------------------------------------------------------
# reload the whole database of posterior samples
MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)

# all parameter names in Posterior
MegaLMM_state$Posterior$posteriorSample_params
MegaLMM_state$Posterior$posteriorMean_params  # these ones only have the posterior mean saved, not individual posterior samples

# instead, load only a specific parameter
Lambda = load_posterior_param(MegaLMM_state,'Lambda')

# boxplots are good ways to visualize Posterior distributions on sets of related parameters
boxplot(MegaLMM_state$Posterior$F_h2[,1,])

# get posterior distribution on a function of parameters
# This is how to calculate the G-matrix for random effect #1 (ie animal above.)
G_samples = get_posterior_FUN(MegaLMM_state,t(Lambda) %*% diag(F_h2['animal',]) %*% Lambda + diag(resid_h2['animal',]/tot_Eta_prec[1,]))

# get posterior mean of a parameter
G = get_posterior_mean(G_samples)

# get Highest Posterior Density intervals for paramters
F_h2_HPD = get_posterior_HPDinterval(MegaLMM_state,F_h2)

# make a boxplot to summarize a subset of parameters.
boxplot(MegaLMM_state$Posterior$B1[,1,],outline=F);abline(h=0)

