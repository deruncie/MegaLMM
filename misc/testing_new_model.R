library(MegaLMM)
# 
# a = list()
# for(i in 1:20) {
#   a[[i]] = list(A = rstdnorm_mat(1000,20))
# }
# 
# a2 = do.call(cbind,lapply(a,function(x) x$A))
# a3 = rstdnorm_mat(nrow(a2),ncol(a2))
# cbind_list_withName(a3,a,'A')
# all(a2==a3)
# 
# library(microbenchmark)
# microbenchmark(do.call(cbind,lapply(a,function(x) x$A)),cbind_list_withName(a3,a,'A'),times=100)



formula = ~ 0 + (1|Line)
h2s_matrix = make_h2s_matrix(formula,line_data,20)

MegaLMM_state_new = MegaLMM_F(formula, sample_data, ID = "Line", relmat = list(Line=K),
  K=20, h2s_matrix = h2s_matrix,
  tot_F_var_prior = list(V = 18/20, nu = 20),
  F_h2_prior = matrix(1,1,20),
  output_folder = 'MegaLMM_model_new',
  X2 = NULL, U2 = NULL, V2 = NULL,
  B2_F_prior = list(
    sampler = sample_B2_prec_horseshoe,
    prop_0 = 0.1
  ))

data_partitions = partition_data_missing(Y_train,max_NA_groups = ncol(Y_train), starting_groups = NULL, verbose=FALSE)
data_partition = data_partitions$Missing_data_map_list[[3]]
data_partition[[5]] = data_partition[[2]]
data_partition[[2]]$Y_cols = data_partition[[2]]$Y_cols[1:3]
data_partition[[5]]$Y_cols = data_partition[[5]]$Y_cols[-c(1:3)]

for(i in 1:length(data_partition)) {
  if(length(data_partition[[i]]$Y_cols)==0) next
  MegaLMM_state_new = MegaLMM_add_trait_matrix(MegaLMM_state_new,Y_train[data_partition[[i]]$Y_obs,data_partition[[i]]$Y_cols,drop=FALSE],
                                                          fixed_formula = ~1,data = sample_data[data_partition[[i]]$Y_obs,,drop=FALSE],
                                                          tot_Y_var_prior = list(V = 0.5,   nu = 3),
                                                          resid_h2_prior = rep(1,20),
                                                          use_X2 = FALSE,
                                                          B2_R_prior = list(
                                                            sampler = sample_B2_prec_horseshoe,
                                                            prop_0 = 0.1
                                                          ),
                                                          center = T,scale = T) 
}

MegaLMM_state_new = MegaLMM_add_Lambda_model(MegaLMM_state_new,list(),
                                             Lambda_prior = list(
                                               sampler = Lambda_prec_ARD_sampler,
                                               Lambda_df = 3,
                                               delta_1 = list(shape = 20,rate=1),
                                               delta_2 = list(shape = 3,rate=1),
                                               delta_iterations_factor = 100
                                             ),
                                             Lambda_beta_prior = list(V=1,nu=3))

MegaLMM_state_new = MegaLMM_preliminary_calculations(MegaLMM_state_new,verbose = T)

MegaLMM_state_new = MegaLMM_initialize_variables(MegaLMM_state_new)


MegaLMM_state = setup_model_MegaLMM(
  Y = Y_train,  
  # The n x p trait matrix
  formula = ~ Population + ibc + (1|Line),  
  # This is syntax like lme4 for mixed effect models. 
  # We specify a fixed effect of population and a random effect for genotype (Line)
  data = sample_data,         
  # the data.frame with information for constructing the model matrices
  relmat = list(Line = K), 
  # A list of covariance matrices to link to the random effects in formula.
  # each grouping variable in formula can be linked to a covariance matrix.
  # If so, every level of the grouping variable must be in the rownames of K.
  # additional rows of K not present in data will still be predicted 
  # (and therefore will use memory and computational time!)
  run_parameters=run_parameters,
  # This list of control parameters created above
  run_ID = sprintf('MegaLMM_fold_%02d',fold_ID)
  # A run identifier. The function will create a folder with this name 
  # and store lots of useful data inside it
)