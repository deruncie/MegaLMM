make_h2s_matrix = function(formula,data, h2_divisions) {
  # get the names of the random effects, make a table with rows as the random effects and columns as the grid
  model_setup = make_model_setup(formula,data,list())

  RE_setup = model_setup$RE_setup

  # add names to RE_setup if needed
  n_RE = length(RE_setup)
  for(i in 1:n_RE){
    if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
      names(RE_setup)[i] = paste0('RE.',i)
    }
  }
  RE_names = names(RE_setup)

  if(length(h2_divisions) < n_RE){
    if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
    h2_divisions = rep(h2_divisions,n_RE)
  }
  if(is.null(names(h2_divisions))) {
    names(h2_divisions) = RE_names
  }
  h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  colnames(h2s_matrix) = NULL

  return(h2s_matrix)
}

# 1. F model - specify formula (with ID), ID, full data with all individuals, K, var_F prior, F_h2 prior
# 2. Add traits - for each, specify formula (with ID), reduced data with ID, var_Y prior, resid_h2 prior. Could be different trait models for each set of traits
#    - check that no Lambda model yet. If so, clear it
#    - The IDs here define the missing data map for columns
#    - must have same random effects as F model, but fixed model can be different
# 3. Add Lambda model. Easiest would be a list of formulas all with same data, each with a prior, plus Lambda prior function and parameters
# 4. Initialize sampling matrices. May be able to combine across trait sets and F
# 5. Initialize variables
# 6. Define posterior functions

MegaLMM_F = function(formula, data, ID = "ID", relmat = NULL,
                     K, h2s_matrix,
                     tot_F_var_prior = list(V = 18/20, nu = 20),
                     F_h2_prior,
                     output_folder = 'MegaLMM_model',
                     X2 = NULL, U2 = NULL, V2 = NULL,
                     B2_F_prior = list(
                       sampler = sample_B2_prec_horseshoe,
                       prop_0 = 0.1
                     )) {
  try(dir.create(output_folder,recursive = T, showWarnings = FALSE),silent=T)

  # first check ID
  if(!ID %in% colnames(data)) stop(sprintf('ID: "%s% missing from data',ID))
  if(!length(data[[ID]]) == length(unique(data[[ID]]))) stop(sprintf('ID: "%s% not unique in data',ID))

  model_setup = make_model_setup(formula,data,relmat)

  # X2 is the "fixed effects" for F
  X2_F = model_setup$lmod$X
  b2_F = ncol(X2_F)
  n = nrow(data)

  # alternatively, (or in addition) these can be specified directly with X2 or U2*V2

  if(!is.null(X2)) {
    if(nrow(X2) != n) stop("nrow(X2) != nrow(data)")
    X2_F = cbind(X2_F,X2)
    b2_F = ncol(X2_F)
  }

  U2_F = NULL
  V2_F = NULL
  if(!is.null(U2) || !is.null(V2)) {
    if(is.null(U2) || is.null(V2)) stop("must specify both U2 and V2 if either is specified")
    if(nrow(U2) != n)  stop("nrow(U2) != nrow(data)")
    if(ncol(U2) != nrow(V2)) stop("sizes of U2 and V2 are not compatible. U2 %*% V2 is calculated")
    if(ncol(X2_F)>0) {
      # expand V2
      V2_F = as.matrix(Matrix::bdiag(diag(1,ncol(X2_F)),V2))
      X2_F = cbind(X2_F,U2)
    }
  }

  if((!is.null(X2_F) || !is.null(U2_F)) && is.null(B2_F_prior)) stop("B2_F_prior is needed")

  # -------- Random effects ---------- #

  # RE_setup is a list like from BGLR that describes the random effects (for both residuals and factors)
  RE_setup = model_setup$RE_setup

  # add names to RE_setup if needed
  n_RE = length(RE_setup)
  for(i in 1:n_RE){
    if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
      names(RE_setup)[i] = paste0('RE.',i)
    }
  }
  RE_names = names(RE_setup)

  # combine Z matrices
  Z = do.call(cbind,lapply(RE_setup,function(x) x$Z))
  Z = as(Z,'dgCMatrix')


  # find RE indices
  RE_lengths = sapply(RE_setup,function(x) ncol(x$Z))
  RE_starts = cumsum(c(0,RE_lengths)[1:n_RE])
  names(RE_starts) = RE_names
  RE_indices = lapply(RE_names,function(re) RE_starts[re] + 1:RE_lengths[re])
  names(RE_indices) = RE_names

  h2_divisions = run_parameters$h2_divisions
  if(length(h2_divisions) < n_RE){
    if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
    h2_divisions = rep(h2_divisions,n_RE)
  }
  if(is.null(names(h2_divisions))) {
    names(h2_divisions) = RE_names
  }
  h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
  colnames(h2s_matrix) = RE_names
  h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
  colnames(h2s_matrix) = NULL

  # # # -------------------------------------------#
  # # # identify groups of traits with same pattern of missingness
  # # # ideally, want to be able to restrict the number of sets. Should be possible to merge sets of columngs together.
  #
  # Missing_data_map = list(list(
  #   Y_obs = 1:n,
  #   Y_cols = 1:p
  # ))
  # Missing_row_data_map = list(list(
  #   Y_obs = 1:n,
  #   Y_cols = 1:p
  # ))


  # # # -------------------------------------------#
  # # # Fixed factors
  # if(is.null(Lambda_fixed)) {
  #   fixed_factors = rep(F,run_parameters$K)
  # } else{
  #   Lambda_fixed = as.matrix(Lambda_fixed)
  #   if(ncol(Lambda_fixed) != p) stop("wrong dimensions of Lambda_fixed")
  #   if(nrow(Lambda_fixed) >= run_parameters$K) stop("nrow(Lambda_fixed) >= K")
  #   fixed_factors = rep(F,run_parameters$K)
  #   fixed_factors[1:nrow(Lambda_fixed)] = T
  # }
  # Kr = run_parameters$K - sum(fixed_factors)

  priors = list()
  priors$tot_Y_var = list(V = c(),nu = c())
  priors$tot_F_var = tot_F_var_prior
  if(length(priors$tot_F_var$V) == 1) {
    priors$tot_F_var$V = rep(priors$tot_F_var$V,K)
    priors$tot_F_var$nu = rep(priors$tot_F_var$nu,K)
  }
  priors$h2_priors_resids = matrix(0,nrow = length(F_h2_prior),ncol = 0)
  priors$h2_priors_factors = F_h2_prior
  # recover()
  if(nrow(priors$h2_priors_factors) != n_RE || !all(rownames(F_h2_prior == RE_names))) stop("F_h2_prior not compatible with random effect formula")

  trait_sets = list(
    list(
      setID = 1,
      Y = NULL,
      rows = 1:nrow(data),
      cols = NULL,
      colnames = NULL,
      X1_cols = NULL,
      use_X2 = is.null(X2_F) && is.null(U2_F),
      Y_model = NULL
    )
  )

  run_variables = list(
    K = K,
    t = 0,
    n      = n,
    # r_RE   = r_RE,
    RE_names = RE_names,
    b2 = b2_F,
    F_formula = formula,
    F_data       = data,
    ID = ID,
    trait_sets = trait_sets
    # Kr = Kr
  )

  data_matrices = list(
    X1 = matrix(0,n,0),
    X2 = X2_F,
    U2 = U2_F,
    V2 = V2_F,
    Z           = Z,
    # ZL          = ZL,
    # ZL_list     = list(ZL),
    RE_setup    = RE_setup,
    # RE_L_list     = list(RE_L),  # List of matrices necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
    RE_indices  = RE_indices,
    h2s_matrix  = h2s_matrix
  )

  current_state = list(
    Eta = matrix(NA,nrow = nrow(data),ncol = 0,dimnames = list(data[[ID]]))
  )
  #
  # run_parameters$observation_model = observation_model
  # run_parameters$observation_model_parameters = observation_model_parameters
  # run_parameters$traitnames = traitnames


  # ----------------------------- #
  # -- create MegaLMM_state object - #
  # ----------------------------- #

  MegaLMM_state = list(
    current_state  = current_state,
    output_folder         = output_folder,
    data_matrices  = data_matrices,
    priors         = priors,
    run_parameters = run_parameters,
    run_variables  = run_variables
  )
  class(MegaLMM_state) = append('MegaLMM_state',class(MegaLMM_state))


  MegaLMM_state$Posterior = list(
    # posteriorSample_params = posteriorSample_params,
    # posteriorMean_params = posteriorMean_params,
    # posteriorFunctions = posteriorFunctions,
    total_samples = 0,
    folder = file.path(output_folder,'Posterior'),
    files = c()
  )

  MegaLMM_state
}

# 2. Add traits - for each, specify formula (with ID), reduced data with ID, var_Y prior, resid_h2 prior. Could be different trait models for each set of traits
#    - check that no Lambda model yet. If so, clear it
#    - The IDs here define the missing data map for columns
#    - must have same random effects as F model, but fixed model can be different

MegaLMM_add_trait_matrix = function(MegaLMM_state,Y,fixed_formula,data,
                                    tot_Y_var_prior = list(V = 0.5,   nu = 3),
                                    resid_h2_prior,
                                    use_X2 = FALSE,
                                    B2_R_prior = list(
                                      sampler = sample_B2_prec_horseshoe,
                                      prop_0 = 0.1
                                    ),
                                    center = T,scale = T) {

  # Set ID
  setID = length(MegaLMM_state$run_variables$trait_sets) + 1

  # check model
  if(length(lme4::findbars(fixed_formula))>0) stop('do not include random effects in fixed_formula. Random effects are inherited from the model for F')
  # recover()
  # check ID
  ID = MegaLMM_state$run_variables$ID
  F_data = MegaLMM_state$run_variables$F_data
  if(!ID %in% colnames(data)) stop(sprintf('ID: "%s% missing from data',ID))
  if(!length(data[[ID]]) == length(unique(data[[ID]]))) stop(sprintf('ID: "%s% not unique in data',ID))
  if(any(data[[ID]] %in% F_data[[ID]] == F)) stop("all IDs in data must be present in F_data")

  # clear Lambda_prior (if it exists) as this needs to be re-set with new traits
  MegaLMM_state$data_matrices$Lambda_X = NULL
  MegaLMM_state$priors$Lambda_prior = NULL

  # row/column names of Y
  if(is.null(rownames(Y))) {
    rownames(Y) = data[[ID]]
  } else if(!all(rownames(Y) == data[[ID]])) {
    stop(sprintf('rownames(Y) must equal data[["%s"]]',ID))
  }
  if(is.null(colnames(Y))) {
    # add column names of the form SetID.1-t
    colnames(Y) = sprintf('Set%03d.%02d',setID,1:ncol(Y))
  }

  # expand Y to full size of F_data
  Y_full = matrix(NA,nrow = nrow(F_data),ncol = ncol(Y),dimnames = list(F_data[[ID]],colnames(Y)))
  Y_full[rownames(Y),] = Y

  # tot_Y_var priors
  priors = MegaLMM_state$priors
  if(length(tot_Y_var_prior$V) == 1) {
    tot_Y_var_prior$V = rep(tot_Y_var_prior$V,ncol(Y))
    tot_Y_var_prior$nu = rep(tot_Y_var_prior$nu,ncol(Y))
  }
  # priors$tot_Y_var$V = c(priors$tot_Y_var$V,tot_Y_var_prior$V)
  # priors$tot_Y_var$nu = c(priors$tot_Y_var$nu,tot_Y_var_prior$nu)

  # resid_h2_priors
  # priors$h2_priors_resids = cbind(resid_h2_prior,matrix(resid_h2_prior,nrow = length(resid_h2_prior),ncol = ncol(Y)))

  # recover()
  # design matrices
  mf = model.frame(fixed_formula,data)
  if(any(is.na(mf))) stop('NAs in variables in fixed_formula in data')
  X1_Y = model.matrix(fixed_formula,mf)
  X1_full = matrix(NA,nrow = nrow(F_data),ncol = ncol(X1_Y),dimnames = list(F_data[[ID]],colnames(X1_Y)))
  X1_full[rownames(Y),] = X1_Y

  # collect trait info
  traitSet_info = list(
    setID = setID,
    Y = Y,
    rows = which(rownames(Y_full) %in% rownames(Y)),
    cols = ncol(MegaLMM_state$current_state$Eta) + 1:ncol(Y),
    colnames = colnames(Y),
    X1_cols = ncol(MegaLMM_state$data_matrices$X1) + 1:ncol(X1_full),
    use_X2 = use_X2,
    tot_Eta_prec_rate   = with(tot_Y_var_prior,V * nu),
    tot_Eta_prec_shape  = with(tot_Y_var_prior,nu - 1),
    h2_priors_resids = matrix(resid_h2_prior,nrow = length(resid_h2_prior),ncol = ncol(Y)),
    B2_R_prior = B2_R_prior
  )
  traitSet_info$sample_Eta = partial_matrix_model_sampler
  setup = partial_matrix_model_setup(Y_full,center=center,scale=scale)
  Eta = setup$Eta
  traitSet_info$sampler_parameters = setup$parameters

  # add to MegaLMM_state
  MegaLMM_state$run_variables$trait_sets = append(MegaLMM_state$run_variables$trait_sets,list(traitSet_info))
  MegaLMM_state$current_state$Eta = cbind(MegaLMM_state$current_state$Eta,Eta)
  MegaLMM_state$run_variables$p = ncol(MegaLMM_state$current_state$Eta)
  MegaLMM_state$data_matrices$X1 = cbind(MegaLMM_state$data_matrices$X1,X1_full)
  # MegaLMM_state$priors = priors

  MegaLMM_state
}

MegaLMM_add_Lambda_model = function(MegaLMM_state,design_matrices,
                                    Lambda_prior,
                                    Lambda_beta_prior) {

  # recover()
  n_matrices = length(design_matrices)
  if(n_matrices == 0) {
    Lambda_X = matrix(0,MegaLMM_state$run_parameters$K,ncol = 0)
    Lambda_X_groups = c()
  } else {
    if(!all(sapply(design_matrices,function(x) all(rownames(x) == colnames(MegaLMM_state$current_state$Eta))))) {
      stop('all design matrices must have rownames matching the column names of MegaLMM_state$current_state$Eta')
    }
    Lambda_X = do.call(cbind,design_matrices)
    Lambda_X_groups = do.call(c,lapply(1:n_matrices,function(i) rep(i,ncol(design_matrices[[i]]))))
  }

  if(ncol(Lambda_X) == 0 || length(Lambda_beta_prior$V) == 1) {
    Lambda_beta_prior$V = rep(Lambda_beta_prior$V,n_matrices)
    Lambda_beta_prior$nu = rep(Lambda_beta_prior$nu,n_matrices)
  }
  Lambda_beta_rate = Lambda_beta_prior$V * Lambda_beta_prior$nu
  Lambda_beta_shape = Lambda_beta_prior$nu - 1

  MegaLMM_state$priors$Lambda_prior = Lambda_prior

  MegaLMM_state$priors$Lambda_prior$Lambda_beta_rate = Lambda_beta_rate
  MegaLMM_state$priors$Lambda_prior$Lambda_beta_shape = Lambda_beta_shape
  MegaLMM_state$priors$Lambda_prior$Lambda_X_groups = Lambda_X_groups

  MegaLMM_state$data_matrices$Lambda_X = Lambda_X

  MegaLMM_state

}

MegaLMM_preliminary_calculations = function(MegaLMM_state,ncores = 1,verbose = FALSE ) {

  run_parameters = MegaLMM_state$run_parameters
  data_matrices = MegaLMM_state$data_matrices
  priors = MegaLMM_state$priors
  h2s_matrix = MegaLMM_state$data_matrices$h2s_matrix
  RE_setup = data_matrices$RE_setup
  F_data = MegaLMM_state$run_variables$F_data
  ID = MegaLMM_state$run_variables$ID

  trait_sets      = MegaLMM_state$run_variables$trait_sets
  # trait_sets[[1]] is the factors
  # rows should be any rows with non-missing data
  trait_sets[[1]]$rows = unique(do.call(c,lapply(trait_sets,function(x) x$rows)))
  # Missing_row_data_map  = MegaLMM_state$run_variables$Missing_row_data_map


  n = MegaLMM_state$run_variables$n
  p = MegaLMM_state$run_variables$p
  n_RE = length(RE_setup)
  RE_names = names(RE_setup)
  n_matrices = 2*ncol(h2s_matrix)

  pb = 0
  if(verbose) {
    print(sprintf("Pre-calculating random effect inverse matrices for %d groups of traits and %d sets of random effect weights", length(trait_sets), ncol(h2s_matrix)))
    pb = txtProgressBar(min=0,max = 2 + length(trait_sets) * n_matrices,style=3)
  }

  X1   = data_matrices$X1
  X2 = data_matrices$X2
  U2 = data_matrices$U2
  Z   = data_matrices$Z

  # function to ensure that covariance matrices are sparse and symmetric
  fix_K = function(x) forceSymmetric(drop0(x,tol = run_parameters$drop0_tol))

  # construct RE_L and ZL for each random effect
  for(i in 1:length(RE_setup)){
    re_name = names(RE_setup)[i]
    RE_setup[[i]] = within(RE_setup[[i]],{
      # recover()
      if(!'ZL' %in% ls()){
        if('K' %in% ls() && !is.null(K)){
          id_names = rownames(K)
          if(is(K,'Matrix') & isDiagonal(K)) {
            L = as(diag(1,nrow(K)),'dgCMatrix')
            K_inv = as(diag(1/diag(K)),'dgCMatrix')
          } else {
            ldl_k = LDLt(K)
            large_d = ldl_k$d > run_parameters$K_eigen_tol
            r_eff = sum(large_d)
            # if need to use reduced rank model, then use D of K in place of K and merge L into Z
            # otherwise, use original K, set L = Diagonal(1,r)
            if(r_eff < length(ldl_k$d)) {
              K = K_inv = as(diag(1,r_eff),'dgCMatrix')
              L = t(ldl_k$P) %*% t(sqrt(ldl_k$d[large_d])*t(ldl_k$L[,large_d]))
              if(is(L,'dgeMatrix')) L = as.matrix(L)
            } else{
              L = as(diag(1,nrow(K)),'dgCMatrix')
              K_inv = as(with(ldl_k,t(P) %*% crossprod(diag(1/sqrt(d)) %*% solve(L)) %*% P),'dgCMatrix')
            }
            rm(list=c('ldl_k','large_d','r_eff'))
          }
          if(is.null(rownames(K))) rownames(K) = 1:nrow(K)
          rownames(K_inv) = rownames(K)
        } else if ('K_inv' %in% ls() && !is.null(K_inv)){
          id_names = rownames(K_inv)
          if(is.null(rownames(K_inv))) rownames(K_inv) = 1:nrow(K_inv)
          K = solve(K_inv)
          rownames(K) = rownames(K_inv)
          L = as(diag(1,nrow(K)),'dgCMatrix')
        } else{
          K = as(diag(1,ncol(Z)),'dgCMatrix')
          rownames(K) = colnames(Z)
          id_names = rownames(K)
          K_inv = K
          L = as(diag(1,nrow(K)),'dgCMatrix')
        }
        # if(is.null(id_names)) id_names = 1:length(id_names)
        rownames(L) = paste(id_names,re_name,sep='::')
        K = fix_K(K)
        ZL = Z %*% L
      }
    })
  }
  if(verbose) setTxtProgressBar(pb,getTxtProgressBar(pb)+1)

  ZL = do.call(cbind,lapply(RE_setup,function(re) re$ZL))
  if(nnzero(ZL)/length(ZL) < .25) {
    ZL = as(ZL,'dgCMatrix')
  } else{
    ZL = as.matrix(ZL)
  }


  if(length(RE_setup) > 1) {
    RE_L = do.call(bdiag,lapply(RE_setup,function(re) re$L))
    rownames(RE_L) = do.call(c,lapply(RE_setup,function(re) rownames(re$L)))
  } else{
    RE_L = RE_setup[[1]]$L
  }
  if(nnzero(RE_L)/length(RE_L) < 0.25) {
    RE_L = as(RE_L,'dgCMatrix')
  } else{
    RE_L = as.matrix(RE_L)
  }
  r_RE = sapply(RE_setup,function(re) ncol(re$ZL))


  # cholesky decompositions (RtR) of each K_inverse matrix
  chol_Ki_mats = lapply(RE_setup,function(re) {
    K_inv = re$K_inv
    if(is(K_inv,'Matrix')) {
      if(isDiagonal(K_inv)) {
        sparseMatrix(i=1:nrow(K_inv),j=1:nrow(K_inv),x=sqrt(diag(K_inv)))
      } else{
        as(chol(forceSymmetric(K_inv)),'dgCMatrix')
      }
    } else {
      as(chol(K_inv),'dgCMatrix')
    }
  })
  if(verbose) setTxtProgressBar(pb,getTxtProgressBar(pb)+1)

  Qt_list_rows = matrix(FALSE,nrow = length(F_data[[ID]]),ncol = 0)
  Qt_list = list()
  QtX1_list = list()
  QtX2_R_list = list()
  ZL_list = list()
  RE_L_list = list()

  chol_V_list_list = list()
  chol_ZtZ_Kinv_list_list = list()

  svd_K1 = NULL
  for(set in seq_along(trait_sets)){
    if(verbose>1) print(sprintf('Set %d',set))
    rows = trait_sets[[set]]$rows
    cols = trait_sets[[set]]$cols
    X1_cols = trait_sets[[set]]$X1_cols
    if(sum(rows) == 0) next

    # check if this set matches with any previous set - we can the re-use the matrices
    Qt_list_ID = 1
    if(set > 1) {
      while(Qt_list_ID <= length(Qt_list)) {
        if(all(Qt_list_rows[rows,Qt_list_ID]) && all(!Qt_list_rows[-rows,Qt_list_ID])) {
          # we found a match
          break
        }
        Qt_list_ID = Qt_list_ID + 1
      }
    }
    trait_sets[[set]]$Qt_list_ID = Qt_list_ID
    if(Qt_list_ID <= length(Qt_list)) {
      Qt = Qt_list[[Qt_list_ID]]
      if(trait_sets[[set]]$use_X2 & (length(QtX2_R_list) < Qt_list_ID || is.null(QtX2_R_list[[Qt_list_ID]]))) {
        QtX2_R_set = Qt %**% X2_R[rows,,drop=FALSE]
        QtX2_R_list[[Qt_list_ID]]  = QtX2_R_set
      }
      # update progress bar
      if(verbose) setTxtProgressBar(pb,getTxtProgressBar(pb)+n_matrices)
      next
    }

    # add new column to Qt_list_rows
    Qt_list_rows = cbind(Qt_list_rows,FALSE)
    Qt_list_rows[rows,Qt_list_ID] = TRUE

    # find Qt = svd(ZLKLtZt)$u
    if(ncol(ZL) < nrow(ZL[rows,])*.9) {
      # a faster way of taking the SVD of ZLKZLt, particularly if ncol(ZL) < nrow(ZL). Probably no benefit if ncol(K) > nrow(ZL)
      if(is.null(svd_K1)){
        svd_K1 = svd(RE_setup[[1]]$K)
      }
      qr_ZU = qr(RE_setup[[1]]$ZL[rows,,drop=FALSE] %*% svd_K1$u)
      R_ZU = drop0(qr.R(qr_ZU,complete=F),tol=run_parameters$drop0_tol)
      Q_ZU = drop0(qr.Q(qr_ZU,complete=T),tol=run_parameters$drop0_tol)
      RKRt = R_ZU %*% diag(svd_K1$d) %*% t(R_ZU)
      svd_RKRt = svd(RKRt)
      RKRt_U = svd_RKRt$u
      if(ncol(Q_ZU) > ncol(RKRt_U)) RKRt_U = bdiag(RKRt_U,diag(1,ncol(Q_ZU)-ncol(RKRt_U)))
      Qt = t(Q_ZU %**% RKRt_U)
    } else{
      ZKZt = with(RE_setup[[1]],ZL[rows,,drop=FALSE] %*% K %*% t(ZL[rows,,drop=FALSE]))
      if(is(ZKZt,'Matrix') & isDiagonal(ZKZt)) {
        nn = nrow(ZKZt)
        result = list(d = diag(ZKZt),u = sparseMatrix(i=1:nn,j=1:nn,x=1))
      } else{
        result = svd(ZKZt)
      }
      Qt = t(result$u)
    }
    Qt = as(drop0(as(Qt,'dgCMatrix'),tol = run_parameters$drop0_tol),'dgCMatrix')
    if(nnzero(Qt)/length(Qt) > 0.5) Qt = as.matrix(Qt)  # only store as sparse if it is sparse

    QtZL_matrices_set = lapply(RE_setup,function(re) Qt %*% re$ZL[rows,,drop=FALSE])
    # QtZL_set = do.call(cbind,QtZL_matrices_set[RE_names])
    # if(nnzero(QtZL_set)/length(QtZL_set) > 0.5)  QtZL_set = as(QtZL_set,'dgCMatrix')
    QtX1_set = Qt %**% X1[rows,X1_cols,drop=FALSE]
    QtX1_keepColumns = logical(0)
    if(ncol(QtX1_set)>0)  QtQtX1_set_keepColumns = c(1:ncol(QtX1_set)) %in% caret::findLinearCombos(QtX1_set)$remove == F # assess which columns of X1 are identifiable in this data set
    QtX1_set = QtX1_set[,QtX1_keepColumns,drop=FALSE]  # drop extra columns

    Qt_list[[Qt_list_ID]]   = Qt
    QtX1_list[[Qt_list_ID]] = list( X1 = QtX1_set,
                                    keepColumns = QtX1_keepColumns
                                    )

    ZKZts_set = list()
    for(i in 1:n_RE){
      ZKZts_set[[i]] = forceSymmetric(drop0(QtZL_matrices_set[[i]] %*% RE_setup[[i]]$K %*% t(QtZL_matrices_set[[i]]),tol = run_parameters$drop0_tol))
      ZKZts_set[[i]] = as(as(ZKZts_set[[i]],'CsparseMatrix'),'dgCMatrix')
      if(nnzero(ZKZts_set[[i]])/length(ZKZts_set[[i]]) > 0.5) {
        ZKZts_set[[i]] = as.matrix(ZKZts_set[[i]])
      }
    }

    chol_V_list_list[[Qt_list_ID]] = make_chol_V_list(ZKZts_set,h2s_matrix,run_parameters$drop0_tol,
                                               verbose,pb,setTxtProgressBar,getTxtProgressBar,ncores)
    # convert any to dense if possible
    for(i in 1:length(chol_V_list_list[[Qt_list_ID]])){
      chol_V_list_list[[Qt_list_ID]][[i]] = drop0(chol_V_list_list[[Qt_list_ID]][[i]],tol = run_parameters$drop0_tol)
      if(nnzero(chol_V_list_list[[Qt_list_ID]][[i]])/length(chol_V_list_list[[Qt_list_ID]][[i]]) > 0.25){
        chol_V_list_list[[Qt_list_ID]][[i]] = as.matrix(chol_V_list_list[[Qt_list_ID]][[i]])
      }
    }

    if(verbose>1) print(sprintf('Set %d S',set))
    ZL_list[[Qt_list_ID]] = ZL
    chol_Ki_mats_set = chol_Ki_mats
    ZtZ_set = as(forceSymmetric(drop0(crossprod(ZL_list[[Qt_list_ID]][rows,]),tol = run_parameters$drop0_tol)),'dgCMatrix')
    if(length(RE_setup) == 1) {
      S = simultaneous_diagonalize(ZtZ_set,solve(chol_Ki_mats[[1]]))$S
      if(nnzero(S)/length(S) > 0.5) {
        S = as.matrix(S)  # only store as sparse if it is sparse
      } else {
        S = as(S,'dgCMatrix')
      }
      ZL_list[[Qt_list_ID]] = ZL_list[[Qt_list_ID]] %**% S
      ZtZ_set = Diagonal(ncol(ZL),colSums(ZL_list[[Qt_list_ID]][rows,]^2))
      # ZtZ_set = as(forceSymmetric(drop0(crossprod(ZL_list[[Qt_list_ID]][x,]),tol = run_parameters$drop0_tol)),'dgCMatrix')  # Not needed because must be symmetric
      RE_L_list[[Qt_list_ID]] = RE_setup[[1]]$L %**% S
      chol_Ki_mats_set[[1]] = as(diag(1,nrow(ZtZ_set)),'dgCMatrix')
    }

    if(verbose>1) print(sprintf('Set %d chol_ZtZ_Kinv_list_list',set))
    if(length(RE_setup) == 1 & isDiagonal(ZtZ_set) & all(sapply(chol_Ki_mats_set,isDiagonal))) {
      # in the case of 1 random effect with diagonal covariance matrix, we can skip the expensive calculations
      chol_ZtZ_Kinv_list_list[[Qt_list_ID]] = sapply(1:length(h2s_matrix),function(i) {
        nn = nrow(ZtZ_set)
        sparseMatrix(i=1:nn,j=1:nn,x = sqrt(1/(1-h2s_matrix[1,i])*diag(ZtZ_set) + 1/h2s_matrix[1,i]*diag(chol_Ki_mats_set[[1]])^2))
      })
      if(verbose) setTxtProgressBar(pb,getTxtProgressBar(pb)+length(h2s_matrix))
    } else if(set == 1 || sum(priors$h2_priors_resids[-1]) > 0) {
      # Note: set >1 only applies to the resids (factors always use set==1).
      chol_ZtZ_Kinv_list_list[[Qt_list_ID]] = make_chol_ZtZ_Kinv_list(chol_Ki_mats_set,h2s_matrix,ZtZ_set,
                                                               run_parameters$drop0_tol,
                                                               verbose,pb,setTxtProgressBar,getTxtProgressBar,ncores)
    } else{
      # if the prior is concentrated at h2s==0, then we don't need to sample U. All values will be 0.
      chol_ZtZ_Kinv_list_list[[Qt_list_ID]] = lapply(seq_along(priors$h2_priors_resids),function(x) matrix(1,0,0))
      if(verbose) setTxtProgressBar(pb,getTxtProgressBar(pb)+length(h2s_matrix))
    }
  }
  if(verbose) close(pb)


  # Qt matrices for factors are only used with row set 1
  rows = trait_sets[[1]]$rows
  Qt1_U2_F = Qt_list[[1]] %**% U2[rows,,drop=FALSE]
  Qt1_X2_F = Qt_list[[1]] %**% X2[rows,,drop=FALSE]


  MegaLMM_state$run_variables = c(MegaLMM_state$run_variables,
                                  list(
                                    Qt_list    = Qt_list,
                                    QtX1_list   = QtX1_list,
                                    QtX2_R_list = QtX2_R_list,
                                    Qt1_U2_F = Qt1_U2_F,
                                    Qt1_X2_F = Qt1_X2_F,
                                    chol_V_list_list          = chol_V_list_list,
                                    chol_ZtZ_Kinv_list_list = chol_ZtZ_Kinv_list_list
                                  ))
  MegaLMM_state$run_variables$trait_sets = trait_sets

  MegaLMM_state$data_matrices$ZL = ZL
  MegaLMM_state$data_matrices$RE_setup = RE_setup
  MegaLMM_state$data_matrices$ZL_list = ZL_list
  MegaLMM_state$data_matrices$RE_L_list = RE_L_list
  MegaLMM_state$run_variables$r_RE = r_RE

  return(MegaLMM_state)
}

MegaLMM_initialize_variables = function(MegaLMM_state,...){
  run_parameters = MegaLMM_state$run_parameters
  run_variables = MegaLMM_state$run_variables
  data_matrices = MegaLMM_state$data_matrices
  priors = MegaLMM_state$priors
  
  trait_sets = run_variables$trait_sets
  
  if(!'Qt_list' %in% names(run_variables)) stop('run MegaLMM_preliminary_calculations() first')
  
  
  
  current_state = with(c(run_parameters,run_variables,data_matrices,priors),{
  # recover()
  # MegaLMM_state$current_state = list()
  # MegaLMM_state$current_state$trait_sets_current_state = list()
    Eta = MegaLMM_state$current_state$Eta
    
    # Factor discrete variances
    # K-matrix of n_RE x K with
    F_h2_index = sample(c(1:ncol(h2s_matrix))[h2_priors_factors>0],K,replace=T)
    F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]
    tot_F_prec = matrix(1,nrow=1,ncol=K)
    
    
    U_F = matrix(rnorm(sum(r_RE) * K, 0, sqrt(F_h2[1,] / tot_F_prec)),ncol = K, byrow = T)
    rownames(U_F) = colnames(ZL)
    
    # Factor fixed effects
    B2_F = 0*matrix(rnorm(b2 * K),b2,K)
    rownames(B2_F) = colnames(X2)
    
    # recover()
    F = X2 %*% B2_F + ZL %*% U_F + matrix(rnorm(n * K, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = K, byrow = T)
    F = as.matrix(F)
    
    Lambda_mean = matrix(0,K,p)
    colnames(Lambda_mean) = colnames(Eta)
    
    
    trait_sets_current_state = list()
  
    for(set in seq_along(trait_sets)) {
      if(length(trait_sets[[set]]$cols) == 0) next
      trait_sets_current_state[[set]] = with(trait_sets[[set]],{
        
      # recover()
      # trait_set = trait_sets[[set]]
      # MegaLMM_state$current_state$trait_sets_current_state[[set]] = list()
      # if(length(trait_sets[[set]]$cols) == 0) next
      
      # p = length(trait_set$cols)
      # traitnames = trait_set$colnames
        p = length(cols)
        traitnames = colnames
      
        Lambda_prec = matrix(1,K,p)
        Lambda = matrix(rnorm(K*p,0,sqrt(1/Lambda_prec)),nr = K,nc = p)
        colnames(Lambda) = traitnames
        
        tot_Eta_prec = matrix(rgamma(p,shape = tot_Eta_prec_shape,rate = tot_Eta_prec_rate),nrow = 1)
        colnames(tot_Eta_prec) = traitnames
        
        resid_h2_index = matrix(apply(h2_priors_resids,2,function(x) sample(c(1:ncol(h2s_matrix))[x>0],1)),nrow=1)
        resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]
        
        U_R = matrix(rnorm(sum(r_RE) * p, 0, sqrt(resid_h2[1,] / tot_Eta_prec)),ncol = p, byrow = T)
        colnames(U_R) = traitnames
        rownames(U_R) = colnames(ZL)
        
        # Fixed effects
        b1 = ncol(X1)
        B1 = matrix(rnorm(b1*p), ncol = p)
        colnames(B1) = traitnames
        
        b2_R = ifelse(use_X2,b2,0)
        B2_R = 0*matrix(rnorm(b2_R*p),b2_R,ncol = p)
        colnames(B2_R) = traitnames
        rownames(B2_R) = colnames(X2)
       
        XB = X1[,X1_cols,drop=FALSE] %**% B1[X1_cols,,drop=FALSE]
        if(use_X2) XB = XB + X2 %*% B2_R
        
        var_Eta = rep(1,p)
        
        Eta_mean = XB + F %**% Lambda + ZL_list[[Qt_list_ID]] %**% U_R
        # recover()
        list(
          cols = cols,
          Eta = Eta[,cols],
          Eta_mean = Eta_mean,
          Lambda = Lambda,
          tot_Eta_prec = tot_Eta_prec,
          resid_h2_index = resid_h2_index,
          resid_h2 = resid_h2,
          U_R = U_R,
          B1 = B1,
          B2_R = B2_R,
          XB = XB,
          var_Eta = var_Eta
        )
      })
      trait_sets_current_state[[set]] = trait_sets[[set]]$sample_Eta(trait_sets_current_state[[set]],trait_sets[[set]]$sampler_parameters)
    #   trait_sets_current_state[[set]] = trait_sets[[set]]$B2_R_prior$sampler(list(current_state = trait_sets_current_state[[set]],
    #                                                                             data_matrices = data_matrices,
    #                                                                              = trait_sets[[set]]$B2_R_prior))
    }
    list(
      K              = K,
      Kr = K,
      tot_F_prec     = tot_F_prec,
      F_h2_index     = F_h2_index,
      F_h2           = F_h2,
      U_F            = U_F,
      F              = F,
      B2_F            = B2_F,
      Lambda_mean = Lambda_mean,
      trait_sets_current_state = trait_sets_current_state,
      nrun           = 0,
      total_time     = 0
    )
  })
  
  MegaLMM_state$current_state = assemble_current_state(current_state)
  # recover()
  # Initialize parameters for Lambda_prior, B_prior, and QTL_prior (may be model-specific)
  # MegaLMM_state$current_state = MegaLMM_state$priors$Lambda_prior$sampler(MegaLMM_state)
  
  # save the initial RNG state
  MegaLMM_state$current_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )

  # ----------------------------- #
  # --- Initialize MegaLMM_state --- #
  # ----------------------------- #

  # Posterior = reset_Posterior(MegaLMM_state$Posterior,MegaLMM_state)
  # MegaLMM_state$Posterior = Posterior
  # 
  # MegaLMM_state$Posterior$posteriorSample_params = unique(c(MegaLMM_state$Posterior$posteriorSample_params,observation_model_state$posteriorSample_params))
  # MegaLMM_state$Posterior$posteriorMean_params = unique(c(MegaLMM_state$Posterior$posteriorMean_params,observation_model_state$posteriorMean_params))

  return(MegaLMM_state)
}

assemble_matrix = function(current_state,trait_sets,matrix,by = 'cols') {
  if(by=='cols') {
    new_matrix = do.call(cbind,lapply(current_state$trait_sets_current_state[trait_sets],function(x) x[[matrix]]))
  } else{
    new_matrix = do.call(rbind,lapply(current_state$trait_sets_current_state[trait_sets],function(x) x[[matrix]]))
  }
  new_matrix
}
assemble_current_state = function(current_state,
                                  sets = 1:length(current_state$trait_sets_current_state),
                                  matrices = c('Lambda','B1','B2_R','U_R','tot_Eta_prec','resid_h2_index','resid_h2','XB','Eta_mean','Eta')) {
  for(i in seq_along(matrices)) {
    current_state[[matrices[i]]] = assemble_matrix(current_state,sets,matrices[i],'cols')
  }
  current_state
}


partition_data_missing = function(Y,max_NA_groups = ncol(Y), starting_groups = NULL, verbose=FALSE) {
  Y_missing = is.na(Y)
  Y_missing_mat = as.matrix(Y_missing)
  non_missing_rows = unname(which(rowSums(!Y_missing)>0))

  # if(is.null(max_NA_groups)) max_NA_groups = MegaLMM_state$run_parameters$max_NA_groups
  # if(is.null(verbose)) verbose = MegaLMM_state$run_parameters$verbose


  # create a base map
  Missing_data_map = list(list(
    Y_obs = non_missing_rows,
    Y_cols = 1:ncol(Y_missing)
  ))

  # if no groups allowed, just return base map
  if(max_NA_groups <= 1) {
    return(list(Missing_data_map = Missing_data_map, Missing_data_map_list = NULL,map_results=c()))
  }
  if(is.infinite(max_NA_groups)) {
    Y_col_obs = lapply(1:ncol(Y_missing_mat),function(x) {
      obs = which(!Y_missing_mat[,x],useNames=F)
      names(obs) = NULL
      obs
    })
    non_missing_rows = unname(which(rowSums(!Y_missing_mat)>0))
    unique_Y_col_obs = unique(c(list(non_missing_rows),Y_col_obs))
    unique_Y_col_obs_str = lapply(unique_Y_col_obs,paste,collapse='')
    Y_col_obs_index = sapply(Y_col_obs,function(x) which(unique_Y_col_obs_str == paste(x,collapse='')))

    Missing_data_map = lapply(1:max(unique(Y_col_obs_index)),function(i) {
      Y_cols = which(Y_col_obs_index==i)
      if(length(Y_cols) == 0) {
        return(list(
          Y_obs = non_missing_rows,
          Y_cols = Y_cols
        ))
      } else{
        return(list(
          Y_obs = unname(which(rowSums(!Y_missing[,Y_cols,drop=FALSE])>0)),  # now find the rows with any non-missing data in this set of columns
          Y_cols = Y_cols
        ))
      }
    })
    return(list(
      Missing_data_map = Missing_data_map,
      Missing_data_map_list = NULL,
      map_results = c()
    ))
  }

  # if a starting map is provided, initialize this
  if(!is.null(starting_groups)) {
    if(length(starting_groups) != ncol(Y_missing)) stop("starting_groups must be vector with length ncol(Y).")
    starting_groups = as.numeric(as.factor(starting_groups))
    Missing_data_map = lapply(unique(starting_groups),function(i) {
      Y_cols = which(starting_groups == i)
      Y_obs = unname(which(rowSums(!Y_missing[,Y_cols,drop=FALSE])>0))
      list(
        Y_obs = Y_obs,
        Y_cols = Y_cols
      )
    })
    if(max(sapply(Missing_data_map,function(x) length(x$Y_obs))) < length(non_missing_rows)) {
      Missing_data_map = c(list(list(
        Y_obs = non_missing_rows,
        Y_cols = c()
      )), Missing_data_map)
    }
  }

  # otherwise, use kmeans to sequentially split map with the most non-excluded NAs
  # save the sequence of maps so the user can select the best

  Missing_data_map_list = list(
    Missing_data_map
  )

  iter = 1
  map_results = data.frame(
    map = iter,
    N_groups = length(Missing_data_map),
    max_group_size = length(non_missing_rows),
    total_kept_NAs = sum(sapply(Missing_data_map,function(x) {
      if(length(x$Y_cols) == 0) return(0)
      sum(Y_missing_mat[x$Y_obs,x$Y_cols])
    }))
  )
  if(verbose) print(map_results)

  while(length(Missing_data_map) < max_NA_groups && iter < 2*max_NA_groups){
    iter = iter+1
    # count the number of non-excluded NAs in each group
    group_NAs = sapply(Missing_data_map,function(x) {
      if(length(x$Y_cols)<= 1) return(0) # don't split a group with only one column
      sum(Y_missing_mat[x$Y_obs,x$Y_cols])
    })
    if(max(group_NAs) == 0) break

    target = which.max(group_NAs)  # group with the most NAs
    group = Missing_data_map[[target]]

    # split this cluster into two
    if(length(group$Y_cols) == 2) {
      clusters = list(cluster = c(1,2))
    } else{
      clusters = with(group,kmeans(t(Y_missing[Y_obs,Y_cols]),centers=2))
    }

    new_groups = lapply(seq_len(max(clusters$cluster)),function(i) {
      Y_cols = group$Y_cols[which(clusters$cluster==i)]
      Y_obs = unname(which(rowSums(!Y_missing_mat[,Y_cols,drop=FALSE])>0))  # now find the rows with any non-missing data in this set of columns
      # if((length(non_missing_rows)-length(Y_obs))/length(non_missing_rows) <= min_perc_drop) {
      #   Y_obs = non_missing_rows  # if Y_obs is too big, not worth keeping this cluster
      #   # this might make an infinite loop!
      # }
      list(
        Y_obs = Y_obs,
        Y_cols = Y_cols
      )
    })
    new_groups = new_groups[order(sapply(new_groups,function(x) length(x$Y_obs)),decreasing=T)]

    if(target == 1) {
      if(length(new_groups[[1]]$Y_obs) == length(non_missing_rows)) {
        Missing_data_map[1] = new_groups[1]
        new_groups = new_groups[-1]
      } else{
        Missing_data_map[[1]]$Y_cols = numeric()
      }
      Missing_data_map = c(Missing_data_map,new_groups)
    } else{
      Missing_data_map = c(Missing_data_map[-target],new_groups)
    }

    Missing_data_map_list[[iter]] = Missing_data_map

    results_i = data.frame(
      map = iter,
      N_groups = length(Missing_data_map),
      max_group_size = max(sapply(Missing_data_map,function(x) length(x$Y_obs))[-1]),
      total_kept_NAs = sum(sapply(Missing_data_map,function(x) {
        if(length(x$Y_cols) == 0) return(0)
        sum(Y_missing_mat[x$Y_obs,x$Y_cols])
      }))
    )

    if(verbose) print(results_i)

    map_results = rbind(map_results,results_i)
  }

  return(list(
    Missing_data_map = Missing_data_map_list[[length(Missing_data_map_list)]],
    Missing_data_map_list = Missing_data_map_list,
    map_results = map_results
  ))
}

