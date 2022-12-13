my_detectCores = function() {
  ncores = try(suppressWarnings(as.numeric(system('printenv SLURM_CPUS_PER_TASK',intern=T))),silent=T)
  if(length(ncores) == 0 || is.na(ncores)) ncores = parallel::detectCores()
  ncores
}

optimize_n_threads = function(MegaLMM_state,n_threads,times = 10) {
  exprs = lapply(n_threads,function(threads) bquote({set_MegaLMM_nthreads(.(threads));sample_latent_traits(MegaLMM_state)}))
  names(exprs) = n_threads
  res=microbenchmark::microbenchmark(list = exprs,times = times)
  summary_res = summary(res)
  list(optim = as.numeric(as.character(summary_res$expr[order(summary_res$mean)[1]])),results = summary_res)
}


make_model_setup = function(formula,data,relmat = NULL) {
  # ensure data has rownames
  if(is.null(rownames(data))) rownames(data) = 1:nrow(data)

  # check that RE's have approporiate levels in data
  RE_levels = list() # a list of levels for each of the random effects
  if(is.null(relmat)) relmat = list()
  for(re in names(relmat)) {
    # check that K is a matrix, then convert to Matrix
    if(is.list(relmat[[re]]) && !is.data.frame(relmat[[re]])){
      if(is.data.frame(relmat[[re]]$K)) relmat[[re]]$K = as.matrix(relmat[[re]]$K)
      if(is.matrix(relmat[[re]]$K)) relmat[[re]]$K = Matrix(relmat[[re]]$K,sparse=T)
      if(is.null(rownames(relmat[[re]]$K))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]]$K)
    } else{
      if(is.data.frame(relmat[[re]])) relmat[[re]] = as.matrix(relmat[[re]])
      if(is.matrix(relmat[[re]])) relmat[[re]] = Matrix(relmat[[re]],sparse=T)
      if(is.null(rownames(relmat[[re]]))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]])
    }
  }
  for(re in names(RE_levels)){
    if(!re %in% colnames(data)) stop(sprintf('Column "%s" required in data',re))
    data[[re]] = as.factor(data[[re]]) # ensure 'data[[re]]' is a factor
    if(!all(data[[re]] %in% RE_levels[[re]])) stop(sprintf('Levels of random effect %s missing.',re))
    data[[re]] = factor(data[[re]],levels = RE_levels[[re]]) # add levels to data[[re]]
  }

  # Use lme4 to evaluate formula in data
  # ensure there is a response
  response = 'y'
  while(response %in% all.vars(formula)){
    response = paste0(response,response)
  }
  formula = sprintf('%s~%s',response,as.character(formula)[2])
  data[[response]] = 1
  lmod <- lme4::lFormula(formula,data=data,#weights=weights,
                         control = lme4::lmerControl(check.nobs.vs.nlev = 'ignore',check.nobs.vs.nRE = 'ignore'))

  # compute RE_setup
  RE_terms = lmod$reTrms

  # construct the RE_setup list
  # contains:
  # Z: n x r design matrix
  # K: r x r PSD covariance matrix
  RE_setup = list()
  for(i in 1:length(RE_terms$cnms)){
    term = names(RE_terms$cnms)[i]
    n_factors = length(RE_terms$cnms[[i]])  # number of factors for this grouping factor

    # extract combined Z matrix
    combined_Zt = RE_terms$Ztlist[[i]]
    Zs_term = tapply(1:nrow(combined_Zt),gl(n_factors,1,nrow(combined_Zt),labels = RE_terms$cnms[[i]]),function(x) Matrix::t(combined_Zt[x,,drop=FALSE]))

    # extract K from relmat. If missing, assign to NULL
    K = NULL
    if(term %in% names(relmat)) {
      K = relmat[[term]]
      if(is(K,'dsCMatrix')) K = as(K,'dgCMatrix')
      if(nnzero(K)/length(K) > .5) K = as.matrix(K)  # if not really sparse
    }

    if(!is.null(K)) {
      if(!all(colnames(Zs_term[[1]]) %in% rownames(K))) stop('rownames of K not lining up with Z')
      colnames(K) = rownames(K)
      K = K[colnames(Zs_term[[1]]),colnames(Zs_term[[1]])]
    } else {
      K = as(diag(1,ncol(Zs_term[[1]])),'dgCMatrix')
      rownames(K) = colnames(K) = colnames(Zs_term[[1]])
    }


    # make an entry in RE_setup for each random effect
    for(j in 1:n_factors){
      # name of variance component
      name = term
      if(n_factors > 1) name = paste(name,RE_terms$cnms[[i]][[j]],sep='.')
      while(name %in% names(RE_setup)) name = paste0(name,'.1') # hack for when same RE used multiple times

      # Z matrix
      Z = as(Zs_term[[j]],'dgCMatrix')


      RE_setup[[name]] = list(
        term = term,
        Z = Z,
        K = K
      )
    }

    # add names to RE_setup if needed
    n_RE = length(RE_setup)
    for(i in 1:n_RE){
      if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
        names(RE_setup)[i] = paste0('RE.',i)
      }
    }

  }
  return(list(lmod = lmod, RE_setup = RE_setup))
}

make_Missing_data_map = function(MegaLMM_state,max_NA_groups = NULL, starting_groups = NULL, verbose=NULL) {
  Y_missing = MegaLMM_state$run_parameters$observation_model_parameters$observation_setup$Y_missing
  Y_missing_mat = as.matrix(Y_missing)
  non_missing_rows = unname(which(rowSums(!Y_missing)>0))

  if(is.null(max_NA_groups)) max_NA_groups = MegaLMM_state$run_parameters$max_NA_groups
  if(is.null(verbose)) verbose = MegaLMM_state$run_parameters$verbose


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


set_Missing_data_map = function(MegaLMM_state,Missing_data_map) {

  MegaLMM_state$run_variables$Missing_data_map = Missing_data_map

  # find the values that are excluded by the groups in Missing_data_map
  Y_missing_mat = matrix(1,MegaLMM_state$run_variables$n,MegaLMM_state$run_variables$p)
  for(map in Missing_data_map) {
    Y_missing_mat[map$Y_obs,map$Y_cols] = 0
  }

  # create an object that indexes the non-missing values by row
  Y_row_obs = lapply(1:nrow(Y_missing_mat),function(x) {
    obs = which(!Y_missing_mat[x,],useNames=F)
    names(obs) = NULL
    obs
  })
  non_missing_cols = unname(which(colSums(!Y_missing_mat)>0))
  unique_Y_row_obs = unique(c(list(non_missing_cols),Y_row_obs))
  unique_Y_row_obs_str = lapply(unique_Y_row_obs,paste,collapse=',')
  Y_row_obs_index = sapply(Y_row_obs,function(x) which(unique_Y_row_obs_str == paste(x,collapse=',')))

  Missing_row_data_map = lapply(seq_along(unique_Y_row_obs),function(i) {
    x = unique_Y_row_obs[[i]]
    return(list(
      Y_cols = x,
      Y_obs = which(Y_row_obs_index == i)
    ))
  })
  
  # expand ZL_list and RE_L_list for missing data map
  
  MegaLMM_state$data_matrices$ZL_list = lapply(seq_along(Missing_data_map),function(x) MegaLMM_state$data_matrices$ZL_list[[1]])
  MegaLMM_state$data_matrices$RE_L_list = lapply(seq_along(Missing_data_map),function(x) MegaLMM_state$data_matrices$RE_L_list[[1]])
  

  MegaLMM_state$run_variables$Missing_row_data_map = Missing_row_data_map

  return(MegaLMM_state)

}

# finds a matrix S that simultaneously diagonalizes A and B
# following algorithm here: https://math.stackexchange.com/questions/1079627/simultaneously-diagonalization-of-two-matrices
simultaneous_diagonalize = function(A,Binvsq) {
  sBAB = svd(t(Binvsq) %*% A %*% Binvsq)
  O = sBAB$u
  S = Binvsq %*% O
  return(list(S=S,d=sBAB$d))
}
