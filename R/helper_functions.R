#' Estimate the memory requirment for a fully initialized MegaLMM model
#' 
#' the \code{initialize_MegaLMM()} function can take a while to run, especially
#' when there are multiple random effects, large numbers of observations, and multiple
#' groups of traits (with different missing data patterns, stored in the \code{Missing_data_map}).
#' Approximately, the initialization time and memory requirements will be linear in the 
#' number of groups of traits, and scale with h2_divisions^(# random effects).
#' 
#' This function will initialize the MegaLMM model for a single h2 vector and extrapolate
#' to the full grid of h2 vectors, enabeling you to estimate if you have enough
#' memory allocated to call \code{initialize_MegaLMM()}
#'
#' @param MegaLMM_state The model after calling \code{initialize_variables_MegaLMM}
#'
#' @return The estimated memory size in bytes
#' @seealso \code{\link{estimate_memory_posterior}}
#' @export
#'
#' @examples
#' estimate_memory_initialization_MegaLMM(MegaLMM_state)
estimate_memory_initialization_MegaLMM = function(MegaLMM_state) {
  object_size = pryr::object_size
  print(sprintf('Random effects: %s',paste(names(MegaLMM_state$data_matrices$RE_setup),collapse=', ')))
  print(sprintf('%d groups of traits and %d h2 grid cells',
                length(MegaLMM_state$run_variables$Missing_data_map),
                ncol(MegaLMM_state$data_matrices$h2s_matrix)))
  
  base_size = object_size(MegaLMM_state)
  base_sizes = lapply(MegaLMM_state,object_size)
  
  # create new MegaLMM_state object with only 1 h2-grid cell
  # initialize this model and measure its size
  # then extrapolate what the full model would be
  MegaLMM_state_grid1 = MegaLMM_state
  h2s_matrix = MegaLMM_state_grid1$data_matrices$h2s_matrix
  MegaLMM_state_grid1$data_matrices$h2s_matrix = h2s_matrix[,ncol(h2s_matrix),drop=FALSE]
  
  MegaLMM_state_grid1 = initialize_MegaLMM(MegaLMM_state_grid1,verbose = FALSE)
  grid1_sizes = lapply(MegaLMM_state_grid1,object_size)
  run_variables_base = MegaLMM_state_grid1$run_variables
  run_variables_base = run_variables_base[names(run_variables_base) %in% c('chol_V_list_list','chol_ZtZ_Kinv_list_list') == F]
  
  estimated_run_variable_size = object_size(run_variables_base) + ncol(h2s_matrix) * object_size(MegaLMM_state_grid1$run_variables[c('chol_V_list_list','chol_ZtZ_Kinv_list_list')])
  
  print(sprintf('Estimated initialized size: %s Gb',
                format(signif((estimated_run_variable_size+base_size)/(1024^3),2))))
  rm('MegaLMM_state_grid1')
  gc()
  invisible(estimated_run_variable_size+base_size)
}

#' Estimates the memory required to store a set of posterior samples collected by \code{sample_MegaLMM}
#' 
#' A call to \code{sample_MegaLMM(MegaLMM_state,n_iter)} will run \code{n_iter}
#' of the Gibbs sampler. If \code{nrun > burn}, then a posterior sample of all variables
#' stored in \code{MegaLMM_state$Posterior} every \code{thin} iteration. If you are doing
#' a long run, and storing a large number of parameters, this will take a lot of memory.
#' This function will estimate the memory requirements.
#' 
#' Note 1: The estimated value will assume all iterations are post-burnin
#' 
#' Note 2: \code{sample_MegaLMM()} will instantiate all arrays to hold the posterior samples
#' prior to running the iterations, so memory requirements will not increase much during the sampling.
#' 
#' Note 3: It is generally not needed to run \code{sample_MegaLMM(MegaLMM_state,n_iter)} 
#' with a large \code{n_iter}. Instead, run the function many times, each with a small \code{n_iter},
#' calling \code{\link{save_posterior_chunk}} between each run. This gives you the ability
#' to diagnose problems during the run, and keeps the memory requirments low. You can always
#' reload the posterior samples from the database on the disk using \code{\link{reload_Posterior}} or
#' \code{\link{load_posterior_param}}.
#'
#' @param MegaLMM_state The model after calling \code{clear_Posterior}
#' @param n_iter number of iterations of the Gibbs sampler
#'
#' @return The estimated memory size in bytes
#' @seealso \code{\link{estimate_memory_initialization_MegaLMM}}, 
#' \code{\link{save_posterior_chunk}}, \code{\link{reload_Posterior}},
#' \code{\link{load_posterior_param}}
#' @export
#'
#' @examples
#' estimate_memory_posterior(MegaLMM_state,100)
estimate_memory_posterior = function(MegaLMM_state,n_iter) {
  Posterior = MegaLMM_state$Posterior
  n_samples = n_iter/MegaLMM_state$run_parameters$thin
  
  estimated_size = 0
  for(param in c(Posterior$posteriorSample_params,names(Posterior$posteriorFunctions))) {
    # 8 bytes per value
    estimated_size = estimated_size + 8*prod(dim(Posterior[[param]])[-1])*n_samples
  }
  for(param in Posterior$posteriorMean_params) {
    estimated_size = estimated_size + 8*prod(dim(Posterior[[param]]))
  }
  
  sprintf('Estimated posterior size for n_samples: %s Gb',
          signif(estimated_size/(1024^3)))
}

#' Make a trace plot of a set of related parameters from MegaLMM
#' 
#' MegaLMM parameters are all matrices. Parameters related to the factor loadings matrix
#' \code{Lambda} all have \code{K} rows. Factors related to the factor scores
#' all have \code{K} columns. Some parameters only have 1 row or column, but are still stored as matrices.
#' Posterior samples of these parameters are arrays with dimensions: \code{Ixnxm} for 
#' \code{I} samples of a parameter with dimensions \code{nxm}.
#' 
#' This function makes trace plots of parameters using a faceting scheme to try to 
#' display as many parameter traces as possible in an organized way. In all plots,
#' the sample number will be the x-axis. The plot(s) will be divided into facets by either the
#' second (default) or 3rd dimension of the \code{sample_array}. Within a facet there will be
#' lines showing the traces of individual values within that dimension of the matrix.
#' In cases where the matrix is large in that dimension, only the \code{n_per_facet}
#' traces with the largest (absolute) posterior mean will be shown, as these are probably the
#' most important.
#' 
#' Sometimes (particularly in a cross-validation framework) we may have a matrix with 
#' many parameters, but only some of them are interesting to inspect traces. You can supply
#' an \code{nxm} logical matrix \code{mask} where \code{TRUE} means the value will be masked
#' from the plot, and \code{FALSE} means the value will be plotted.
#'
#' @param sample_array \code{Ixnxm} array of posterior samples for a matrix of dimension \code{nxm}
#' @param facet_dim either 2 or 3, giving the dimension of \code{sample_array} that should be used as the facets
#' @param name string giving a name to assign to the plot. The file will be a pdf booklet stored in the `run_ID` folder with this name.
#' @param n_per_facet maximum number of values to make traces per facet
#' @param mask optional \code{nxm} logical matrix giving values of the parameter matrix that should NOT be plotted
#' @param include_zero should each facet include the value zero?
#'
#' @return None
#' @export
#'
#' @examples
#' traceplot_array(load_posterior_param(MegaLMM_state,'Lambda'),2,'Lambda')
traceplot_array = function(sample_array,facet_dim = 2,name = 'param.pdf',
                           n_per_facet = 5, mask = NULL, include_zero = TRUE) {
  
  if(length(dim(sample_array)) != 3) stop('sample_array is not an array with 3 dimensions. Note: dimension 1 should be the sampleID')
  if(dim(sample_array)[1] == 0) stop('sample_array is empty')
  if(facet_dim %in% c(2:3) == F) stop('facet_dim should be 2 or 3')
  if(facet_dim == 3) {
    sample_array = aperm(sample_array,c(1,3,2))
    if(!is.null(mask)) mask = t(mask)
  }
  if(!is.null(mask) & !all(dim(mask) == dim(sample_array)[-1])) 
    stop(sprintf('mask has wrong dimension. Should be %s',paste(dim(sample_array)[-1],collapse='x')))
  if(!is.null(mask)) {
    mask[is.na(mask)] = T
    for(i in 1:dim(sample_array)[1]) {
      sample_array[i,,] = sample_array[i,,] * !mask
    }
  }
  # file = sprintf('%s/traceplot_%s.pdf',MegaLMM_state$run_ID,name)
  file = name
  tempfile = sprintf('%s.temp',name)
  pdf(tempfile)
  
  n_facets = dim(sample_array)[2]
  rows = min(4,ceiling(sqrt(n_facets)))
  cols = min(4,ceiling(n_facets / rows))
  # recover()
  par(mfrow=c(rows,cols))
  
  ylim = range(sample_array)
  if(include_zero) {
    if(ylim[1] > 0) ylim[1] = 0
    if(ylim[2] < 0) ylim[2] = 0
  }
  
  for(k in 1:dim(sample_array)[2]){
    sample_matrix_k = matrix(sample_array[,k,],ncol = dim(sample_array)[3])
    # if(!is.null(mask)) sample_matrix_k = sweep(sample_matrix_k,2,!mask[k,],'*')
    o = order(-abs(colMeans(sample_matrix_k)))
    o = o[1:min(n_per_facet,length(o))]
    traces = sample_matrix_k[,o,drop=FALSE]
    trace_plot(traces,main = sprintf('%s %d',name,k),ylim = ylim)
    abline(h=0)
  }
  dev.off()
  system(sprintf('mv %s %s',tempfile,file))
}

#' Create data matrices for MegaLMM from a "tall" data.frame
#' 
#' Often multi-trait data will be provided in tall format with a single value
#' per row, and individuals represented on multiple rows. MegaLMM expects the 
#' data to be in wide format with all observations on the same individual in a row,
#' and only one row per individual.
#' 
#' Specifically, the input to \code{\link{setup_model_MegaLMM}} expect \code{Y}, a 
#' \code{nxp} matrix for \code{n} individuals and \code{p} traits (Missing values coded as \code{NA}
#' are fine in \code{Y}), and \code{data}, a data.frame with the variables used in \code{formula}.
#' This function takes the tall format data and pivots it using the \code{\link{tidyr::pivot_wider}}
#' function, returning the \code{Y} and \code{data} objects necessary to set up a
#' MegaLMM model
#'
#' @param tall_data data.frame in the tall format
#' @param id_cols vector of strings representing the set of columns of \code{tall_data} 
#'  necessary to uniquely specify a single individual. See \code{\link{tidyr::pivot_wider}}.
#' @param names_from vector of strings representing the set of columns of \code{tall_data} 
#'  necessary to uniquely specify a single trait See \code{\link{tidyr::pivot_wider}}.
#' @param values_from string giving the column of \code{tall_data} with the observations.
#' @param ... additional parameters passed to \code{\link{tidyr::pivot_wider}}.
#'
#' @return list with elements:
#' * \code{data} a data.frame with the elements of \code{id_cols}
#' * \code{Y} a numeric matrix of the traits
#' * \code{trait_info} a data.frame with 1 row per trait and columns from \code{names_from}
#' @export
#'
#' @examples
create_data_matrices = function(tall_data,id_cols,names_from,values_from,...) {
  wide_data = tidyr::pivot_wider(tall_data,id_cols = all_of(id_cols), names_from = all_of(names_from),values_from = all_of(values_from))
  
  sample_info = data.frame(wide_data[,id_cols])
  Y = as.matrix(wide_data[,colnames(wide_data) %in% id_cols == F])
  # add ID column to sample_info
  if(max(table(sample_info[[1]])) == 1) {
    rownames(Y) = sample_info[[1]]
  } else {
    id_name = 'ID'
    while(T) {
      if(!id_name %in% colnames(sample_info)) {
        sample_info = cbind(1:nrow(sample_info),sample_info)
        colnames(sample_info)[1] = id_name
        rownames(Y) = sample_info[[id_name]]
        break
      }
      id_name = paste0(id_name,'.')
    }
  }
    
  # now build a table of the trait info
  id_name = 'unique_ID'
  if(id_name %in% colnames(tall_data)) id_name = paste0(id_name,sample(1:1e4,1))
  tall_data[[id_name]] = 1:nrow(tall_data)
  wide_data_id = tidyr::pivot_wider(tall_data,id_cols = all_of(id_cols), names_from = all_of(names_from),values_from = all_of(id_name))
  trait_info = data.frame(Trait = colnames(Y),tall_data[apply(wide_data_id[,colnames(Y)],2,min,na.rm=T),names_from])
  
  return(list(data = sample_info, Y = Y, trait_info = trait_info))
}
