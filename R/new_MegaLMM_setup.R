

#' Initializes a MegaLMM model
#' 
#' This function initializes a basic MegaLMM model with a single matrix of traits.
#' The return in an object of class \code{"MegaLMM_state"} with the structures
#' to hold all aspects of the model.
#' 
#' The user will be able to add additional trait matrices (with reduced sets of individuals)
#' and create a regression model for Lambda using additional functions
#'
#' @param Y trait matrix (n x t). Can be null in which case all traits are added in 
#' additional trait blocks afterwards (e.g. if no individuals are observed for all traits)
#' @param formula formula specifying model for each column of \code{"Y"} 
#' and each column of \code{"F"}. Must include at least 1 random effect.
#' Multiple correlated random effects are not allowed (e.g. \code{"(1+z|x)"}). 
#' Please use \code{"(1+z||x)}"} instead. 
#' If no random effects are desired, include \code{"(1|ID)"} where \code{"ID"}
#' is a column of \code{"data"} with a unique level for each row of "Y"
#' @param data data.frame with variables used to construct design matrices. 
#' Must have the same number of rows as \code{"Y"} if \code{"Y"} is non-null.
#' @param rowID column in \code{"data"} specifying a unique identifier for each \
#' row in \code{"Y"}. Only required for adding additional trait matrices
#' @param relmat 
#' @param output_folder 
#'
#' @return
#' @export
#'
#' @examples
MegaLMM_model = function(Y,formula,data,rowID,relmat,output_folder = 'MegaLMM_model',
                         X2 = NULL, U2 = NULL, V2 = NULL) {
  try(dir.create(output_folder,recursive = T, showWarnings = FALSE),silent=T)
 
  model_setup = make_model_setup(formula,data,relmat)

  # X1 is the "fixed effects", un-shrunk covariates that only affect residuals
  X1 = model_setup$lmod$X
  b1 = ncol(X1)

  # -------- regressions ---------- #
  # if(is.null(X2)) X2 = matrix(0,n,0)
  if(!is.null(X2)) {
    # project out intercept
    if(nrow(X2) != n) stop("nrow(X2) != nrow(data)")
    M = diag(1,n) - matrix(1/n,n,n)
    X2 = M %*% X2
    X2 = X2[,colSums(abs(X2))>1e-10,drop=FALSE]
    b2 = ncol(X2)
    if(!is.null(U2) || !is.null(V2)) stop("only one of X2 or (U2/V2) should be provided")
  }

  if(!is.null(U2)) {
    if(is.null(V2)) stop("V2 must be provided if U2 is provided")
    if(nrow(U2) != n) stop("nrow(U2) != nrow(data)")
    m = nrow(V2)
    if(ncol(U2) != m) stop("U2 and V2 are not compatible")
    # project out intercept
    M = diag(1,m) - matrix(1/m,m,m)
    V2 = M %*% V2
    V2 = V2[,colSums(abs(V2))>1e-10,drop=FALSE]
    b2 = ncol(V2)
  }


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
              K = as(diag(ldl_k$d[large_d]),'dgCMatrix')
              K_inv = as(diag(1/ldl_k$d[large_d]),'dgCMatrix')
              L = t(ldl_k$P) %*% ldl_k$L[,large_d]
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


  data_matrices = list(
    X1          = X1,
    X2          = X2,
    U2          = U2,
    V2          = V2,
    Z           = Z,
    ZL          = ZL,
    ZL_list     = list(ZL),
    RE_setup    = RE_setup,
    RE_L_list     = list(RE_L),  # List of matrices necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
    RE_indices  = RE_indices,
    h2s_matrix  = h2s_matrix,
    Lambda_fixed = Lambda_fixed,
    data       = data
  )

  run_parameters = list(
    K = K,
    n = n,
    r_RE   = r_RE,
    RE_names = RE_names,
    b1       = b1,
    b2_R = b2_R,
    b2_F = b2_F
  )

  return(list(
    data_matrices = data_matrices,
    run_parameters = run_parameters
  ))

}
  
  
   
#   MegaLMM_state = list(
#     current_state = list(),
#     run_parameters = list(),
#     design_matrices = list(),
#     transformed_matrices = list(),
#     priors = list(),
#     Posterior = list(),
#   )
#   class(MegaLMM_state) = append('MegaLMM_state',class(MegaLMM_state))
#   return(MegaLMM_state)
# }


#' 
#' #' Create F model
#' #'
#' #' Species the model for the factors. 
#' #' \code{"F"} will be learned for all observations in \code{"data"}.
#' #' \code{"formula"} must contain at least 1 random effect term (e.g. \code{"(1|x)"}).
#' #' Multiple correlated random effects are not allowed (e.g. \code{"(1+z|x)"}) 
#' #' and will be converted to independent random effects (e.g. \code{"(1+z||x)}"}).
#' #' 
#' #' @param formula 
#' #' @param data 
#' #' @param relmat 
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' F_model = function(K, formula,data,relmat = NULL, h2_prior = NULL, X2 = NULL, U2 = NULL, V2 = NULL, B2_prior = NULL, F_var_prior = NULL) {
#'   
#'   # ---------------------------------- #
#'   # -------- model matrices ---------- #
#'   # ---------------------------------- #
#'   
#'   model_setup = make_model_setup(formula,data,relmat)
#'   
#'   # X1 is the "fixed effects", un-shrunk covariates that only affect residuals
#'   X1 = model_setup$lmod$X
#'   b1 = ncol(X1)
#'   
#'   # -------- regressions ---------- #
#'   # if(is.null(X2)) X2 = matrix(0,n,0)
#'   if(!is.null(X2)) {
#'     # project out intercept
#'     if(nrow(X2) != n) stop("nrow(X2) != nrow(data)")
#'     M = diag(1,n) - matrix(1/n,n,n)
#'     X2 = M %*% X2
#'     X2 = X2[,colSums(abs(X2))>1e-10,drop=FALSE]
#'     b2 = ncol(X2)
#'     if(!is.null(U2) || !is.null(V2)) stop("only one of X2 or (U2/V2) should be provided")
#'   }
#'   
#'   if(!is.null(U2)) {
#'     if(is.null(V2)) stop("V2 must be provided if U2 is provided")
#'     if(nrow(U2) != n) stop("nrow(U2) != nrow(data)")
#'     m = nrow(V2)
#'     if(ncol(U2) != m) stop("U2 and V2 are not compatible")
#'     # project out intercept
#'     M = diag(1,m) - matrix(1/m,m,m)
#'     V2 = M %*% V2
#'     V2 = V2[,colSums(abs(V2))>1e-10,drop=FALSE]
#'     b2 = ncol(V2)
#'   }
#'   
#'   
#'   # -------- Random effects ---------- #
#'   
#'   # RE_setup is a list like from BGLR that describes the random effects (for both residuals and factors)
#'   RE_setup = model_setup$RE_setup
#'   
#'   # add names to RE_setup if needed
#'   n_RE = length(RE_setup)
#'   for(i in 1:n_RE){
#'     if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
#'       names(RE_setup)[i] = paste0('RE.',i)
#'     }
#'   }
#'   RE_names = names(RE_setup)
#'   
#'   # combine Z matrices
#'   Z = do.call(cbind,lapply(RE_setup,function(x) x$Z))
#'   Z = as(Z,'dgCMatrix')
#'   
#'   
#'   # find RE indices
#'   RE_lengths = sapply(RE_setup,function(x) ncol(x$Z))
#'   RE_starts = cumsum(c(0,RE_lengths)[1:n_RE])
#'   names(RE_starts) = RE_names
#'   RE_indices = lapply(RE_names,function(re) RE_starts[re] + 1:RE_lengths[re])
#'   names(RE_indices) = RE_names
#'   
#'   # function to ensure that covariance matrices are sparse and symmetric
#'   fix_K = function(x) forceSymmetric(drop0(x,tol = run_parameters$drop0_tol))
#'   
#'   # construct RE_L and ZL for each random effect
#'   for(i in 1:length(RE_setup)){
#'     re_name = names(RE_setup)[i]
#'     RE_setup[[i]] = within(RE_setup[[i]],{
#'       # recover()
#'       if(!'ZL' %in% ls()){
#'         if('K' %in% ls() && !is.null(K)){
#'           id_names = rownames(K)
#'           if(is(K,'Matrix') & isDiagonal(K)) {
#'             L = as(diag(1,nrow(K)),'dgCMatrix')
#'             K_inv = as(diag(1/diag(K)),'dgCMatrix')
#'           } else {
#'             ldl_k = LDLt(K)
#'             large_d = ldl_k$d > run_parameters$K_eigen_tol
#'             r_eff = sum(large_d)
#'             # if need to use reduced rank model, then use D of K in place of K and merge L into Z
#'             # otherwise, use original K, set L = Diagonal(1,r)
#'             if(r_eff < length(ldl_k$d)) {
#'               K = as(diag(ldl_k$d[large_d]),'dgCMatrix')
#'               K_inv = as(diag(1/ldl_k$d[large_d]),'dgCMatrix')
#'               L = t(ldl_k$P) %*% ldl_k$L[,large_d]
#'               if(is(L,'dgeMatrix')) L = as.matrix(L)
#'             } else{
#'               L = as(diag(1,nrow(K)),'dgCMatrix')
#'               K_inv = as(with(ldl_k,t(P) %*% crossprod(diag(1/sqrt(d)) %*% solve(L)) %*% P),'dgCMatrix')
#'             }
#'             rm(list=c('ldl_k','large_d','r_eff'))
#'           }
#'           if(is.null(rownames(K))) rownames(K) = 1:nrow(K)
#'           rownames(K_inv) = rownames(K)
#'         } else if ('K_inv' %in% ls() && !is.null(K_inv)){
#'           id_names = rownames(K_inv)
#'           if(is.null(rownames(K_inv))) rownames(K_inv) = 1:nrow(K_inv)
#'           K = solve(K_inv)
#'           rownames(K) = rownames(K_inv)
#'           L = as(diag(1,nrow(K)),'dgCMatrix')
#'         } else{
#'           K = as(diag(1,ncol(Z)),'dgCMatrix')
#'           rownames(K) = colnames(Z)
#'           id_names = rownames(K)
#'           K_inv = K
#'           L = as(diag(1,nrow(K)),'dgCMatrix')
#'         }
#'         # if(is.null(id_names)) id_names = 1:length(id_names)
#'         rownames(L) = paste(id_names,re_name,sep='::')
#'         K = fix_K(K)
#'         ZL = Z %*% L
#'       }
#'     })
#'   }
#'   
#'   
#'   ZL = do.call(cbind,lapply(RE_setup,function(re) re$ZL))
#'   if(nnzero(ZL)/length(ZL) < .25) {
#'     ZL = as(ZL,'dgCMatrix')
#'   } else{
#'     ZL = as.matrix(ZL)
#'   }
#'   
#'   
#'   if(length(RE_setup) > 1) {
#'     RE_L = do.call(bdiag,lapply(RE_setup,function(re) re$L))
#'     rownames(RE_L) = do.call(c,lapply(RE_setup,function(re) rownames(re$L)))
#'   } else{
#'     RE_L = RE_setup[[1]]$L
#'   }
#'   if(nnzero(RE_L)/length(RE_L) < 0.25) {
#'     RE_L = as(RE_L,'dgCMatrix')
#'   } else{
#'     RE_L = as.matrix(RE_L)
#'   }
#'   r_RE = sapply(RE_setup,function(re) ncol(re$ZL))
#'   
#'   
#'   h2_divisions = run_parameters$h2_divisions
#'   if(length(h2_divisions) < n_RE){
#'     if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
#'     h2_divisions = rep(h2_divisions,n_RE)
#'   }
#'   if(is.null(names(h2_divisions))) {
#'     names(h2_divisions) = RE_names
#'   }
#'   h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
#'   colnames(h2s_matrix) = RE_names
#'   h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
#'   colnames(h2s_matrix) = NULL
#'   
#'   
#'   data_matrices = list(
#'     X1          = X1,
#'     X2          = X2,
#'     U2          = U2,
#'     V2          = V2,
#'     Z           = Z,
#'     ZL          = ZL,
#'     ZL_list     = list(ZL),
#'     RE_setup    = RE_setup,
#'     RE_L_list     = list(RE_L),  # List of matrices necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
#'     RE_indices  = RE_indices,
#'     h2s_matrix  = h2s_matrix,
#'     Lambda_fixed = Lambda_fixed,
#'     data       = data
#'   )
#'   
#'   run_parameters = list(
#'     K = K,
#'     n = n,
#'     r_RE   = r_RE,
#'     RE_names = RE_names,
#'     b1       = b1,
#'     b2_R = b2_R,
#'     b2_F = b2_F
#'   )
#'   
#'   return(list(
#'     data_matrices = data_matrices,
#'     run_parameters = run_parameters
#'   ))
#'   
#' }
#' 
#' #' Create Eta model
#' #' 
#' #' Species the model for a set of traits. 
#' #' Any missing values in \code{"Eta"} will be imputed as missing parameters
#' #' Rows of \code{"data"} must correspond to rows of Eta.
#' #' \code{"formula"} must contain at least 1 random effect term (e.g. \code{"(1|x)"}).
#' #' The special term \code{"FLambda"} will be added if missing and corresponds to the 
#' #' latent factors that are components of \code{"Eta"}. 
#' #' Multiple correlated random effects are not allowed (e.g. \code{"(1+z|x)"}) 
#' #' and will be converted to independent random effects (e.g. \code{"(1+z||x)}"}).
#' #' 
#' #'
#' #' @param Eta 
#' #' @param formula 
#' #' @param data 
#' #'
#' #' @return List with design matrices for the model for this \code{"Eta"}. 
#' #' These can be combined to merge multiple \code{"Eta"}'s into a joint model.
#' #' 
#' #' @export
#' #'
#' #' @examples
#' Eta_model = function(Eta,formula,data,relmat = NULL, h2_prior = NULL, observation_model = 'identity', X2 = NULL, U2 = NULL, V2 = NULL, B2_prior = NULL, Eta_var_prior = NULL) {
#'   
#'   
#'   n = nrow(Eta)
#'   p = ncol(Eta)
#'   
#'   if(nrow(data) != n) stop("nrow(data) != nrow(Eta)")
#'   
#'   model_list = F_model(p,formula,data,relmat,h2_prior,X2,U2,V2,B2_prior,Eta_var_prior)
#'   
#' }
#' 
#' combine_Eta_models = function(Eta_model_list) {
#'   
#' }
#' 
#' 
#' Lambda_model = function(formula,data,Lambda_var_prior, relmat = NULL, h2_prior = NULL, X2 = NULL, B2_prior = NULL) {
#'   
#' }