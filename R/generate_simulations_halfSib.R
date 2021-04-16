# Copyright 2020 Daniel Runcie
# Use of this source code is governed by the PolyForm Noncommercial License 1.0.0
# that can be found in the LICENSE file and available at
# https://polyformproject.org/licenses/noncommercial/1.0.0/

z_transform = function(X) {
  X <- sweep(X,2,colMeans(X),'-')
  X <- sweep(X,2,apply(X,2,sd),'/')
}
family_simulation = function(name = 'simulation_1', family_sizes,within_group_cor = .25,p, b, factor_h2s, prop_Factor_p = 0.5,total_h2 = 0.5,Vb = 0, numeff = NULL){
  nTot = sum(family_sizes)
  Sire = as.factor(do.call(c,lapply(1:length(family_sizes),function(x) rep(x,family_sizes[x]))))
  K = within_group_cor*tcrossprod(Matrix(model.matrix(~0+Sire))) + diag(1-within_group_cor,length(Sire))
  # K = tcrossprod(matrix(rnorm(nTot^2,nTot))/nTot/2)+diag(1,nTot)
  rownames(K) = 1:nrow(K)

  # for more informative K
  # K[K>0 & K<1] = 1

  # for low-rank K
  # rownames_K = rownames(K)
  # K = tcrossprod(matrix(rnorm(nrow(K)*100),ncol=100))
  # sd_K = sqrt(diag(K))
  # K = t(K / sd_K)/sd_K
  # rownames(K) = rownames_K

  # in case K is low-rank
  K_chol = chol(K+diag(1e-6,nrow(K)))


  # Lambda matrix
  # factor_h2s = rep(0,k)
  # factor_h2s[2+1:k_G] = runif(k_G)
  k = length(factor_h2s)
  Lambda = matrix(0,k,p)
  if(is.null(numeff)) numeff = sort(sample((p/30):(p/4),k,replace=T),decreasing = TRUE)
  for(h in 1:k){
    Lambda[h,sample(1:p,numeff[h])] = 1 #rnorm(numeff[h])
  }
  # cols=1:k
  # g_cols = factor_h2s>0
  # if(sum(g_cols) == 0) g_cols = 1:k
  # Lambda = Lambda[do.call("order", unname(split(-abs(Lambda[cols[g_cols],drop=FALSE]), col(Lambda[cols[g_cols],drop=FALSE])))),]

  V_tot = rep(1,p)
  V_factor = runif(p)
  Lambda_scale = colSums(Lambda^2) / V_factor
  Lambda_scale[Lambda_scale == 0] = 1
  Lambda = sweep(Lambda,2,sqrt(Lambda_scale),'/')
  V_resid = V_tot - colSums(Lambda^2)

  # fixed effect design
  n = nrow(K)
  X = cbind(1,z_transform(matrix(rnorm(n*(b-1)),nrow=n)))
  colnames_X = 'intercept'
  if(b > 1) colnames_X = c(colnames_X,paste0('Fixed',1:max(1,(b-1))))
  colnames(X) = colnames_X
  X_F = X[,-1,drop=FALSE]

  B_F = rstdnorm_mat(ncol(X_F),k) * sqrt(Vb)
  B   = matrix(0,ncol(X),p)

  # RE design
  data = droplevels(data.frame(X,Sire=Sire,animal = rownames(K)))
  Z = diag(1,n)


  U_F = sweep(t(K_chol) %*% z_transform(rstdnorm_mat(nTot,k)),2,sqrt(factor_h2s),'*')
  E_F = sweep(z_transform(rstdnorm_mat(nTot,k)),2,sqrt(1-factor_h2s),'*')
  F = as.matrix(X_F %*% B_F + Z %*% U_F + E_F)
  F = z_transform(F)

  Vg_factor = colSums(sweep(Lambda^2,1,factor_h2s,'*'))
  Vg_resid = pmax(0,rep(total_h2,p) - Vg_factor) * V_resid

  U_R = z_transform(t(K_chol) %*% rstdnorm_mat(nTot,p))
  E_R = z_transform(rstdnorm_mat(nTot,p))
  U_R = sweep(U_R,2,sqrt(Vg_resid),'*')
  E_R = sweep(E_R,2,sqrt(V_resid-Vg_resid),'*')
  Y = X %*% B + F %*% Lambda + U_R + E_R

  G = crossprod(sqrt(factor_h2s)*Lambda) + diag(Vg_resid)
  R = crossprod(sqrt(1-factor_h2s)*Lambda) + diag(V_resid-Vg_resid)

  # recover()
  Y = as.matrix(Y)
  colnames(Y) = paste0('gene',1:p)


  setup = list(
    Y = Y,
    data = data,
    K = K,
    B = B,
    B_F = B_F,
    Lambda = Lambda,
    h2 = diag(G)/(diag(G)+diag(R)),
    factor_h2s = factor_h2s,
    G = G,
    R = R,
    X = X,
    X_F = X_F,
    U_F = U_F,
    U_R = U_R,
    F = F,
    # E_F = E_F,
    # E_R = E_R,
    name = name
  )
  return(setup)
  # save(setup,file='setup.RData')
}

new_halfSib_simulation = function(name, nSire,nRep,p, b, factor_h2s, Va = 0.2, Ve = 0.2,Vb = 0, numeff = NULL){
  nTot = mean(nRep) * nSire
  nRep = rep(nRep,nSire/length(nRep))
  Sire = as.factor(do.call(c,lapply(1:nSire,function(x) rep(x,nRep[x]))))
  # K = .25*tcrossprod(Matrix(model.matrix(~0+Sire))) + diag(.75,length(Sire))
  K = tcrossprod(matrix(rnorm(nTot^2,nTot))/nTot/2)+diag(1,nTot)
  rownames(K) = 1:nrow(K)

  # for more informative K
  # K[K>0 & K<1] = 1

  # for low-rank K
  # rownames_K = rownames(K)
  # K = tcrossprod(matrix(rnorm(nrow(K)*100),ncol=100))
  # sd_K = sqrt(diag(K))
  # K = t(K / sd_K)/sd_K
  # rownames(K) = rownames_K

  # in case K is low-rank
  K_chol = chol(K+diag(1e-6,nrow(K)))

  # Lambda matrix
  # factor_h2s = rep(0,k)
  # factor_h2s[2+1:k_G] = runif(k_G)
  k = length(factor_h2s)
  Lambda = matrix(0,k,p)
  if(is.null(numeff)) numeff = sample((p/30):(p/4),k,replace=T)
  for(h in 1:k){
    Lambda[h,sample(1:p,numeff[h])] = rnorm(numeff[h])
  }
  Lambda = Lambda[order(-diag(Lambda %*% t(Lambda))),,drop=FALSE]
  cols=1:k
  g_cols = factor_h2s>0
  if(sum(g_cols) == 0) g_cols = 1:k
  Lambda = Lambda[,do.call("order", unname(split(-abs(Lambda[cols[g_cols],,drop=FALSE]), row(Lambda[cols[g_cols],,drop=FALSE])))),drop=FALSE]

  # resid variances
  # tot_Y_prec = rep(0.5,p)
  # resid_h2 = runif(p,0,0.6)
  resid_Va = rep(Va,p)
  resid_Ve = rep(Ve,p)
  resid_h2 = resid_Va/(resid_Va+resid_Ve)
  tot_Y_prec = 1/(resid_Va+resid_Ve)

  # G, R matrices
  G = crossprod(sweep(Lambda,1,sqrt(factor_h2s),'*')) + diag(resid_h2/tot_Y_prec)
  R = crossprod(sweep(Lambda,1,sqrt(1-factor_h2s),'*')) + diag((1-resid_h2)/tot_Y_prec)

  # fixed effect design
  n = nrow(K)
  X = cbind(1,matrix(rnorm(n*(b-1)),nrow=n))
  colnames_X = 'intercept'
  if(b > 1) colnames_X = c(colnames_X,paste0('Fixed',1:max(1,(b-1))))
  colnames(X) = colnames_X
  X_F = X[,-1]

  # RE design
  data = droplevels(data.frame(X,Sire=Sire,animal = rownames(K)))
  Z = diag(1,n)
  # Z = model.matrix(~0+Sire,data)
  # r = ncol(Z)

  B_F = matrix(rnorm((b-1)*k),nrow = (b-1),ncol=k) * sqrt(Vb)
  B = rbind(rnorm(p),matrix(0,(b-1),p)) * sqrt(Vb)

  U_F = t(K_chol) %*% matrix(rnorm(n*k,0,sqrt(factor_h2s)),n,k,byrow=T)
  E_F = matrix(rnorm(n*k,0,sqrt(1-factor_h2s)),n,k,byrow=T)
  F = X_F %*% B_F + Z %*% U_F + E_F

  U_R = t(K_chol) %*% matrix(rnorm(n*p,0,sqrt(resid_h2/tot_Y_prec)),n,p,byrow=T)
  E_R = matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)
  Y = X %*% B + F %*% Lambda + U_R + E_R

  # recover()
  Y = as.matrix(Y)
  colnames(Y) = paste0('gene',1:p)


  setup = list(
    Y = Y,
    data = data,
    K = K,
    B = B,
    B_F = B_F,
    Lambda = Lambda,
    h2 = resid_h2,
    factor_h2s = factor_h2s,
    G = G,
    R = R,
    X = X,
    X_F = X_F,
    U_F = U_F,
    U_R = U_R,
    F = F,
    # E_F = E_F,
    # E_R = E_R,
    name = name
  )
  save(setup,file='setup.RData')
}

