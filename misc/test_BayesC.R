library(microbenchmark)
n = 1000
p = 500

X = matrix(sample(c(0,1),n*p,replace=T),n,p)

beta_act = c(1,1,1,rep(0,p-3))
y = X %*% beta_act + rnorm(n)


## Gibbs sampler
df = 4
scale = 1
df_alpha = 4
scale_alpha = 1

nIter = 1000
thin = 100
beta_samples = matrix(0,nIter,p)
var_samples = matrix(0,nIter,3)
# invVarRes_samples = rep(NA,nIter)
#init
beta = beta_act#rep(0,p)
delta = rep(1,p)
alpha = beta
invVarRes = 1
varEffects = rep(1,p)
pi = 1-3/1000
yCorr = y - X %*% beta
diag_XtX = colSums(X^2)
for(i in 1:nIter){
  for(j in 1:thin) {
    invVarRes = regression_sampler_v4(X,diag_XtX,yCorr,alpha,beta,delta,invVarRes,varEffects,pi,df,scale)
    varEffects[] = (sum(alpha^2) + df_alpha*scale_alpha)/rchisq(1,nLoci+df_alpha)
    pi = rbeta(1,p-nLoci+1,nLoci+1)
  }
  nLoci = sum(delta)
  print(nLoci)
  var_samples[i,] = c(invVarRes,varEffects[1],pi)
  beta_samples[i,] = beta
  # invVarRes_samples[i] = invVarRes
}
