set_device = function(device) {
  # create new devices up to device
  devices = dev.list()
  devices = devices[names(devices) != 'RStudioGD']
  while(length(devices) < device){
    dev.new(noRStudioGD = TRUE)
    devices = dev.list()
    devices = devices[names(devices) != 'RStudioGD']
  }
  dev.set(devices[device])
}

trace_plot = function(data,main = NULL,ylim = NULL){
  if(is.null(ylim)) {
    range_y = range(data)
    max_range = max(abs(range_y))
    ylim = c(-max_range,max_range)
  }
  plot(NA,NA,xlim = c(0,nrow(data)),ylim = ylim,main= main,xlab = "iteration")
  for(i in 1:ncol(data)) lines(data[,i],col=i)
}

trace_plot_h2s = function(F_h2_samples, n_factors = 8, device = NULL){
  if(!is.null(device)) {
    set_device(device)
  }

  rows = ceiling(n_factors / 2)
  cols = 2*ceiling(n_factors / rows)
  par(mfrow=c(rows,cols))
  n_h2 = dim(F_h2_samples)[2]
  h2_names = dimnames(F_h2_samples)[[2]]

  # for(k in 1:min(n_factors,dim(F_h2_samples)[3])){
  for(k in 1:dim(F_h2_samples)[3]){
    trace_plot(matrix(F_h2_samples[,,k],ncol = n_h2),main = sprintf('Factor %d h2s',k),ylim = c(0,1))
    hist(F_h2_samples[,1,k],breaks=seq(0,1,length=100),xlim = c(-0.1,1),main = sprintf('Factor %d',k))
    if(n_h2 > 1){
      for(i in 2:n_h2){
        hist(F_h2_samples[,i,k],breaks=seq(0,1,length=100),col = i,border=i,add=T)
      }
    }
  }
}

trace_plot_Lambda = function(Lambda, n_factors = 16, device = NULL,main = 'Lambda'){
  if(!is.null(device)) {
    set_device(device)
  }

  rows = ceiling(sqrt(n_factors))
  cols = ceiling(n_factors / rows)
  par(mfrow=c(rows,cols))

  # for(k in 1:min(n_factors,dim(Lambda)[3])){
  for(k in 1:dim(Lambda)[2]){
    o = order(-abs(colMeans(Lambda[,k,])))
    o = o[1:min(5,length(o))]
    traces = Lambda[,k,o]
    trace_plot(traces,main = sprintf('Factor %d %s',k,main))
    abline(h=0)
  }
}

boxplot_Bs = function(B,main = ''){
  k = dim(B)[3]
  b = dim(B)[2]
  par(mfrow=c(3,3))
  for(i in 1:min(27,k)) {
    data = B[,,i]
    ylim = range(c(0,data))
    boxplot(data,ylim=ylim,main = sprintf('%s %d',main,i),outpch=NA,boxlty=0,whisklty=1,staplelty=0)
    abline(h=0)
  }
}

plot_factor_correlations = function(Lambda, sim_Lambda){
  cors = abs(cor(t(Lambda),t(sim_Lambda)))
  cors = rbind(cors,apply(cors,2,max))
  cors = cbind(cors,apply(cors,1,max))
  image(t(cors)[,nrow(cors):1],xlab = 'Actual factors', ylab = 'Estimated factors',xaxt='n',yaxt='n', main = 'Correlation of fitted and\nsimulated factor loadings'
  )
  axis(1,at=seq(0,1,length=ncol(cors)),labels = c(1:(ncol(cors)-1),'best'),las=2)
  axis(2,at=seq(0,1,length=nrow(cors)),labels = c('best',(nrow(cors)-1):1),las=2)
}

plot_element_wise_covariances = function(actual_G, estimated_G, main = NULL){
  actual_G = actual_G[upper.tri(actual_G)]
  estimated_G = estimated_G[upper.tri(estimated_G)]
  xlim = ylim = range(c(actual_G,estimated_G))
  i = 1:length(actual_G)
  i = sample(i,min(10000,length(i)))
  # i = order(-c(abs(actual_G)))[1:1000]
  plot(actual_G[i],estimated_G[i],xlim = xlim,ylim=ylim,main = main,xlab = 'actual',ylab = 'estimated')
  abline(0,1)
}


plot_diagonal_covariances = function(actual_G, estimated_G, main = NULL){
  actual_G = diag(actual_G)
  estimated_G = diag(estimated_G)
  xlim = ylim = range(c(actual_G,estimated_G))
  i = 1:length(actual_G)
  i = sample(i,min(10000,length(i)))
  # i = order(-c(abs(actual_G)))[1:1000]
  plot(actual_G[i],estimated_G[i],xlim = xlim,ylim=ylim,main = main,xlab = 'actual',ylab = 'estimated')
  abline(0,1)
}

plot_current_state_simulation = function(MegaLMM_state, device = NULL){
  if(!is.null(device)) {
    set_device(device)
  }
  par(mfrow=c(3,3))

  setup = MegaLMM_state$setup
  run_parameters = MegaLMM_state$run_parameters
  run_variables = MegaLMM_state$run_variables

  current_state = within(MegaLMM_state$current_state,{
    U_R = as.matrix(MegaLMM_state$data_matrices$RE_L %*% U_R)
    U_F = as.matrix(MegaLMM_state$data_matrices$RE_L %*% U_F)
    # transform variables so that the variance of each column of F is 1.
    F_var = 1/tot_F_prec
    U_F = sweep(U_F,2,sqrt(F_var),'/')
    B2_F = sweep(B2_F,2,sqrt(F_var),'/')
    F = sweep(F,2,sqrt(F_var),'/')
    Lambda = sweep(Lambda,1,sqrt(F_var),'*')

    if(!'var_Eta' %in% ls()) var_Eta = rep(1,ncol(Lambda))
    U_R[] = sweep(U_R,2,sqrt(var_Eta),'*')
    B1[] = sweep(B1,2,sqrt(var_Eta),'*')
    B2_R[] = sweep(B2_R,2,sqrt(var_Eta),'*')
    Lambda[] = sweep(Lambda,2,sqrt(var_Eta),'*')
    tot_Eta_prec[] = tot_Eta_prec / var_Eta
  })

  Lambda = current_state$Lambda
  F_h2 = current_state$F_h2
  if(is.null(dim(F_h2))) F_h2 = matrix(F_h2,nrow=1)
  U_R_prec = current_state$tot_Eta_prec / current_state$resid_h2
  resid_Eta_prec = current_state$tot_Eta_prec / (1-current_state$resid_h2)
  p = run_variables$p


  # correlation of factors
  sim_Lambda = setup$Lambda
  plot_factor_correlations(Lambda, sim_Lambda)

  # element-wise correlations
  G_plots = lapply(1:nrow(F_h2),function(re) {
    G_est = crossprod(sweep(Lambda,1,sqrt(F_h2[re,]),'*')) + diag(c(current_state$resid_h2[re,]/current_state$tot_Eta_prec))
    G_act = setup$G
    RE_name = rownames(F_h2)[re]
    if(is.null(RE_name)) RE_name = 'Va'
    plot_element_wise_covariances(G_act, G_est, main = sprintf('G: %s elements',RE_name))
    plot_diagonal_covariances(G_act, G_est, main = sprintf('G: %s diagonal',RE_name))
  })
  U_R_est = crossprod(sweep(Lambda,1,sqrt(1-colSums(F_h2)),'*')) + diag(c((1-colSums(current_state$resid_h2))/current_state$tot_Eta_prec))
  U_R_act = setup$R
  plot_element_wise_covariances(U_R_act, U_R_est,main = 'E elements')
  plot_diagonal_covariances(U_R_act, U_R_est,main = 'E diagonal')

  # factor h2s
  plot_factor_h2s(F_h2)

  # B's
  if(dim(setup$B)[1] > 1) {
    B2_R = MegaLMM_state$current_state$B2_R
    B2_F = MegaLMM_state$current_state$B2_F %*% MegaLMM_state$current_state$Lambda
    try({plot(c(B2_R),c(setup$B))},silent = TRUE)
    abline(0,1)
    if(dim(B2_F)[1] > 1) {
      try({plot(c(B2_F),c(setup$B_F %*% setup$Lambda))},silent = TRUE)
      abline(0,1)
      xlim = ylim = range(c(B2_R,B2_F))
      try({plot(c(B2_F),c(B2_R),xlim = xlim,ylim=ylim);abline(0,1)},silent = TRUE)
    }
  }

  # cis_effects
  if(length(setup$cis_effects) > 0 &&  length(setup$cis_effects) == length(MegaLMM_state$current_state$cis_effects)){
    try({plot(setup$cis_effects,MegaLMM_state$current_state$cis_effects);abline(0,1)},silent = TRUE)
  }
}

plot_fixed_effects = function(B_act, B_resid,B_factor,B_total){
  xlim = ylim = range(unlist(c(B_act, B_resid,B_factor,B_total)))
  plot(c(B_act),c(B_resid),xlim=xlim,ylim=ylim,col=1,main = 'Fixed effects')
  points(c(B_act),c(B_factor),col=2)
  points(c(B_act),c(B_total),pch=21,col=0,bg=2)
  abline(0,1)
}

calc_posterior_mean_cov = function(Posterior,random_effect){
  with(Posterior,{
    p = dim(Lambda)[3]
    G = matrix(0,p,p)
    for(i in 1:total_samples){
      if(random_effect == 'Ve'){
        factor_h2s_i = 1-colSums(F_h2[i,,,drop=FALSE])
        resid_h2s_i = 1-colSums(resid_h2[i,,,drop=FALSE])
      } else{
        factor_h2s_i = F_h2[i,random_effect,]
        resid_h2s_i = resid_h2[i,random_effect,]
      }
      G_i = crossprod(sweep(Lambda[i,,],1,sqrt(factor_h2s_i),'*')) + diag(c(resid_h2s_i/tot_Eta_prec[i,1,]))
      G = G + G_i/total_samples
    }
    G
  })
}

calc_posterior_mean_Lambda = function(Posterior){
  return(with(Posterior,apply(Lambda,c(2,3),mean)))
}

plot_factor_h2s = function(F_h2) {
  if(is.null(rownames(F_h2))) rownames(F_h2) = 'Va'
  F_h2 = rbind(F_h2,1-colSums(F_h2))
  rownames(F_h2)[nrow(F_h2)] = 'Ve'
  colnames(F_h2) = 1:ncol(F_h2)
  F_h2_data = melt(F_h2)
  barplot(F_h2,main = 'Factor h2s',legend.text = rownames(F_h2))
}


plot_posterior_simulation = function(MegaLMM_state, device = NULL){
  if(!is.null(device)) {
    set_device(device)
  }
  par(mfrow=c(3,3))

  setup = MegaLMM_state$setup
  Posterior = MegaLMM_state$Posterior
  p = dim(Posterior$Lambda)[3]
  estimated_R = calc_posterior_mean_cov(Posterior,'Ve')
  plot_element_wise_covariances(setup$R, estimated_R,main = 'E elements')
  plot_diagonal_covariances(setup$R, estimated_R,main = 'E diagonal')
  for(random_effect in 1:dim(Posterior$F_h2)[2]) {
    estimated_G = calc_posterior_mean_cov(Posterior,random_effect)
    plot_element_wise_covariances(setup$G, estimated_G,main = sprintf('G: %s elements',random_effect))
    plot_diagonal_covariances(setup$G, estimated_G,main = sprintf('G: %s diagonal',random_effect))
  }

  plot_factor_correlations(calc_posterior_mean_Lambda(Posterior),setup$Lambda)

  plot_factor_h2s(apply(Posterior$F_h2,c(2,3),mean))

  if(dim(setup$B)[1] > 1) {
    B_mean = apply(Posterior$B1,c(2,3),mean)
    try({plot(c(B_mean),c(setup$B))},silent = TRUE)
    abline(0,1)
    if(!is.null(setup$B_F) & ncol(MegaLMM_state$data_matrices$X2_F) == nrow(setup$B_F)){
      B_factor_mean = with(c(Posterior,MegaLMM_state$data_matrices), {
        if(ncol(X2_F) == 0) return(rep(0,dim(Lambda)[1]))
        matrix(rowMeans(sapply(1:total_samples,function(i) B2_F[i,,] %*% Lambda[i,,])),nrow = ncol(X2_F))
      })
      try({plot(c(B_factor_mean),c(setup$B_F %*% setup$Lambda));abline(0,1)},silent = TRUE)
      xlim = ylim = range(c(B_mean[-1,],B_factor_mean))
      try({plot(c(B_factor_mean),c(B_mean),xlim = xlim,ylim=ylim);abline(0,1)},silent = TRUE)
      B_f_HPD = HPDinterval(mcmc(Posterior$B2_F[,1,]))
      B_f_mean = colMeans(Posterior$B2_F[,1,])

      try({
        plot(1:length(B_f_mean),B_f_mean,xlim = c(1,length(B_f_mean)),ylim = range(B_f_HPD),xlab = '',main = 'Posterior B2_F')
        arrows(seq_along(B_f_mean),B_f_HPD[,1],seq_along(B_f_mean),B_f_HPD[,2],length=0)
        abline(h=0)
      },silent = TRUE)
    }
  }
}


summarize_posterior = function(X) {
  X = coda::mcmc(X)
  Xdata = data.frame(ID = 1:ncol(X),Mean = colMeans(X),Median = apply(X,2,median))
  Xi = coda::HPDinterval(X,prob = 0.95)
  Xdata$low_95 = Xi[,1]
  Xdata$high_95 = Xi[,2]
  Xi = coda::HPDinterval(X,prob = 0.8)
  Xdata$low_80 = Xi[,1]
  Xdata$high_80 = Xi[,2]
  return(Xdata)
}
posterior_plot = function(X,xlab='',ylab='',colorSig = T,ylim = NULL,colorGroup = NULL) {
  require(cowplot)
  Xdata = summarize_posterior(X)
  if(!is.null(colorGroup)) {
    Xdata$color = colorGroup
    if(colorSig) {
      Xdata$color[sign(Xdata$low_95) == sign(Xdata$high_95)] = 'Significant'
      Xdata$color = factor(Xdata$color,levels = c('Significant',unique(colorGroup)))
    }
  } else{
    Xdata$color = 'NS'
    if(colorSig) {
      Xdata$color[sign(Xdata$low_95) == sign(Xdata$high_95)] = 'Sig'
    }
  }
  if(is.null(ylim)) ylim = range(Xdata[,4:5])

  p = ggplot(Xdata,aes(x=ID)) + geom_hline(yintercept = 0) +
    xlab(xlab) + ylab(ylab) + ylim(ylim)+
    geom_segment(aes(xend = ID,y = low_95,yend=high_95,color=color),size=.5) +
    geom_segment(aes(xend = ID,y = low_80,yend=high_80,color=color),size = .9) +
    geom_point(aes(y=Median,color = color)) +
    theme(legend.position = 'none')
  p
}

plot_diagnostics_simulation = function(MegaLMM_state){
  MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
  plot_current_state_simulation(MegaLMM_state)
  if(MegaLMM_state$Posterior$total_samples > 0) {
    plot_posterior_simulation(MegaLMM_state)
    trace_plot_h2s(MegaLMM_state$Posterior$F_h2)
    trace_plot_Lambda(MegaLMM_state$Posterior$Lambda)
  }
}


#' Plots diagnostic plots for the fit of a MegaLMM model
#'
#' Currently, only traceplots of F_h2 and elements of Lambda are shown
#'
#' @param MegaLMM_state a MegaLMM_state object
plot_diagnostics = function(MegaLMM_state){
  if(MegaLMM_state$Posterior$total_samples > 0) {
    trace_plot_h2s(load_posterior_param(MegaLMM_state,'F_h2'))
    trace_plot_Lambda(load_posterior_param(MegaLMM_state,'Lambda'))
    try({trace_plot_Lambda(load_posterior_param(MegaLMM_state,'B2_F'),main='B2_F')},silent=T)
    try({
      B2_R = load_posterior_param(MegaLMM_state,'B2_R')[,-1,]
      trace_plot_Lambda(B2_R,main='B2_R')},silent=T)
    # try({boxplot_Bs(load_posterior_param(MegaLMM_state,'B'),'B')},silent=T)
    try({
      B1 = load_posterior_param(MegaLMM_state,'B1')[,-1,]
      trace_plot_Lambda(B1,main='B1')},silent=T)

}
}



