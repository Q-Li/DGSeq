# two functions of DGseq


# Here we don't assume equal shape parameter for the two conditions. This function doesn't rely on edgeR but requires large sample size. 

GamCAT.II = function (X,n1,n2,M=2000) { # X is the data of one gene, n1 is the number of replicates in condition 1, n2 is the number of replicates in condition2, M is the number of bootstrap samples
  #mle1 = fitdistr(X[1:n1],'Negative Binomial')$estimate # estimate the two parameters of NB from condition 1
  #shape1 = mle1[1]
  #mu1 = mle1[2]
  #mle2 = fitdistr(X[(n1+1):(n1+n2)],'Negative Binomial')$estimate #estimate the two parameters of NB from condition 2
  #shape2=mle2[1]
  #mu2 = mle2[2]
  
  nlogL = function (s,m){
    -sum(dnbinom(X[1:n1],size=s,mu = m,log=T))
  }
  fit1 = mle(minuslogl=nlogL, start = list(s=1, m=floor(mean(X))))
  shape1 = coef(fit1)[1]
  mu1 = coef(fit1)[2]
  
  nlogL = function (s,m){
    -sum(dnbinom(X[(n1+1):(n1+n2)],size=s,mu = m,log=T))
  }
  fit2 = mle(minuslogl=nlogL, start = list(s=1, m=floor(mean(X))))
  shape2 = coef(fit2)[1]
  mu2 = coef(fit2)[2]
  ##### calculate a statistic eta, which indicates the variance of two means, for the original data #########
  beta1 = log(mu1)
  beta2 = log(mu2)
  beta.mean = (beta1+beta2)/2
  eta = (beta1-beta.mean)^2+(beta2-beta.mean)^2
  ############################################################
  
  
  #### assuming Mean_1 = Mean_2, generate M samples under this assumption via parametric bootstrap. #########
  mu = mean(X)
  Gamma1 = matrix(rgamma(n1*M,shape=shape1,rate=shape1/mu),nrow=M,ncol=n1,byrow=T)
  Gamma2 = matrix(rgamma(n2*M,shape=shape2,rate=shape2/mu),nrow=M,ncol=n2,byrow=T)
  ############################################################################################################
  
  
  #### calculate eta for every bootstrap sample##############################################################
  #sim.mu1=apply(Gamma1,1,function (x) fitdistr(x,'Gamma')$estimate[1]/fitdistr(x,'Gamma')$estimate[2])
  #sim.mu2=apply(Gamma2,1,function (x) fitdistr(x,'Gamma')$estimate[1]/fitdistr(x,'Gamma')$estimate[2])
  sim.mu1=rowMeans(Gamma1)
  sim.mu2=rowMeans(Gamma2)
  sim.beta1 = log(sim.mu1)
  sim.beta2 = log(sim.mu2)
  sim.beta.mean=(sim.beta1+sim.beta2)/2
  sim.eta = (sim.beta1-sim.beta.mean)^2+(sim.beta2-sim.beta.mean)^2
  ############################################################################################################
  
  ### find the percentile of original eta
  return(sum(sim.eta>=eta)/M)
}



#### When sample size is small, we reply on edge R to estimate the overdispersion parameter


GamCAT.I.edgeR = function (X,n1,n2,shape,M=2000) { # X is the data of one gene, n1 is the number of replicates in condition 1, n2 is the number of replicates in condition2, M is the number of bootstrap samples. Shape is 
  
  
  mu1 = mean(X[1:n1])
  mu2 = X[(n1+1):(n1+n2)]
  
  ##### calculate a statistic eta, which indicates the variance of two means, for the original data #########
  beta1 = log(mu1)
  beta2 = log(mu2)
  beta.mean = (beta1+beta2)/2
  eta = (beta1-beta.mean)^2+(beta2-beta.mean)^2
  ############################################################
  
  
  #### assuming Mean_1 = Mean_2, generate M samples under this assumption via parametric bootstrap. #########
  mu = mean(X)
  Gamma1 = matrix(rgamma(n1*M,shape=shape,rate=shape/mu),nrow=M,ncol=n1,byrow=T)
  Gamma2 = matrix(rgamma(n2*M,shape=shape,rate=shape/mu),nrow=M,ncol=n2,byrow=T)
  ############################################################################################################
  
  
  #### calculate eta for every bootstrap sample##############################################################
  #sim.mu1=apply(Gamma1,1,function (x) fitdistr(x,'Gamma')$estimate[1]/fitdistr(x,'Gamma')$estimate[2])
  #sim.mu2=apply(Gamma2,1,function (x) fitdistr(x,'Gamma')$estimate[1]/fitdistr(x,'Gamma')$estimate[2])
  sim.mu1=rowMeans(Gamma1)
  sim.mu2=rowMeans(Gamma2)
  sim.beta1 = log(sim.mu1)
  sim.beta2 = log(sim.mu2)
  sim.beta.mean=(sim.beta1+sim.beta2)/2
  sim.eta = (sim.beta1-sim.beta.mean)^2+(sim.beta2-sim.beta.mean)^2
  ############################################################################################################
  
  ### find the percentile of original eta
  return(sum(sim.eta>=eta)/M)
}
