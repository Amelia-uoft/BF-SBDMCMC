BF.weight.adj=function(m.index,X.geno,y,a0,weight,Z.mat){
 # print(m.index)
  X=matrix(X.geno[,m.index],nrow=n,ncol=length(m.index))
  KK=ncol(X)
  y0 =exp(X %*% weight[m.index])/(1+exp(X %*% weight[m.index]))#estimate of y0 only depends on genotypes rather than confounders
  
  JJ=ncol(Z.mat)
  X.full=cbind(rep(1,n),Z.mat,X)
  X.cov=cbind(rep(1,n),Z.mat)
  
  # prior
  Data.prior = as.data.frame(cbind(y0,Z.mat,X))
  
  names(Data.prior)=c("pheno",paste("z",1:JJ,sep=""),paste("v",1:KK,sep=""))
  
  MLE.fit = summary(glm(as.numeric(pheno) ~ .,Data.prior,family=binomial))
  beta.hat = MLE.fit$coefficients[,1]
  eta.prior =  X.full %*% beta.hat
  V.vec.prior=rep(NA,length(eta.prior))
  for(vv in 1:length(eta.prior)){
    V.vec.prior[vv] = exp(eta.prior[vv])/((1+exp(eta.prior[vv]))^2)
  }
  v = diag(V.vec.prior)
  
  T.hat = t(X.full) %*% v %*% X.full
  beta.var.prior = a0^(-1)*solve(T.hat)
  
  # posterior
  y.pos = (a0*y0+y)/(a0+1)
  Data.pos = as.data.frame(cbind(y.pos,Z.mat,X))
  names(Data.pos)=c("pheno.pos",paste("z",1:JJ,sep=""),paste("v",1:KK,sep=""))
  
  #MLE under large model
  MLE.fit.pos = summary(glm(as.numeric(pheno.pos) ~ .,Data.pos,family=binomial))
  beta.hat.pos = MLE.fit.pos$coefficients[,1]
  eta.pos = X.full %*% beta.hat.pos
  V.vec.pos=rep(NA,length(eta.pos))
  for(vv in 1:length(eta.pos)){
    V.vec.pos[vv] = exp(eta.pos[vv])/((1+exp(eta.pos[vv]))^2)
  }
  v.pos = diag(V.vec.pos)
  
  T.hat.pos = t(X.full) %*% v.pos %*% X.full
  beta.var.pos = (a0+1)^(-1)*solve(T.hat.pos)
  
  y0 = matrix(y0,nrow=n,ncol=1)
  logL_prior_function = function(beta.vec){
    eta.vec.full = X.full %*% beta.vec
    eta.vec.null = X.cov %*% beta.vec[1:(JJ+1)]
    
    ll.full = 0
    ll.null = 0
    
    ll.full = t(y0) %*% eta.vec.full-matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.full))
    ll.null = t(y0) %*% eta.vec.null-matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.null))
    #ll=exp(a0*(ll.null-ll.full))
    ll=a0*(ll.null-ll.full)
    
    return(ll)
  }
  
  y.pos = matrix(y.pos,nrow=n,ncol=1)
  logL_post_function = function(beta.vec){
    
    eta.vec.full = X.full %*% beta.vec
    eta.vec.null = X.cov %*% beta.vec[1:(JJ+1)]
    
    ll.full = 0
    ll.null = 0
   
    ll.full = t((a0+1)*y.pos) %*% eta.vec.full-(a0+1)*matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.full))
    ll.null = t((a0+1)*y.pos) %*% eta.vec.null-(a0+1)*matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.null))
    
    #ll=exp(ll.null-ll.full)
    ll=ll.null-ll.full
    return(ll)
  }
  
  logw_prior_function = function(beta.vec){
    return(dmvnorm(beta.vec,beta.hat,beta.var.prior,log=T)-dmvnorm(beta.vec[1:(JJ+1)],beta.hat[1:(JJ+1)],beta.var.prior[1:(JJ+1),1:(JJ+1)],log=T))
  }
  
  logw_post_function = function(beta.vec){
    return(dmvnorm(beta.vec,beta.hat.pos,beta.var.pos,log=T)-dmvnorm(beta.vec[1:(JJ+1)],beta.hat.pos[1:(JJ+1)],beta.var.pos[1:(JJ+1),1:(JJ+1)],log=T))
  }
  nMC=100
  set.seed(20)
  
  beta.sim.prior = mvrnorm(nMC,beta.hat,beta.var.prior)
  beta.sim.post = mvrnorm(nMC,beta.hat.pos,beta.var.pos)
  C.ratio.prior = NULL
  
  Fun.ratio.prior = function(beta.vec){
    logL_ratio = NULL
    logw_fun = NULL
    logL_ratio=logL_prior_function(beta.vec)
    logw_fun = logw_prior_function(beta.vec)
    return(logL_ratio+logw_fun)
  }
  
  C.ratio.prior = unlist(lapply(data.frame(t(beta.sim.prior)),Fun.ratio.prior))
  
  C.ratio.post = NULL
  
 
  
  Fun.ratio.post = function(beta.vec){
    logL_ratio = NULL
    logw_fun = NULL
    logL_ratio=logL_post_function(beta.vec)
    logw_fun = logw_post_function(beta.vec)
    return(logL_ratio+logw_fun)
  }
  
  C.ratio.post = unlist(lapply(data.frame(t(beta.sim.post)),Fun.ratio.post))
  
  logBF =mean(C.ratio.prior)-mean(C.ratio.post) #compare alternative model to null model
  return(round(logBF,2))
}
