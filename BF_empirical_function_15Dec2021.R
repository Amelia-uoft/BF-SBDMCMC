#replace for loop with sapply function and improve the computation of covariance matrix
BF.weight=function(m.index,X.raw,y,a0,weight){
  
  X=X.raw[,m.index]
  KK=ncol(X)
  y0 =exp(X %*% weight[m.index])/(1+exp(X %*% weight[m.index]))
  
  
  
  # prior
  Data.prior = as.data.frame(cbind(y0,X))
  names(Data.prior)=c("pheno",paste("v",1:KK,sep=""))
  
  MLE.fit = summary(glm(as.numeric(pheno) ~ .-1,Data.prior,family=binomial))
  beta.hat = MLE.fit$coefficients[,1]
  eta.prior =  X %*% beta.hat
  #V.vec.prior=rep(NA,length(eta.prior))
  #for(vv in 1:length(eta.prior)){
  #  V.vec.prior[vv] = exp(eta.prior[vv])/((1+exp(eta.prior[vv]))^2)
  #}
  #v = diag(V.vec.prior)
  #use sapply rather than for loop
  V1.fun = function(eta.val){
     VV=exp(eta.val)/((1+exp(eta.val))^2)
     return(VV)
  }
  V.vec.prior = sapply(eta.prior,V1.fun)
  
  
 # T.hat = t(X) %*% v %*% X
 
  #compute covariance matrix T.hat in an efficient way
  V2.fun = function(x.vec,w){
        x.vec.new = as.matrix(x.vec,nrow=length(x.vec),ncol=1)
        return(x.vec.new %*% w %*% t(x.vec.new))
  }
  listofmatrix.prior = mapply(V2.fun,as.list(data.frame(t(X))),V.vec.prior,SIMPLIFY=F)
  T.hat=Reduce('+',listofmatrix.prior)
  beta.var.prior = a0^(-1)*solve(T.hat)
  
  # posterior
  y.pos = (a0*y0+y)/(a0+1)
  Data.pos = as.data.frame(cbind(y.pos,X))
  names(Data.pos)=c("pheno.pos",paste("v",1:KK,sep=""))
  
  #MLE under large model
  MLE.fit.pos = summary(glm(as.numeric(pheno.pos) ~ .-1,Data.pos,family=binomial))
  beta.hat.pos = MLE.fit.pos$coefficients[,1]
  eta.pos = X %*% beta.hat.pos
  #V.vec.pos=rep(NA,length(eta.pos))
  #for(vv in 1:length(eta.pos)){
  # V.vec.pos[vv] = exp(eta.pos[vv])/((1+exp(eta.pos[vv]))^2)
  #}
  #v.pos = diag(V.vec.pos)
   V.vec.pos = sapply(eta.pos,V1.fun)
   
   listofmatrix.pos = mapply(V2.fun,as.list(data.frame(t(X))),V.vec.pos,SIMPLIFY = F)
  
  T.hat.pos = Reduce('+',listofmatrix.pos)

  
  #T.hat.pos = t(X) %*% v.pos %*% X
  beta.var.pos = (a0+1)^(-1)*solve(T.hat.pos)
  
  y0 = matrix(y0,nrow=n,ncol=1)
  logL_prior_function = function(beta.vec){
    eta.vec.full = X %*% beta.vec
    eta.vec.null = as.matrix(Data.prior$v1) %*% beta.vec[1]
    
    ll.full = 0
    ll.null = 0
    
    # for(i in 1:n){
    #   ll.full = ll.full+y0[i]*eta.vec.full[i]-log(1+exp(eta.vec.full[i]))
    #   ll.null = ll.null+y0[i]*eta.vec.null[i]-log(1+exp(eta.vec.null[i]))
    # 
    # }
    ll.full = t(y0) %*% eta.vec.full-matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.full))
    ll.null = t(y0) %*% eta.vec.null-matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.null))
    #ll=exp(a0*(ll.null-ll.full))
    ll=a0*(ll.null-ll.full)
    
    return(ll)
  }
  
  y.pos = matrix(y.pos,nrow=n,ncol=1)
  logL_post_function = function(beta.vec){
    
    eta.vec.full = X %*% beta.vec
    eta.vec.null = as.matrix(Data.prior$v1) %*% beta.vec[1]
    
    ll.full = 0
    ll.null = 0
    # for(i in 1:n){
    #   ll.full = ll.full+(y[i]+a0*y0[i])*eta.vec.full[i]-(a0+1)*log(1+exp(eta.vec.full[i]))
    #   ll.null = ll.null+(y[i]+a0*y0[i])*eta.vec.null[i]-(a0+1)*log(1+exp(eta.vec.null[i]))
    # }
    
    ll.full = t((a0+1)*y.pos) %*% eta.vec.full-(a0+1)*matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.full))
    ll.null = t((a0+1)*y.pos) %*% eta.vec.null-(a0+1)*matrix(1,nrow=1,ncol=n) %*% log(1+exp(eta.vec.null))
    
    #ll=exp(ll.null-ll.full)
    ll=ll.null-ll.full
    return(ll)
  }
  
  logw_prior_function = function(beta.vec){
    return(dmvnorm(beta.vec,beta.hat,beta.var.prior,log=T)-dnorm(beta.vec[1],beta.hat[1],beta.var.prior[1,1]^0.5,log=T))
  }
  
  logw_post_function = function(beta.vec){
    return(dmvnorm(beta.vec,beta.hat.pos,beta.var.pos,log=T)-dnorm(beta.vec[1],beta.hat.pos[1],beta.var.pos[1,1]^0.5,log=T))
  }
  nMC=100
  set.seed(20)
  
  beta.sim.prior = mvrnorm(nMC,beta.hat,beta.var.prior)
  beta.sim.post = mvrnorm(nMC,beta.hat.pos,beta.var.pos)
  C.ratio.prior = NULL
  
  
  # for(r in 1:nMC){
  #   L_ratio = NULL
  #   w_fun = NULL
  #   L_ratio=L_prior_function(beta.sim.prior[r,])
  #   w_fun = w_prior_function(beta.sim.prior[r,])
  #   check.prior=rbind(check.prior,c(L_ratio,w_fun))
  #   
  #  C.ratio.prior = c(C.ratio.prior,L_ratio*w_fun) # reference model is alternative model
  # }
  
  Fun.ratio.prior = function(beta.vec){
    logL_ratio = NULL
    logw_fun = NULL
    logL_ratio=logL_prior_function(beta.vec)
    logw_fun = logw_prior_function(beta.vec)
    return(logL_ratio+logw_fun)
  }
  
  C.ratio.prior = unlist(lapply(data.frame(t(beta.sim.prior)),Fun.ratio.prior))
  
  C.ratio.post = NULL
  
  # for(r in 1:nMC){
  #   L_ratio = NULL
  #   w_fun = NULL
  #   L_ratio=L_post_function(beta.sim.post[r,])
  #   w_fun = w_post_function(beta.sim.post[r,])
  #   check.post=rbind(check.post,c(L_ratio,w_fun))
  #   C.ratio.post = c(C.ratio.post,L_ratio*w_fun)  # reference model is alternative model
  #  }
  
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
