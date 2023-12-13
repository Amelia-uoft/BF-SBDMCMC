# Parameter:
# X: genotype matrix
# weight: coefficients for all RVs within the region
# Z: matrix for confounding factors

bdmcmcbin_empirical  <-  function(X,y,N.max,a0,p0,weight,Z=NULL)   { #p0 is the probability cut-off point for variable selection
  glm_family = binomial(link = "logit")
  n<-dim(X)[1] # sample size
  p<-dim(X)[2] # number of variables, excluding intercept
  X.raw=cbind(rep(1,n),X) # include intercept and covariates
  X.geno = X
  #add error report if genotype and weight dimensions are mismatch (Aug 14, 2023)
  if(ncol(X)!=length(weight)){
    stop("length of weight is different from number of rare variants")
  }
  
  #estimate coefficients of intercept and confounders excluding RVs in the regression
  # Data.original = as.data.frame(cbind(y,Z))
  # names(Data.original)=c("pheno",paste("z",1:ncol(Z),sep=""))
  # 
  # 
  # MLE.original = summary(glm(y ~ .,Data.original,family=binomial))
  # beta.confounder.hat = MLE.original$coefficients[,1]
  # eta.confounder =  matrix(cbind(rep(1,n),Z)) %*% beta.confounder.hat
  
  nei<-matrix(0,N.max+1,p) # matrix to store the results of MCMC iterations
  
  w<-matrix(0,N.max+1,1)  # Posterior probability of the selected model in each MCMC iteration.
  
  temp<-nei[1,]
  
  index <- which(temp>0)   # variable "index" return the varialbes in the starting model, here it is empty as we start from empty model.
  
  if (length(index)==0) {
    
    bf0<-0
    
  } else{
    
    #X1 <- as.matrix(X[,index])
    
    #fit  <-  glm(y~X1,glm_family)  
    #fit0  <-  glm(y~1,glm_family) 
    #bf0<-bayesfactor_models(fit,denominator = fit0)[1,2]
  #  BF.index=c(1,index+1)
  # bf0<-BF.weight(BF.index,X.raw,y,a0,weight)
    if(is.null(Z)==T){
      BF.index=c(1,index+1)
      bf0<-BF.weight(BF.index,X.raw,y,a0,c(0,weight))
    }else{
      bf0<-BF.weight.adj(index,X.geno,y,a0,weight,Z)
    }
  
  }
  
  
  len0 = length(index)
  
  # set=1
  # k1 = (set-1)*set.range+1
  # k2 = set*set.range
  
  k=1
  repeat{
    
    
    # cat("Iteration set", set , "\n")
    # bf.set=NULL
    
      
      cat("Iteration number", k , "\n")
      
      pg <- matrix(0,1,p)
      len <- rep(0,p)
      
      #temp2  <-  foreach (i=1:p, .combine = cbind) %dopar%{
      temp2 = sapply( 1:p,  function(i) {  
        
        temp  =  nei[k,]
        temp[i]  =  abs(temp[i]-1)
        index  =  which(temp > 0)
        
        if (  length(index) == 0 ) {
          
          # fit = lm( y ~ 1 )  
          bf<-0
          
        } else {
          
          #X1 <- as.matrix(X[,index])
          
          #fit  <-  glm(y~X1,glm_family)  
          #fit0  <-  glm(y~1,glm_family) 
          #bf<-bayesfactor_models(fit,denominator = fit0)[1,2]
          #BF.index=c(1,index+1)
          if(is.null(Z)==T){
            BF.index=c(1,index+1)
            bf<-BF.weight(BF.index,X.raw,y,a0,c(0,weight))
          }else{
            bf<-BF.weight.adj(index,X.geno,y,a0,weight,Z)
            
          }
         }
        
        
        #cat("loglik = ", model_logLik,"\n")
        
        
        #pg<-BIC(fit)
        len = length(index)
        
        rbind(bf,len)
      }
      )
      
      bf<-temp2[1,]
      
      #if (max(pg)>5000) {
      # pg=pg/100
      #  pg0=pg0/100
      #}
      
      len<-temp2[2,]
      #pg1<-exp(0.5*(pg0-pg))
      bf1=bf-bf0
      pg1<-exp(bf1) # BF of neighbor models vs. current model
      w[k]<-1/sum(pg1) # BF of current model vs. neighbor models
      pg1<-pg1/sum(pg1)
      sn<-sample(x=c(1:p),size=1,prob=pg1)
      #sn=which(pg1==max(pg1))
      nei[k+1,]<-nei[k,]
      nei[k+1,sn]<-abs(nei[k+1,sn]-1)
      bf0<-bf[sn]
      len0<-length(which(nei[k+1,]>0))
      cat("Bayes factor of selected model: ", sprintf("%10.3f",bf0), "\n" )
      print(which(nei[k+1,]>0))
      
      # bf.set = c(bf.set,bf0)
    
      if(k>1){
        if(k==N.max | length(which((nei[k+1,]==nei[k-1,])==F))==0){break}
        else{
          k=k+1 
        }
      }else{
        k=k+1
      }
   
      
   

    # bf.last.set=bf.set
    
    # set = set+1
    # if(set==N.max/set.range+1) break
    # k1 = (set-1)*set.range+1
    # k2 = set*set.range
  }
  result<-cbind(nei[1:k,] ,w[1:k])
  
  #result<-cbind(nei[(Nburn+1):(N+Nburn),] ,w[(Nburn+1):(N+Nburn)])
  
  result[,p+1]<-result[,p+1]/sum(result[,p+1])
  prob<-NULL
  for(i in 1:p){
    index<-which(result[,i]>0)
    prob[i]<-sum(result[index,p+1])
  }
  
 # final.index = c(1,which(prob>p0)+1)
  final.index = which(prob>p0)
  if(is.null(Z)==T){
    final.index=c(1,final.index+1)
    return(list(VariantProb = prob, VariantSet = final.index,final.BF = BF.weight(final.index,X.raw,y,a0,c(0,weight))))
    
  }else{
    return(list(VariantProb = prob, VariantSet = final.index,final.BF = BF.weight.adj(final.index,X.geno,y,a0,weight,Z)))
    
  }
}
