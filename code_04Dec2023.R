setwd(getwd())

############################
#Step1: Note: use this step only if you want to create your own dataset, otherwise use the existing datafiles in the dataset folder.
## The RV genotype data are simulated based on the 1000 Genomes data for genes given by the user, i.e., needs chromosome (chr) and gene index (gene.index)
#############################
chr = 20 #chromosome used in the simulated data, choose between 1:22
gene.index = 18 #gene index usde in the simulated data
n = 1000 #sample size
c.value = 1 #controls the effect size of the rare variants on the outcome: C=0.5 small; C=0.75 moderate; C=1 strong
w.type = "annotation" #Indicate the type of weights for the RVs: options are "MAF","equal", and "annotation"

install.packages("sim1000G_2.0.5_v2.tar.gz")

library(sim1000G)
packageVersion("sim1000G")
#[1] ‘2.0.5’

readGeneList = function() {
  url = "https://raw.githubusercontent.com/adimitromanolakis/sim1000G-gene-regions/main/ensembl-genes-grch37.csv"
  genes = read.csv(url, as=T)
  genes = genes[ genes$Transcript.type == "protein_coding" , ]
  genes = genes[ genes$Chromosome.scaffold.name %in% as.character(1:22) , ]
  genes
}

if( ! exists("genes") )
  genes = readGeneList()#include genes across whole genome, including column "Gene.name" and "Chromosome.scaffold.name"
genes = genes[which(duplicated(genes$Gene.stable.ID)==F),]
# number of genes on every chromosome
# 1   10   11   12   13   14   15   16   17   18   19    2   20   21   22    3 
# 2067  766 1306 1061  328  645  612  867 1189  285 1444 1276  554  242  448 1070 
# 4    5    6    7    8    9 
# 767  890 1039  915  700  799 

genes = subset(genes,Chromosome.scaffold.name==chr)
geneVCFGitHubLocation = function(gene) {
  
  s = genes[genes$Gene.name == gene,]
  chrom = s$Chromosome.scaffold.name[1]
  
  f = "https://raw.githubusercontent.com/adimitromanolakis/sim1000G-gene-regions/main/ceu-tsi-gbr/chr%s/genes-chr%s-%s.vcf.gz"
  f = sprintf(f,chrom,chrom,gene)
  
  cat("-> Gene location for download is",f, "\n")
  f
}

genename = genes$Gene.name[gene.index]
vcf.file = readVCF(geneVCFGitHubLocation(genename),maxNumberOfVariants = 1e6 , min_maf = 1e-6 ,max_maf = 0.01)

set.seed(100)
startSimulation(vcf.file, totalNumberOfIndividuals = n)
ids <- generateUnrelatedIndividuals(n)


## retrieve all the genotypes from individuals
gt = retrieveGenotypes( ids )
#apply(gt,2,sum,na.rm=T)
reverse.gt = function(xx){
  if(sum(xx==2)>n/2){
    xx=2-xx
  }
  return(xx)
}
gt = apply(gt,2,reverse.gt)#some SNP, we need to switch minor and major alleles
convert.gt = function(xx){
  xx = ifelse(xx==2,1,xx)
  return(xx)
}
X.covariate = matrix(unlist(lapply(gt,convert.gt)),n,ncol(gt))

# remove highly correlated variants, r2>0.99
var.maf0 = which(apply(X.covariate,2,sum,na.rm=T)==0)

if(length(var.maf0)>0){
  X.covariate = X.covariate[,-var.maf0]
}

cor.mat = cor(X.covariate)
rm.list = NULL
for(rr in 1:ncol(X.covariate)){
  if(rr %in% rm.list) next
  
  for(cc in setdiff(1:ncol(X.covariate),1:rr)){
    
    if(cc %in% rm.list) next
    
    if(cor.mat[rr,cc]>0.99){
      print(c(rr,cc))
      rm.list = c(rm.list,cc)
    }
    
  }
  
}
if(is.null(rm.list)==F){
  X.covariate = X.covariate[,-rm.list]
}

########### choose causal variants that are independent ##############
cor.mat.new = cor(X.covariate)
diag(cor.mat.new) = rep(NA,nrow(cor.mat.new))
SNP.ind.set = which(apply(cor.mat.new,1,max,na.rm=T)<0.2)

set.seed(10)
causal = sample(SNP.ind.set,5) #choose 5 causal RVs which are independent with other RVs (r2<0.2)

if(c.value>0){
  
  freq = apply(X.covariate,2,sum)/(2*nrow(X.covariate))
  p1 = length(causal)
  
  beta.sign = rep(1,p1)
  beta.abs = c.value*abs(log10(freq[causal]))
  beta.val = beta.sign*beta.abs
  
  eta = beta.val %*% t(X.covariate[,causal])
  prob = exp(eta)/(1+exp(eta))
}else{
  p1=0
}  

#simulate phenotype
set.seed(100)
y=rep(NA,n)
if(c.value>0){
  for(i in 1:n){
    y[i] = rbinom(1,1,prob[1,i])
  }
}else{
  for(i in 1:n){
    y[i] = rbinom(1,1,0.5)
  }
}

#simulate weights
if(w.type=="equal"){
  weight.est = rep(0.1,ncol(X.covariate))
  
}else if(w.type=="MAF"){
  maf.vec=apply(X.covariate,2,sum,na.rm=T)/(2*nrow(X.covariate))
  weight.est = -0.05*log10(maf.vec)
  
}else{
  weight.est = rgamma(ncol(X.covariate),1,3)#
  weight.est = ifelse(weight.est>2,2,weight.est)#
  weight.est[causal] = rnorm(length(causal),1.5,0.1)#
  
}
#simulate covariates "cov"
set.seed(100)
cov=rep(NA,n)
for(i in 1:n){
  if(y[i]==0){
    cov[i]=rbinom(1,1,0.2)
  }else{
    cov[i]=rbinom(1,1,0.3)
  }
}


write.table(X.covariate,"./Dataset/Genotype.txt",row.names=F,col.names=F)
write.table(as.matrix(y),"./Dataset/Phenotype.txt",row.names=F,col.names=F)
write.table(as.matrix(weight.est),"./Dataset/Weights.txt",row.names=F,col.names=F)
write.table(as.matrix(cov),"./Dataset/Confounders.txt",row.names=F,col.names=F)

######################################
### Step 2: run analyses with BF #####
######################################
library(mvtnorm)
library(MASS)
setwd(getwd())
#Source code for the unadjusted BF-BDMCMC
source("BF_empirical_function_15Dec2021.R")

#Source code for the adjusted BF-BDMCMC
source("BF_empirical_function_adj_14Aug2023.R")

#Source code to run the SBDMCMC algorithm
source("bdmcmc_emprical_BF_adj_01Dec2023.R")

a0=1 #Indicates the extent in prior belief; a0=1 means the prior has the same weight has the data; a0=0.5 means the prior has half the weight of the data (less prior belief)

X.covariate=read.table("./Dataset/Genotype.txt",head=F)
y=read.table("./Dataset/Phenotype.txt",head=F)[,1]
weight.est=read.table("./Dataset/Weights.txt",head=F)[,1]
cov=read.table("./Dataset/Confounders.txt",head=F)
n=nrow(X.covariate) #sample size

# bdmcmcbin_empirical  <-  function(X,y,N.max,a0,p0,weight,Z=NULL)  
# Arguments for the bdmcmc function
# X: genotype matrix
# y: phenotype, 0 or 1
# N.max: number of iterations for the SBDMCMC algorithm
# a0: Indicates the extent in prior belief; a0=1 means the prior has the same weight has the data; a0=0.5 means the prior has half the weight of the data (less prior belief)
# p0: Posterior Inclusion Probability cutoff; only RVs with PIP>p0 are included into the BF test statistic
# weight: Weights for all RVs (e.g., annotation score)
# Z: matrix for confounding factors


##############################
# Model1: no adjustment of confounder
###########################
fit.unadjust = bdmcmcbin_empirical(as.matrix(X.covariate),y,30,a0,0.5,weight.est) #BF-SBDMCMC not adjusted for confounding factors
variant.select = fit.unadjust$VariantSet
print(variant.select[-1]-1)#selected RVs
#print(causal)#causal RVs if using step 1 to generate data
logBF.unadjust = fit.unadjust$final.BF
logBF.unadjust

#PIP plot
plot(fit.unadjust$VariantProb,  pch=18,frame = FALSE,
     col = "orange", xlab = "RV index", ylab = "posterior inclusion probability",ylim=c(0,1))

###############
# We assume that three quartiles are used to define risk levels of RVs
# Modifier, if weight.est<q1
# Low risk, if q1<=weight.est<q2
# Moderate risk, if q2<=weight.est<q3
# high risk, if weight.est>=q3
###############
q1=quantile(weight.est,0.25)
q2=quantile(weight.est,0.5)
q3=quantile(weight.est,0.75)

select.index=which(weight.est>=q3)
points(select.index,fit.unadjust$VariantProb[select.index],col="red",pch=18)
select.index=which(weight.est>=q1 & weight.est<q2)
points(select.index,fit.unadjust$VariantProb[select.index],col="purple",pch=18)
select.index=which(weight.est<q1)
points(select.index,fit.unadjust$VariantProb[select.index],col="blue",pch=18)

abline(h=0.5,lty=2)
legend("topleft", legend=c("High risk","Moderate risk", "Low risk","Modifier"),pch=rep(18,4),
       col=c("red","orange", "purple","blue"),  cex=0.4)

######################
#Model 2: adjusted for confounder
######################
fit.adjust = bdmcmcbin_empirical(as.matrix(X.covariate),y,30,a0,0.5,weight.est,as.matrix(cov))#BF-SBDMCMC adjusted for confounding factors

variant.select = fit.adjust$VariantSet 
print(variant.select[-1]-1)#selected RVs
#print(causal)#causal RVs if using step 1 to generate data
logBF.adjust = fit.adjust$final.BF
logBF.adjust

#PIP plot for adjusted model
plot(fit.adjust$VariantProb,  pch=18,frame = FALSE,
     col = "orange", xlab = "RV index", ylab = "posterior inclusion probability",ylim=c(0,1))

select.index=which(weight.est>=q3)
points(select.index,fit.adjust$VariantProb[select.index],col="red",pch=18)
select.index=which(weight.est>=q1 & weight.est<q2)
points(select.index,fit.adjust$VariantProb[select.index],col="purple",pch=18)
select.index=which(weight.est<q1)
points(select.index,fit.adjust$VariantProb[select.index],col="blue",pch=18)

abline(h=0.5,lty=2)
legend("topleft", legend=c("High risk","Moderate risk", "Low risk","Modifier"),pch=rep(18,4),
       col=c("red","orange", "purple","blue"),  cex=0.4)


###############################################################
#Step3: Assess significance of the BF using permutation testing
###############################################################

replicates=20
library(snow)
library(iterators)
library(foreach)
library(doSNOW)
cl <- makeCluster(8, type="SOCK") #define how many jobs you want the computer to run at the same time, define number of cores to use

clusterEvalQ(cl, library(mvtnorm))
clusterEvalQ(cl,library(MASS))
registerDoSNOW(cl) #use cluster cl

res.BF <- foreach(ii=1:replicates, .packages=NULL,.combine=rbind,.inorder=FALSE) %dopar% { 
  set.seed(ii)
  
  #generate phenotype y under H0
  y=rep(NA,n)
  for(i in 1:n){
    y[i] = rbinom(1,1,0.5)
  }
  
  fit.unadjust = bdmcmcbin_empirical(as.matrix(X.covariate),y,30,a0,0.5,weight.est)
  fit.adjust = bdmcmcbin_empirical(as.matrix(X.covariate),y,30,a0,0.5,weight.est,as.matrix(cov))
  
  c(iteration=ii,
    logBF.unadjust = fit.unadjust$final.BF,
    logBF.adjust = fit.adjust$final.BF
  )
}
head(res.BF)
#compute permutation p-value
p.unadjust = length(which(res.BF[,2]>logBF.unadjust))/nrow(res.BF)
p.adjust = length(which(res.BF[,3]>logBF.adjust))/nrow(res.BF)
