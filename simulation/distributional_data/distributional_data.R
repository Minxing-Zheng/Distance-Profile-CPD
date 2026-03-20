# rm(list=ls())
require(gtools)
library(boot)
library(magrittr)
library(MASS)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)
library(mvtnorm)
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##### new package
require("gSeg")
require(Rcpp)
require("kerSeg") 
library(optparse)
source("../../functions/ecp_distmat_input.R")
sourceCpp('../../functions/energyChangePoint.cpp')
sourceCpp("../../functions/depth_CPDcpp.cpp")
sourceCpp("../../functions/depth_CPDcpp_ALL.cpp")
source("../../functions/gen_data.R")
Sigma_list <- make_Sigma_list(c(30, 90, 180))
names(Sigma_list) <- c("S1", "S2", "S3")
S1 <- Sigma_list$S1
S2 <- Sigma_list$S2
S3 <- Sigma_list$S3
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(min(cores[1],11)-1) #not to overload your computer
registerDoParallel(cl=cl,cores=min(cores[1],11)-1) 
detectCores()

option_list <- list(
  make_option(c("--num_permut"), type = "integer", default = 1, 
              help = "Number of permutations for the permutation test [default: %default]"),
  make_option(c("--delta"), type = "numeric", default = 0.5, 
              help = "Signal strength for change points [default: %default]."),
  make_option(c("--monte_carlo"), type = "numeric", default = 1, 
              help = "Index of Monte Carlo simulations [default: %default]."),
  make_option(c("--random_seed"), type = "numeric", default = 1, 
              help = "Random seed Monte Carlo simulations [default: %default]."),
  make_option(c("--change_type"), type = "numeric", default = 1, 
              help = "type of change points: 1,2,3,4 for mean, var, mixture, tail changes [default: %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))
type<-opt$change_type
num_permut<-opt$num_permut
delta<-opt$delta
monte_carlo<-opt$monte_carlo
rand_seed<-opt$random_seed
set.seed(rand_seed)

n1=100;n2=200;n=n1+n2
I<-diag(1,2)

if (type==1){
  # mean shift: delta in seq(0,1,0.1), 
  z1<-mvrnorm(n=n1,mu=c(0,0),Sigma = 0.25*I)
  z2<-mvrnorm(n=n2,mu=c(delta,0),Sigma = 0.25*I)
  
}else if(type==2){
  # scale shift: delta in seq(0,0.4,0.04), 
  z1<-mvrnorm(n=n1,mu=c(0,0),Sigma = 0.16*I)
  z2<-mvrnorm(n=n2,mu=c(0,0),Sigma = diag(c((0.4+delta)**2,0.4**2)))
}

df<-expand.grid(seq(-3,3,len=100),seq(-3,3,len=100))
df_matrix<-data.matrix(df)
l<-lapply(seq_len(nrow(df_matrix)), function(i) df_matrix[i,])
Data<-c()
for (i in seq_len(nrow(rbind(z1,z2)))){
  mu<-rbind(z1,z2)[i,]
  result<-sapply(l,pmvnorm,lower=-Inf,mean=mu,sigma=I)
  Data<-rbind(Data,result)
}
inc<-seq(-3,3,len=100)[2]-seq(-3,3,len=100)[1]
distmat<-as.matrix(dist(Data,method = "manhattan")*inc)

result_dist_profile<-foreach (i = (0:num_permut),.noexport=c('depth_CPD_cpp','depth_CPDcpp_ALL')) %dopar%{
  require(Rcpp)

  sourceCpp("../../functions/depth_CPDcpp.cpp")
  sourceCpp("../../functions/depth_CPDcpp_ALL.cpp")
  if(i==0){
    distmat_temp = distmat
  }else{
    ind<-sample(nrow(distmat))
    distmat_temp<-distmat[ind,ind]
  }
  res_all = depth_CPDcpp_ALL(distmat_temp,num_permut = 0) # dF, AD, W
  res_list <- lapply(1:3, function(i) {
    list(loc = res_all$loc[i], observed_test_statistics = res_all$observed_stat[i])
  })
  res1=depth_CPD_cpp(distmat_temp,num_permut = 0)
  res2 <- res_list[[1]]
  res3 <- res_list[[2]]
  res4 <- res_list[[3]]
  l = list(res1,res2,res3,res4)
  names(l)<-c('dist_cpd_uniform','dist_cpd','dist_cpd_AD','dist_cpd_W')
  l  
}
#### Graph CPD
E1 = mstree(dist( Data, method = 'euclidean' ),ngmax = 5)
result_graph = gseg1(nrow(distmat),E1, statistics="all",B = num_permut,pval.perm=TRUE)
#### Energy CPD
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=num_permut,k=NULL,min.size=50,alpha=1)
#### MMD test
result_MMD<-c()
for (j in 0:num_permut){
  if (j!=0){
    # set.seed(j)
    ind<-sample(nrow(distmat))
    D<-distmat[ind,ind]
  }else{
    D<-distmat
  }
  result_mmd<-MMD_test(D)
  loc<-which.max(result_mmd)
  ob_stat<-result_mmd[loc]
  result_MMD<-rbind(result_MMD,c(loc,ob_stat))
}
r<-list(result_dist_profile,result_graph,result_ecp,result_MMD)
names(r)<-c("dist_profile",'graph','ecp','mmd')
result<-list()
result[[1]]=r
type_name<-c('mean','var')
path<-paste("distr_",type_name[type],'_delta_',delta,'_run_',monte_carlo,'.Rdata',sep="")
save(result, file=path)


