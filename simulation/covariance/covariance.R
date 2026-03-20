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
library(ecp)
library(optparse)
library(pracma)
source("../../functions/objTest_fctns.R")
source("../../functions/depth_CPD_func.R") 
source("../../functions/ecp_distmat_input.R")
source("../../functions/kcp_distmat_input.R")
sourceCpp('../../functions/energyChangePoint.cpp')
sourceCpp("../../functions/getTcpp.cpp")
sourceCpp("../../functions/depth_CPDcpp.cpp")
sourceCpp("../../functions/depth_CPDcpp_ALL.cpp")

library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(min(cores[1],11)-1) #not to overload your computer
registerDoParallel(cl=cl,cores=min(cores[1],11)-1) 
detectCores()

gen_cov<-function(r,d1,d2){
  topleft = r*matrix(1,d1,d1)+(1-r)*diag(d1)
  top=cbind(topleft,matrix(0,d1,d2))
  bottom=cbind(matrix(0,d2,d1),diag(d2))
  return(rbind(top,bottom))
}

option_list <- list(
  make_option(c("--num_permut"), type = "integer", default = 1, 
              help = "Number of permutations for the permutation test [default: %default]"),
  make_option(c("--delta"), type = "numeric", default = 0.5, 
              help = "Signal strength for change points [default: %default]."),
  make_option(c("--rho"), type = "numeric", default = 0.5, 
              help = "rho [default: %default]."),
  make_option(c("--monte_carlo"), type = "numeric", default = 1, 
              help = "Index of Monte Carlo simulations [default: %default]."),
  make_option(c("--random_seed"), type = "numeric", default = 1, 
              help = "Random seed Monte Carlo simulations [default: %default]."),
  make_option(c("--change_type"), type = "numeric", default = 1, 
              help = "type of change points: 1,2,3,4 for mean, var, mixture, tail changes [default: %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))
type<-opt$change_type
B<-opt$num_permut
delta<-opt$delta # 0 in setting 1, 0 to 2 in setting 2. 
rho<-opt$rho # 0 to 1.5 in setting 1, rho = 0.5 in setting 2 
monte_carlo<-opt$monte_carlo
rand_seed<-opt$random_seed
set.seed(rand_seed)

#setting 1 - difference in covariance structures - identity to spiked model (try both exponentiated and non-exponentiated)
if (type==1){
  #  set delta = 0, vary rho from 0 to 1.5 
  U=randortho(20,type = "orthonormal")
  S1=U%*%t(U)
  S2=U%*%diag(c(rep(1+rho,5),rep(1,15)))%*%t(U)
  mu1=rep(0,20)
  mu2=c(rep(delta,5),rep(0,15))
}else if(type==2){
  # delta from 0 to 2
  # rho=0.5
  U=randortho(20,type = "orthonormal")
  S1=U%*%t(U)
  S2=U%*%diag(c(1.5,rep(1+rho,4),rep(1,15)))%*%t(U)
  mu1=U[,1]
  mu2=delta*U[,1]
}

X1=mvrnorm(n=50,mu=mu1,Sigma = S1)
X2=mvrnorm(n=100,mu=mu2,Sigma = S2 )
X=rbind(X1,X2)
distmat=as.matrix(dist(X))

result_dist_profile<-foreach (i = (0:num_permut),.noexport=c('depth_CPD_cpp','depth_CPDcpp_ALL')) %dopar%{
  require(Rcpp)
  sourceCpp("../../functions/getTcpp.cpp")
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
result_graph = gseg1(nrow(distmat),E1, statistics="g",B = num_permut,pval.perm=TRUE)
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
result<-r

type_name<-c('diag','mean')
path<-paste('covariance_',type_name[type],'_delta_',delta,'_run_',monte_carlo,'.Rdata',sep="")
save(result, file=path)
