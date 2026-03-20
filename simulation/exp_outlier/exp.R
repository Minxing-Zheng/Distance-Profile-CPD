### expo

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
source("../../functions/ecp_distmat_input.R")
source("../../functions/kcp_distmat_input.R")
sourceCpp('../../functions/energyChangePoint.cpp')
sourceCpp("../../functions/getTcpp.cpp")
sourceCpp("../../functions/depth_CPDcpp.cpp")
sourceCpp("../../functions/getTcpp_dF.cpp")
sourceCpp("../../functions/depth_CPDcpp_dF.cpp")
sourceCpp("../../functions/getTcpp_dF_AD.cpp")
sourceCpp("../../functions/depth_CPDcpp_dF_AD.cpp")
sourceCpp("../../functions/getTcpp_dF_W.cpp")
sourceCpp("../../functions/depth_CPDcpp_dF_W.cpp")

library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(min(cores[1],11)-1) #not to overload your computer
registerDoParallel(cl=cl,cores=min(cores[1],11)-1) 

MMD_test<-function(D){
  
  bw<-median(D)
  gram_matrix<-exp(D^2*(-1/(2*bw^2))) 
  T<-nrow(D)
  
  obj_vals<-matrix(0,T)
  for (t in seq(2,T-1)){
    term1 = (T-t)/(t*T)*sum(gram_matrix[1:t, 1:t])
    term2 = -2/T*sum(gram_matrix[1:t, t:T])
    term3 = t/((T-t)*T)*sum(gram_matrix[t:T, t:T])
    obj_vals[t] = term1 + term2 + term3
  }
  return(obj_vals)
}


# gen_data<-function(total_length,cp_loc,type){
#   
# }

frac_exp<-as.numeric(commandArgs(TRUE)[1])
delta<-as.numeric(commandArgs(TRUE)[2])
monte_carlo<-as.numeric(commandArgs(TRUE)[3])

S1<-get(load("../../Sigma_30.Rdata"))
S2<-get(load("../../Sigma_90.Rdata"))
S3<-get(load("../../Sigma_180.Rdata"))
n1=100;n2=200;n=n1+n2
num_permut<-300
i<-2
cat("seed",as.integer(10000*monte_carlo+100*frac_exp))


p=c(30,90,180)[i] # for mean_diff, scale_diff, mixture gaussian cases
#p=c(5,15,60)[i] # for heavy tail case
I<-diag(x = 1, p, p)

#### 1. mean diff
Sigma=list(S1,S2,S3)[[i]]
n_exponentiate = floor(frac_exp * p)
signal_1<-rep(0,p)
signal_2<-rep(delta,p)
Data<-rbind(mvrnorm(n1,mu=signal_1,Sigma=Sigma),mvrnorm(n2,mu=signal_2,Sigma=Sigma))

if (n_exponentiate==0){
}else{
  exp_coordinate<-seq(1,n_exponentiate)
  cat("exp_cor:",exp_coordinate,'\n')
  Data[,exp_coordinate]<-exp(Data[,exp_coordinate])
}
cat(n_exponentiate)
distmat<-as.matrix(dist( Data, method = 'euclidean' ))

result_dist_profile<-foreach (i = (0:num_permut),.noexport=c('depth_CPD_cpp','depth_CPD_cpp_dF','depth_CPD_cpp_dF_AD','depth_CPD_cpp_dF_W')) %dopar%{
  library(foreach)
  library(doParallel)
  require(Rcpp)
  sourceCpp("../../functions/getTcpp.cpp")
  sourceCpp("../../functions/depth_CPDcpp.cpp")
  sourceCpp("../../functions/getTcpp_dF.cpp")
  sourceCpp("../../functions/depth_CPDcpp_dF.cpp")
  sourceCpp("../../functions/getTcpp_dF_AD.cpp")
  sourceCpp("../../functions/depth_CPDcpp_dF_AD.cpp")
  sourceCpp("../../functions/getTcpp_dF_W.cpp")
  sourceCpp("../../functions/depth_CPDcpp_dF_W.cpp")
  if(i==0){
    distmat_temp = distmat
  }else{
    set.seed(i)
    ind<-sample(nrow(distmat))
    distmat_temp<-distmat[ind,ind]
  }
  
  res1=depth_CPD_cpp(distmat_temp,num_permut = 0)
  res2=depth_CPD_cpp_dF(distmat_temp,num_permut = 0)
  res3=depth_CPD_cpp_dF_AD(distmat_temp,num_permut = 0)
  res4=depth_CPD_cpp_dF_W(distmat_temp,num_permut = 0)
  l = list(res1,res2,res3,res4)
  names(l)<-c('dist_cpd_uniform','dist_cpd','dist_cpd_AD','dist_cpd_W')
  l  }
#### Graph CPD
E1 = mstree(dist( Data, method = 'euclidean' ),ngmax = 5)
result_graph = gseg1(nrow(distmat),E1, statistics="all",B = num_permut)
#### Energy CPD
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=num_permut,k=NULL,min.size=50,alpha=1)

#### MMD test
result_MMD<-c()
for (j in 0:num_permut){
  if (j!=0){
    set.seed(j)
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

result<-list(result_dist_profile,result_graph,result_ecp,result_MMD)
names(result)<-c("dist_profile",'graph','ecp','mmd')
path<-paste("exponentiate",'_delta_',delta,'_run_',monte_carlo,'.Rdata',sep="")
save(result, file=path)
