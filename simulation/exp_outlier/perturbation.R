require(gtools)
library(boot)
library(magrittr)
library(MASS)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)
library(mvtnorm)
library(optparse)
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##### new package
require("gSeg")
require(Rcpp)
require("kerSeg") 
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

# install.packages(c('Hotelling', 'gTests', 'matrixStats',"ade4","mvtnorm","gSeg","Rcpp","kerSeg"))
# install.packages(c("foreach","doParallel"))

library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(min(cores[1],11)-1)
registerDoParallel(cl=cl,cores=min(cores[1],11)-1) 
detectCores()

# Define command-line options
option_list <- list(
  make_option(c("--num_permut"), type = "integer", default = 5, 
              help = "Number of permutations for the permutation test [default: %default]"),
  make_option(c("--perturbate_type"), type = "character", default = 'outlier',
              help = "Type of perturbation: 'outlier' or 'exp'."),
  make_option(c("--perturbation_level"), type = "numeric", default = 0.1, 
              help = "Level of perturbation applied [default: %default]."),
  make_option(c("--delta"), type = "numeric", default = 0.5, 
              help = "Signal strength for change points [default: %default]."),
  make_option(c("--monte_carlo"), type = "numeric", default = 1, 
              help = "Index of Monte Carlo simulations [default: %default]."),
  make_option(c("--feature_outlier_perturbate"), type = "numeric", default = 0.1, 
              help = "Fraction of features to add outliers to [default: %default]."),
  make_option(c("--random_seed"), type = "numeric", default = 1, 
              help = "Random seed Monte Carlo simulations [default: %default]."),
  make_option(c("--data_type"), type = "numeric", default = 1, 
              help = "type of change points: 1,2,3,4 for mean, var, mixture, tail changes [default: %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))
n1=100;n2=200;n=n1+n2
data_type<-opt$data_type
num_permut<-opt$num_permut
perturbate_type<-opt$perturbate_type
perturbation_level<-opt$perturbation_level
delta<-opt$delta
monte_carlo<-opt$monte_carlo
feature_outlier_perturbate<-opt$feature_outlier_perturbate
rand_seed<-opt$random_seed
cat(perturbate_type,perturbation_level,delta,monte_carlo,feature_outlier_perturbate)
cat(class(perturbate_type),class(perturbation_level),class(delta),class(monte_carlo),class(feature_outlier_perturbate))
# perturbate_type<-'outlier';perturbation_level<-0.1;delta<-1;monte_carlo<-1

set.seed(rand_seed)
result=list()
for (dim_type in c(1,2,3)){
  dim=c(30,90,180)[dim_type] 
  Sigma=list(S1,S2,S3)[[dim_type]];I<-diag(x = 1, dim, dim)
  ##### generate base data with change point and exponentiate columns
  # Data<-rbind(mvrnorm(n1,mu=rep(0,dim),Sigma=Sigma),mvrnorm(n2,mu=c(rep(delta,dim)),Sigma=Sigma))
  if (data_type==1){
    # 1. mean diff
    Data<-rbind(mvrnorm(n1,mu=rep(0,dim),Sigma=Sigma),mvrnorm(n2,mu=c(rep(delta,dim)),Sigma=Sigma))
  }else if(data_type==2){
    # 2. scale diff
    Data<-rbind(mvrnorm(n1,mu=rep(0,dim),Sigma=0.8*I),mvrnorm(n2,mu=c(rep(0,dim)),Sigma=(0.8-delta)*I))
  }else if(data_type==3){
    # 3. mixture gaussian distributions
    A<-rbinom(n2,1,0.5)
    mu<-c(rep(delta,0.1*dim),rep(0,0.9*dim))
    Z1<-mvrnorm(n2,-mu,I);Z2<-mvrnorm(n2,mu,I)
    Data<-rbind(mvrnorm(n1,mu=rep(0,dim),Sigma=I),A*Z1+(1-A)*Z2)
  }else if(data_type==4){
    # 4. heavy tailed distribution: different degrees of freedom
    Data<-rbind(mvrnorm(n1,mu=rep(0,dim),Sigma=I),matrix(rt(n2*dim, df=delta),nrow=n2))
  }
  
  
  if (perturbate_type=='exp'){
    n_exponentiate = ceiling(perturbation_level * dim) # perturbation_level in seq(0,1,0.1)
    if (n_exponentiate>0){
      exp_coordinate<-seq(1,n_exponentiate)
      cat("number of exponentiate coordinate (features):",n_exponentiate,'\n')
      Data[,exp_coordinate]<-exp(Data[,exp_coordinate])
    }
  }else if(perturbate_type=='outlier'){
    cat('outlier perturbation')
    t_df = 10
    U = runif(n) # row level outlier (observation-level)
    noise<-10*matrix(rt(n*dim, t_df), nrow = n, ncol = dim) # 10 * t-distributed random noises
    mask_matrix<-matrix(rbinom(n * dim, size = 1, prob = feature_outlier_perturbate), nrow = n, ncol = dim)
    masked_noise<-noise * mask_matrix # mask entries with B_ij = 0, where B_ij \sim Bernoulli(feature_outlier_perturbate) 
    
    outliers <- matrix(0, nrow = n, ncol = dim)
    rows_to_perturb <- which(U <= perturbation_level)
    outliers[rows_to_perturb, ] <- masked_noise[rows_to_perturb, ]
    Data<-Data+outliers
  }
  distmat<-as.matrix(dist( Data, method = 'euclidean' ))
  #### distance profile variants CPD 
  result_dist_profile<-foreach (i = (0:num_permut),.noexport=c('depth_CPD_cpp','depth_CPD_cpp_dF','depth_CPD_cpp_dF_AD','depth_CPD_cpp_dF_W')) %dopar%{
    # library(foreach)
    # library(doParallel)
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
    l  }
  #### Graph CPD
  E1 = mstree(dist( Data, method = 'euclidean' ),ngmax = 5)
  result_graph = gseg1(nrow(distmat),E1, statistics="all",B = num_permut,pval.perm = TRUE)
  #### Energy CPD
  result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=num_permut,k=NULL,min.size=50,alpha=1)
  
  #### Kernel MMD test based CPD 
  result_MMD<-c()
  for (j in 0:num_permut){
    if (j!=0){
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
  result[[dim_type]]=r
}
change_type <- switch(as.character(data_type),
                      "1" = "mean",
                      "2" = "var",
                      "3" = "mixture",
                      "4" = "tail",
                      stop("Invalid data_type"))
path<-paste(change_type,'_',perturbate_type,'_perturbation_level_',perturbation_level,'_feature_perturbate_',feature_outlier_perturbate,'_delta_',delta,'_run_',monte_carlo,'.Rdata',sep="")
save(result, file=path)
