require(gtools)
library(boot)
library(magrittr)
library(MASS)
#library(energy)
library(Hotelling)
library(gTests)
library(matrixStats)
library(ade4)
library(mvtnorm)
require(usedist)
require(reshape2)
library(ggplot2)
library(igraph)
library(jsonlite)
### set parent dir as working dir
# setwd(dirname(normalizePath(sys.frame(1)$ofile)))
### set current dir as working dir
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("/home/mzheng/research/Prof_Dubey/code_final_new/real_data")
# getwd()
require("gSeg")
require(Rcpp)
require("kerSeg") 
library(optparse)
source("../functions/ecp_distmat_input.R")
sourceCpp('../functions/energyChangePoint.cpp')
sourceCpp("../functions/depth_CPDcpp.cpp")
sourceCpp("../functions/depth_CPDcpp_ALL.cpp")
# sourceCpp("../functions/depth_one_perm_cpp.cpp")
####### first load the dataset and compress it to daily data
####### since original dataset has time interval of 4-hour
####### then calculate the laplacian matrix

Data<-get(load("mix_human_machine_text/embeddings_word_boundary_300words/ex_0001.Rdata"))
json_list <- fromJSON("mix_human_machine_text/embeddings_word_boundary_300words/ex_0001.json")
cpd_loc<-json_list$cp_word
cpd_loc
dim(Data)



n_cores <- max(1, detectCores() - 1)

#### compute cosine similarity matrix of word embeddings
X <- as.matrix(Data)
# cosine similarity: S_ij = <x_i, x_j> / (||x_i|| ||x_j||)
norms <- sqrt(rowSums(X^2))
Xn <- X / pmax(norms, .Machine$double.eps)   # avoid divide-by-zero
cos_sim <- Xn %*% t(Xn)                      # n x n similarity matrix

# if you want a cosine *distance* matrix (like dist()):
cos_dist <- 1 - cos_sim
diag(cos_dist) <- 0

####### run depth_CPD function, and get the estimated change point location
num_permut<-500
res1=depth_CPD_cpp(cos_dist,num_permut = num_permut)
res_all = depth_CPDcpp_ALL(cos_dist,num_permut = 0) # dF, AD, W




#### Graph CPD
E1 = mstree(as.dist(cos_dist),ngmax = 5)
result_graph = gseg1(nrow(cos_dist),E1, statistics="g",B = num_permut,pval.perm=TRUE)
result_graph$pval.perm$generalized$pval
#### Energy CPD
result_ecp<-e.divisive_distmat(D=cos_dist,sig.lvl=.05,R=num_permut,k=1,min.size=50,alpha=1)
result_ecp$estimates[2]
result_ecp$p.values[1]
#### MMD test
library(parallel)


run_one <- function(j, cos_dist) {
  if (j != 0) {
    ind <- sample(nrow(cos_dist))
    D <- cos_dist[ind, ind]
  } else {
    D <- cos_dist
  }
  
  result_mmd <- MMD_test(D)   # vector of MMD statistics over candidate locations
  loc <- which.max(result_mmd)
  ob_stat <- result_mmd[loc]
  
  c(loc = loc, stat = ob_stat)
}
n_cores <- max(1, detectCores() - 1)

result_list <- mclapply(
  X = 0:n_perm,
  FUN = run_one,
  cos_dist = cos_dist,
  mc.cores = n_cores
)

result_MMD <- do.call(rbind, result_list)
result_MMD <- as.data.frame(result_MMD)

# observed result
obs_loc  <- result_MMD$loc[1]
obs_stat <- result_MMD$stat[1]
perm_stats <- result_MMD$stat[-1]
p_value <- (1 + sum(perm_stats >= obs_stat)) / (n_perm + 1)

list(
  obs_loc = obs_loc,
  obs_stat = obs_stat,
  p_value = p_value
)
