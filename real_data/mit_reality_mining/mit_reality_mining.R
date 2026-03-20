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
##### new package
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
require("gSeg")
require(Rcpp)
require("kerSeg") 
library(optparse)
source("../functions/objTest_fctns.R")
source("../functions/depth_CPD_func.R") 
source("../functions/ecp_distmat_input.R")
source("../functions/kcp_distmat_input.R")
sourceCpp('../functions/energyChangePoint.cpp')
sourceCpp("../functions/getTcpp.cpp")
sourceCpp("../functions/depth_CPDcpp.cpp")
sourceCpp("../functions/depth_CPDcpp_ALL.cpp")

####### first load the dataset and compress it to daily data
####### since original dataset has time interval of 4-hour
####### then calculate the laplacian matrix
MIT<-get(load("reality_mining_1392.RData"))
G_list<-list()
for (i in 1:(dim(MIT)[3]/6)){
  G_list[[i]]<-rowSums(MIT[,,(i*6-5):(i*6)],dims=2)
}
Data<-c()
for (i in 1:(dim(MIT)[3]/6)){
  g<-as.matrix(laplacian_matrix(graph_from_adjacency_matrix(G_list[[i]])))
  Data<-rbind(Data,c(2*g[upper.tri(g)],diag(g)))
}

####### calculate the distance matrix 
distmat<-as.matrix(dist(Data , method = 'manhattan' ))

####### run depth_CPD function, and get the estimated change point location
num_permut<-0
res1=depth_CPD_cpp(distmat,num_permut = num_permut)
res_all = depth_CPDcpp_ALL(distmat,num_permut = num_permut) # dF, AD, W
res1$loc # 94
res_all$loc # 94 93 94


### Graph CPD

E1 = mstree(as.dist(distmat),ngmax = 5)
result_graph = gseg1(nrow(distmat),E1, statistics="g",B = 1000,pval.perm=TRUE)
result_graph$pval.perm$generalized$pval
result_graph$scanZ$generalized$tauhat # 214
#### Energy CPD
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=100,k=1,min.size=50,alpha=1)
result_ecp$estimates[2]# 94
### MMD
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

n_perm <- 300
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
obs_loc  <- result_MMD$loc[1] # 91
obs_stat <- result_MMD$stat[1]
perm_stats <- result_MMD$stat[-1]
p_value <- (1 + sum(perm_stats >= obs_stat)) / (n_perm + 1)



#calculate the scan statistics sequence and plot it 
n=nrow(distmat);c<-0.1
test<-c();test_all<-c()
for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
  testStat <- getTcpp( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
  testStat_all <- getTestStats_All(distmat = distmat, indices = 1:n, n = cp, m = n-cp)
  test<-rbind(test,unlist(testStat))
  test_all<-rbind(test_all,unlist(testStat_all))
}
est_value<-test_all[,3]
est_loc<-res_all$loc[3]
scan_stat<-c(rep(0,ceiling(n*c)-1),est_value,rep(0,ceiling(n*c)))
date_seq<-seq(as.Date("2004/9/14"), as.Date("2004-09-14")+dim(distmat)[1]-1, "days")

df<-data.frame(scan_stat=scan_stat,day=date_seq)
ggplot(df,aes(x=day,y=scan_stat))+geom_line()+xlab("Day") + ylab("Scan statistic")+
  geom_vline(aes(xintercept = as.numeric(date_seq[est_loc])-1,col='change point location \n 2004-12-15'))+
  scale_color_manual(name = "Dist-CP", values = c("change point location \n 2004-12-15" = "red"))+
  scale_x_date(date_breaks = "22 days", date_labels =  "%d %b %Y")+
  theme_bw()+theme(legend.position=c(0.85, 0.9))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))
ggsave("mit_scan_stat.pdf",width = 7, height = 5,units="in")
 
