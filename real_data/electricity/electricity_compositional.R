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
##### new package
require("gSeg")
require(Rcpp)
require("kerSeg") 

# setwd(dirname(rstudioapi::getSourceEditorContext()$path))
require("gSeg")
require(Rcpp)
require("kerSeg") 
library(optparse)
source("../functions/ecp_distmat_input.R")
sourceCpp('../functions/energyChangePoint.cpp')
sourceCpp("../functions/depth_CPDcpp.cpp")
sourceCpp("../functions/depth_CPDcpp_ALL.cpp")

##### define the geodesic distance for compositional data 
geo_dist<-function(x,y){
  acos(sum(sqrt(x)*sqrt(y)))
}

##### load the preprocessed dataset
elec_data<-read.csv("elec_processed.csv")
distmat<-as.matrix(dist_make(elec_data,geo_dist))
##### run depth_CPD function, to get the estimated change point location
num_permut<-500
res1=depth_CPD_cpp(distmat,num_permut = num_permut)
res_all = depth_CPDcpp_ALL(distmat,num_permut = num_permut) # dF, AD, W
res1$p_val
res1$loc
res_all$p_val
res_all$loc
result<-list()
result[['dist-CP-U']]<-res1
result[['dist-CP-ALL']]<-res_all
save(result,file='elec_single_cp_result.Rdata')
result<-get(load('elec_single_cp_result.Rdata'))
names(result)
result$`dist-CP-U` # 171
result$`dist-CP-ALL`$loc #132 132 130



### Graph CPD

E1 = mstree(as.dist(distmat),ngmax = 5)
result_graph = gseg1(nrow(distmat),E1, statistics="g",B = 1000,pval.perm=TRUE)
result_graph$pval.perm$generalized$pval
result_graph$scanZ$generalized$tauhat # 129
#### Energy CPD
result_ecp<-e.divisive_distmat(D=distmat,sig.lvl=.05,R=100,k=1,min.size=50,alpha=1)
result_ecp$estimates[2]# 171
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

n=nrow(distmat);c<-0.1
test<-c();test_all<-c() 
for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
  testStat <- getTcpp( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
  testStat_all <- getTestStats_All(distmat = distmat, indices = 1:n, n = cp, m = n-cp)
  test<-rbind(test,unlist(testStat))
  test_all<-rbind(test_all,unlist(testStat_all))
}
##### we combine the result and formulate into a long dataframe form
df_elec<-as.data.frame(elec_data)
names(df_elec)<-c("coal","petroleum","gas","nuclear",'conventional hydro',"renewable",'solar')
df_elec$time<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
df <- melt(df_elec ,  id.vars = 'time', variable.name = 'Resources')

##### plot the preprocess data with each composition of each resource,
##### along with a red vertical line indicating the estiamted change point location
ggplot(df, aes(time,value)) + 
  geom_line(aes(colour = Resources,lty=Resources))+
  xlab("Month") + ylab("Composition")+geom_vline(xintercept = as.numeric(df$time[170]), color = "black",lty="dashed")+
  theme_bw()+
  scale_x_date(date_breaks = "40 month", date_labels =  "%b %Y")+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))+
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave("elec_cpd.pdf",width = 7, height = 5,units="in")


##### calculate the scan statistics on interval I_c and plot it
n=nrow(distmat);c<-0.1
test<-c()
for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
  testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
  test<-cbind(test,testStat)
}

scan_stat<-c(rep(0,ceiling(n*c)-1),test,rep(0,ceiling(n*c)))
date_seq<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
df<-data.frame(scan_stat=scan_stat,day=date_seq)
ggplot(df,aes(x=day,y=scan_stat))+geom_line()+xlab("Month") + ylab("Scan statistic")+
  geom_vline(aes(xintercept = as.numeric(date_seq[depth_result$loc-1]),col='change point location \n Feb 2015'))+
  scale_color_manual(name = "Dist-CP", values = c("change point location \n Feb 2015" = "red"))+
  scale_x_date(date_breaks = "22 month", date_labels =  "%b %Y")+
  theme_bw()+theme(legend.position=c(0.85, 0.9))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  theme(panel.grid.minor = element_line(linetype="dashed"))+theme(panel.grid.major = element_line(linetype="dashed"))
ggsave("elec_scan_stat.pdf",width = 7, height = 5,units="in")


# 
# Date<-c()
# date_seq<-seq(as.Date("2001/1/1"), as.Date("2022/12/1"), "months")
# for (x in 1:length(date_seq)){
#   Date<-c(Date,format(date_seq[x], "%Y-%m"))
# }
# 
# library(lubridate)
# (d <- ymd("2001-01-01"))
# d %m+% months(170)
# 
