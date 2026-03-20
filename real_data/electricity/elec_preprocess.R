#https://www.eia.gov/electricity/data/browser/#/topic/0?agg=2,0,1&fuel=vtvv&geo=g&sec=g&linechart=ELEC.GEN.ALL-US-99.M~ELEC.GEN.COW-US-99.M~ELEC.GEN.NG-US-99.M~ELEC.GEN.NUC-US-99.M~ELEC.GEN.HYC-US-99.M~ELEC.GEN.WND-US-99.M~ELEC.GEN.TSN-US-99.M&columnchart=ELEC.GEN.ALL-US-99.M~ELEC.GEN.COW-US-99.M~ELEC.GEN.NG-US-99.M~ELEC.GEN.NUC-US-99.M~ELEC.GEN.HYC-US-99.M~ELEC.GEN.WND-US-99.M&map=ELEC.GEN.ALL-US-99.M&freq=M&start=200101&end=202212&ctype=linechart&ltype=pin&rtype=s&maptype=0&rse=0&pin=
#suppose the monthly net generation data is downloaded from above link.
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
citation("gSeg")

#load the raw data
headers = read.csv("elec_raw.csv", skip = 4, header = F, nrows = 1, as.is = T)
df = read.csv("elec_raw.csv", skip = 8, header = F)
colnames(df)= headers

#remove hydro-electric pumped storage, which is incorrect data
elec<-df[-c(15),]

#remove several unuseful rows and cols
elec_sub<-elec[-c(8,10,13,14,16),-c(1,2,3)] 

#impute the missing values
elec_sub[elec_sub=="--"]<-0

elec_matrix<-matrix(as.numeric(unlist(elec_sub)),nrow=nrow(elec_sub))
composition_matrix<-sweep(elec_matrix,2,colSums(elec_matrix),FUN="/")
#apply(composition_matrix[,c(1,2,3)],2,FUN=sum) 

#eacch row of the elec_data is a compositional vector
elec_data<-t(composition_matrix)

#combine several resources together
colnames(elec_data)<-c("coal","petroleum liquids","petroleum coke","natural gas","other gases","nuclear",'conventional hydro',"wind","geothermal",'biomass','other','sloar pho','scale solar')
elec_data[,2]<-elec_data[,2]+elec_data[,3]# combine petroleum
elec_data[,4]<-elec_data[,4]+elec_data[,5]# combine gases
elec_data[,8]<-elec_data[,8]+elec_data[,9]+elec_data[,10]+elec_data[,11]# renewable
elec_data[,12]<-elec_data[,12]+elec_data[,13] #solar
elec_data<-elec_data[,c(-3,-5,-9,-10,-11,-13)]

write.csv(elec_data,file='elec_processed.csv', row.names=FALSE)

geo_dist<-function(x,y){
  acos(sum(sqrt(x)*sqrt(y)))
}
X<-elec_data #264*7 dimension, 264 months and 7 resources
distmat<-as.matrix(dist_make(X,geo_dist)) 


