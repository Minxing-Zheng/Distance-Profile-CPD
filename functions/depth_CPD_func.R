#source("objTest_fctns.R")
#require(MASS)

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
########## This is the function to run the depth change point detection algorith
#' @param distmat The pairwise distance matrix 
#' @param c the cutoff for the interval, I_c
#' @param num_permut number of permutation
depth_CPD<-function(distmat,c=0.1,num_permut=500){
  n<-nrow(distmat)
  
  max_stat<-c()
  for (j in 0:num_permut){
    #data permutation
    cat(j,"th iteration","\n")
    if (j!=0){
      set.seed(j)
      ind<-sample(nrow(distmat))
      distmat_temp<-distmat[ind,ind]
    }else{
      distmat_temp<-distmat
    }
    ######## calculate the scan statistics
    test<-c()
    for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
      if (cp==1 || cp==n-1){next}
      testStat <- getT( distmat = distmat_temp, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
      test<-cbind(test,testStat)
    }
    ######## check the maximum scan statistics and its index
    max_stat<-cbind(max_stat,max(test))
    if (j==0){
      max_stat_index<-which.max(test)+ceiling(n*c)
    }
  }
  
  if (num_permut==0) {p_val<-NULL}else{
    p_val<-(1+sum(max_stat[,1]<max_stat[,2:ncol(max_stat)]))/(1+num_permut)
  }
  loc<-max_stat_index
  
  result=list()
  result[['p_val']]=p_val
  result[['loc']]=loc
  result[['observed_test_statistics']]=max_stat[1]
  #result[['permuted_test_statistics']]=max_stat[-1]
  result
}



depth_CPD_permutated<-function(distmat,c=0.1,num_permut=500){
  n<-nrow(distmat)
  
  max_stat<-c()
  for (j in 0:num_permut){
    #data permutation
    cat(j,"th iteration","\n")
    if (j!=0){
      set.seed(j)
      ind<-sample(nrow(distmat))
      distmat<-distmat[ind,ind]
    }
    ######## calculate the scan statistics
    test<-c()
    for (cp in seq(ceiling(n*c),n-ceiling(n*c),1)){
      if (cp==1 || cp==n-1){next}
      testStat <- getT( distmat = distmat, indices = 1:n, n = cp, m = n-cp, cut_off = 0 )
      test<-cbind(test,testStat)
    }
    ######## check the maximum scan statistics and its index
    max_stat<-cbind(max_stat,max(test))
    if (j==0){
      max_stat_index<-which.max(test)+ceiling(n*c)
    }
  }
  if (num_permut==0) {p_val<-NULL}else{
    p_val<-(1+sum(max_stat[,1]<max_stat[,2:ncol(max_stat)]))/(1+num_permut)
  }
  loc<-max_stat_index
  
  result=list()
  result[['p_val']]=p_val
  result[['loc']]=loc
  result[['observed_test_statistics']]=max_stat[1]
  result[['permuted_test_statistics']]=max_stat[-1]
  result
}



scan_depth_CPD<-function(distmat,scan_loc){
  n<-nrow(distmat)

  testStat <- getT( distmat = distmat, indices = 1:n, n = scan_loc, m = n-scan_loc, cut_off = 0 )
  testStat
}

