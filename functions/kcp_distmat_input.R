srcGetV = '
 #define SUM(A) std::accumulate(A.begin(), A.end(), 0.0)
 using namespace Rcpp;
 NumericMatrix K(K_);
 int N = K.nrow();
 NumericMatrix V(N,N);

 for(int i=0;i<N;++i)
 	for(int j=i;j<N;++j)
 		V(i,j) = V(j,i) = SUM(diag(K(Range(i,j),Range(i,j))))-SUM(K(Range(i,j),Range(i,j)))/(j-i+1);
 return wrap(V);
'
getVCpp = inline::cxxfunction(signature(K_='matrix'),srcGetV,plugin='Rcpp')


srcGetBandwidth = '
 	using namespace Rcpp;
 	NumericMatrix X(X_);
 	NumericVector rows(rows_);
 	int N = rows.size();
 	NumericVector u(N*N,(double)0);
 	for(int i=0;i<N;++i){
 		for(int j=0;j<N;++j){
 			u[i*N+j] = sum((X(rows[i],_)-X(rows[j],_))*(X(rows[i],_)-X(rows[j],_)));
 		}
 	}
 	return wrap(u);
 '
getBandwidthCpp = inline::cxxfunction(signature(X_='matrix',rows_='numeric'),srcGetBandwidth,plugin='Rcpp')

srcKcpa = '
 	using namespace Rcpp;
 	NumericMatrix II(II_);
 	NumericMatrix V(V_);
 	IntegerMatrix H(H_);
 	int N = V.nrow();
 	int L = H.nrow();

 	for(int k=1;k<L;++k){
 		for(int i=k;i<N;++i){
 			for(int j=k-1;j<i;++j){
 				double tmp = II(k-1,j) + V(j+1,i);
 				if(tmp < II(k,i)){
 					II(k,i) = tmp;
 					H(k,i) = j+1;//fix indexing differences between R and C++
 				}
 			}
 		}
 	}
 	return List::create(II,H);
 '
kcpaCpp = inline::cxxfunction(signature(II_='matrix',V_='matrix',H_='matrix'),srcKcpa,plugin='Rcpp')

#input euclidean distance matrix
getBandwidth_distmat = function(D){
  if(nrow(D) > 250)
    rows = sample(0:(nrow(D)-1),250,replace=FALSE)
  else
    rows = 0:(nrow(D)-1)
  return(median(as.vector(D^2)))
}

#input original data matrix
getV_distmat = function(D){
  DD = D**2
  h = getBandwidth_distmat(D)
  K = exp(-DD/h)
  V = getVCpp(K)
  return(V)
}


kcpa_distmat = function(D,L){
  #The maximum number of change points
  L = L+1 #user input L=max # of cps, algo uses L=max # of segments
  
  V = getV_distmat(D)
  N = nrow(D)
  II = matrix(Inf,L,N)
  H = matrix(0,L,N)
  II[1,] = V[1,]
  answer = kcpaCpp(II,V,H)
  
  II = answer[[1]]
  H = answer[[2]]
  S = II[,N]
  
  
  x1<-(function (i) 1/N * log(choose(N-1,i-1))) (1:L)[(0.6*L):L]
  x2<-(function (i) i/N) (1:L)[(0.6*L):L]
  lm_re<-lm(S[(0.6*L):L]~x1+x2)
  c1<-lm_re$coefficients[2]*(-2)
  c2<-lm_re$coefficients[3]*(-2)
  penalty<-(function (i) 1/N*(c1*log(choose(N-1,i-1))+c2*i)) (1:L)
  
  for(i in 1:L){
    #S[i] = S[i] + (C*i/N)*(1+log(N/i))
    #S[i] = S[i] + (C*v/N)*(i+log(choose(N-1,i-1)))
    S[i] = S[i] +penalty[i]
  }
  k = which.min(S)
  cps = N
  cp = cps
  while(k>0){
    cp = H[k,cp]
    cps = c(cp,cps)
    k = k-1
  }
  cps = cps+1
  return(cps)
}
kcpa_distmat_1cpd = function(D,L){
  #The maximum number of change points
  L = L+1 #user input L=max # of cps, algo uses L=max # of segments
  
  V = getV_distmat(D)
  N = nrow(D)
  II = matrix(Inf,L,N)
  H = matrix(0,L,N)
  II[1,] = V[1,]
  answer = kcpaCpp(II,V,H)
  
  II = answer[[1]]
  H = answer[[2]]
  S = II[,N]
  
  
  x1<-(function (i) 1/N * log(choose(N-1,i-1))) (1:L)[(0.6*L):L]
  x2<-(function (i) i/N) (1:L)[(0.6*L):L]
  lm_re<-lm(S[(0.6*L):L]~x1+x2)
  c1<-lm_re$coefficients[2]*(-2)
  c2<-lm_re$coefficients[3]*(-2)
  penalty<-(function (i) 1/N*(c1*log(choose(N-1,i-1))+c2*i)) (1:L)
  
  for(i in 1:L){
    #S[i] = S[i] + (C*i/N)*(1+log(N/i))
    #S[i] = S[i] + (C*v/N)*(i+log(choose(N-1,i-1)))
    S[i] = S[i] +penalty[i]
  }
  k = 2
  cps = N
  cp = cps
  while(k>0){
    cp = H[k,cp]
    cps = c(cp,cps)
    k = k-1
  }
  cps = cps+1
  return(cps)
}
