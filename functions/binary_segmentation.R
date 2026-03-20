
#decay ranges from 1 to 2

seeded.intervals <- function(n, decay = sqrt(2), unique.int = F){
  n	<- as.integer(n)
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  #n_layer=ceiling(log((min_len-1)/n+1,base=sqrt(0.5)))
  
  
  M	<- sum(2^(1:depth)-1)
  
  boundary_mtx           <- matrix(NA, ncol = 2)
  colnames(boundary_mtx) <- c("st", "end")
  boundary_mtx[1, ]      <- c(1, n)
  
  depth	<- log(n, base = decay)
  depth	<- ceiling(depth)
  
  
  for(i in 2:depth){
    int_length	<- n * (1/decay)^(i-1)
    
    n_int		<- ceiling(round(n/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
    
    boundary_mtx	<- rbind(boundary_mtx,
                          cbind(floor(seq(1, n-int_length, length.out = (n_int))), 
                                ceiling(seq(int_length, n, length.out = (n_int)))))
  }
  
  if(unique.int){return(unique(boundary_mtx))}
  boundary_mtx
}
