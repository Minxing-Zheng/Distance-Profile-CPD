require(pracma)
make_Sigma_list <- function(p_list) {
  Sigma_s <- vector("list", length(p_list))
  names(Sigma_s) <- paste0("p_", p_list)
  
  for (i in seq_along(p_list)) {
    p <- p_list[i]
    
    M <- matrix(1, 1, p) / sqrt(p)
    N <- nullspace(M)
    U <- cbind(t(M), N)
    
    v <- c()
    for (k in 1:p) {
      v <- c(v, cos(k * pi / p) + 1.5)
    }
    
    Delta <- diag(v)
    Sigma <- U %*% Delta %*% t(U)
    
    path <- paste0("Sigma_", p, ".RData")
    save(Sigma, file = path)
    
    Sigma_s[[i]] <- Sigma
  }
  
  return(Sigma_s)
}