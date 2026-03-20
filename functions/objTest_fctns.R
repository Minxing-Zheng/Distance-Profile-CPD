# library(pracma)
# library(MASS)
# library(energy)
# library(gTests)
# library(igraph)
# library(ade4)
# library(matrixStats)
# library(magrittr)

# function computing Hodge rank ----
#'@param d A matrix holding the pairwise distances between individuals.
getHodgeDepth<-function(d){
  qtgrid <- seq(0,1,1/50)
  sortd <- sort(as.vector(d))
  outgrid <- seq(sortd[1]-0.005,sortd[length(d)]+0.005,(sortd[length(d)]+0.005-sortd[1]+0.005)/100)
  cdf <- matrix(unlist(lapply(1:nrow(d),function(i){sort(unlist(lapply(1:length(outgrid),function(j){length(which(d[i,-i]<=outgrid[j]))/length(d[i,-i])})))})),ncol = length(outgrid), byrow = TRUE)
  qtOnMidPts <- matrix(unlist(lapply(1:nrow(cdf),function(i){approx(cdf[i,],outgrid,qtgrid,ties="min")$y})),ncol = length(qtgrid), byrow = TRUE)
  n <- nrow(qtOnMidPts)
  HodgeRank <- sapply( seq_len(n), function(i) {
    qtDiff <- qtOnMidPts[-i,] - matrix(qtOnMidPts[i,], nrow = n-1, ncol = ncol(qtOnMidPts), byrow = TRUE)
    sum( sign( rowMeans( qtDiff ) ) * ( rowMeans( abs( qtDiff ) ) ) )
  }) / n
  return(HodgeRank)
}

# function computing empirical quantile functions ----
#' Empirical quantile function, i.e., left continuous inverse of empirical distribution function.
#' @param x A vector holding the data.
#' @param qSup Support grid of the output empirical quantile function, all the elements of which should be between 0 and 1.
#' @return A vector of the same length as \code{qSup} holding the empirical quantile function evaluated on \code{qSup}.
empQf <- function ( x, qSup ) {
  ind <- empQfInd( n = length(x), qSup = qSup )
  if ( is.unsorted(x) ) {
    xSorted <- sort(x)
  } else {
    xSorted <- x
  }
  xSorted[ind]
}
#' @param n Sample size.
empQfInd <- function ( n, qSup ) {
  ind <- findInterval( x = qSup, vec = seq_len(n) / n, left.open = TRUE ) + 1
  ind[ind > n] <- n
  ind
}

# function computing empirical distribution functions ----
#' Right continuous empirical distribution function. 
#' @param x A vector holding the data.
#' @param dsup Support grid of the output empirical distribution function.
#' @return A vector of the same length as \code{dSup} holding the empirical distribution function evaluated on \code{dSup}.
empDf <- function ( x, dSup ) {
  ind <- findInterval( x = dSup, vec = sort(x) ) + 1
  n <- length(x)
  val <- c( 0, seq_len(n) / n )
  val[ind]
}

# function computing the test statistic ----
#' @param distmat An (\code{n+m})-by-(\code{n+m}) distance matrix.
#' @param indices A vector of length \code{n+m} holding the indices of rows and columns of \code{distmat} to be used.
#' @param n,m Sample sizes of the two samples.
#' @param cut_off A vector holding the values of a thresholding parameter. 
#' The test statistic only include observations of which the rank depths are not less than 
#' the \code{cut_off}-quantile of the rank depths within the corresponding sample. Default: 0.
#' @param ndSup,nqSup The number of equidistant points in the support grids of 
#' distribution functions and quantile functions, respectively.
getT <- function( distmat, indices, n, m, cut_off = 0, nqSup = 1000, ndSup = 1000 ){
  
  DummyD <- distmat[indices,indices]
  dXX <- DummyD[1:n,1:n]
  dYY <- DummyD[(n+1):(n+m),(n+1):(n+m)]
  dXY <- DummyD[1:n,(n+1):(n+m)]
  dYX <- t(dXY)
  
  
  qSup <- ( seq_len(nqSup) - 0.5 ) / nqSup
  dSup <- range(DummyD)
  dSup <- seq( dSup[1], dSup[2], length.out = ndSup )
  
  offdiag <- function (A) {
    dimA <- dim(A)
    if ( dimA[1] != dimA[2] ) {
      stop ("A is not a square matrix.")
    }
    matrix( A[ upper.tri(A) | lower.tri(A) ], nrow = dimA[1]-1, ncol = dimA[2] )
  }
  
  FXX <- apply( offdiag( dXX ), 2, empDf, dSup = dSup )
 
  FYY <- apply( offdiag( dYY ), 2, empDf, dSup = dSup )

  FXY <- apply( dYX, 2, empDf, dSup = dSup )

  FYX <- apply( dXY, 2, empDf, dSup = dSup )
  
  TFX <- colSums( ( FXX[,1:ncol(FXX)] - FXY[,1:ncol(FXX)] ) ^ 2 ) * ( dSup[ndSup] - dSup[1] ) / ( ndSup - 1 )
 
  TFY <- colSums( ( FYY[,1:ncol(FYY)] - FYX[,1:ncol(FYY)] ) ^ 2 ) * ( dSup[ndSup] - dSup[1] ) / ( ndSup - 1 )

  teststat <- (ncol(FXX)*ncol(FYY)/(ncol(FXX)+ncol(FYY)))*(mean( TFX ) + mean( TFY ))
  teststat2 <-(ncol(FXX)*ncol(FYY)/(ncol(FXX)+ncol(FYY)))*((sum( TFX ) + sum( TFY ))/(length(TFX)+length(TFY)))
  teststat3 <-sum( TFX ) + sum( TFY )
  #return( c(teststat,teststat2,teststat3))
  return(c(teststat2))
}

# function computing the test statistic and p-value for our test ----
#'@param n,m Sample sizes for the two samples.
#'@param distmat An \eqn{(n+m)\times(n+m)} dimensional distance matrix.
#'@param cut_off A vector of cut-offs.
#'@param B Number of bootstrap runs (default: 1000).
#'@param ncores Number of cores for parallel computing (if greater than 1). Default: 1.
twoSampleTest <- function( n, m, distmat, cut_off, B = 1000, ncores = 1, replace = FALSE ){
  
  testStat <- getT( distmat = distmat, indices = 1:(n+m), n = n, m = m, cut_off = cut_off )
  
  if ( ncores == 1 ){
    testStat.boot <- lapply( 1:B, function(i) {
      getT(distmat, sample(1:(n+m), replace = replace), n, m, cut_off)
    })
  } else {
    require(doParallel); cl = makeCluster(ncores); registerDoParallel(cl)
    testStat.boot <- foreach( i = 1:B, .export = c("getHodgeDepth","getT","empDf","empQfInd","empQf") ) %dopar% {
      getT( distmat, sample(1:(n+m), replace = replace), n, m, cut_off )
    }
    stopCluster(cl)
  }
  
  methods <- names(testStat)
  methods.boot <- names( testStat.boot[[1]] )
  
  pval <- sapply( seq_along(testStat), function (j) {
    res <- testStat[[j]]
    if ( methods[j] %in% methods.boot ) {
      # permutation tests
      resboot <- sapply( testStat.boot, function (res) res[[j]] )
      if ( length( res ) > 1 ) {
        # multiple cut-offs
        return ( min( rowSums( ( resboot - res ) > 0 ) ) / B )
      } else {
        return ( sum( resboot > res ) / B )
      }
    } else {
      # asymptotically N(0,1)
      return ( pnorm( res, lower.tail = ( res < 0 ) ) * 2 )
    }
  })
  names(pval) <- methods
  
  return ( pval )
}

# function computing the power of energy test, graph-based test and our test ----
#'@param n,m Sample sizes for the two samples.
#'@param D An object of class \code{dist} holding the pairwise distance between observations.
#'@param cut_off A vector of cut-offs.
#'@param R_energy Number of bootstrap runs for energy tests (default: 100).
#'@param R_graph Number of permutation runs for graph-based tests. The default value is 0, i.e., 
#'the permutation is not performed and only approximate p-value based on asymptotic theory is provided.
#'@param R_new Number of bootstrap runs for new tests (default: 100).
#'@param ncores Number of cores for parallel computing (if greater than 1). Default: 1.
tests <- function (n, m, D, cut_off = 0, R_energy = 1000, R_graph = 0, R_new = 1000, ncores = 1){
  # sam=rbind(sam1,sam2)
  # n=nrow(sam1)
  # m=nrow(sam2)
  # D<-dist(sam)
  
  # graph-based test
  # g <- ade4::neig2mat(ade4::mstree(D,ngmax=5))
  # E <- as_edgelist(as.undirected(graph.adjacency(g)))
  # E <- cbind(as.numeric(E[,1]),as.numeric(E[,2]))
  # res <- g.tests(E, 1:n, (n+1):(n+m), test.type="all", perm=R_graph)
  # g.pval <- res$generalized$pval.approx
  
  # energy test
  # e.pval = eqdist.etest(D, c(n,m), distance = TRUE, method="original", R = R_energy)$p.val 
  
  # our test
  new.pval <- twoSampleTest( n = n, m = m, distmat = as.matrix(D), cut_off = cut_off, B = R_new, ncores = ncores, replace = FALSE )
  
  return( list( e.pval = e.pval, g.pval=g.pval, new.pval = new.pval ) )
}

