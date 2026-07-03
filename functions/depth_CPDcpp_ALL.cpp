#include <Rcpp.h>
#include <algorithm>
#include <random>

using namespace Rcpp;

struct TestStats {
  double stat_dF;
  double stat_AD;
  double stat_W;
};

////////////////////////////////////////////////////////////
//////////////////// FAST HELPERS ///////////////////////////
////////////////////////////////////////////////////////////

inline
NumericVector empDfCpp_All(const NumericVector& sorted_x,
                                const NumericMatrix& D,
                                int col){
  
  const int n = sorted_x.size();
  const int N = D.nrow();
  
  NumericVector result(N);
  
  for(int i = 0; i < N; ++i){
    
    int ind =
      std::upper_bound(sorted_x.begin(),
                       sorted_x.end(),
                       D(i, col)) - sorted_x.begin();
    
    result[i] =
      static_cast<double>(ind) /
        static_cast<double>(n);
  }
  
  return result;
}


//----------------------------------------------------------
// Column sums
//----------------------------------------------------------

inline
NumericVector colSumsCpp_All(const NumericMatrix& mat){
  
  const int nr = mat.nrow();
  const int nc = mat.ncol();
  
  NumericVector sums(nc);
  
  for(int j = 0; j < nc; ++j)
    for(int i = 0; i < nr; ++i)
      sums[j] += mat(i, j);
  
  return sums;
}


//----------------------------------------------------------
// Squared differences
//----------------------------------------------------------

inline
NumericMatrix squaredDiff_All(const NumericMatrix& A,
                                   const NumericMatrix& B){
  
  if(A.nrow() != B.nrow() || A.ncol() != B.ncol())
    stop("Matrices A and B must have the same dimensions.");
  
  const int nr = A.nrow();
  const int nc = A.ncol();
  
  NumericMatrix result(nr, nc);
  
  for(int i = 0; i < nr; ++i){
    for(int j = 0; j < nc; ++j){
      
      double d = A(i, j) - B(i, j);
      result(i, j) = d * d;
    }
  }
  
  return result;
}

// Unified test function for dF, AD, and W
inline
TestStats getTestStats_All(const NumericMatrix& distmat,
                                const IntegerVector& indices,
                                int n,
                                int m){
  
  const int N = indices.size();
  
  //--------------------------------------------------------
  // Construct sub-distance matrix
  //--------------------------------------------------------
  
  NumericMatrix DummyD(N, N);
  
  for(int i = 0; i < N; ++i){
    
    const int ii = indices[i] - 1;
    
    for(int j = 0; j < N; ++j)
      DummyD(i, j) = distmat(ii, indices[j] - 1);
  }
  
  //--------------------------------------------------------
  // Allocate CDF matrices
  //--------------------------------------------------------
  
  NumericMatrix FXX(N, n);
  NumericMatrix FYY(N, m);
  NumericMatrix FXY(N, n);
  NumericMatrix FYX(N, m);
  
  //--------------------------------------------------------
  // Pre-sort rows once
  //--------------------------------------------------------
  
  std::vector<NumericVector> sortedXX(n);
  std::vector<NumericVector> sortedXY(n);
  std::vector<NumericVector> sortedYY(m);
  std::vector<NumericVector> sortedYX(m);
  std::vector<NumericVector> sortedPool(N);
  
  for(int j = 0; j < n; ++j){
    
    sortedXX[j] = NumericVector(n);
    for(int k = 0; k < n; ++k)
      sortedXX[j][k] = DummyD(j, k);
    std::sort(sortedXX[j].begin(), sortedXX[j].end());
    
    sortedXY[j] = NumericVector(m);
    for(int k = 0; k < m; ++k)
      sortedXY[j][k] = DummyD(j, n + k);
    std::sort(sortedXY[j].begin(), sortedXY[j].end());
  }
  
  for(int j = 0; j < m; ++j){
    
    sortedYY[j] = NumericVector(m);
    for(int k = 0; k < m; ++k)
      sortedYY[j][k] = DummyD(n + j, n + k);
    std::sort(sortedYY[j].begin(), sortedYY[j].end());
    
    sortedYX[j] = NumericVector(n);
    for(int k = 0; k < n; ++k)
      sortedYX[j][k] = DummyD(n + j, k);
    std::sort(sortedYX[j].begin(), sortedYX[j].end());
  }
  
  for(int j = 0; j < N; ++j){
    
    sortedPool[j] = NumericVector(N);
    
    for(int i = 0; i < N; ++i)
      sortedPool[j][i] = DummyD(i, j);
    
    std::sort(sortedPool[j].begin(), sortedPool[j].end());
  }
  
  //--------------------------------------------------------
  // Populate FXX
  //--------------------------------------------------------
  
  for(int j = 0; j < n; ++j)
    FXX(_, j) =
      empDfCpp_All(sortedXX[j], DummyD, j);
  
  //--------------------------------------------------------
  // Populate FYY
  //--------------------------------------------------------
  
  for(int j = 0; j < m; ++j)
    FYY(_, j) =
      empDfCpp_All(sortedYY[j], DummyD, j + n);
  
  //--------------------------------------------------------
  // Populate FXY
  //--------------------------------------------------------
  
  for(int j = 0; j < n; ++j)
    FXY(_, j) =
      empDfCpp_All(sortedXY[j], DummyD, j);
  
  //--------------------------------------------------------
  // Populate FYX
  //--------------------------------------------------------
  
  for(int j = 0; j < m; ++j)
    FYX(_, j) =
      empDfCpp_All(sortedYX[j], DummyD, j + n);
  
  //--------------------------------------------------------
  // Pooled empirical CDF
  //--------------------------------------------------------
  
  NumericMatrix Fpooled(N, N);
  
  for(int j = 0; j < N; ++j)
    Fpooled(_, j) =
      empDfCpp_All(sortedPool[j], DummyD, j);
  
  //--------------------------------------------------------
  // Calculate all three statistics
  //--------------------------------------------------------
  
  double sum_dF = 0.0;
  double sum_AD = 0.0;
  double sum_W  = 0.0;
  
  for(int j = 0; j < n; ++j){
    
    for(int i = 0; i < N; ++i){
      
      double diff = FXX(i,j) - FXY(i,j);
      diff *= diff;
      
      sum_dF += diff;
      
      const double f = Fpooled(i,j);
      
      if(f > 0.0){
        
        sum_W += diff / f;
        
        if(f < 1.0)
          sum_AD += diff / (f * (1.0 - f));
      }
    }
  }
  
  for(int j = 0; j < m; ++j){
    
    for(int i = 0; i < N; ++i){
      
      double diff = FYY(i,j) - FYX(i,j);
      diff *= diff;
      
      sum_dF += diff;
      
      const double f = Fpooled(i,j+n);
      
      if(f > 0.0){
        
        sum_W += diff / f;
        
        if(f < 1.0)
          sum_AD += diff / (f * (1.0 - f));
      }
    }
  }
  
  const double factor =
    (static_cast<double>(n) * static_cast<double>(m)) /
      static_cast<double>(N);
  
  const double norm =
    static_cast<double>(N) * static_cast<double>(N);
  
  TestStats out;
  
  out.stat_dF = factor * (sum_dF / norm);
  out.stat_AD = factor * (sum_AD / norm);
  out.stat_W  = factor * (sum_W  / norm);
  
  return out;
}

// Output check
// [[Rcpp::export]]
List getTestStats_All_debug(const NumericMatrix& distmat,
                                 const IntegerVector& indices,
                                 int n,
                                 int m){
  
  TestStats s =
    getTestStats_All(distmat,
                          indices,
                          n,
                          m);
  
  return List::create(
    Named("stat_dF") = s.stat_dF,
    Named("stat_AD") = s.stat_AD,
    Named("stat_W")  = s.stat_W
  );
}

// Permutation-based p-values
// [[Rcpp::export]]
List depth_CPDcpp_ALL(const NumericMatrix& distmat,
                           double c = 0.1,
                           int num_permut = 500){
  
  const int n = distmat.nrow();
  
  const int first_cp = static_cast<int>(std::ceil(n * c));
  const int last_cp  = n - first_cp;
  const int num_cp   = last_cp - first_cp + 1;
  
  NumericVector observed_stat(3);
  NumericMatrix max_stat(3, num_permut + 1);
  IntegerMatrix max_stat_index(3, num_permut + 1);
  NumericVector p_val(3);
  IntegerVector loc(3);
  
  //--------------------------------------------------------
  // Reusable objects
  //--------------------------------------------------------
  
  NumericMatrix current_distmat(n, n);
  
  IntegerVector indices(n);
  for(int i = 0; i < n; ++i)
    indices[i] = i + 1;
  
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);
  
  std::random_device rd;
  std::default_random_engine eng(rd());
  
  //--------------------------------------------------------
  // Permutations
  //--------------------------------------------------------
  
  for(int b = 0; b <= num_permut; ++b){
    
    // Rcout << b << "th iteration" << std::endl;
    
    if(b == 0){
      
      std::copy(distmat.begin(),
                distmat.end(),
                current_distmat.begin());
      
    }else{
      
      std::iota(perm.begin(), perm.end(), 0);
      std::shuffle(perm.begin(), perm.end(), eng);
      
      for(int i = 0; i < n; ++i){
        
        const int ii = perm[i];
        
        for(int j = 0; j < n; ++j)
          current_distmat(i,j) = distmat(ii, perm[j]);
      }
    }
    
    double best_dF = R_NegInf;
    double best_AD = R_NegInf;
    double best_W  = R_NegInf;
    
    int best_dF_cp = first_cp;
    int best_AD_cp = first_cp;
    int best_W_cp  = first_cp;
    
    for(int cp = first_cp; cp <= last_cp; ++cp){
      
      TestStats stats =
        getTestStats_All(current_distmat,
                              indices,
                              cp,
                              n - cp);
      
      if(stats.stat_dF > best_dF){
        best_dF = stats.stat_dF;
        best_dF_cp = cp;
      }
      
      if(stats.stat_AD > best_AD){
        best_AD = stats.stat_AD;
        best_AD_cp = cp;
      }
      
      if(stats.stat_W > best_W){
        best_W = stats.stat_W;
        best_W_cp = cp;
      }
    }
    
    max_stat(0,b) = best_dF;
    max_stat(1,b) = best_AD;
    max_stat(2,b) = best_W;
    
    max_stat_index(0,b) = best_dF_cp;
    max_stat_index(1,b) = best_AD_cp;
    max_stat_index(2,b) = best_W_cp;
  }
  
  //--------------------------------------------------------
  // Observed statistics
  //--------------------------------------------------------
  
  for(int k = 0; k < 3; ++k){
    observed_stat[k] = max_stat(k,0);
    loc[k] = max_stat_index(k,0);
  }
  
  //--------------------------------------------------------
  // p-values
  //--------------------------------------------------------
  
  for(int k = 0; k < 3; ++k){
    
    int cnt = 0;
    
    for(int b = 1; b <= num_permut; ++b)
      cnt += (max_stat(k,b) > observed_stat[k]);
    
    p_val[k] =
      (1.0 + cnt) /
        (1.0 + num_permut);
  }
  
  return List::create(
    Named("p_val") = p_val,
    Named("loc") = loc,
    Named("observed_stat") = observed_stat,
    Named("max_stat") = max_stat,
    Named("max_stat_index") = max_stat_index
  );
}