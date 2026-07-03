#include <Rcpp.h>
#include <algorithm> // for std::sort and std::upper_bound
#include <random>    // for std::shuffle and random number generation
using namespace Rcpp;

inline double empDf_scalar_fast(const NumericVector& sorted_x,
                                double t){
  int n = sorted_x.size();
  int ind = std::upper_bound(sorted_x.begin(),
                             sorted_x.end(),
                             t) - sorted_x.begin();
  return static_cast<double>(ind) / static_cast<double>(n);
}

// Function to compute test statistics for uniform integrating measure
// [[Rcpp::export]]
double getTcpp(const NumericMatrix& distmat,
                    const IntegerVector& indices,
                    int n,
                    int m,
                    double cut_off = 0,
                    int nqSup = 1000,
                    int ndSup = 1000){
  
  const int N = indices.size();
  
  NumericMatrix D(N, N);
  
  for(int i = 0; i < N; ++i){
    const int ii = indices[i] - 1;
    for(int j = 0; j < N; ++j)
      D(i,j) = distmat(ii, indices[j] - 1);
  }
  
  double minVal = D(0,0);
  double maxVal = D(0,0);
  
  for(int i = 0; i < N; ++i){
    for(int j = 0; j < N; ++j){
      if(D(i,j) < minVal) minVal = D(i,j);
      if(D(i,j) > maxVal) maxVal = D(i,j);
    }
  }
  
  const double step =
    (maxVal - minVal) / static_cast<double>(ndSup - 1);
  
  //--------------------------------------------------------
  // Pre-sort required distance vectors
  //--------------------------------------------------------
  
  std::vector<NumericVector> sortedXX(n);
  std::vector<NumericVector> sortedXY(n);
  std::vector<NumericVector> sortedYY(m);
  std::vector<NumericVector> sortedYX(m);
  
  for(int j = 0; j < n; ++j){
    
    sortedXX[j] = NumericVector(n - 1);
    int idx = 0;
    for(int k = 0; k < n; ++k){
      if(k != j)
        sortedXX[j][idx++] = D(k, j);
    }
    std::sort(sortedXX[j].begin(), sortedXX[j].end());
    
    sortedXY[j] = NumericVector(m);
    for(int k = 0; k < m; ++k)
      sortedXY[j][k] = D(n + k, j);
    std::sort(sortedXY[j].begin(), sortedXY[j].end());
  }
  
  for(int j = 0; j < m; ++j){
    
    sortedYY[j] = NumericVector(m - 1);
    int idx = 0;
    for(int k = 0; k < m; ++k){
      if(k != j)
        sortedYY[j][idx++] = D(n + k, n + j);
    }
    std::sort(sortedYY[j].begin(), sortedYY[j].end());
    
    sortedYX[j] = NumericVector(n);
    for(int k = 0; k < n; ++k)
      sortedYX[j][k] = D(k, n + j);
    std::sort(sortedYX[j].begin(), sortedYX[j].end());
  }
  
  //--------------------------------------------------------
  // Accumulate statistic directly
  //--------------------------------------------------------
  
  double sumX = 0.0;
  double sumY = 0.0;
  
  for(int j = 0; j < n; ++j){
    for(int r = 0; r < ndSup; ++r){
      
      double t = minVal + r * step;
      
      double fx  = empDf_scalar_fast(sortedXX[j], t);
      double fxy = empDf_scalar_fast(sortedXY[j], t);
      
      double diff = fx - fxy;
      sumX += diff * diff;
    }
  }
  
  for(int j = 0; j < m; ++j){
    for(int r = 0; r < ndSup; ++r){
      
      double t = minVal + r * step;
      
      double fy  = empDf_scalar_fast(sortedYY[j], t);
      double fyx = empDf_scalar_fast(sortedYX[j], t);
      
      double diff = fy - fyx;
      sumY += diff * diff;
    }
  }
  
  const double factor =
    (static_cast<double>(n) * static_cast<double>(m)) /
      static_cast<double>(n + m);
  
  const double integral_weight = step;
  
  double teststat =
    factor *
    ((sumX + sumY) * integral_weight /
      static_cast<double>(n + m));
  
  return teststat;
}

// [[Rcpp::export]]
List depth_CPD_cpp(NumericMatrix distmat,
                        double c = 0.1,
                        int num_permut = 500){
  
  const int n = distmat.nrow();
  
  const int first_cp = static_cast<int>(std::ceil(n * c));
  const int last_cp  = n - first_cp;
  
  NumericMatrix max_stat(num_permut + 1, 1);
  IntegerVector max_stat_index(num_permut + 1);
  
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
    
    //Rcout << b << "th iteration" << std::endl;
    
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
    
    double best = R_NegInf;
    int best_cp = first_cp;
    
    for(int cp = first_cp; cp <= last_cp; ++cp){
      
      double stat =
        getTcpp(current_distmat,
                     indices,
                     cp,
                     n - cp);
      
      if(stat > best){
        
        best = stat;
        best_cp = cp;
      }
    }
    
    max_stat(b,0) = best;
    
    if(b == 0)
      max_stat_index[b] = best_cp;
  }
  
  //--------------------------------------------------------
  // p-value
  //--------------------------------------------------------
  
  const double observed_stat = max_stat(0,0);
  
  int count = 0;
  
  for(int b = 1; b <= num_permut; ++b)
    count += (max_stat(b,0) > observed_stat);
  
  const double p_val =
    (1.0 + count) /
      (1.0 + num_permut);
  
  return List::create(
    Named("p_val") = p_val,
    Named("loc") = max_stat_index[0],
    Named("observed_test_statistics") = observed_stat,
    Named("max_stat") = max_stat,
    Named("max_stat_index") = max_stat_index
  );
}
