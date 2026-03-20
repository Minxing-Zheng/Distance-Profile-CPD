#include <Rcpp.h>
#include <algorithm> // for std::sort and std::upper_bound
#include <random>    // for std::shuffle and random number generation
using namespace Rcpp;

// Function to calculate empirical distribution function
// [[Rcpp::export]]
NumericVector empDfCpp_dF_AD(const NumericVector& x, const NumericVector& dSup) {
  // Sort the vector 'x'
  NumericVector sorted_x = Rcpp::clone(x);
  std::sort(sorted_x.begin(), sorted_x.end());
  
  // Create a vector to store the indices
  IntegerVector ind(dSup.size());
  
  // Find intervals (similar to R's findInterval)
  for (int i = 0; i < dSup.size(); ++i) {
    ind[i] = std::upper_bound(sorted_x.begin(), sorted_x.end(), dSup[i]) - sorted_x.begin();
  }
  
  // Calculate the values as seq_len(n) / n
  int n = x.size();
  NumericVector val(n + 1);
  for (int i = 0; i < n; ++i) {
    val[i + 1] = (i + 1) / static_cast<double>(n);
  }
  
  // Return the corresponding values for ind
  NumericVector result(ind.size());
  for (int i = 0; i < ind.size(); ++i) {
    result[i] = val[ind[i]];
  }
  
  return result;
}

// Function to compute the column sums
NumericVector colSumsCustom_dF_AD(const NumericMatrix& mat) {
  int ncol = mat.ncol();           // Get number of columns
  NumericVector sums(ncol);        // Initialize vector to store sums
  
  for (int j = 0; j < ncol; ++j) { // Loop over each column
    for (int i = 0; i < mat.nrow(); ++i) { // Loop over each row
      sums[j] += mat(i, j);        // Sum the elements of each column
    }
  }
  
  return sums;                     // Return the resulting sums
}

// Function to compute element-wise squared differences between two matrices
NumericMatrix squaredDiff_dF_AD(const NumericMatrix& A, const NumericMatrix& B) {
  // Check if the matrices have the same dimensions
  if (A.nrow() != B.nrow() || A.ncol() != B.ncol()) {
    stop("Matrices A and B must have the same dimensions.");
  }
  
  int nrow = A.nrow();  // Number of rows
  int ncol = A.ncol();  // Number of columns
  
  NumericMatrix result(nrow, ncol);  // Initialize result matrix
  
  // Calculate squared differences
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      result(i, j) = pow(A(i, j) - B(i, j), 2);  // Compute squared difference
    }
  }
  
  return result;  // Return the result matrix
}

// [[Rcpp::export]]
double getTcpp_dF_AD(const NumericMatrix& distmat, const IntegerVector& indices, int n, int m, 
                     double cut_off = 0) {
  
  // Extract submatrices based on indices
  NumericMatrix DummyD(indices.size(), indices.size());
  for (int i = 0; i < indices.size(); ++i) {
    for (int j = 0; j < indices.size(); ++j) {
      DummyD(i, j) = distmat(indices[i] - 1, indices[j] - 1); // Adjust for 1-based indexing
    }
  }
  
  NumericMatrix dXX = DummyD(Range(0, n - 1), Range(0, n - 1));
  NumericMatrix dYY = DummyD(Range(n, n + m - 1), Range(n, n + m - 1));
  NumericMatrix dXY = DummyD(Range(0, n - 1), Range(n, n + m - 1));
  NumericMatrix dYX = transpose(dXY);
  
  // Create Fpooled
  NumericMatrix Fpooled(DummyD.nrow(), DummyD.ncol());
  for (int j = 0; j < DummyD.ncol(); ++j) {
    Fpooled(_, j) = empDfCpp_dF_AD(DummyD(_, j), DummyD(_, j)); // Column-wise computation
  }
  
  // Create W = 1 / (Fpooled * (1 - Fpooled)) if > 0, else 0
  NumericMatrix W(DummyD.nrow(), DummyD.ncol());
  for (int i = 0; i < Fpooled.nrow(); ++i) {
    for (int j = 0; j < Fpooled.ncol(); ++j) {
      double value = Fpooled(i, j) * (1.0 - Fpooled(i, j)); // Compute Fpooled * (1 - Fpooled)
      if (value > 0.0) {                                   // Check for positive value
        W(i, j) = 1.0 / value;                             // Compute W
      } else {
        W(i, j) = 0.0;                                     // Assign 0 otherwise
      }
    }
  }
  
  // Initialize FXX, FYY, FXY, FYX with appropriate sizes
  NumericMatrix FXX(DummyD.nrow(), dXX.nrow());
  NumericMatrix FYY(DummyD.nrow(), dYY.nrow());
  NumericMatrix FXY(DummyD.nrow(), dXY.nrow());
  NumericMatrix FYX(DummyD.nrow(), dYX.nrow());
  
  // Populate FXX, FYY, FXY, FYX
  for (int j = 0; j < dXX.nrow(); ++j) {
    FXX(_, j) = empDfCpp_dF_AD(dXX(j, _), DummyD(_, j));
  }
  for (int j = 0; j < dYY.nrow(); ++j) {
    FYY(_, j) = empDfCpp_dF_AD(dYY(j, _), DummyD(_, j + n));
  }
  for (int j = 0; j < dXY.nrow(); ++j) {
    FXY(_, j) = empDfCpp_dF_AD(dXY(j, _), DummyD(_, j)); // Use DummyD(_, j) as dSup
  }
  for (int j = 0; j < dYX.nrow(); ++j) {
    FYX(_, j) = empDfCpp_dF_AD(dYX(j, _), DummyD(_, j + n)); // Use DummyD(_, j + n) as dSup
  }
  
  // Calculate squared differences
  NumericMatrix diffFXX_FXY = squaredDiff_dF_AD(FXX, FXY);
  NumericMatrix diffFYY_FYX = squaredDiff_dF_AD(FYY, FYX);
  
  // Modify TFX: Element-wise multiplication with the first `n` columns of W, then take column sums
  NumericMatrix weighted_diffFXX_FXY(diffFXX_FXY.nrow(), diffFXX_FXY.ncol());
  for (int i = 0; i < diffFXX_FXY.nrow(); ++i) {
    for (int j = 0; j < diffFXX_FXY.ncol(); ++j) {
      weighted_diffFXX_FXY(i, j) = diffFXX_FXY(i, j) * W(i, j); // Element-wise multiplication
    }
  }
  NumericVector TFX = colSumsCustom_dF_AD(weighted_diffFXX_FXY) * (1.0 / static_cast<double>(n + m)); // Adjusted TFX
  
  // Modify TFY: Element-wise multiplication with the last `m` columns of W, then take column sums
  NumericMatrix weighted_diffFYY_FYX(diffFYY_FYX.nrow(), diffFYY_FYX.ncol());
  for (int i = 0; i < diffFYY_FYX.nrow(); ++i) {
    for (int j = 0; j < diffFYY_FYX.ncol(); ++j) {
      weighted_diffFYY_FYX(i, j) = diffFYY_FYX(i, j) * W(i, j + n); // Element-wise multiplication
    }
  }
  NumericVector TFY = colSumsCustom_dF_AD(weighted_diffFYY_FYX) * (1.0 / static_cast<double>(n + m)); // Adjusted TFY
  
  // Calculate the test statistic
  double teststat = (static_cast<double>(FXX.ncol()) * static_cast<double>(FYY.ncol()) / 
                     (static_cast<double>(FXX.ncol()) + static_cast<double>(FYY.ncol()))) * 
                     ((sum(TFX) + sum(TFY)) / static_cast<double>(TFX.size() + TFY.size()));
  
  return teststat;
}





// [[Rcpp::export]]
List depth_CPD_cpp_dF_AD(NumericMatrix distmat, double c = 0.1, int num_permut = 500) {
  int n = distmat.nrow(); // Number of rows in the distance matrix
  
  // Initialize matrix to store maximum statistics and indices
  NumericMatrix max_stat(num_permut + 1, 1);
  IntegerVector max_stat_index(num_permut + 1);
  
  for (int j = 0; j <= num_permut; j++) {
    Rcout << j << "th iteration" << std::endl;
    
    // If not the first iteration, permute the distance matrix
    NumericMatrix current_distmat = clone(distmat);
    if (j != 0) {
      IntegerVector ind = seq_len(n) - 1; // Generate sequence of indices
      std::random_device rd;             // Random number generator
      std::default_random_engine eng(rd());
      std::shuffle(ind.begin(), ind.end(), eng); // Shuffle indices
      
      // Apply the permutation to the rows and columns of the matrix
      for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
          current_distmat(i, k) = distmat(ind[i], ind[k]);
        }
      }
    }
    
    // Initialize vector to store test statistics for all possible change points
    NumericVector test;
    for (int cp = ceil(n * c); cp <= n - ceil(n * c); cp++) {
      // Call the updated getTcpp_dF function
      double testStat = getTcpp_dF_AD(current_distmat, seq_len(n), cp, n - cp);
      test.push_back(testStat);
    }
    
    // Store the maximum test statistic and its index
    double max_value = max(test);
    max_stat(j, 0) = max_value;
    if (j == 0) {
      max_stat_index(j) = which_max(test) + static_cast<int>(ceil(n * c))+1; // Ensure numeric division
    }
  }
  
  // Calculate the p-value based on the permutation results
  double observed_stat = max_stat(0, 0); // Observed maximum statistic
  int count = 0;
  for (int j = 1; j <= num_permut; j++) {
    if (max_stat(j, 0) > observed_stat) {
      count++;
    }
  }
  double p_val = (1.0 + static_cast<double>(count)) / (1.0 + static_cast<double>(num_permut)); // Ensure numeric division
  int loc = max_stat_index[0]; // Location of the detected change point
  
  // Return the results
  return List::create(Named("p_val") = p_val,
                      Named("loc") = loc,
                      Named("observed_test_statistics") = observed_stat);
}
