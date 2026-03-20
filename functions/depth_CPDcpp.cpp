#include <Rcpp.h>
#include <algorithm> // for std::sort and std::upper_bound
#include <random>    // for std::shuffle and random number generation
using namespace Rcpp;

// Function to calculate empirical distribution function
// [[Rcpp::export]]
NumericVector empDfCpp(const NumericVector& x, const NumericVector& dSup) {
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

// Function to extract off-diagonal elements
// [[Rcpp::export]]
NumericMatrix offdiagCpp(const NumericMatrix& A) {
  // Check if the matrix is square
  if (A.nrow() != A.ncol()) {
    stop("A is not a square matrix.");
  }
  
  int n = A.nrow();
  
  // Create a NumericMatrix to hold the off-diagonal elements
  NumericMatrix offdiagMat(n - 1, n); // Size (n-1) x n
  
  // Fill the offdiagMat with the off-diagonal elements
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      if (i != j) { // Skip diagonal elements
        int rowIndex = (i < j) ? i : (i - 1); // Adjust row index
        offdiagMat(rowIndex, j) = A(i, j);
      }
    }
  }
  
  return offdiagMat;
}

// Function to compute the column sums
NumericVector colSumsCustom(const NumericMatrix& mat) {
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
NumericMatrix squaredDiff(const NumericMatrix& A, const NumericMatrix& B) {
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
double getTcpp(const NumericMatrix& distmat, const IntegerVector& indices, int n, int m, 
               double cut_off = 0, int nqSup = 1000, int ndSup = 1000) {
  
  // Extract submatrices based on indices
  NumericMatrix DummyD(indices.size(), indices.size());
  for (int i = 0; i < indices.size(); ++i) {
    for (int j = 0; j < indices.size(); ++j) {
      DummyD(i, j) = distmat(indices[i] - 1, indices[j] - 1); // Adjust for 1-based indexing
    }
  }
  
  NumericMatrix dXX = DummyD(Range(0, n-1), Range(0, n-1));
  NumericMatrix dYY = DummyD(Range(n, n+m-1), Range(n, n+m-1));
  NumericMatrix dXY = DummyD(Range(0, n-1), Range(n, n+m-1));
  NumericMatrix dYX = transpose(dXY);
  
  // Create qSup manually
  NumericVector qSup(nqSup);
  for (int i = 0; i < nqSup; ++i) {
    qSup[i] = (i + 1 - 0.5) / nqSup; // Fill qSup with desired values
  }
  
  // Create dSup manually
  NumericVector dSup(ndSup);
  double minVal = min(DummyD);
  double maxVal = max(DummyD);
  double step = (maxVal - minVal) / static_cast<double>(ndSup - 1);
  
  for (int i = 0; i < ndSup; ++i) {
    dSup[i] = minVal + i * step; // Fill dSup with values from minVal to maxVal
  }
  
  // Calculate FXX, FYY, FXY, FYX using offdiagCpp and empDfCpp
  NumericMatrix offdiagFXX = offdiagCpp(dXX);
  NumericMatrix FXX(dSup.size(), offdiagFXX.ncol());
  for (int j = 0; j < offdiagFXX.ncol(); ++j) {
    // Check the length of dSup and empDfCpp result
    FXX(_, j) = empDfCpp(offdiagFXX(_, j), dSup);
  }
  
  NumericMatrix offdiagFYY = offdiagCpp(dYY);
  NumericMatrix FYY(dSup.size(), offdiagFYY.ncol());
  for (int j = 0; j < offdiagFYY.ncol(); ++j) {
    FYY(_, j) = empDfCpp(offdiagFYY(_, j), dSup);
  }
  
  NumericMatrix FXY(dSup.size(), dYX.ncol());
  for (int j = 0; j < dYX.ncol(); ++j) {
    FXY(_, j) = empDfCpp(dYX(_, j), dSup);
  }
  
  NumericMatrix FYX(dSup.size(), dXY.ncol());
  for (int j = 0; j < dXY.ncol(); ++j) {
    FYX(_, j) = empDfCpp(dXY(_, j), dSup);
  }
  
  // Calculate squared differences using the custom function
  NumericMatrix diffFXX_FXY = squaredDiff(FXX, FXY);
  
  // Calculate TFX using the custom column sums function
  NumericVector TFX = colSumsCustom(diffFXX_FXY) * (dSup[ndSup - 1] - dSup[0]) / static_cast<double>(ndSup - 1);
  
  // Calculate squared differences for TFY
  NumericMatrix diffFYY_FYX = squaredDiff(FYY, FYX);
  
  // Calculate TFY
  NumericVector TFY = colSumsCustom(diffFYY_FYX) * (dSup[ndSup - 1] - dSup[0]) / static_cast<double>(ndSup - 1);
  
  // Calculate the test statistic
  double teststat = (static_cast<double>(FXX.ncol()) * FYY.ncol() / 
                     (static_cast<double>(FXX.ncol()) + FYY.ncol())) * 
                     ((sum(TFX) + sum(TFY)) / (TFX.size() + TFY.size()));
  
  return teststat;
}


// [[Rcpp::export]]
List depth_CPD_cpp(NumericMatrix distmat, double c = 0.1, int num_permut = 500) {
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
      double testStat = getTcpp(current_distmat, seq_len(n), cp, n - cp);
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


// [[Rcpp::export]]
NumericVector depth_one_perm_single_cpp(NumericMatrix distmat,
                                        double c = 0.1,
                                        bool do_permute = true) {
  int n = distmat.nrow();

  NumericMatrix current_distmat = clone(distmat);

  if (do_permute) {
    IntegerVector ind = seq_len(n) - 1;
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::shuffle(ind.begin(), ind.end(), eng);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        current_distmat(i, k) = distmat(ind[i], ind[k]);
      }
    }
  }

  NumericVector test;
  for (int cp = std::ceil(n * c); cp <= n - std::ceil(n * c); cp++) {
    double testStat = getTcpp(current_distmat, seq_len(n), cp, n - cp);
    test.push_back(testStat);
  }

  double max_value = max(test);
  int loc = which_max(test) + static_cast<int>(std::ceil(n * c)) + 1;

  return NumericVector::create(
    Named("stat") = max_value,
    Named("loc") = loc
  );
}