#include <Rcpp.h>
#include <algorithm> // for std::sort and std::upper_bound

using namespace Rcpp;

// Function to calculate empirical distribution function
// [[Rcpp::export]]
NumericVector empDfCpp_dF(const NumericVector& x, const NumericVector& dSup) {
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
NumericVector colSumsCustom_dF(const NumericMatrix& mat) {
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
NumericMatrix squaredDiff_dF(const NumericMatrix& A, const NumericMatrix& B) {
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
double getTcpp_dF(const NumericMatrix& distmat, const IntegerVector& indices, int n, int m, 
                  double cut_off = 0) {
  
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
  
  // Initialize FXX, FYY, FXY, FYX with appropriate sizes
  NumericMatrix FXX(DummyD.nrow(), dXX.nrow());
  NumericMatrix FYY(DummyD.nrow(), dYY.nrow());
  NumericMatrix FXY(DummyD.nrow(), dXY.nrow());
  NumericMatrix FYX(DummyD.nrow(), dYX.nrow());
  
  // Populate FXX, FYY, FXY, FYX by operating on rows of dXX, dYY, dXY, and dYX
  for (int j = 0; j < dXX.nrow(); ++j) {
    FXX(_, j) = empDfCpp_dF(dXX(j, _), DummyD(_, j)); // Use j-th column of DummyD as dSup
  }
  for (int j = 0; j < dYY.nrow(); ++j) {
    FYY(_, j) = empDfCpp_dF(dYY(j, _), DummyD(_, j + n)); // Use (j+n)-th column of DummyD as dSup
  }
  for (int j = 0; j < dXY.nrow(); ++j) {
    FXY(_, j) = empDfCpp_dF(dXY(j, _), DummyD(_, j)); // Use (j+n)-th column of DummyD as dSup
  }
  for (int j = 0; j < dYX.nrow(); ++j) {
    FYX(_, j) = empDfCpp_dF(dYX(j, _), DummyD(_, j + n)); // Use j-th column of DummyD as dSup
  }
  
  // Calculate squared differences using the custom function
  NumericMatrix diffFXX_FXY = squaredDiff_dF(FXX, FXY);
  NumericMatrix diffFYY_FYX = squaredDiff_dF(FYY, FYX);
  
  // Calculate TFX and TFY with column sums multiplied by 1/(n+m)
  NumericVector TFX = colSumsCustom_dF(diffFXX_FXY) * (1.0 / (n + m));
  NumericVector TFY = colSumsCustom_dF(diffFYY_FYX) * (1.0 / (n + m));
  
  // Calculate the test statistic
  double teststat = (static_cast<double>(FXX.ncol()) * FYY.ncol() / 
                     (static_cast<double>(FXX.ncol()) + FYY.ncol())) * 
                     ((sum(TFX) + sum(TFY)) / (TFX.size() + TFY.size()));
  
  return teststat;
}
