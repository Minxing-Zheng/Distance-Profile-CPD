#include <Rcpp.h>
#include <algorithm> // for std::sort and std::upper_bound

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
