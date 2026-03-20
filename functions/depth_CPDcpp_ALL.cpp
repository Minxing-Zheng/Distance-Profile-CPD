#include <Rcpp.h>
#include <algorithm>
#include <random>
#include <limits>
using namespace Rcpp;

// Function to calculate empirical distribution function
// [[Rcpp::export]]
NumericVector empDfCpp_All(const NumericVector& x, const NumericVector& dSup) {
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
NumericVector colSumsCpp_All(const NumericMatrix& mat) {
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
NumericMatrix squaredDiff_All(const NumericMatrix& A, const NumericMatrix& B) {
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

// Unified test function for dF, AD, and W
// [[Rcpp::export]]
List getTestStats_All(const NumericMatrix& distmat, const IntegerVector& indices, int n, int m) {
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
    FXX(_, j) = empDfCpp_All(dXX(j, _), DummyD(_, j)); // Use j-th column of DummyD as dSup
  }
  for (int j = 0; j < dYY.nrow(); ++j) {
    FYY(_, j) = empDfCpp_All(dYY(j, _), DummyD(_, j + n)); // Use (j+n)-th column of DummyD as dSup
  }
  for (int j = 0; j < dXY.nrow(); ++j) {
    FXY(_, j) = empDfCpp_All(dXY(j, _), DummyD(_, j)); // Use (j+n)-th column of DummyD as dSup
  }
  for (int j = 0; j < dYX.nrow(); ++j) {
    FYX(_, j) = empDfCpp_All(dYX(j, _), DummyD(_, j + n)); // Use j-th column of DummyD as dSup
  }
  
  // Calculate squared differences using the custom function
  NumericMatrix diffFXX_FXY = squaredDiff_All(FXX, FXY);
  NumericMatrix diffFYY_FYX = squaredDiff_All(FYY, FYX);
  
  // Calculate TFX and TFY with column sums multiplied by 1/(n+m)
  NumericVector TFX = colSumsCpp_All(diffFXX_FXY) * (1.0 / (n + m));
  NumericVector TFY = colSumsCpp_All(diffFYY_FYX) * (1.0 / (n + m));
  double stat_dF = (static_cast<double>(FXX.ncol()) * FYY.ncol() / 
                    (static_cast<double>(FXX.ncol()) + FYY.ncol())) * 
                    ((sum(TFX) + sum(TFY)) / (TFX.size() + TFY.size()));
  
  // Adaptive weights for stat_AD
  // Create Fpooled
  NumericMatrix Fpooled(DummyD.nrow(), DummyD.ncol());
  for (int j = 0; j < DummyD.ncol(); ++j) {
    Fpooled(_, j) = empDfCpp_All(DummyD(_, j), DummyD(_, j)); // Column-wise computation
  }
  
  // Create W = 1 / (Fpooled * (1 - Fpooled)) if > 0, else 0
  NumericMatrix W_AD(DummyD.nrow(), DummyD.ncol());
  for (int i = 0; i < Fpooled.nrow(); ++i) {
    for (int j = 0; j < Fpooled.ncol(); ++j) {
      double value = Fpooled(i, j) * (1.0 - Fpooled(i, j)); // Compute Fpooled * (1 - Fpooled)
      if (value > 0.0) {                                   // Check for positive value
        W_AD(i, j) = 1.0 / value;                             // Compute W
      } else {
        W_AD(i, j) = 0.0;                                     // Assign 0 otherwise
      }
    }
  }
  
  // Modify TFX: Element-wise multiplication with the first `n` columns of W, then take column sums
  NumericMatrix weighted_diffFXX_FXY(diffFXX_FXY.nrow(), diffFXX_FXY.ncol());
  for (int i = 0; i < diffFXX_FXY.nrow(); ++i) {
    for (int j = 0; j < diffFXX_FXY.ncol(); ++j) {
      weighted_diffFXX_FXY(i, j) = diffFXX_FXY(i, j) * W_AD(i, j); // Element-wise multiplication
    }
  }
  NumericVector TFX_AD = colSumsCpp_All(weighted_diffFXX_FXY) * (1.0 / static_cast<double>(n + m)); // Adjusted TFX
  
  // Modify TFY: Element-wise multiplication with the last `m` columns of W, then take column sums
  NumericMatrix weighted_diffFYY_FYX(diffFYY_FYX.nrow(), diffFYY_FYX.ncol());
  for (int i = 0; i < diffFYY_FYX.nrow(); ++i) {
    for (int j = 0; j < diffFYY_FYX.ncol(); ++j) {
      weighted_diffFYY_FYX(i, j) = diffFYY_FYX(i, j) * W_AD(i, j + n); // Element-wise multiplication
    }
  }
  NumericVector TFY_AD = colSumsCpp_All(weighted_diffFYY_FYX) * (1.0 / static_cast<double>(n + m)); // Adjusted TFY
  
  double stat_AD = (static_cast<double>(FXX.ncol()) * static_cast<double>(FYY.ncol()) / 
                    (static_cast<double>(FXX.ncol()) + static_cast<double>(FYY.ncol()))) * 
                    ((sum(TFX_AD) + sum(TFY_AD)) / static_cast<double>(TFX_AD.size() + TFY_AD.size()));
  
  // New adaptive weights for stat_W
  // Create W = 1 / Fpooled if > 0, else 0
  NumericMatrix W_W(DummyD.nrow(), DummyD.ncol());
  for (int i = 0; i < Fpooled.nrow(); ++i) {
    for (int j = 0; j < Fpooled.ncol(); ++j) {
      if (Fpooled(i, j) > 0.0) {                         // Check for positive Fpooled
        W_W(i, j) = 1.0 / Fpooled(i, j);                   // Compute W
      } else {
        W_W(i, j) = 0.0;                                   // Assign 0 otherwise
      }
    }
  }
  
  // Modify TFX: Element-wise multiplication with the first `n` columns of W, then take column sums
  NumericMatrix weighted_diffFXX_FXY_W(diffFXX_FXY.nrow(), diffFXX_FXY.ncol());
  for (int i = 0; i < diffFXX_FXY.nrow(); ++i) {
    for (int j = 0; j < diffFXX_FXY.ncol(); ++j) {
      weighted_diffFXX_FXY_W(i, j) = diffFXX_FXY(i, j) * W_W(i, j); // Element-wise multiplication
    }
  }
  NumericVector TFX_W = colSumsCpp_All(weighted_diffFXX_FXY_W) * (1.0 / static_cast<double>(n + m)); // Adjusted TFX
  
  // Modify TFY: Element-wise multiplication with the last `m` columns of W, then take column sums
  NumericMatrix weighted_diffFYY_FYX_W(diffFYY_FYX.nrow(), diffFYY_FYX.ncol());
  for (int i = 0; i < diffFYY_FYX.nrow(); ++i) {
    for (int j = 0; j < diffFYY_FYX.ncol(); ++j) {
      weighted_diffFYY_FYX_W(i, j) = diffFYY_FYX(i, j) * W_W(i, j + n); // Element-wise multiplication
    }
  }
  NumericVector TFY_W = colSumsCpp_All(weighted_diffFYY_FYX_W) * (1.0 / static_cast<double>(n + m)); // Adjusted TFY
  
  double stat_W = (static_cast<double>(FXX.ncol()) * static_cast<double>(FYY.ncol()) / 
                   (static_cast<double>(FXX.ncol()) + static_cast<double>(FYY.ncol()))) * 
                   ((sum(TFX_W) + sum(TFY_W)) / static_cast<double>(TFX_W.size() + TFY_W.size()));
  
  return List::create(
    Named("stat_dF") = stat_dF,
    Named("stat_AD") = stat_AD,
    Named("stat_W") = stat_W
  );
}

// Permutation-based p-values for dF, AD, and W
// [[Rcpp::export]]
List depth_CPDcpp_ALL(NumericMatrix distmat, double c = 0.1, int num_permut = 500) {
  int n = distmat.nrow(); // Number of rows in the distance matrix
  
  // Compute the number of change points
  int num_cp = n - 2 * ceil(n * c) + 1; // Total number of change points
  
  // Initialize matrices and arrays to store results
  NumericVector observed_stat(3); // Observed test statistics (stat_dF, stat_AD, stat_W)
  NumericMatrix max_stat(3, num_permut + 1); // Row-wise maxima for each permutation (3 rows: stat_dF, stat_AD, stat_W)
  IntegerMatrix max_stat_index(3, num_permut + 1); // Indices of the row-wise maxima for each permutation
  IntegerVector count(3); // Counts for each statistic (stat_dF, stat_AD, stat_W)
  NumericVector p_val(3); // P-values for each statistic
  IntegerVector loc(3);   // Locations of detected change points (for observed statistics)
  
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
    
    // Initialize the test matrix for this permutation
    NumericMatrix test(3, num_cp); // 3 rows (stat_dF, stat_AD, stat_W), num_cp columns
    
    int col_index = 0; // Column index for `test`
    for (int cp = ceil(n * c); cp <= n - ceil(n * c); cp++) {
      // Call the updated getTestStats_All function
      List testStats = getTestStats_All(current_distmat, seq_len(n), cp, n - cp); // Returns a list of 3 values
      test(0, col_index) = as<double>(testStats["stat_dF"]); // First row for stat_dF
      test(1, col_index) = as<double>(testStats["stat_AD"]); // Second row for stat_AD
      test(2, col_index) = as<double>(testStats["stat_W"]);  // Third row for stat_W
      col_index++;
    }
    
    // Compute row-wise maxima and indices for the current permutation
    for (int row = 0; row < 3; row++) {
      max_stat(row, j) = max(test(row, _)); // Maximum value for this row
      max_stat_index(row, j) = which_max(test(row, _)) + static_cast<int>(ceil(n * c)); // Adjust index for change points
    }
  }
  
  // Process observed statistics (j = 0)
  for (int row = 0; row < 3; row++) {
    observed_stat[row] = max_stat(row, 0); // Observed statistic is the first column
    loc[row] = max_stat_index(row, 0);    // Location of the observed change point
  }
  
  // Calculate p-values for each statistic
  for (int row = 0; row < 3; row++) {
    count[row] = 0;
    for (int j = 1; j <= num_permut; j++) {
      if (max_stat(row, j) > observed_stat[row]) {
        count[row]++;
      }
    }
    p_val[row] = (1.0 + static_cast<double>(count[row])) / (1.0 + static_cast<double>(num_permut)); // Ensure numeric division
  }
  
  // Return the results
  return List::create(
    Named("p_val") = p_val,               // P-values for stat_dF, stat_AD, and stat_W
    Named("loc") = loc,                   // Locations of detected change points for each statistic
    Named("observed_stat") = observed_stat, // Observed test statistics
    Named("max_stat") = max_stat,         // Row-wise maxima for each statistic
    Named("max_stat_index") = max_stat_index // Indices of row-wise maxima
  );
}

// assumes getTestStats_All(...) is already available in another sourced file,
// OR you can paste its definition into this file as well.

// [[Rcpp::export]]
NumericVector depth_one_perm_cpp(NumericMatrix distmat,
                                 double c = 0.1,
                                 bool do_permute = true) {
  int n = distmat.nrow();
  int cp_min = std::ceil(n * c);
  int cp_max = n - std::ceil(n * c);

  NumericMatrix current_distmat = clone(distmat);

  if (do_permute) {
    IntegerVector ind = seq(0, n - 1);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::shuffle(ind.begin(), ind.end(), eng);

    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        current_distmat(i, k) = distmat(ind[i], ind[k]);
      }
    }
  }

  double best_dF = -std::numeric_limits<double>::infinity();
  double best_AD = -std::numeric_limits<double>::infinity();
  double best_W  = -std::numeric_limits<double>::infinity();

  int loc_dF = cp_min;
  int loc_AD = cp_min;
  int loc_W  = cp_min;

  IntegerVector idx = seq_len(n);

  for (int cp = cp_min; cp <= cp_max; cp++) {
    List testStats = getTestStats_All(current_distmat, idx, cp, n - cp);

    double s_dF = as<double>(testStats["stat_dF"]);
    double s_AD = as<double>(testStats["stat_AD"]);
    double s_W  = as<double>(testStats["stat_W"]);

    if (s_dF > best_dF) {
      best_dF = s_dF;
      loc_dF = cp;
    }
    if (s_AD > best_AD) {
      best_AD = s_AD;
      loc_AD = cp;
    }
    if (s_W > best_W) {
      best_W = s_W;
      loc_W = cp;
    }
  }

  return NumericVector::create(
    Named("stat_dF") = best_dF,
    Named("stat_AD") = best_AD,
    Named("stat_W")  = best_W,
    Named("loc_dF")  = loc_dF,
    Named("loc_AD")  = loc_AD,
    Named("loc_W")   = loc_W
  );
}