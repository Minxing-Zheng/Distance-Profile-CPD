#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <random>

using namespace Rcpp;

struct CombinedStats {
  double uniform;
  double dF;
  double AD;
  double W;
};

inline double emp_cdf_scalar(const NumericVector& sorted_x, double t) {
  const int n = sorted_x.size();
  const int ind = std::upper_bound(sorted_x.begin(), sorted_x.end(), t) - sorted_x.begin();
  return static_cast<double>(ind) / static_cast<double>(n);
}

inline NumericVector emp_cdf_column(const NumericVector& sorted_x,
                                    const NumericMatrix& D,
                                    int col) {
  const int n = sorted_x.size();
  const int N = D.nrow();
  NumericVector out(N);

  for (int i = 0; i < N; ++i) {
    const int ind = std::upper_bound(sorted_x.begin(), sorted_x.end(), D(i, col)) - sorted_x.begin();
    out[i] = static_cast<double>(ind) / static_cast<double>(n);
  }

  return out;
}

inline NumericVector column_sums(const NumericMatrix& mat) {
  const int nr = mat.nrow();
  const int nc = mat.ncol();
  NumericVector out(nc);

  for (int j = 0; j < nc; ++j) {
    for (int i = 0; i < nr; ++i) {
      out[j] += mat(i, j);
    }
  }

  return out;
}

inline double uniform_stat_from_D(const NumericMatrix& D,
                                  int n_left,
                                  int n_right,
                                  int ndSup) {
  const int N = n_left + n_right;

  double min_val = D(0, 0);
  double max_val = D(0, 0);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (D(i, j) < min_val) min_val = D(i, j);
      if (D(i, j) > max_val) max_val = D(i, j);
    }
  }

  if (ndSup < 2 || max_val <= min_val) {
    return 0.0;
  }

  const double step = (max_val - min_val) / static_cast<double>(ndSup - 1);

  std::vector<NumericVector> sortedXX(n_left);
  std::vector<NumericVector> sortedXY(n_left);
  std::vector<NumericVector> sortedYY(n_right);
  std::vector<NumericVector> sortedYX(n_right);

  for (int j = 0; j < n_left; ++j) {
    sortedXX[j] = NumericVector(n_left - 1);
    int idx = 0;
    for (int k = 0; k < n_left; ++k) {
      if (k != j) {
        sortedXX[j][idx++] = D(k, j);
      }
    }
    std::sort(sortedXX[j].begin(), sortedXX[j].end());

    sortedXY[j] = NumericVector(n_right);
    for (int k = 0; k < n_right; ++k) {
      sortedXY[j][k] = D(n_left + k, j);
    }
    std::sort(sortedXY[j].begin(), sortedXY[j].end());
  }

  for (int j = 0; j < n_right; ++j) {
    sortedYY[j] = NumericVector(n_right - 1);
    int idx = 0;
    for (int k = 0; k < n_right; ++k) {
      if (k != j) {
        sortedYY[j][idx++] = D(n_left + k, n_left + j);
      }
    }
    std::sort(sortedYY[j].begin(), sortedYY[j].end());

    sortedYX[j] = NumericVector(n_left);
    for (int k = 0; k < n_left; ++k) {
      sortedYX[j][k] = D(k, n_left + j);
    }
    std::sort(sortedYX[j].begin(), sortedYX[j].end());
  }

  double sumX = 0.0;
  double sumY = 0.0;

  for (int j = 0; j < n_left; ++j) {
    for (int r = 0; r < ndSup; ++r) {
      const double t = min_val + r * step;
      const double diff = emp_cdf_scalar(sortedXX[j], t) - emp_cdf_scalar(sortedXY[j], t);
      sumX += diff * diff;
    }
  }

  for (int j = 0; j < n_right; ++j) {
    for (int r = 0; r < ndSup; ++r) {
      const double t = min_val + r * step;
      const double diff = emp_cdf_scalar(sortedYY[j], t) - emp_cdf_scalar(sortedYX[j], t);
      sumY += diff * diff;
    }
  }

  const double factor = static_cast<double>(n_left) * static_cast<double>(n_right) /
    static_cast<double>(N);

  return factor * ((sumX + sumY) * step / static_cast<double>(N));
}

inline void all_stats_from_D(const NumericMatrix& D,
                             int n_left,
                             int n_right,
                             bool need_dF,
                             bool need_AD,
                             bool need_W,
                             double& stat_dF,
                             double& stat_AD,
                             double& stat_W) {
  const int N = n_left + n_right;

  std::vector<NumericVector> sortedXX(n_left);
  std::vector<NumericVector> sortedXY(n_left);
  std::vector<NumericVector> sortedYY(n_right);
  std::vector<NumericVector> sortedYX(n_right);
  const bool need_weight = need_AD || need_W;
  std::vector<NumericVector> sortedPool;
  if (need_weight) {
    sortedPool.resize(N);
  }

  for (int j = 0; j < n_left; ++j) {
    sortedXX[j] = NumericVector(n_left);
    for (int k = 0; k < n_left; ++k) {
      sortedXX[j][k] = D(j, k);
    }
    std::sort(sortedXX[j].begin(), sortedXX[j].end());

    sortedXY[j] = NumericVector(n_right);
    for (int k = 0; k < n_right; ++k) {
      sortedXY[j][k] = D(j, n_left + k);
    }
    std::sort(sortedXY[j].begin(), sortedXY[j].end());
  }

  for (int j = 0; j < n_right; ++j) {
    sortedYY[j] = NumericVector(n_right);
    for (int k = 0; k < n_right; ++k) {
      sortedYY[j][k] = D(n_left + j, n_left + k);
    }
    std::sort(sortedYY[j].begin(), sortedYY[j].end());

    sortedYX[j] = NumericVector(n_left);
    for (int k = 0; k < n_left; ++k) {
      sortedYX[j][k] = D(n_left + j, k);
    }
    std::sort(sortedYX[j].begin(), sortedYX[j].end());
  }

  if (need_weight) {
    for (int j = 0; j < N; ++j) {
      sortedPool[j] = NumericVector(N);
      for (int i = 0; i < N; ++i) {
        sortedPool[j][i] = D(i, j);
      }
      std::sort(sortedPool[j].begin(), sortedPool[j].end());
    }
  }

  double sum_dF = 0.0;
  double sum_AD = 0.0;
  double sum_W = 0.0;

  for (int j = 0; j < n_left; ++j) {
    for (int i = 0; i < N; ++i) {
      double diff = emp_cdf_scalar(sortedXX[j], D(i, j)) -
        emp_cdf_scalar(sortedXY[j], D(i, j));
      diff *= diff;
      if (need_dF) {
        sum_dF += diff;
      }

      if (need_weight) {
        const double f = emp_cdf_scalar(sortedPool[j], D(i, j));
        if (f > 0.0) {
          if (need_W) {
            sum_W += diff / f;
          }
          if (need_AD && f < 1.0) {
            sum_AD += diff / (f * (1.0 - f));
          }
        }
      }
    }
  }

  for (int j = 0; j < n_right; ++j) {
    const int col = j + n_left;
    for (int i = 0; i < N; ++i) {
      double diff = emp_cdf_scalar(sortedYY[j], D(i, col)) -
        emp_cdf_scalar(sortedYX[j], D(i, col));
      diff *= diff;
      if (need_dF) {
        sum_dF += diff;
      }

      if (need_weight) {
        const double f = emp_cdf_scalar(sortedPool[col], D(i, col));
        if (f > 0.0) {
          if (need_W) {
            sum_W += diff / f;
          }
          if (need_AD && f < 1.0) {
            sum_AD += diff / (f * (1.0 - f));
          }
        }
      }
    }
  }

  const double factor = static_cast<double>(n_left) * static_cast<double>(n_right) /
    static_cast<double>(N);
  const double norm = static_cast<double>(N) * static_cast<double>(N);

  stat_dF = need_dF ? factor * (sum_dF / norm) : R_NegInf;
  stat_AD = need_AD ? factor * (sum_AD / norm) : R_NegInf;
  stat_W = need_W ? factor * (sum_W / norm) : R_NegInf;
}

inline CombinedStats combined_stats_for_cp(const NumericMatrix& current_distmat,
                                           int cp,
                                           int ndSup,
                                           bool need_uniform,
                                           bool need_dF,
                                           bool need_AD,
                                           bool need_W) {
  const int N = current_distmat.nrow();

  CombinedStats out;
  out.uniform = need_uniform ? uniform_stat_from_D(current_distmat, cp, N - cp, ndSup) : R_NegInf;
  if (need_dF || need_AD || need_W) {
    all_stats_from_D(
      current_distmat, cp, N - cp,
      need_dF, need_AD, need_W,
      out.dF, out.AD, out.W
    );
  } else {
    out.dF = R_NegInf;
    out.AD = R_NegInf;
    out.W = R_NegInf;
  }
  return out;
}

// [[Rcpp::export]]
List depth_CPD_combined(const NumericMatrix& distmat,
                        double c = 0.1,
                        int num_permut = 500,
                        int ndSup = 1000,
                        int seed = -1,
                        Nullable<LogicalVector> variants = R_NilValue) {
  const int n = distmat.nrow();

  if (distmat.ncol() != n) {
    stop("distmat must be square.");
  }
  if (c <= 0.0 || c >= 0.5) {
    stop("c must be in (0, 0.5).");
  }
  if (num_permut < 0) {
    stop("num_permut must be nonnegative.");
  }
  if (ndSup < 2) {
    stop("ndSup must be at least 2.");
  }
  LogicalVector variant_mask = variants.isNull() ?
    LogicalVector::create(true, true, true, true) :
    LogicalVector(variants);

  if (variant_mask.size() != 4) {
    stop("variants must be a logical vector of length 4.");
  }

  const bool need_uniform = LogicalVector::is_na(variant_mask[0]) ? false : static_cast<bool>(variant_mask[0]);
  const bool need_dF = LogicalVector::is_na(variant_mask[1]) ? false : static_cast<bool>(variant_mask[1]);
  const bool need_AD = LogicalVector::is_na(variant_mask[2]) ? false : static_cast<bool>(variant_mask[2]);
  const bool need_W = LogicalVector::is_na(variant_mask[3]) ? false : static_cast<bool>(variant_mask[3]);

  if (!need_uniform && !need_dF && !need_AD && !need_W) {
    stop("At least one variant must be requested.");
  }

  const int first_cp = static_cast<int>(std::ceil(n * c));
  const int last_cp = n - first_cp;

  NumericMatrix max_stat(4, num_permut + 1);
  IntegerMatrix max_stat_index(4, num_permut + 1);
  NumericVector observed_stat(4);
  NumericVector p_val(4);
  IntegerVector loc(4);

  NumericMatrix current_distmat(n, n);
  std::vector<int> perm(n);
  std::iota(perm.begin(), perm.end(), 0);

  std::mt19937 eng;
  if (seed >= 0) {
    eng.seed(static_cast<std::mt19937::result_type>(seed));
  } else {
    std::random_device rd;
    eng.seed(rd());
  }

  for (int b = 0; b <= num_permut; ++b) {
    if (b == 0) {
      std::copy(distmat.begin(), distmat.end(), current_distmat.begin());
    } else {
      std::iota(perm.begin(), perm.end(), 0);
      std::shuffle(perm.begin(), perm.end(), eng);

      for (int i = 0; i < n; ++i) {
        const int ii = perm[i];
        for (int j = 0; j < n; ++j) {
          current_distmat(i, j) = distmat(ii, perm[j]);
        }
      }
    }

    double best_uniform = R_NegInf;
    double best_dF = R_NegInf;
    double best_AD = R_NegInf;
    double best_W = R_NegInf;
    int best_uniform_cp = first_cp;
    int best_dF_cp = first_cp;
    int best_AD_cp = first_cp;
    int best_W_cp = first_cp;

    for (int cp = first_cp; cp <= last_cp; ++cp) {
      const CombinedStats stats = combined_stats_for_cp(
        current_distmat, cp, ndSup,
        need_uniform, need_dF, need_AD, need_W
      );

      if (stats.uniform > best_uniform) {
        best_uniform = stats.uniform;
        best_uniform_cp = cp;
      }
      if (stats.dF > best_dF) {
        best_dF = stats.dF;
        best_dF_cp = cp;
      }
      if (stats.AD > best_AD) {
        best_AD = stats.AD;
        best_AD_cp = cp;
      }
      if (stats.W > best_W) {
        best_W = stats.W;
        best_W_cp = cp;
      }
    }

    max_stat(0, b) = best_uniform;
    max_stat(1, b) = best_dF;
    max_stat(2, b) = best_AD;
    max_stat(3, b) = best_W;
    max_stat_index(0, b) = best_uniform_cp;
    max_stat_index(1, b) = best_dF_cp;
    max_stat_index(2, b) = best_AD_cp;
    max_stat_index(3, b) = best_W_cp;
  }

  for (int k = 0; k < 4; ++k) {
    observed_stat[k] = max_stat(k, 0);
    loc[k] = max_stat_index(k, 0);

    int count = 0;
    for (int b = 1; b <= num_permut; ++b) {
      count += (max_stat(k, b) >= observed_stat[k]);
    }
    p_val[k] = (1.0 + static_cast<double>(count)) /
      (1.0 + static_cast<double>(num_permut));
  }

  CharacterVector method_names = CharacterVector::create(
    "dist_cpd_uniform", "dist_cpd", "dist_cpd_AD", "dist_cpd_W"
  );

  return List::create(
    Named("method") = method_names,
    Named("p_val") = p_val,
    Named("loc") = loc,
    Named("observed_stat") = observed_stat,
    Named("max_stat") = max_stat,
    Named("max_stat_index") = max_stat_index,
    Named("c") = c,
    Named("num_permut") = num_permut,
    Named("ndSup") = ndSup
  );
}
