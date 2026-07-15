// Combined distance-profile change-point detection kernel.
//
// Instead of rebuilding sortedXX/sortedXY/sortedYY/sortedYX from scratch at every
// candidate change point (full O(n log n) sort per column, per candidate), this
// version maintains them incrementally: as the candidate change point advances by
// one position, exactly one point moves from the "right" group to the "left"
// group, so each affected column's sorted comparison array only needs a single
// insertion or removal (O(n) with a std::vector, due to shifting) instead of a
// full resort (O(n log n)).
//
// This removes the log(n) factor from the array-construction cost only. The
// CDF-evaluation step still costs O(n^2 log n) per candidate cp, so this is a
// bounded speedup rather than an asymptotic change.

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

inline double emp_cdf_scalar(const std::vector<double>& sorted_x, double t) {
  const int n = sorted_x.size();
  const int ind = std::upper_bound(sorted_x.begin(), sorted_x.end(), t) - sorted_x.begin();
  return static_cast<double>(ind) / static_cast<double>(n);
}

inline void matrix_min_max(const NumericMatrix& D, double& min_val, double& max_val) {
  const int N = D.nrow();
  min_val = D(0, 0);
  max_val = D(0, 0);

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (D(i, j) < min_val) min_val = D(i, j);
      if (D(i, j) > max_val) max_val = D(i, j);
    }
  }
}

inline void sorted_insert(std::vector<double>& v, double val) {
  v.insert(std::upper_bound(v.begin(), v.end(), val), val);
}

// Erases exactly one occurrence equal to val. Caller guarantees val is present
// (it was inserted earlier by sorted_insert, or was part of the initial build).
inline void sorted_erase_value(std::vector<double>& v, double val) {
  std::vector<double>::iterator it = std::lower_bound(v.begin(), v.end(), val);
  v.erase(it);
}

inline std::vector<double> build_excl_self(const NumericMatrix& D,
                                           int idx,
                                           const std::vector<int>& group) {
  std::vector<double> out;
  out.reserve(group.size() > 0 ? group.size() - 1 : 0);
  for (std::size_t t = 0; t < group.size(); ++t) {
    const int k = group[t];
    if (k != idx) {
      out.push_back(D(idx, k));
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

inline std::vector<double> build_cross(const NumericMatrix& D,
                                       int idx,
                                       const std::vector<int>& other_group) {
  std::vector<double> out(other_group.size());
  for (std::size_t t = 0; t < other_group.size(); ++t) {
    out[t] = D(idx, other_group[t]);
  }
  std::sort(out.begin(), out.end());
  return out;
}

// dF/AD/W's "own-group" arrays include the self distance (D(idx,idx) == 0 for a
// valid distance matrix), while the uniform variant's convention excludes it.
// Rather than maintaining two incremental copies, we derive the self-included
// version on demand by inserting a single 0.0 into the (maintained) self-excluded
// array -- cheap (O(n)) relative to the O(n log n) this design already avoids.
inline std::vector<double> with_zero_inserted(const std::vector<double>& v) {
  std::vector<double> out(v);
  sorted_insert(out, 0.0);
  return out;
}

// Incremental state for one permutation's cp-scan. Cross-group arrays (XY, YX)
// are identical between the uniform and dF/AD/W conventions (both are plain
// cross-group distance sets, symmetric in D), so only one copy of each is kept.
struct CPState {
  int N;
  int cp; // current n_left
  std::vector<std::vector<double> > coreXX; // valid for idx in left: self-excluded own group
  std::vector<std::vector<double> > coreYY; // valid for idx in right: self-excluded own group
  std::vector<std::vector<double> > coreXY; // valid for idx in left: cross to right
  std::vector<std::vector<double> > coreYX; // valid for idx in right: cross to left
};

inline void init_state(const NumericMatrix& D, int first_cp, CPState& st) {
  const int N = D.nrow();
  st.N = N;
  st.cp = first_cp;
  st.coreXX.assign(N, std::vector<double>());
  st.coreYY.assign(N, std::vector<double>());
  st.coreXY.assign(N, std::vector<double>());
  st.coreYX.assign(N, std::vector<double>());

  std::vector<int> left(first_cp);
  std::vector<int> right(N - first_cp);
  for (int i = 0; i < first_cp; ++i) left[i] = i;
  for (int i = first_cp; i < N; ++i) right[i - first_cp] = i;

  for (std::size_t t = 0; t < left.size(); ++t) {
    const int idx = left[t];
    st.coreXX[idx] = build_excl_self(D, idx, left);
    st.coreXY[idx] = build_cross(D, idx, right);
  }
  for (std::size_t t = 0; t < right.size(); ++t) {
    const int idx = right[t];
    st.coreYY[idx] = build_excl_self(D, idx, right);
    st.coreYX[idx] = build_cross(D, idx, left);
  }
}

// Moves point m = st.cp from right to left; st.cp becomes st.cp + 1.
inline void advance_state(const NumericMatrix& D, CPState& st) {
  const int m = st.cp;
  const int N = st.N;
  const int old_n_left = st.cp;

  for (int j = 0; j < old_n_left; ++j) {
    const double d = D(j, m);
    sorted_insert(st.coreXX[j], d);
    sorted_erase_value(st.coreXY[j], d);
  }
  for (int r = m + 1; r < N; ++r) {
    const double d = D(r, m);
    sorted_erase_value(st.coreYY[r], d);
    sorted_insert(st.coreYX[r], d);
  }

  std::vector<int> new_left(old_n_left + 1);
  for (int i = 0; i <= old_n_left; ++i) new_left[i] = i;
  std::vector<int> new_right(N - old_n_left - 1);
  for (int i = m + 1; i < N; ++i) new_right[i - m - 1] = i;

  st.coreXX[m] = build_excl_self(D, m, new_left);
  st.coreXY[m] = build_cross(D, m, new_right);
  // st.coreYY[m] / st.coreYX[m] are now stale and unused (m left the right group).

  st.cp = old_n_left + 1;
}

inline double uniform_stat_from_state(const CPState& st,
                                      int ndSup,
                                      double min_val,
                                      double max_val) {
  const int n_left = st.cp;
  const int N = st.N;
  const int n_right = N - n_left;

  if (ndSup < 2 || max_val <= min_val) {
    return 0.0;
  }

  const double step = (max_val - min_val) / static_cast<double>(ndSup - 1);
  double sumX = 0.0;
  double sumY = 0.0;

  for (int j = 0; j < n_left; ++j) {
    const std::vector<double>& xx = st.coreXX[j];
    const std::vector<double>& xy = st.coreXY[j];
    for (int r = 0; r < ndSup; ++r) {
      const double t = min_val + r * step;
      const double diff = emp_cdf_scalar(xx, t) - emp_cdf_scalar(xy, t);
      sumX += diff * diff;
    }
  }

  for (int j = n_left; j < N; ++j) {
    const std::vector<double>& yy = st.coreYY[j];
    const std::vector<double>& yx = st.coreYX[j];
    for (int r = 0; r < ndSup; ++r) {
      const double t = min_val + r * step;
      const double diff = emp_cdf_scalar(yy, t) - emp_cdf_scalar(yx, t);
      sumY += diff * diff;
    }
  }

  const double factor = static_cast<double>(n_left) * static_cast<double>(n_right) /
    static_cast<double>(N);

  return factor * ((sumX + sumY) * step / static_cast<double>(N));
}

inline void all_stats_from_state(const NumericMatrix& D,
                                 const CPState& st,
                                 bool need_dF,
                                 bool need_AD,
                                 bool need_W,
                                 const std::vector<std::vector<double> >& sortedPool,
                                 double& stat_dF,
                                 double& stat_AD,
                                 double& stat_W) {
  const int n_left = st.cp;
  const int N = st.N;
  const bool need_weight = need_AD || need_W;

  double sum_dF = 0.0;
  double sum_AD = 0.0;
  double sum_W = 0.0;

  for (int j = 0; j < n_left; ++j) {
    const std::vector<double> xx_incl = with_zero_inserted(st.coreXX[j]);
    const std::vector<double>& xy = st.coreXY[j];
    for (int i = 0; i < N; ++i) {
      const double Dij = D(i, j);
      double diff = emp_cdf_scalar(xx_incl, Dij) - emp_cdf_scalar(xy, Dij);
      diff *= diff;
      if (need_dF) {
        sum_dF += diff;
      }
      if (need_weight) {
        const double f = emp_cdf_scalar(sortedPool[j], Dij);
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

  for (int j = n_left; j < N; ++j) {
    const std::vector<double> yy_incl = with_zero_inserted(st.coreYY[j]);
    const std::vector<double>& yx = st.coreYX[j];
    for (int i = 0; i < N; ++i) {
      const double Dij = D(i, j);
      double diff = emp_cdf_scalar(yy_incl, Dij) - emp_cdf_scalar(yx, Dij);
      diff *= diff;
      if (need_dF) {
        sum_dF += diff;
      }
      if (need_weight) {
        const double f = emp_cdf_scalar(sortedPool[j], Dij);
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

  const double factor = static_cast<double>(n_left) * static_cast<double>(N - n_left) /
    static_cast<double>(N);
  const double norm = static_cast<double>(N) * static_cast<double>(N);

  stat_dF = need_dF ? factor * (sum_dF / norm) : R_NegInf;
  stat_AD = need_AD ? factor * (sum_AD / norm) : R_NegInf;
  stat_W = need_W ? factor * (sum_W / norm) : R_NegInf;
}

// [[Rcpp::export]]
List distCPD_combined(const NumericMatrix& distmat,
                      double c = 0.1,
                      int num_permut = 500,
                      int ndSup = 200,
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

  double min_val = 0.0;
  double max_val = 0.0;
  if (need_uniform) {
    matrix_min_max(distmat, min_val, max_val);
  }

  const bool need_weight = need_AD || need_W;

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

  std::vector<std::vector<double> > sortedPool;
  if (need_weight) {
    sortedPool.resize(n);
    for (int j = 0; j < n; ++j) {
      sortedPool[j].resize(n);
    }
  }

  CPState st;

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

    if (need_weight) {
      for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
          sortedPool[j][i] = current_distmat(i, j);
        }
        std::sort(sortedPool[j].begin(), sortedPool[j].end());
      }
    }

    init_state(current_distmat, first_cp, st);

    double best_uniform = R_NegInf;
    double best_dF = R_NegInf;
    double best_AD = R_NegInf;
    double best_W = R_NegInf;
    int best_uniform_cp = first_cp;
    int best_dF_cp = first_cp;
    int best_AD_cp = first_cp;
    int best_W_cp = first_cp;

    for (int cp = first_cp; cp <= last_cp; ++cp) {
      CombinedStats stats;
      stats.uniform = need_uniform ?
        uniform_stat_from_state(st, ndSup, min_val, max_val) : R_NegInf;
      if (need_dF || need_AD || need_W) {
        all_stats_from_state(
          current_distmat, st,
          need_dF, need_AD, need_W,
          sortedPool,
          stats.dF, stats.AD, stats.W
        );
      } else {
        stats.dF = R_NegInf;
        stats.AD = R_NegInf;
        stats.W = R_NegInf;
      }

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

      if (cp < last_cp) {
        advance_state(current_distmat, st);
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
