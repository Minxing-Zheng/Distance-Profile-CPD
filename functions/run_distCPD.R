# Shared wrappers for distance-profile change-point scans.
#
# Main entry point:
#   run_distCPD(distmat, c = 0.1, num_permut = 500)
#
# Input:
#   distmat: n by n pairwise distance matrix.
#   c: boundary cutoff. Candidate locations are in [ceiling(n*c), n-ceiling(n*c)].
#   num_permut: number of permutation replicates.
#   ndSup: number of grid points for the uniform-integral statistic.
#   variants: "all" or any subset of dist_cpd_uniform, dist_cpd, dist_cpd_AD,
#             and dist_cpd_W.
#
# Output:
#   A named list for dist_cpd_uniform, dist_cpd, dist_cpd_AD, dist_cpd_W.
#   Each method contains test_stat, permuted_test_stat, critical_value, p_val,
#   reject, loc, candidate_loc, raw scan statistics, settings, and runtime_sec.

.distCPD_runner_env <- new.env(parent = emptyenv())
.distCPD_runner_env$cpp_loaded <- FALSE

.find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, "functions"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find project root containing a functions/ directory.")
    }
    path <- parent
  }
}

.validate_distmat <- function(distmat) {
  if (!is.matrix(distmat)) {
    distmat <- as.matrix(distmat)
  }
  if (nrow(distmat) != ncol(distmat)) {
    stop("distmat must be a square matrix.")
  }
  if (anyNA(distmat)) {
    stop("distmat contains NA values.")
  }
  storage.mode(distmat) <- "double"
  distmat
}

.normalize_dist_profile_variants <- function(variants = "all") {
  method_names <- c("dist_cpd_uniform", "dist_cpd", "dist_cpd_AD", "dist_cpd_W")
  aliases <- c(
    uniform = "dist_cpd_uniform",
    basic = "dist_cpd",
    df = "dist_cpd",
    dF = "dist_cpd",
    AD = "dist_cpd_AD",
    ad = "dist_cpd_AD",
    W = "dist_cpd_W",
    w = "dist_cpd_W"
  )

  if (length(variants) == 1 && identical(variants, "all")) {
    return(method_names)
  }
  if (length(variants) == 0) {
    stop("variants must contain at least one method.")
  }

  out <- character(length(variants))
  for (i in seq_along(variants)) {
    variant <- variants[[i]]
    if (variant %in% method_names) {
      out[[i]] <- variant
    } else if (variant %in% names(aliases)) {
      out[[i]] <- unname(aliases[[variant]])
    } else {
      stop(
        "Unknown dist-profile variant: ", variant,
        ". Use 'all' or one of: ", paste(method_names, collapse = ", ")
      )
    }
  }

  unique(out)
}

load_distCPD_cpp <- function(project_root = NULL, force = FALSE) {
  if (isTRUE(.distCPD_runner_env$cpp_loaded) && !force) {
    return(invisible(TRUE))
  }
  if (is.null(project_root)) {
    project_root <- .find_project_root()
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' is required.")
  }

  function_dir <- file.path(project_root, "functions")
  Rcpp::sourceCpp(file.path(function_dir, "distCPD_combined.cpp"), rebuild = force)
  .distCPD_runner_env$cpp_loaded <- TRUE
  invisible(TRUE)
}

.make_method_result <- function(method, stats_matrix, alpha, c, num_permut, ndSup, runtime_sec) {
  stat_name <- method
  loc_name <- paste0("loc_", method)
  test_stat <- unname(stats_matrix[1, stat_name])
  permuted_test_stat <- if (nrow(stats_matrix) > 1) {
    unname(stats_matrix[-1, stat_name])
  } else {
    numeric(0)
  }
  p_val <- if (length(permuted_test_stat) == 0) {
    NA_real_
  } else {
    (1 + sum(permuted_test_stat >= test_stat)) / (1 + length(permuted_test_stat))
  }
  critical_value <- if (length(permuted_test_stat) == 0) {
    NA_real_
  } else {
    unname(stats::quantile(permuted_test_stat, probs = 1 - alpha, type = 1, names = FALSE))
  }
  reject <- if (is.na(p_val)) NA else p_val <= alpha
  candidate_loc <- unname(stats_matrix[1, loc_name])
  loc <- if (isTRUE(reject)) candidate_loc else NA_integer_

  list(
    method = method,
    test_stat = test_stat,
    permuted_test_stat = permuted_test_stat,
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = loc,
    candidate_loc = candidate_loc,
    observed_stat = test_stat,
    permuted_stats = permuted_test_stat,
    all_stats = unname(stats_matrix[, stat_name]),
    detected = reject,
    c = c,
    num_permut = num_permut,
    ndSup = ndSup,
    alpha = alpha,
    runtime_sec = runtime_sec
  )
}

run_distCPD <- function(distmat,
                        c = 0.1,
                        num_permut = 500,
                        alpha = 0.05,
                        seed = NULL,
                        project_root = NULL,
                        progress = TRUE,
                        force_compile = FALSE,
                        ndSup = 1000,
                        variants = "all") {
  distmat <- .validate_distmat(distmat)
  variants <- .normalize_dist_profile_variants(variants)
  n <- nrow(distmat)
  if (c <= 0 || c >= 0.5) {
    stop("c must be in (0, 0.5).")
  }
  if (num_permut < 0 || num_permut != as.integer(num_permut)) {
    stop("num_permut must be a nonnegative integer.")
  }
  if (ndSup < 2 || ndSup != as.integer(ndSup)) {
    stop("ndSup must be an integer greater than or equal to 2.")
  }

  load_distCPD_cpp(project_root = project_root, force = force_compile)

  runtime <- system.time({
    if (isTRUE(progress)) {
      message(
        "running combined dist-profile permutations: 0/", num_permut,
        " variants=", paste(variants, collapse = ",")
      )
    }

    method_names <- c("dist_cpd_uniform", "dist_cpd", "dist_cpd_AD", "dist_cpd_W")
    variant_mask <- method_names %in% variants
    res_combined <- distCPD_combined(
      distmat = distmat,
      c = c,
      num_permut = num_permut,
      ndSup = ndSup,
      seed = if (is.null(seed)) -1L else as.integer(seed),
      variants = variant_mask
    )

    stats_matrix <- matrix(
      NA_real_,
      nrow = num_permut + 1,
      ncol = 8,
      dimnames = list(
        NULL,
        c(
          "dist_cpd_uniform", "dist_cpd", "dist_cpd_AD", "dist_cpd_W",
          "loc_dist_cpd_uniform", "loc_dist_cpd", "loc_dist_cpd_AD", "loc_dist_cpd_W"
        )
      )
    )

    stats_matrix[, "dist_cpd_uniform"] <- as.numeric(res_combined$max_stat[1, ])
    stats_matrix[, "dist_cpd"] <- as.numeric(res_combined$max_stat[2, ])
    stats_matrix[, "dist_cpd_AD"] <- as.numeric(res_combined$max_stat[3, ])
    stats_matrix[, "dist_cpd_W"] <- as.numeric(res_combined$max_stat[4, ])
    stats_matrix[1, "loc_dist_cpd_uniform"] <- res_combined$loc[1]
    stats_matrix[1, "loc_dist_cpd"] <- res_combined$loc[2]
    stats_matrix[1, "loc_dist_cpd_AD"] <- res_combined$loc[3]
    stats_matrix[1, "loc_dist_cpd_W"] <- res_combined$loc[4]
  })[["elapsed"]]

  methods <- variants
  result <- lapply(
    methods,
    .make_method_result,
    stats_matrix = stats_matrix,
    alpha = alpha,
    c = c,
    num_permut = num_permut,
    ndSup = ndSup,
    runtime_sec = runtime
  )
  names(result) <- methods
  attr(result, "runtime_sec") <- runtime
  attr(result, "n") <- n
  attr(result, "c") <- c
  attr(result, "num_permut") <- num_permut
  attr(result, "ndSup") <- ndSup
  attr(result, "variants") <- variants
  result
}

load_dist_profile_cpp <- load_distCPD_cpp
run_dist_profile <- run_distCPD
