# Shared wrapper to run dist-CP and baseline methods from common inputs.
#
# Main entry point:
#   run_all_methods(data = X, distmat = NULL, num_permut = 500)
#
# At least one of data or distmat must be supplied. If distmat is missing, it is
# computed from data using stats::dist(data, method = distance_method).

.run_all_runner_env <- new.env(parent = emptyenv())
.run_all_runner_env$energy_loaded <- FALSE

.run_all_find_project_root <- function(start = getwd()) {
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

.run_all_validate_distmat <- function(distmat) {
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

source_run_dist_profile <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- .run_all_find_project_root()
  }
  source(file.path(project_root, "functions", "run_dist_profile.R"))
  invisible(TRUE)
}

.ensure_run_dist_profile_loaded <- function(project_root = NULL) {
  if (!exists("run_dist_profile", mode = "function")) {
    source_run_dist_profile(project_root)
  }
}

.load_energy_cp <- function(project_root = NULL) {
  if (isTRUE(.run_all_runner_env$energy_loaded)) {
    return(invisible(TRUE))
  }
  if (is.null(project_root)) {
    project_root <- .run_all_find_project_root()
  }
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Package 'Rcpp' is required for energy-CP.")
  }
  function_dir <- file.path(project_root, "functions")
  source(file.path(function_dir, "ecp_distmat_input.R"))
  Rcpp::sourceCpp(file.path(function_dir, "energyChangePoint.cpp"))
  .run_all_runner_env$energy_loaded <- TRUE
  invisible(TRUE)
}

.compute_distmat <- function(data, distance_method = "euclidean") {
  if (is.null(data)) {
    stop("data is required when distmat is not supplied.")
  }
  as.matrix(stats::dist(data, method = distance_method))
}

.mmd_scan_stat <- function(distmat, c = 0.1) {
  distmat <- .run_all_validate_distmat(distmat)
  n <- nrow(distmat)
  cp_min <- ceiling(n * c)
  cp_max <- n - ceiling(n * c)
  if (cp_min < 2) cp_min <- 2
  if (cp_max > n - 1) cp_max <- n - 1

  off_diag <- distmat[upper.tri(distmat)]
  bw <- stats::median(off_diag[off_diag > 0], na.rm = TRUE)
  if (!is.finite(bw) || bw <= 0) {
    bw <- 1
  }
  gram_matrix <- exp(-(distmat^2) / (2 * bw^2))

  obj_vals <- rep(NA_real_, n)
  for (t in cp_min:cp_max) {
    left <- seq_len(t)
    right <- (t + 1):n
    term1 <- length(right) / (length(left) * n) * sum(gram_matrix[left, left, drop = FALSE])
    term2 <- -2 / n * sum(gram_matrix[left, right, drop = FALSE])
    term3 <- length(left) / (length(right) * n) * sum(gram_matrix[right, right, drop = FALSE])
    obj_vals[t] <- term1 + term2 + term3
  }

  loc <- which.max(obj_vals)
  c(stat = unname(obj_vals[loc]), loc = loc)
}

run_kernel_mmd <- function(distmat,
                           c = 0.1,
                           num_permut = 500,
                           alpha = 0.05,
                           seed = NULL,
                           progress = FALSE) {
  distmat <- .run_all_validate_distmat(distmat)
  n <- nrow(distmat)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  runtime <- system.time({
    stats <- matrix(NA_real_, nrow = num_permut + 1, ncol = 2)
    colnames(stats) <- c("stat", "loc")
    for (b in 0:num_permut) {
      if (isTRUE(progress)) {
        message("kernel permutation ", b, "/", num_permut)
      }
      current_distmat <- if (b == 0) {
        distmat
      } else {
        ind <- sample.int(n)
        distmat[ind, ind, drop = FALSE]
      }
      stats[b + 1, ] <- .mmd_scan_stat(current_distmat, c = c)
    }
  })[["elapsed"]]

  test_stat <- unname(stats[1, "stat"])
  permuted_test_stat <- if (num_permut > 0) unname(stats[-1, "stat"]) else numeric(0)
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
  reject <- if (is.na(critical_value)) NA else test_stat > critical_value
  candidate_loc <- unname(stats[1, "loc"])
  loc <- if (isTRUE(reject)) candidate_loc else NA_integer_

  list(
    method = "kernel_mmd",
    test_stat = test_stat,
    permuted_test_stat = permuted_test_stat,
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = loc,
    candidate_loc = candidate_loc,
    observed_stat = test_stat,
    permuted_stats = permuted_test_stat,
    all_stats = unname(stats[, "stat"]),
    detected = reject,
    c = c,
    num_permut = num_permut,
    alpha = alpha,
    runtime_sec = runtime
  )
}

run_energy_cp <- function(distmat,
                          num_permut = 500,
                          alpha = 0.05,
                          min_size = 50,
                          project_root = NULL) {
  distmat <- .run_all_validate_distmat(distmat)
  .load_energy_cp(project_root)

  runtime <- system.time({
    fit <- e.divisive_distmat(
      D = distmat,
      sig.lvl = alpha,
      R = num_permut,
      k = NULL,
      min.size = min_size,
      alpha = 1
    )
  })[["elapsed"]]

  test_stat <- if (length(fit$test.statistics) > 0) fit$test.statistics[1] else NA_real_
  permuted_test_stat <- if (length(fit$permuted.statistics) > 0) {
    fit$permuted.statistics[[1]]
  } else {
    numeric(0)
  }
  critical_value <- if (length(permuted_test_stat) > 0) {
    unname(stats::quantile(permuted_test_stat, probs = 1 - alpha, type = 1, names = FALSE))
  } else {
    NA_real_
  }
  candidate_loc <- if (length(fit$considered.last) > 0) fit$considered.last[1] else NA_integer_
  reject <- length(fit$estimates) > 2
  loc <- if (isTRUE(reject)) fit$order.found[3] else NA_integer_
  p_val <- if (length(fit$p.values) > 0) fit$p.values[1] else NA_real_

  list(
    method = "energy_cp",
    test_stat = test_stat,
    permuted_test_stat = unname(permuted_test_stat),
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = loc,
    candidate_loc = candidate_loc,
    significant_locs = if (length(fit$estimates) > 2) fit$estimates[-c(1, length(fit$estimates))] else integer(0),
    fit = fit,
    detected = reject,
    num_permut = num_permut,
    alpha = alpha,
    min_size = min_size,
    runtime_sec = runtime
  )
}

run_graph_cp <- function(data = NULL,
                         distmat = NULL,
                         num_permut = 500,
                         alpha = 0.05,
                         graph_ngmax = 5,
                         statistics = "all") {
  if (!requireNamespace("ade4", quietly = TRUE)) {
    warning("Package 'ade4' is not installed; skipping graph-CP.")
    return(NULL)
  }
  if (!requireNamespace("gSeg", quietly = TRUE)) {
    warning("Package 'gSeg' is not installed; skipping graph-CP.")
    return(NULL)
  }
  if (is.null(distmat)) {
    if (is.null(data)) {
      stop("Either data or distmat is required for graph-CP.")
    }
    dist_obj <- stats::dist(data)
  } else {
    dist_obj <- stats::as.dist(.run_all_validate_distmat(distmat))
  }

  start_time <- proc.time()[["elapsed"]]
  fit <- tryCatch({
    edge <- ade4::mstree(dist_obj, ngmax = graph_ngmax)
    graph_fit <- NULL
    utils::capture.output(
      graph_fit <- gSeg::gseg1(attr(dist_obj, "Size"), edge, statistics = statistics, B = num_permut, pval.perm = TRUE)
    )
    graph_fit
  }, error = function(e) e)
  runtime <- proc.time()[["elapsed"]] - start_time

  if (inherits(fit, "error")) {
    warning("graph-CP failed: ", conditionMessage(fit))
    return(list(
      method = "graph_cp",
      test_stat = NA_real_,
      permuted_test_stat = numeric(0),
      critical_value = NA_real_,
      p_val = NA_real_,
      reject = NA,
      loc = NA_integer_,
      candidate_loc = NA_integer_,
      error = conditionMessage(fit),
      detected = NA,
      num_permut = num_permut,
      alpha = alpha,
      runtime_sec = runtime
    ))
  }

  graph_result <- fit$scanZ$generalized
  graph_perm <- fit$pval.perm$generalized
  p_val <- NA_real_
  test_stat <- NA_real_
  permuted_test_stat <- numeric(0)
  critical_value <- NA_real_
  candidate_loc <- NA_integer_
  reject <- NA

  if (!is.null(fit$pval.perm$generalized$pval)) {
    p_val <- graph_perm$pval
  } else if (!is.null(fit$pval.appr$generalized)) {
    p_val <- fit$pval.appr$generalized
  }
  if (!is.null(graph_result$Zmax)) {
    test_stat <- graph_result$Zmax
  }
  if (!is.null(graph_result$tauhat)) {
    candidate_loc <- graph_result$tauhat
  }
  if (!is.null(graph_perm$maxZs)) {
    permuted_test_stat <- graph_perm$maxZs
    critical_value <- unname(stats::quantile(permuted_test_stat, probs = 1 - alpha, type = 1, names = FALSE))
  }
  reject <- if (!is.na(p_val)) {
    p_val <= alpha
  } else if (!is.na(critical_value)) {
    test_stat > critical_value
  } else {
    NA
  }
  loc <- if (isTRUE(reject)) candidate_loc else NA_integer_

  list(
    method = "graph_cp",
    test_stat = test_stat,
    permuted_test_stat = unname(permuted_test_stat),
    critical_value = critical_value,
    p_val = p_val,
    reject = reject,
    loc = loc,
    candidate_loc = candidate_loc,
    fit = fit,
    detected = reject,
    num_permut = num_permut,
    alpha = alpha,
    runtime_sec = runtime
  )
}

run_all_methods <- function(data = NULL,
                            distmat = NULL,
                            distance_method = "euclidean",
                            c = 0.1,
                            num_permut = 500,
                            alpha = 0.05,
                            min_size = 50,
                            methods = c("dist_profile", "energy", "graph", "kernel"),
                            dist_profile_variants = "all",
                            dist_profile_ndSup = 1000,
                            seed = NULL,
                            project_root = NULL,
                            progress = TRUE) {
  if (is.null(distmat)) {
    distmat <- .compute_distmat(data, distance_method = distance_method)
  } else {
    distmat <- .run_all_validate_distmat(distmat)
  }

  out <- list()

  if ("dist_profile" %in% methods) {
    .ensure_run_dist_profile_loaded(project_root)
    out$dist_profile <- run_dist_profile(
      distmat = distmat,
      c = c,
      num_permut = num_permut,
      alpha = alpha,
      seed = seed,
      project_root = project_root,
      progress = progress,
      variants = dist_profile_variants,
      ndSup = dist_profile_ndSup
    )
  }

  if ("energy" %in% methods) {
    out$energy <- run_energy_cp(
      distmat = distmat,
      num_permut = num_permut,
      alpha = alpha,
      min_size = min_size,
      project_root = project_root
    )
  }

  if ("graph" %in% methods) {
    out$graph <- run_graph_cp(
      data = data,
      distmat = distmat,
      num_permut = num_permut,
      alpha = alpha
    )
  }

  if ("kernel" %in% methods) {
    out$kernel <- run_kernel_mmd(
      distmat = distmat,
      c = c,
      num_permut = num_permut,
      alpha = alpha,
      seed = if (is.null(seed)) NULL else seed + 100000L,
      progress = progress
    )
  }

  attr(out, "n") <- nrow(distmat)
  attr(out, "c") <- c
  attr(out, "num_permut") <- num_permut
  attr(out, "alpha") <- alpha
  attr(out, "dist_profile_variants") <- dist_profile_variants
  attr(out, "dist_profile_ndSup") <- dist_profile_ndSup
  out
}
