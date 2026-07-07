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
  baseline_dir <- file.path(project_root, "functions", "baselines")
  source(file.path(baseline_dir, "ecp_distmat_input.R"))
  Rcpp::sourceCpp(file.path(baseline_dir, "energyChangePoint.cpp"))
  .run_all_runner_env$energy_loaded <- TRUE
  invisible(TRUE)
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
