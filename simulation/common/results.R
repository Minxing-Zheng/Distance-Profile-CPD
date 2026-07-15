.result_scalar <- function(x, name, default = NA) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0L) return(default)
  value[[1L]]
}

.config_scalar <- function(cfg, name, default = NA) {
  value <- cfg[[name]]
  if (is.null(value) || length(value) == 0L || is.na(value) || identical(value, "")) return(default)
  value[[1L]]
}

flatten_method_result <- function(result, method, cfg, replicate, data_seed,
                                  test_seed, true_cp, distance_runtime_sec,
                                  error = NULL) {
  scenario_id <- as.character(cfg$scenario_id)
  if (!is.null(error) || is.null(result)) {
    return(data.frame(
      scenario_id = scenario_id, replicate = replicate, data_seed = data_seed,
      test_seed = test_seed, generator = as.character(cfg$generator), method = method,
      n = as.integer(cfg$n), tau = as.numeric(cfg$tau), true_cp = true_cp,
      p = as.integer(cfg$p), signal = as.numeric(cfg$signal), null = as.logical(cfg$null),
      c = as.numeric(cfg$c), rho = as.numeric(.config_scalar(cfg, "rho", NA_real_)),
      covariance = as.character(.config_scalar(cfg, "covariance", NA_character_)),
      grid_size = as.integer(.config_scalar(cfg, "grid_size", NA_integer_)),
      object_size = as.integer(.config_scalar(cfg, "object_size", NA_integer_)),
      num_permut = as.integer(cfg$num_permut), alpha = as.numeric(cfg$alpha),
      test_stat = NA_real_, critical_value = NA_real_, p_value = NA_real_,
      reject = NA, location = NA_integer_, candidate_location = NA_integer_,
      absolute_error = NA_real_, candidate_absolute_error = NA_real_,
      distance_runtime_sec = distance_runtime_sec, method_runtime_sec = NA_real_,
      error = if (is.null(error)) "method returned NULL" else conditionMessage(error),
      stringsAsFactors = FALSE
    ))
  }

  reject <- .result_scalar(result, "reject", NA)
  location <- as.integer(.result_scalar(result, "loc", NA_integer_))
  candidate <- as.integer(.result_scalar(result, "candidate_loc", NA_integer_))
  has_truth <- !is.na(true_cp)
  data.frame(
    scenario_id = scenario_id, replicate = replicate, data_seed = data_seed,
    test_seed = test_seed, generator = as.character(cfg$generator), method = method,
    n = as.integer(cfg$n), tau = as.numeric(cfg$tau), true_cp = true_cp,
    p = as.integer(cfg$p), signal = as.numeric(cfg$signal), null = as.logical(cfg$null),
    c = as.numeric(cfg$c), rho = as.numeric(.config_scalar(cfg, "rho", NA_real_)),
    covariance = as.character(.config_scalar(cfg, "covariance", NA_character_)),
    grid_size = as.integer(.config_scalar(cfg, "grid_size", NA_integer_)),
    object_size = as.integer(.config_scalar(cfg, "object_size", NA_integer_)),
    num_permut = as.integer(cfg$num_permut), alpha = as.numeric(cfg$alpha),
    test_stat = as.numeric(.result_scalar(result, "test_stat", NA_real_)),
    critical_value = as.numeric(.result_scalar(result, "critical_value", NA_real_)),
    p_value = as.numeric(.result_scalar(result, "p_val", NA_real_)),
    reject = as.logical(reject), location = location, candidate_location = candidate,
    absolute_error = if (has_truth && isTRUE(reject) && !is.na(location)) abs(location - true_cp) else NA_real_,
    candidate_absolute_error = if (has_truth && !is.na(candidate)) abs(candidate - true_cp) else NA_real_,
    distance_runtime_sec = distance_runtime_sec,
    method_runtime_sec = as.numeric(.result_scalar(result, "runtime_sec", NA_real_)),
    error = "", stringsAsFactors = FALSE
  )
}

summarize_tidy_results <- function(detail) {
  if (nrow(detail) == 0L) return(data.frame())
  keys <- interaction(detail$scenario_id, detail$method, drop = TRUE, lex.order = TRUE)
  pieces <- split(detail, keys)
  rows <- lapply(pieces, function(x) {
    valid <- is.na(x$error) | x$error == ""
    detected <- valid & !is.na(x$reject) & x$reject & !is.na(x$absolute_error)
    candidate_valid <- valid & !is.na(x$candidate_absolute_error)
    data.frame(
      scenario_id = x$scenario_id[1L], method = x$method[1L],
      replicates = nrow(x), successful_replicates = sum(valid), failures = sum(!valid),
      rejection_rate = if (any(valid)) mean(x$reject[valid], na.rm = TRUE) else NA_real_,
      conditional_mae = if (any(detected)) mean(x$absolute_error[detected], na.rm = TRUE) else NA_real_,
      conditional_median_ae = if (any(detected)) stats::median(x$absolute_error[detected], na.rm = TRUE) else NA_real_,
      candidate_mae = if (any(candidate_valid)) mean(x$candidate_absolute_error[candidate_valid]) else NA_real_,
      mean_method_runtime_sec = if (any(valid)) mean(x$method_runtime_sec[valid], na.rm = TRUE) else NA_real_,
      mean_distance_runtime_sec = if (any(valid)) mean(x$distance_runtime_sec[valid], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$scenario_id, out$method), , drop = FALSE]
}
