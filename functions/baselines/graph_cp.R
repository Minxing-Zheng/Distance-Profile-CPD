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
    # gSeg can report p = 0 when the observed maximum exceeds every
    # permutation. Use the same finite-sample plus-one correction as the
    # other permutation wrappers so p-values and rejection flags agree.
    if (!is.na(test_stat)) {
      p_val <- (1 + sum(permuted_test_stat >= test_stat)) /
        (1 + length(permuted_test_stat))
    }
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
