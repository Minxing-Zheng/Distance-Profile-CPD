parse_args <- function(args) {
  out <- list(
    mc = 30L,
    num_permut = 300L,
    p = 10L,
    n1 = 100L,
    n2 = 200L,
    delta = 0.5,
    alpha = 0.05,
    dist_profile_variants = "all",
    dist_profile_ndSup = 1000L,
    cores = 1L,
    seed = 1L,
    out_dir = file.path("simulation", "mul_vec", "results_mean_p10")
  )

  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(out)) {
      stop("Unknown argument: ", args[i])
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", args[i])
    }
  old_value <- out[[key]]
  value <- args[i + 1L]
    out[[key]] <- if (identical(key, "dist_profile_variants")) {
      strsplit(value, ",", fixed = TRUE)[[1]]
    } else if (is.integer(old_value)) {
      as.integer(value)
    } else if (is.numeric(old_value)) {
      as.numeric(value)
    } else {
      value
    }
    i <- i + 2L
  }
  out
}

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, "functions")) &&
        dir.exists(file.path(path, "simulation"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find project root.")
    }
    path <- parent
  }
}

make_sigma <- function(p) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is required.")
  }
  m <- matrix(1, 1, p) / sqrt(p)
  n <- pracma::nullspace(m)
  u <- cbind(t(m), n)
  v <- cos(seq_len(p) * pi / p) + 1.5
  u %*% diag(v) %*% t(u)
}

method_rows <- function(res, run_id, true_tau) {
  rows <- list()

  for (name in names(res$dist_profile)) {
    x <- res$dist_profile[[name]]
    rows[[length(rows) + 1L]] <- data.frame(
      run = run_id,
      method = name,
      test_stat = x$test_stat,
      critical_value = x$critical_value,
      p_val = x$p_val,
      reject = isTRUE(x$reject),
      loc = ifelse(is.na(x$loc), NA_integer_, x$loc),
      candidate_loc = x$candidate_loc,
      abs_error = ifelse(isTRUE(x$reject), abs(x$loc - true_tau), NA_real_),
      candidate_abs_error = abs(x$candidate_loc - true_tau),
      runtime_sec = x$runtime_sec
    )
  }

  for (name in c("energy", "graph", "kernel")) {
    x <- res[[name]]
    rows[[length(rows) + 1L]] <- data.frame(
      run = run_id,
      method = name,
      test_stat = x$test_stat,
      critical_value = x$critical_value,
      p_val = x$p_val,
      reject = isTRUE(x$reject),
      loc = ifelse(is.na(x$loc), NA_integer_, x$loc),
      candidate_loc = x$candidate_loc,
      abs_error = ifelse(isTRUE(x$reject), abs(x$loc - true_tau), NA_real_),
      candidate_abs_error = abs(x$candidate_loc - true_tau),
      runtime_sec = x$runtime_sec
    )
  }

  do.call(rbind, rows)
}

summarize_results <- function(detail) {
  methods <- unique(detail$method)
  rows <- lapply(methods, function(method) {
    x <- detail[detail$method == method, , drop = FALSE]
    data.frame(
      method = method,
      monte_carlo = nrow(x),
      rejection_rate = mean(x$reject, na.rm = TRUE),
      conditional_mae = mean(x$abs_error, na.rm = TRUE),
      conditional_median_ae = stats::median(x$abs_error, na.rm = TRUE),
      candidate_mae = mean(x$candidate_abs_error, na.rm = TRUE),
      mean_runtime_sec = mean(x$runtime_sec, na.rm = TRUE)
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$method), , drop = FALSE]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- find_project_root()
setwd(project_root)
source(file.path("functions", "run_all.R"))

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(args$seed)
sigma <- make_sigma(args$p)
true_tau <- args$n1

run_one <- function(run_id) {
  run_seed <- args$seed + run_id - 1L
  set.seed(run_seed)
  data <- rbind(
    MASS::mvrnorm(args$n1, mu = rep(0, args$p), Sigma = sigma),
    MASS::mvrnorm(args$n2, mu = rep(args$delta, args$p), Sigma = sigma)
  )
  distmat <- as.matrix(stats::dist(data, method = "euclidean"))

  message(
    "run ", run_id, "/", args$mc,
    " seed=", run_seed,
    " n=", args$n1 + args$n2,
    " p=", args$p,
    " B=", args$num_permut
  )

  res <- run_all_methods(
    data = data,
    distmat = distmat,
    c = 0.1,
    num_permut = args$num_permut,
    alpha = args$alpha,
    min_size = 50,
    methods = c("dist_profile", "energy", "graph", "kernel"),
    dist_profile_variants = args$dist_profile_variants,
    dist_profile_ndSup = args$dist_profile_ndSup,
    seed = run_seed,
    progress = FALSE
  )
  method_rows(res, run_id = run_id, true_tau = true_tau)
}

if (args$cores > 1L) {
  if (.Platform$OS.type == "windows") {
    stop("--cores > 1 uses parallel::mclapply and is only supported on Unix/macOS.")
  }
  source_run_dist_profile(project_root)
  load_dist_profile_cpp(project_root = project_root)
  .load_energy_cp(project_root = project_root)
  if (!requireNamespace("ade4", quietly = TRUE)) {
    stop("Package 'ade4' is required for graph-CP.")
  }
  if (!requireNamespace("gSeg", quietly = TRUE)) {
    stop("Package 'gSeg' is required for graph-CP.")
  }
  message("running Monte Carlo in parallel with cores=", args$cores)
  detail <- parallel::mclapply(
    seq_len(args$mc),
    run_one,
    mc.cores = args$cores,
    mc.preschedule = FALSE
  )
  errors <- vapply(detail, inherits, logical(1), what = "error")
  if (any(errors)) {
    stop("Monte Carlo worker failed: ", conditionMessage(detail[[which(errors)[1]]]))
  }
} else {
  detail <- lapply(seq_len(args$mc), run_one)
}

detail <- do.call(rbind, detail)
summary <- summarize_results(detail)

stamp <- paste0(
  "mean_p", args$p,
  "_n", args$n1 + args$n2,
  "_delta", args$delta,
  "_mc", args$mc,
  "_B", args$num_permut
)
detail_path <- file.path(args$out_dir, paste0(stamp, "_detail.csv"))
summary_path <- file.path(args$out_dir, paste0(stamp, "_summary.csv"))
rds_path <- file.path(args$out_dir, paste0(stamp, ".rds"))

utils::write.csv(detail, detail_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)
saveRDS(list(args = args, detail = detail, summary = summary), rds_path)

print(summary)
message("wrote detail: ", detail_path)
message("wrote summary: ", summary_path)
message("wrote rds: ", rds_path)
