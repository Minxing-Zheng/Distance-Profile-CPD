# Shared wrapper to run benchmark change-point methods from common inputs.
#
# Main entry point:
#   run_benchmark_methods(data = X, distmat = NULL, num_permut = 500)
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

source_run_distCPD <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- .run_all_find_project_root()
  }
  source(file.path(project_root, "functions", "run_distCPD.R"))
  invisible(TRUE)
}

.ensure_run_distCPD_loaded <- function(project_root = NULL) {
  if (!exists("run_distCPD", mode = "function")) {
    source_run_distCPD(project_root)
  }
}

.load_baseline_helpers <- function(project_root = NULL) {
  if (is.null(project_root)) {
    project_root <- .run_all_find_project_root()
  }
  baseline_dir <- file.path(project_root, "functions", "baselines")
  source(file.path(baseline_dir, "energy_cp.R"))
  source(file.path(baseline_dir, "kernel_mmd.R"))
  source(file.path(baseline_dir, "graph_cp.R"))
  source(file.path(baseline_dir, "wbs_common.R"))
  source(file.path(baseline_dir, "sn1_cp.R"))
  source(file.path(baseline_dir, "sn2_cp.R"))
  source(file.path(baseline_dir, "ss_sn_cp.R"))
  invisible(TRUE)
}

.compute_distmat <- function(data, distance_method = "euclidean") {
  if (is.null(data)) {
    stop("data is required when distmat is not supplied.")
  }
  as.matrix(stats::dist(data, method = distance_method))
}

run_benchmark_methods <- function(data = NULL,
                                  distmat = NULL,
                                  distance_method = "euclidean",
                                  c = 0.1,
                                  num_permut = 500,
                                  alpha = 0.05,
                                  min_size = 50,
                            methods = c("dist_profile", "energy", "graph", "kernel"),
                                  dist_profile_variants = "all",
                                  dist_profile_ndSup = 200,
                                  dist_profile_permutation_scheme = "iid",
                                  dist_profile_block_length = NULL,
                                  seed = NULL,
                                  project_root = NULL,
                                  progress = TRUE) {
  if (is.null(distmat)) {
    distmat <- .compute_distmat(data, distance_method = distance_method)
  } else {
    distmat <- .run_all_validate_distmat(distmat)
  }
  if (any(c("sn1", "sn2", "ss_sn") %in% methods) && is.null(data)) {
    stop("data is required for the SN1/SN2/SS-SN baseline methods.")
  }
  .load_baseline_helpers(project_root)

  out <- list()

  if ("dist_profile" %in% methods) {
    .ensure_run_distCPD_loaded(project_root)
    out$dist_profile <- run_distCPD(
      distmat = distmat,
      c = c,
      num_permut = num_permut,
      alpha = alpha,
      seed = seed,
      project_root = project_root,
      progress = progress,
      variants = dist_profile_variants,
      ndSup = dist_profile_ndSup,
      permutation_scheme = dist_profile_permutation_scheme,
      block_length = dist_profile_block_length
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

  if ("sn1" %in% methods) {
    out$sn1 <- run_sn1_cp(
      data = data,
      alpha = alpha,
      num_permut = num_permut,
      seed = if (is.null(seed)) 2L else seed + 400000L
    )
  }

  if ("sn2" %in% methods) {
    out$sn2 <- run_sn2_cp(
      data = data,
      alpha = alpha,
      num_permut = num_permut,
      seed = if (is.null(seed)) 1L else seed + 300000L
    )
  }

  if ("ss_sn" %in% methods) {
    out$ss_sn <- run_ss_sn_cp(
      data = data,
      alpha = alpha,
      num_permut = num_permut,
      seed = if (is.null(seed)) 100L else seed + 200000L
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

source_run_dist_profile <- source_run_distCPD
.ensure_run_dist_profile_loaded <- .ensure_run_distCPD_loaded
run_all_methods <- run_benchmark_methods
run_distCPD_ALL <- run_benchmark_methods

.load_baseline_helpers()
