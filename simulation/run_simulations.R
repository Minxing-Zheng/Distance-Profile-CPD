#!/usr/bin/env Rscript

# Configuration-driven, resumable simulation runner.
# One --array_index maps to exactly one (configuration row, replicate) pair.

parse_cli <- function(args) {
  out <- list(
    config = file.path("simulation", "configs", "smoke.csv"),
    output_dir = file.path("simulation", "results"),
    array_index = NA_integer_, config_row = NA_integer_, replicate = NA_integer_,
    overwrite = FALSE, save_raw = TRUE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument: ", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    value <- args[[i + 1L]]
    if (key %in% c("array_index", "config_row", "replicate")) {
      value <- as.integer(value)
    } else if (key %in% c("overwrite", "save_raw")) {
      value <- tolower(value) %in% c("true", "t", "1", "yes", "y")
    }
    out[[key]] <- value
    i <- i + 2L
  }
  out
}

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "functions", "run_benchmark_methods.R")) &&
        dir.exists(file.path(path, "simulation"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find the Distance-Profile-CPD project root.")
    path <- parent
  }
}

split_names <- function(value, default) {
  if (is.null(value) || is.na(value) || trimws(value) == "") return(default)
  unique(trimws(strsplit(as.character(value), "[;,]")[[1L]]))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

with_default_columns <- function(config) {
  defaults <- list(
    scenario_id = NA_character_, generator = "gaussian_mean", mc = 1L,
    n = 300L, tau = 1 / 3, p = 10L, signal = 0, null = FALSE,
    c = 0.1, num_permut = 99L, alpha = 0.05, min_size = 50L,
    methods = "dist_profile;energy;graph;kernel",
    variants = "dist_cpd;dist_cpd_AD;dist_cpd_W", ndsup = 200L,
    base_seed = 20260712L,
    permutation_scheme = "iid", block_length = ""
  )
  for (name in names(defaults)) {
    if (!name %in% names(config)) config[[name]] <- defaults[[name]]
    missing <- is.na(config[[name]]) | as.character(config[[name]]) == ""
    config[[name]][missing] <- defaults[[name]]
  }
  if (anyNA(config$scenario_id) || any(config$scenario_id == "")) {
    stop("Every configuration row must have a nonempty scenario_id.")
  }
  if (anyDuplicated(config$scenario_id)) stop("scenario_id values must be unique.")
  config
}

map_array_index <- function(config, index) {
  mc <- as.integer(config$mc)
  if (anyNA(mc) || any(mc < 1L)) stop("mc must be a positive integer in every row.")
  total <- sum(mc)
  if (is.na(index) || index < 1L || index > total) {
    stop("array_index must be between 1 and ", total, ".")
  }
  ends <- cumsum(mc)
  row <- which(index <= ends)[1L]
  start <- if (row == 1L) 1L else ends[row - 1L] + 1L
  list(config_row = row, replicate = index - start + 1L, total = total)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Could not move temporary output to ", path)
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".tmp_", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, temporary, compress = "xz")
  if (!file.rename(temporary, path)) stop("Could not move temporary output to ", path)
}

run_one_method <- function(method, generated, distmat, cfg, test_seed, project_root) {
  methods <- c("dist_profile", "energy", "graph", "kernel", "sn1", "sn2", "ss_sn")
  if (!method %in% methods) stop("Unknown method: ", method)
  set.seed(test_seed)
  elapsed <- system.time({
    # sn1/sn2/ss_sn ship with hardcoded asymptotic thresholds (see their
    # source) used only when num_permut=0; forcing that here regardless of
    # cfg$num_permut keeps them on their built-in calibration rather than
    # silently switching to a permutation loop shared with dist_profile.
    method_num_permut <- if (method %in% c("sn1", "sn2", "ss_sn")) {
      0L
    } else {
      as.integer(cfg$num_permut)
    }
    answer <- run_benchmark_methods(
      data = generated$data, distmat = distmat,
      c = as.numeric(cfg$c), num_permut = method_num_permut,
      alpha = as.numeric(cfg$alpha), min_size = as.integer(cfg$min_size),
      methods = method, dist_profile_variants = split_names(cfg$variants, "dist_cpd"),
      dist_profile_ndSup = as.integer(cfg$ndsup),
      dist_profile_permutation_scheme = as.character(cfg$permutation_scheme),
      dist_profile_block_length = if (is.na(cfg$block_length) || cfg$block_length == "") {
        NULL
      } else {
        as.integer(cfg$block_length)
      },
      seed = test_seed,
      project_root = project_root, progress = FALSE
    )
  })[["elapsed"]]

  if (method == "dist_profile") {
    for (name in names(answer$dist_profile)) {
      if (is.null(answer$dist_profile[[name]]$runtime_sec)) {
        answer$dist_profile[[name]]$runtime_sec <- elapsed
      }
    }
    return(answer$dist_profile)
  }
  result <- answer[[method]]
  if (!is.null(result) && is.null(result$runtime_sec)) result$runtime_sec <- elapsed
  setNames(list(result), method)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
project_root <- find_project_root()
setwd(project_root)
source(file.path("functions", "run_benchmark_methods.R"))
source(file.path("simulation", "common", "generators.R"))
source(file.path("simulation", "common", "results.R"))

config_path <- normalizePath(args$config, mustWork = TRUE)
config <- with_default_columns(utils::read.csv(config_path, stringsAsFactors = FALSE, check.names = FALSE))

if (!is.na(args$array_index)) {
  mapped <- map_array_index(config, args$array_index)
  config_row <- mapped$config_row
  replicate <- mapped$replicate
} else {
  config_row <- args$config_row
  replicate <- args$replicate
  if (is.na(config_row) || is.na(replicate)) {
    stop("Supply --array_index, or both --config_row and --replicate.")
  }
}
if (config_row < 1L || config_row > nrow(config)) stop("config_row is out of range.")
cfg <- as.list(config[config_row, , drop = FALSE])
if (replicate < 1L || replicate > as.integer(cfg$mc)) stop("replicate is out of range for this scenario.")

scenario_dir <- file.path(args$output_dir, as.character(cfg$scenario_id))
stem <- sprintf("replicate_%05d", replicate)
detail_path <- file.path(scenario_dir, paste0(stem, ".csv"))
raw_path <- file.path(scenario_dir, paste0(stem, ".rds"))
if (file.exists(detail_path) && !isTRUE(args$overwrite)) {
  message("Skipping completed result: ", detail_path)
  quit(save = "no", status = 0L)
}

base_seed <- as.integer(cfg$base_seed)
data_seed <- base_seed + replicate - 1L
message(
  "scenario=", cfg$scenario_id, " config_row=", config_row,
  " replicate=", replicate, "/", cfg$mc, " data_seed=", data_seed
)

generated <- generate_simulation_data(cfg, seed = data_seed)
distance_runtime <- system.time({
  if (is.null(generated$distmat)) {
    distance_method <- generated$distance_method %||% as.character(.cfg_value(cfg, "distance_method", "euclidean"))
    distance_scale <- generated$distance_scale %||% .cfg_numeric(cfg, "distance_scale", 1)
    distmat <- as.matrix(stats::dist(generated$data, method = distance_method)) * distance_scale
  } else {
    distmat <- generated$distmat
  }
})[["elapsed"]]

requested_methods <- split_names(cfg$methods, c("dist_profile", "energy", "graph", "kernel"))
seed_offsets <- c(
  dist_profile = 100000L, energy = 200000L, graph = 300000L,
  kernel = 400000L, sn1 = 500000L, sn2 = 600000L, ss_sn = 700000L
)
tidy <- list()
raw <- list()
for (method in requested_methods) {
  test_seed <- data_seed + seed_offsets[[method]]
  method_error <- NULL
  method_results <- tryCatch(
    run_one_method(method, generated, distmat, cfg, test_seed, project_root),
    error = function(e) {
      method_error <<- e
      NULL
    }
  )
  if (is.null(method_results)) {
    tidy[[length(tidy) + 1L]] <- flatten_method_result(
      NULL, method, cfg, replicate, data_seed, test_seed, generated$true_cp,
      distance_runtime, error = method_error
    )
  } else {
    for (result_name in names(method_results)) {
      result <- method_results[[result_name]]
      tidy[[length(tidy) + 1L]] <- flatten_method_result(
        result, result_name, cfg, replicate, data_seed, test_seed,
        generated$true_cp, distance_runtime,
        error = if (is.null(result)) simpleError("method returned NULL") else NULL
      )
      raw[[result_name]] <- result
    }
  }
}

tidy <- do.call(rbind, tidy)
atomic_write_csv(tidy, detail_path)
if (isTRUE(args$save_raw)) {
  atomic_save_rds(list(
    config = cfg, config_path = config_path, config_row = config_row,
    replicate = replicate, data_seed = data_seed, true_cp = generated$true_cp,
    generator_metadata = generated$metadata, distance_runtime_sec = distance_runtime,
    methods = raw, session_info = utils::sessionInfo()
  ), raw_path)
}
message("Wrote ", detail_path)
