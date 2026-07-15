#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(
    input_dir = file.path("simulation", "results"),
    output = file.path("simulation", "results", "summary.csv")
  )
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(out)) stop("Unknown argument: ", args[[i]])
    if (i == length(args)) stop("Missing value for argument: ", args[[i]])
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

find_project_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "simulation", "common", "results.R"))) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find project root.")
    path <- parent
  }
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
root <- find_project_root()
setwd(root)
source(file.path("simulation", "common", "results.R"))

files <- list.files(
  args$input_dir, pattern = "^replicate_[0-9]+[.]csv$",
  recursive = TRUE, full.names = TRUE
)
if (length(files) == 0L) stop("No replicate CSV files found under ", args$input_dir)
detail <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
summary <- summarize_tidy_results(detail)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(summary, args$output, row.names = FALSE, na = "")
print(summary)
message("Wrote ", args$output)
