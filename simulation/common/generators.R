# Data generators used by simulation/run_simulations.R.
#
# Every generator returns a list with:
#   data: observations in rows (used by raw-data baselines)
#   distmat: optional precomputed distance matrix
#   true_cp: integer change-point location, or NA for a null scenario
#   metadata: generator-specific values worth saving

.cfg_value <- function(cfg, name, default = NULL) {
  value <- cfg[[name]]
  if (is.null(value) || length(value) == 0L || is.na(value) || identical(value, "")) {
    return(default)
  }
  value
}

.cfg_numeric <- function(cfg, name, default) {
  as.numeric(.cfg_value(cfg, name, default))
}

.cfg_integer <- function(cfg, name, default) {
  as.integer(.cfg_value(cfg, name, default))
}

.cfg_logical <- function(cfg, name, default = FALSE) {
  value <- .cfg_value(cfg, name, default)
  if (is.logical(value)) return(value)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

.spectral_covariance <- function(p) {
  first <- rep(1 / sqrt(p), p)
  basis <- qr.Q(qr(cbind(first, diag(p))))[, seq_len(p), drop = FALSE]
  eigenvalues <- cos(seq_len(p) * pi / p) + 1.5
  basis %*% diag(eigenvalues, p) %*% t(basis)
}

.rmvnorm <- function(n, mean, sigma) {
  p <- length(mean)
  if (n == 0L) return(matrix(numeric(0), nrow = 0L, ncol = p))
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  sweep(z %*% chol(sigma), 2L, mean, "+")
}

.segment_sizes <- function(cfg) {
  n <- .cfg_integer(cfg, "n", 300L)
  tau <- .cfg_numeric(cfg, "tau", 1 / 3)
  is_null <- .cfg_logical(cfg, "null", FALSE)
  if (n < 4L) stop("n must be at least 4.")
  if (!is.finite(tau) || tau <= 0 || tau >= 1) stop("tau must be in (0, 1).")
  n1 <- as.integer(round(n * tau))
  n1 <- max(2L, min(n - 2L, n1))
  list(n = n, n1 = n1, n2 = n - n1, true_cp = if (is_null) NA_integer_ else n1)
}

.base_covariance <- function(cfg, p) {
  covariance <- tolower(as.character(.cfg_value(cfg, "covariance", "spectral")))
  if (covariance == "identity") return(diag(p))
  if (covariance == "spectral") return(.spectral_covariance(p))
  stop("Unknown covariance: ", covariance, ". Use identity or spectral.")
}

.apply_perturbation <- function(data, cfg) {
  transform_fraction <- .cfg_numeric(cfg, "transform_fraction", 0)
  if (transform_fraction > 0) {
    n_transform <- min(ncol(data), ceiling(transform_fraction * ncol(data)))
    data[, seq_len(n_transform)] <- exp(data[, seq_len(n_transform), drop = FALSE])
  }

  observation_fraction <- .cfg_numeric(cfg, "outlier_observation_fraction", 0)
  feature_fraction <- .cfg_numeric(cfg, "outlier_feature_fraction", 0)
  if (observation_fraction > 0 && feature_fraction > 0) {
    scale <- .cfg_numeric(cfg, "outlier_scale", 10)
    df <- .cfg_numeric(cfg, "outlier_df", 10)
    rows <- stats::runif(nrow(data)) <= observation_fraction
    mask <- matrix(
      stats::rbinom(length(data), size = 1L, prob = feature_fraction),
      nrow = nrow(data), ncol = ncol(data)
    )
    noise <- scale * matrix(stats::rt(length(data), df = df), nrow = nrow(data))
    noise[!rows, ] <- 0
    data <- data + mask * noise
  }
  data
}

.generate_gaussian <- function(cfg, kind) {
  sizes <- .segment_sizes(cfg)
  p <- .cfg_integer(cfg, "p", 10L)
  signal <- .cfg_numeric(cfg, "signal", 0)
  is_null <- .cfg_logical(cfg, "null", FALSE)
  effective_signal <- if (is_null) 0 else signal
  identity <- diag(p)

  if (kind == "mean") {
    sigma <- .base_covariance(cfg, p)
    left <- .rmvnorm(sizes$n1, rep(0, p), sigma)
    right <- .rmvnorm(sizes$n2, rep(effective_signal, p), sigma)
  } else if (kind == "scale") {
    base_variance <- .cfg_numeric(cfg, "base_variance", 0.8)
    post_variance <- base_variance - effective_signal
    if (post_variance <= 0) stop("scale signal must be smaller than base_variance.")
    left <- .rmvnorm(sizes$n1, rep(0, p), base_variance * identity)
    right <- .rmvnorm(sizes$n2, rep(0, p), post_variance * identity)
  } else if (kind == "mixture") {
    active <- max(1L, ceiling(.cfg_numeric(cfg, "active_fraction", 0.1) * p))
    mu <- c(rep(effective_signal, active), rep(0, p - active))
    left <- .rmvnorm(sizes$n1, rep(0, p), identity)
    component <- stats::rbinom(sizes$n2, 1L, 0.5)
    right <- .rmvnorm(sizes$n2, rep(0, p), identity)
    right <- right + ifelse(component == 1L, 1, -1) * matrix(mu, sizes$n2, p, byrow = TRUE)
  } else if (kind == "tail") {
    left <- matrix(stats::rnorm(sizes$n1 * p), sizes$n1, p)
    if (is_null || signal <= 0) {
      right <- matrix(stats::rnorm(sizes$n2 * p), sizes$n2, p)
    } else {
      right <- matrix(stats::rt(sizes$n2 * p, df = signal), sizes$n2, p)
    }
  } else if (kind == "covariance") {
    spike_count <- min(p, .cfg_integer(cfg, "spike_count", max(1L, ceiling(p / 4))))
    sigma_left <- identity
    eigenvalues <- c(rep(1 + effective_signal, spike_count), rep(1, p - spike_count))
    sigma_right <- diag(eigenvalues, p)
    left <- .rmvnorm(sizes$n1, rep(0, p), sigma_left)
    right <- .rmvnorm(sizes$n2, rep(0, p), sigma_right)
  } else {
    stop("Unsupported Gaussian generator kind: ", kind)
  }

  data <- .apply_perturbation(rbind(left, right), cfg)
  list(
    data = data,
    distmat = NULL,
    distance_method = "euclidean",
    distance_scale = 1,
    true_cp = sizes$true_cp,
    metadata = list(kind = kind, effective_signal = effective_signal)
  )
}

.generate_ar_mean <- function(cfg) {
  sizes <- .segment_sizes(cfg)
  p <- .cfg_integer(cfg, "p", 10L)
  rho <- .cfg_numeric(cfg, "rho", 0.3)
  signal <- if (.cfg_logical(cfg, "null", FALSE)) 0 else .cfg_numeric(cfg, "signal", 0)
  if (abs(rho) >= 1) stop("rho must have absolute value below 1.")
  innovations <- matrix(stats::rnorm(sizes$n * p, sd = sqrt(1 - rho^2)), sizes$n, p)
  data <- innovations
  data[1L, ] <- stats::rnorm(p)
  for (i in 2:sizes$n) data[i, ] <- rho * data[i - 1L, ] + innovations[i, ]
  data[(sizes$n1 + 1L):sizes$n, ] <- data[(sizes$n1 + 1L):sizes$n, ] + signal
  list(
    data = .apply_perturbation(data, cfg), distmat = NULL,
    distance_method = "euclidean", distance_scale = 1,
    true_cp = sizes$true_cp, metadata = list(kind = "ar_mean", rho = rho)
  )
}

.build_dependent_chain <- function(n, p, dependence, params, innovation_fn) {
  eps <- innovation_fn(n, p)
  x <- matrix(0, n, p)
  x[1L, ] <- eps[1L, ]
  if (n < 2L) return(x)
  if (dependence == "ar1") {
    rho <- params$rho
    if (abs(rho) >= 1) stop("rho must have absolute value below 1.")
    for (i in 2:n) x[i, ] <- rho * x[i - 1L, ] + eps[i, ]
  } else if (dependence == "arma11") {
    phi <- params$phi
    theta <- params$theta
    if (abs(phi) >= 1) stop("phi must have absolute value below 1.")
    for (i in 2:n) x[i, ] <- phi * x[i - 1L, ] + eps[i, ] + theta * eps[i - 1L, ]
  } else {
    stop("Unknown dependence: ", dependence, ". Use ar1 or arma11.")
  }
  x
}

# Serially dependent change-point generator. `dependence` selects the time-
# series structure (ar1: rho; arma11: phi, theta) and `change_kind` selects
# where the break enters (mean: additive shift on one continuous chain;
# scale/tail: independent left/right chains with different innovations, kept
# separate so each segment is a clean stationary realization at its own
# innovation distribution, mirroring the non-dependent scale/tail generators).
.generate_dependent_change <- function(cfg) {
  sizes <- .segment_sizes(cfg)
  p <- .cfg_integer(cfg, "p", 10L)
  change_kind <- tolower(as.character(.cfg_value(cfg, "change_kind", "mean")))
  dependence <- tolower(as.character(.cfg_value(cfg, "dependence", "ar1")))
  is_null <- .cfg_logical(cfg, "null", FALSE)
  signal <- if (is_null) 0 else .cfg_numeric(cfg, "signal", 0)

  params <- if (dependence == "ar1") {
    list(rho = .cfg_numeric(cfg, "rho", 0))
  } else if (dependence == "arma11") {
    list(phi = .cfg_numeric(cfg, "phi", 0), theta = .cfg_numeric(cfg, "theta", 0))
  } else {
    stop("Unknown dependence: ", dependence, ". Use ar1 or arma11.")
  }

  gaussian_innovations <- function(n, p) matrix(stats::rnorm(n * p), n, p)

  if (change_kind == "mean") {
    data <- .build_dependent_chain(sizes$n, p, dependence, params, gaussian_innovations)
    data[(sizes$n1 + 1L):sizes$n, ] <- data[(sizes$n1 + 1L):sizes$n, ] + signal
  } else if (change_kind == "scale") {
    base_variance <- .cfg_numeric(cfg, "base_variance", 0.8)
    post_variance <- base_variance - signal
    if (post_variance <= 0) stop("scale signal must be smaller than base_variance.")
    left <- .build_dependent_chain(
      sizes$n1, p, dependence, params,
      function(n, p) matrix(stats::rnorm(n * p, sd = sqrt(base_variance)), n, p)
    )
    right <- .build_dependent_chain(
      sizes$n2, p, dependence, params,
      function(n, p) matrix(stats::rnorm(n * p, sd = sqrt(post_variance)), n, p)
    )
    data <- rbind(left, right)
  } else if (change_kind == "tail") {
    left <- .build_dependent_chain(sizes$n1, p, dependence, params, gaussian_innovations)
    if (is_null || signal <= 0) {
      right <- .build_dependent_chain(sizes$n2, p, dependence, params, gaussian_innovations)
    } else {
      right <- .build_dependent_chain(
        sizes$n2, p, dependence, params,
        function(n, p) matrix(stats::rt(n * p, df = signal), n, p)
      )
    }
    data <- rbind(left, right)
  } else {
    stop("Unsupported change_kind: ", change_kind, ". Use mean, scale, or tail.")
  }

  list(
    data = .apply_perturbation(data, cfg), distmat = NULL,
    distance_method = "euclidean", distance_scale = 1,
    true_cp = sizes$true_cp,
    metadata = list(kind = paste0("dependent_", change_kind), dependence = dependence, params = params)
  )
}

.generate_distributional <- function(cfg, kind) {
  sizes <- .segment_sizes(cfg)
  signal <- if (.cfg_logical(cfg, "null", FALSE)) 0 else .cfg_numeric(cfg, "signal", 0)
  grid_size <- .cfg_integer(cfg, "grid_size", 40L)
  grid_limit <- .cfg_numeric(cfg, "grid_limit", 3)
  grid <- seq(-grid_limit, grid_limit, length.out = grid_size)
  points <- expand.grid(x = grid, y = grid)

  if (kind == "location") {
    z_left <- cbind(stats::rnorm(sizes$n1, 0, 0.5), stats::rnorm(sizes$n1, 0, 0.5))
    z_right <- cbind(stats::rnorm(sizes$n2, signal, 0.5), stats::rnorm(sizes$n2, 0, 0.5))
  } else if (kind == "scale") {
    post_sd <- 0.4 + signal
    z_left <- cbind(stats::rnorm(sizes$n1, 0, 0.4), stats::rnorm(sizes$n1, 0, 0.4))
    z_right <- cbind(stats::rnorm(sizes$n2, 0, post_sd), stats::rnorm(sizes$n2, 0, 0.4))
  } else {
    stop("Unsupported distributional kind: ", kind)
  }
  centers <- rbind(z_left, z_right)

  # With identity covariance, the bivariate Gaussian CDF factors into two
  # univariate Gaussian CDFs. This is equivalent to the old pmvnorm loop and
  # is much faster.
  data <- t(vapply(seq_len(nrow(centers)), function(i) {
    stats::pnorm(points$x - centers[i, 1L]) * stats::pnorm(points$y - centers[i, 2L])
  }, numeric(nrow(points))))
  step <- if (grid_size > 1L) grid[2L] - grid[1L] else 1
  list(
    data = data, distmat = NULL, distance_method = "manhattan", distance_scale = step,
    true_cp = sizes$true_cp,
    metadata = list(kind = paste0("distributional_", kind), grid_size = grid_size)
  )
}

.generate_network <- function(cfg) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for the network generator.")
  }
  sizes <- .segment_sizes(cfg)
  object_size <- .cfg_integer(cfg, "object_size", 50L)
  signal <- if (.cfg_logical(cfg, "null", FALSE)) 0 else .cfg_numeric(cfg, "signal", 0)
  base_power <- .cfg_numeric(cfg, "base_power", 0)

  one_graph <- function(power) {
    graph <- igraph::sample_pa(object_size, power = power, directed = FALSE)
    adjacency <- as.matrix(igraph::as_adjacency_matrix(graph, sparse = FALSE))
    laplacian <- diag(rowSums(adjacency)) - adjacency
    c(2 * laplacian[upper.tri(laplacian)], diag(laplacian))
  }
  data <- t(vapply(
    c(rep(base_power, sizes$n1), rep(base_power + signal, sizes$n2)),
    one_graph, numeric(object_size * (object_size + 1L) / 2L)
  ))
  list(
    data = data, distmat = NULL, distance_method = "manhattan", distance_scale = 1,
    true_cp = sizes$true_cp,
    metadata = list(kind = "network_pa", object_size = object_size)
  )
}

generate_simulation_data <- function(cfg, seed) {
  set.seed(seed)
  generator <- tolower(as.character(.cfg_value(cfg, "generator", "gaussian_mean")))
  result <- switch(
    generator,
    gaussian_mean = .generate_gaussian(cfg, "mean"),
    gaussian_scale = .generate_gaussian(cfg, "scale"),
    gaussian_mixture = .generate_gaussian(cfg, "mixture"),
    gaussian_tail = .generate_gaussian(cfg, "tail"),
    gaussian_covariance = .generate_gaussian(cfg, "covariance"),
    ar_mean = .generate_ar_mean(cfg),
    dependent_change = .generate_dependent_change(cfg),
    distributional_location = .generate_distributional(cfg, "location"),
    distributional_scale = .generate_distributional(cfg, "scale"),
    network_pa = .generate_network(cfg),
    stop("Unknown generator: ", generator)
  )
  result$generator <- generator
  result
}
