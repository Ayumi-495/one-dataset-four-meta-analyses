script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

# blsmeta 0.1.0 constructs the level-three design matrix in first-appearance
# order, while its likelihood indexes study effects by sorted numeric ID.
# Sorting and renumbering are therefore part of the estimand-preserving input
# contract, not an optional cosmetic step.
prepare_blsmeta_subset <- function() {
  data <- prepare_study_level_data()$matched_data
  data <- data[order(data$study_ID, data$effect_size_ID), , drop = FALSE]
  data$es_blsmeta <- seq_len(nrow(data))
  data$study_blsmeta <- blsmeta::make_study_id(data$study_ID)
  first_rows <- !duplicated(data$study_blsmeta)
  stopifnot(
    nrow(data) == 232L,
    length(unique(data$study_blsmeta)) == 30L,
    identical(as.integer(data$study_blsmeta[first_rows]), seq_len(30L)),
    all(tapply(data$fertilizer, data$study_blsmeta, function(x) length(unique(x))) == 1L)
  )
  data
}

fit_blsmeta_subset <- function(data, response, variance, iter, warmup, seed) {
  fit_data <- data
  fit_data$y_blsmeta <- fit_data[[response]]
  fit_data$v_blsmeta <- fit_data[[variance]]
  set.seed(seed)
  started <- Sys.time()
  fit <- blsmeta::blsmeta(
    yi = y_blsmeta,
    vi = v_blsmeta,
    es_id = es_blsmeta,
    study_id = study_blsmeta,
    mods = ~ fertilizer,
    mods_scale2 = ~ 0 + fertilizer,
    mods_scale3 = ~ 0 + fertilizer,
    iter = iter,
    warmup = warmup,
    chains = 4,
    data = fit_data
  )
  attr(fit, "phase_a_metadata") <- list(
    response = response,
    variance = variance,
    iter = iter,
    warmup = warmup,
    seed = seed,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    package_version = as.character(utils::packageVersion("blsmeta"))
  )
  fit
}

post_warmup_samples <- function(fit) {
  # blsmeta 0.1.0 requests iter + warmup draws but does not discard warmup.
  # Apply the documented intent explicitly before diagnostics or summaries.
  first_iteration <- stats::start(fit$posterior_samples)[[1L]]
  warmup <- attr(fit, "phase_a_metadata")$warmup
  stats::window(fit$posterior_samples, start = first_iteration + warmup)
}

extract_blsmeta_subset <- function(fit) {
  metadata <- attr(fit, "phase_a_metadata")
  samples <- post_warmup_samples(fit)
  parameters <- c("beta[1]", "beta[2]", "gamma[1]", "gamma[2]", "eta[1]", "eta[2]")
  diagnostics <- coda::gelman.diag(samples[, parameters], multivariate = FALSE)$psrf[, 1L]
  effective_size <- coda::effectiveSize(samples[, parameters])
  draws <- do.call(rbind, lapply(samples, as.matrix))

  residual_sd <- exp(draws[, c("gamma[1]", "gamma[2]"), drop = FALSE])
  study_sd <- exp(draws[, c("eta[1]", "eta[2]"), drop = FALSE])
  residual_ratio <- residual_sd[, 2L] / residual_sd[, 1L]
  study_ratio <- study_sd[, 2L] / study_sd[, 1L]
  ratio_chains <- coda::mcmc.list(lapply(samples, function(chain) {
    chain <- as.matrix(chain)
    coda::mcmc(cbind(
      residual_sd_ratio = exp(chain[, "gamma[2]"] - chain[, "gamma[1]"]),
      study_sd_ratio = exp(chain[, "eta[2]"] - chain[, "eta[1]"])
    ))
  }))
  ratio_rhat <- coda::gelman.diag(ratio_chains, multivariate = FALSE)$psrf[, 1L]
  ratio_ess <- coda::effectiveSize(ratio_chains)
  diagnostic_pass <- max(diagnostics) <= 1.01 && min(effective_size) >= 400 &&
    max(ratio_rhat) <= 1.01 && min(ratio_ess) >= 400

  interval_row <- function(component, term, values) {
    interval <- stats::quantile(values, c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(
      response = metadata$response,
      model = "blsmeta_direct_study_scale",
      estimator = "Bayesian",
      component = component,
      term = term,
      estimate = interval[[2L]],
      lower = interval[[1L]],
      upper = interval[[3L]],
      fit_status = "success",
      diagnostic_status = if (diagnostic_pass) "passed" else "failed_convergence",
      max_rhat = max(diagnostics),
      min_ess = min(effective_size),
      ratio_rhat = if (component == "ratio" && grepl("residual", term)) ratio_rhat[[1L]] else if (component == "ratio") ratio_rhat[[2L]] else NA_real_,
      ratio_ess = if (component == "ratio" && grepl("residual", term)) ratio_ess[[1L]] else if (component == "ratio") ratio_ess[[2L]] else NA_real_,
      elapsed_seconds = metadata$elapsed_seconds,
      package_version = metadata$package_version,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, list(
    interval_row("residual_sd", "animal", residual_sd[, 1L]),
    interval_row("residual_sd", "plant", residual_sd[, 2L]),
    interval_row("study_sd", "animal", study_sd[, 1L]),
    interval_row("study_sd", "plant", study_sd[, 2L]),
    interval_row("ratio", "residual_plant_animal", residual_ratio),
    interval_row("ratio", "study_plant_animal", study_ratio)
  ))
}

write_blsmeta_outputs <- function(fits) {
  project_root <- find_project_root()
  cache_dir <- file.path(project_root, "Rdata", "study_level_heterogeneity", "direct")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- do.call(rbind, lapply(fits, extract_blsmeta_subset))
  utils::write.csv(
    summary,
    file.path(cache_dir, "blsmeta_subset_full_summary.csv"),
    row.names = FALSE
  )
  invisible(summary)
}

run_blsmeta_subset <- function(refit = FALSE) {
  if (!requireNamespace("blsmeta", quietly = TRUE)) stop("blsmeta is required.")
  if (!requireNamespace("coda", quietly = TRUE)) stop("coda is required.")
  project_root <- find_project_root()
  cache_dir <- file.path(project_root, "Rdata", "study_level_heterogeneity", "direct")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  settings <- list(
    lnRR = list(variance = "var.lnRR", iter = 20000L, warmup = 5000L, seed = 20260830L),
    lnCVR = list(variance = "var.lnCVR", iter = 5000L, warmup = 1000L, seed = 20260828L)
  )
  data <- prepare_blsmeta_subset()
  fits <- lapply(names(settings), function(response) {
    path <- file.path(cache_dir, paste0("blsmeta_", response, "_full.rds"))
    if (refit || !file.exists(path)) {
      specification <- settings[[response]]
      fit <- fit_blsmeta_subset(
        data, response, specification$variance,
        specification$iter, specification$warmup, specification$seed
      )
      saveRDS(fit, path)
    } else {
      fit <- readRDS(path)
    }
    fit
  })
  names(fits) <- names(settings)
  summary <- write_blsmeta_outputs(fits)
  stopifnot(
    all(c("lnRR", "lnCVR") %in% summary$response),
    all(summary$diagnostic_status == "failed_convergence")
  )
  print(summary[summary$component == "ratio", ])
  invisible(summary)
}

if (sys.nframe() == 0L) {
  run_blsmeta_subset(refit = "--refit" %in% commandArgs(trailingOnly = TRUE))
}
