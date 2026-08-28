# Bayesian location--scale meta-analysis for the study-level heterogeneity arc.
#
# Symbolic alignment (all standard deviations are on the effect-size scale):
#   y_i | beta_c, u_{s,c}, sigma_e,c, v_i ~ Normal(
#     beta_c + u_{s,c}, sqrt(v_i + sigma_e,c^2))
#   u_s ~ Normal(0, sigma_u)                          [baseline]
#   u_{s,c} ~ Normal(0, sigma_u,c), independently by c [extended]
#
# brms implements the known diagonal v_i exactly through
# `y | se(sqrt(vi), sigma = TRUE)`. Thus `sigma` is only additional
# effect-size residual heterogeneity; it is not a substitute for sampling SD.

script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

bayesian_constants <- list(
  smoke_chains = 2L, smoke_iter = 1000L, smoke_warmup = 500L,
  full_chains = 4L, full_iter = 2000L, full_warmup = 1000L,
  smoke_rhat_max = 1.10, full_rhat_max = 1.01, full_ess_min = 400,
  fit_gate_minutes = 90, seed = 20260827L
)

require_bayesian_packages <- function() {
  packages <- c("brms", "cmdstanr", "posterior")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) stop("Missing required package(s): ", paste(missing, collapse = ", "), ".")
  if (!nzchar(cmdstanr::cmdstan_path())) stop("cmdstanr has no configured CmdStan installation.")
}

prepare_bayesian_data <- function(data, response, variance) {
  required <- c(response, variance, "fertilizer", "study_ID")
  if (!all(required %in% names(data))) stop("Bayesian data are missing one or more required columns.")
  result <- data.frame(
    y = as.numeric(data[[response]]), vi = as.numeric(data[[variance]]),
    fertilizer = factor(data$fertilizer, levels = c("animal", "plant")),
    study_ID = factor(data$study_ID)
  )
  if (any(!is.finite(result$y)) || any(!is.finite(result$vi)) || any(result$vi <= 0) ||
      anyNA(result$fertilizer) || anyNA(result$study_ID)) {
    stop("Bayesian data require finite y, positive finite vi, and complete group labels.")
  }
  result
}

bayesian_formula <- function(extended = FALSE) {
  study_term <- if (extended) "(0 + fertilizer || study_ID)" else "(1 | study_ID)"
  brms::bf(
    stats::as.formula(paste0("y | se(sqrt(vi), sigma = TRUE) ~ 0 + fertilizer + ", study_term)),
    sigma ~ 0 + fertilizer
  )
}

bayesian_priors <- function(sensitivity = FALSE) {
  # b_sigma is on log(SD). Exp(2) has prior median 0.35 for study SD;
  # the declared sensitivity weakens this to Exp(1), median 0.69.
  study_sd_prior <- if (sensitivity) "exponential(1)" else "exponential(2)"
  c(
    brms::set_prior("normal(0, 1)", class = "b"),
    brms::set_prior("normal(-1, 1)", class = "b", dpar = "sigma"),
    brms::set_prior(study_sd_prior, class = "sd", group = "study_ID")
  )
}

model_label <- function(response, extended, sensitivity = FALSE) {
  paste(response, if (extended) "extended" else "baseline",
    if (sensitivity) "sd_prior_sensitivity" else "primary", sep = "_")
}

fit_bayesian_model <- function(data, response, variance, extended, settings, sensitivity = FALSE) {
  model_data <- prepare_bayesian_data(data, response, variance)
  fit_seed <- bayesian_constants$seed +
    (if (identical(response, "lnCVR")) 100L else 0L) +
    (if (extended) 10L else 0L) +
    (if (sensitivity) 1L else 0L)
  set.seed(fit_seed)
  started <- Sys.time()
  fit <- brms::brm(
    formula = bayesian_formula(extended), data = model_data, family = gaussian(),
    prior = bayesian_priors(sensitivity), backend = "cmdstanr",
    chains = settings$chains, cores = settings$chains,
    iter = settings$iter, warmup = settings$warmup, seed = fit_seed,
    refresh = 0, silent = 2, control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  attr(fit, "study_heterogeneity_metadata") <- list(
    label = model_label(response, extended, sensitivity), response = response,
    extended = extended, sensitivity = sensitivity,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    formula = if (extended) {
      "y | se(sqrt(vi), sigma = TRUE) ~ 0 + fertilizer + (0 + fertilizer || study_ID); sigma ~ 0 + fertilizer"
    } else {
      "y | se(sqrt(vi), sigma = TRUE) ~ 0 + fertilizer + (1 | study_ID); sigma ~ 0 + fertilizer"
    },
    study_sd_prior = if (sensitivity) "exponential(1)" else "exponential(2)"
  )
  fit
}

parameter_draws <- function(fit) {
  draws <- posterior::as_draws_df(fit)
  parameters <- setdiff(names(draws), c(".chain", ".iteration", ".draw", "lp__"))
  posterior::subset_draws(draws, variable = parameters)
}

fit_diagnostics <- function(fit) {
  draws <- parameter_draws(fit)
  summary <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  nuts <- brms::nuts_params(fit)
  divergences <- sum(nuts$Value[nuts$Parameter == "divergent__"])
  data.frame(variable = summary$variable, rhat = summary$rhat,
    ess_bulk = summary$ess_bulk, ess_tail = summary$ess_tail,
    divergences = divergences, stringsAsFactors = FALSE)
}

validate_fit <- function(fit, mode) {
  diagnostic <- fit_diagnostics(fit)
  rhat_max <- if (identical(mode, "smoke")) bayesian_constants$smoke_rhat_max else bayesian_constants$full_rhat_max
  ess_min <- if (identical(mode, "smoke")) 50 else bayesian_constants$full_ess_min
  valid <- nrow(diagnostic) > 0L && all(is.finite(diagnostic$rhat)) &&
    all(is.finite(diagnostic$ess_bulk)) && all(is.finite(diagnostic$ess_tail)) &&
    max(diagnostic$rhat) <= rhat_max && min(diagnostic$ess_bulk) >= ess_min &&
    min(diagnostic$ess_tail) >= ess_min && identical(unique(diagnostic$divergences), 0)
  if (!valid) stop("Bayesian ", mode, " diagnostics failed: max Rhat=", signif(max(diagnostic$rhat), 4),
    ", min bulk ESS=", signif(min(diagnostic$ess_bulk), 4), ", min tail ESS=",
    signif(min(diagnostic$ess_tail), 4), ", divergences=", unique(diagnostic$divergences)[1L], ".")
  diagnostic
}

draw_column <- function(draws, pattern) {
  matches <- grep(pattern, names(draws), value = TRUE)
  if (length(matches) != 1L) stop("Could not uniquely resolve posterior draw column matching: ", pattern)
  matches
}

posterior_interval <- function(values, response, model, component, term) {
  quantiles <- stats::quantile(values, probs = c(0.025, 0.5, 0.975), names = FALSE)
  data.frame(response = response, model = model, component = component, term = term,
    mean = mean(values), median = quantiles[2L], q2.5 = quantiles[1L], q97.5 = quantiles[3L],
    stringsAsFactors = FALSE)
}

natural_scale_ratio_component <- function(response) {
  switch(response,
    lnRR = "mean_response_ratio",
    lnCVR = "coefficient_of_variation_ratio",
    stop("No natural-scale ratio label is declared for response: ", response)
  )
}

extract_posterior_natural_scales <- function(fit) {
  metadata <- attr(fit, "study_heterogeneity_metadata")
  draws <- posterior::as_draws_df(fit)
  output <- vector("list", 0L)
  for (category in c("animal", "plant")) {
    beta <- draws[[draw_column(draws, paste0("^b_fertilizer", category, "$"))]]
    residual_sd <- exp(draws[[draw_column(draws, paste0("^b_sigma_fertilizer", category, "$"))]])
    study_pattern <- if (metadata$extended) paste0("^sd_study_ID__fertilizer", category, "$") else "^sd_study_ID__Intercept$"
    study_sd <- draws[[draw_column(draws, study_pattern)]]
    output[[length(output) + 1L]] <- posterior_interval(beta, metadata$response, metadata$label, "category_mean", category)
    output[[length(output) + 1L]] <- posterior_interval(exp(beta), metadata$response, metadata$label,
      natural_scale_ratio_component(metadata$response), category)
    output[[length(output) + 1L]] <- posterior_interval(residual_sd, metadata$response, metadata$label, "residual_sd", category)
    output[[length(output) + 1L]] <- posterior_interval(study_sd, metadata$response, metadata$label, "study_sd", category)
  }
  residual_ratio <- exp(draws[[draw_column(draws, "^b_sigma_fertilizerplant$")]] -
    draws[[draw_column(draws, "^b_sigma_fertilizeranimal$")]])
  output[[length(output) + 1L]] <- posterior_interval(residual_ratio, metadata$response, metadata$label, "ratio", "residual_sd_plant_animal")
  if (metadata$extended) {
    study_ratio <- draws[[draw_column(draws, "^sd_study_ID__fertilizerplant$")]] /
      draws[[draw_column(draws, "^sd_study_ID__fertilizeranimal$")]]
    output[[length(output) + 1L]] <- posterior_interval(study_ratio, metadata$response, metadata$label, "ratio", "study_sd_plant_animal")
  }
  do.call(rbind, output)
}

canonical_prediction_draws <- function(beta, study_sd, residual_sd, vi_grid) {
  if (length(beta) != length(study_sd) || length(beta) != length(residual_sd) ||
      any(!is.finite(c(beta, study_sd, residual_sd, vi_grid))) ||
      any(study_sd < 0) || any(residual_sd < 0) || any(vi_grid <= 0)) {
    stop("Prediction draws require equal-length finite posterior vectors and positive vi.")
  }
  if (is.null(names(vi_grid)) || any(names(vi_grid) == "")) {
    stop("Prediction variance scenarios must be named.")
  }
  # Both targets without sampling error are sampled once per
  # outcome/model/category and reused verbatim for every vi scenario. The
  # study mean contains u; the latent effect contains u + e. For k = 1, the
  # mean of one observed effect is the observed effect itself, so both labels
  # reference one draw.
  latent_study_mean <- stats::rnorm(length(beta), beta, study_sd)
  latent_new_effect <- stats::rnorm(
    length(beta), beta, sqrt(study_sd^2 + residual_sd^2)
  )
  result <- list(
    latent_new_study_mean = stats::setNames(
      rep(list(latent_study_mean), length(vi_grid)), names(vi_grid)
    ),
    latent_new_effect = stats::setNames(
      rep(list(latent_new_effect), length(vi_grid)), names(vi_grid)
    ),
    observed_new_effect = vector("list", length(vi_grid)),
    mean_of_1_new_effects = vector("list", length(vi_grid)),
    mean_of_5_new_effects = vector("list", length(vi_grid)),
    mean_of_10_new_effects = vector("list", length(vi_grid))
  )
  names(result$observed_new_effect) <- names(vi_grid)
  names(result$mean_of_1_new_effects) <- names(vi_grid)
  names(result$mean_of_5_new_effects) <- names(vi_grid)
  names(result$mean_of_10_new_effects) <- names(vi_grid)
  for (vi_name in names(vi_grid)) {
    vi <- unname(vi_grid[[vi_name]])
    observed <- stats::rnorm(length(beta), beta, sqrt(study_sd^2 + residual_sd^2 + vi))
    result$observed_new_effect[[vi_name]] <- observed
    result$mean_of_1_new_effects[[vi_name]] <- observed
    result$mean_of_5_new_effects[[vi_name]] <- stats::rnorm(
      length(beta), beta, sqrt(study_sd^2 + (residual_sd^2 + vi) / 5)
    )
    result$mean_of_10_new_effects[[vi_name]] <- stats::rnorm(
      length(beta), beta, sqrt(study_sd^2 + (residual_sd^2 + vi) / 10)
    )
  }
  result
}

new_study_predictions <- function(fit, raw_data) {
  metadata <- attr(fit, "study_heterogeneity_metadata")
  variance <- if (identical(metadata$response, "lnRR")) "var.lnRR" else "var.lnCVR"
  data <- prepare_bayesian_data(raw_data, metadata$response, variance)
  draws <- posterior::as_draws_df(fit)
  output <- vector("list", 0L)
  set.seed(bayesian_constants$seed + if (identical(metadata$response, "lnCVR")) 500L else 400L)
  for (category in c("animal", "plant")) {
    beta <- draws[[draw_column(draws, paste0("^b_fertilizer", category, "$"))]]
    residual_sd <- exp(draws[[draw_column(draws, paste0("^b_sigma_fertilizer", category, "$"))]])
    study_pattern <- if (metadata$extended) paste0("^sd_study_ID__fertilizer", category, "$") else "^sd_study_ID__Intercept$"
    study_sd <- draws[[draw_column(draws, study_pattern)]]
    vi_category <- data$vi[data$fertilizer == category]
    vi_grid <- c(median = stats::median(vi_category),
      iqr_low = unname(stats::quantile(vi_category, .25)),
      iqr_high = unname(stats::quantile(vi_category, .75)))
    canonical_draws <- canonical_prediction_draws(beta, study_sd, residual_sd, vi_grid)
    for (vi_name in names(vi_grid)) {
      vi <- unname(vi_grid[[vi_name]])
      target_draws <- lapply(canonical_draws, `[[`, vi_name)
      for (target in names(target_draws)) {
        result <- posterior_interval(target_draws[[target]], metadata$response, metadata$label, "new_study_prediction", target)
        result$fertilizer <- category
        result$vi_scenario <- vi_name
        result$vi <- vi
        output[[length(output) + 1L]] <- result
      }
    }
  }
  do.call(rbind, output)
}

model_metadata_table <- function(fits) do.call(rbind, lapply(fits, function(fit) {
  metadata <- attr(fit, "study_heterogeneity_metadata")
  data.frame(label = metadata$label, response = metadata$response, extended = metadata$extended,
    sensitivity = metadata$sensitivity, formula = metadata$formula,
    location_prior = "normal(0, 1)", log_residual_sd_prior = "normal(-1, 1)",
    study_sd_prior = metadata$study_sd_prior, elapsed_seconds = metadata$elapsed_seconds,
    stringsAsFactors = FALSE)
}))

bayesian_cache_paths <- function(project_root, mode) {
  if (!identical(mode, "smoke") && !identical(mode, "full")) {
    stop("Bayesian cache mode must be either 'smoke' or 'full'.")
  }
  cache_dir <- file.path(project_root, "Rdata", "study_level_heterogeneity", "bayesian")
  list(
    directory = cache_dir,
    results = file.path(cache_dir, paste0("bayesian_", mode, "_results.rds")),
    models = file.path(cache_dir, paste0("bayesian_", mode, "_models.csv")),
    posterior_summary = file.path(cache_dir, paste0("bayesian_", mode, "_posterior_summary.csv")),
    predictions = file.path(cache_dir, paste0("bayesian_", mode, "_new_study_predictions.csv")),
    diagnostics = file.path(cache_dir, paste0("bayesian_", mode, "_diagnostics.csv"))
  )
}

save_bayesian_artifacts <- function(project_root, fits, summaries, predictions, diagnostics, mode, projected_full_minutes = NA_real_) {
  paths <- bayesian_cache_paths(project_root, mode)
  dir.create(paths$directory, recursive = TRUE, showWarnings = FALSE)
  metadata <- model_metadata_table(fits)
  results <- list(mode = mode, projected_full_minutes = projected_full_minutes, metadata = metadata,
    posterior_summary = summaries, new_study_predictions = predictions, diagnostics = diagnostics, fits = fits)
  saveRDS(results, paths$results)
  utils::write.csv(metadata, paths$models, row.names = FALSE)
  utils::write.csv(summaries, paths$posterior_summary, row.names = FALSE)
  utils::write.csv(predictions, paths$predictions, row.names = FALSE)
  utils::write.csv(diagnostics, paths$diagnostics, row.names = FALSE)
  invisible(results)
}

regenerate_bayesian_predictions <- function(mode) {
  require_bayesian_packages()
  project_root <- find_project_root()
  paths <- bayesian_cache_paths(project_root, mode)
  if (!file.exists(paths$results)) stop("No retained ", mode, " Bayesian cache exists.")
  retained <- readRDS(paths$results)
  if (!identical(retained$mode, mode) || !is.list(retained$fits) || length(retained$fits) == 0L) {
    stop("Retained Bayesian cache is not a valid ", mode, " artifact.")
  }
  prepared <- prepare_study_level_data(project_root)
  summaries <- do.call(rbind, lapply(retained$fits, extract_posterior_natural_scales))
  predictions <- do.call(rbind, lapply(retained$fits, new_study_predictions, raw_data = prepared$data))
  save_bayesian_artifacts(project_root, retained$fits, summaries,
    predictions, retained$diagnostics, mode, retained$projected_full_minutes)
}

run_bayesian_smoke <- function() {
  require_bayesian_packages()
  project_root <- find_project_root()
  prepared <- prepare_study_level_data(project_root)
  settings <- list(chains = bayesian_constants$smoke_chains, iter = bayesian_constants$smoke_iter, warmup = bayesian_constants$smoke_warmup)
  fit <- fit_bayesian_model(prepared$data, "lnRR", "var.lnRR", extended = TRUE, settings = settings)
  diagnostic <- validate_fit(fit, "smoke")
  elapsed <- attr(fit, "study_heterogeneity_metadata")$elapsed_seconds
  projected_full_minutes <- elapsed / 60 * (bayesian_constants$full_chains * bayesian_constants$full_iter) /
    (bayesian_constants$smoke_chains * bayesian_constants$smoke_iter)
  summaries <- extract_posterior_natural_scales(fit)
  predictions <- new_study_predictions(fit, prepared$data)
  save_bayesian_artifacts(project_root, list(fit), summaries, predictions, diagnostic, "smoke", projected_full_minutes)
  if (projected_full_minutes > bayesian_constants$fit_gate_minutes) stop("Projected full fit duration is ",
    signif(projected_full_minutes, 4), " minutes, exceeding the 90-minute gate.")
  cat("BAYESIAN SMOKE PASS\n")
  invisible(list(fit = fit, diagnostics = diagnostic, projected_full_minutes = projected_full_minutes))
}

run_bayesian_full <- function() {
  smoke <- run_bayesian_smoke()
  require_bayesian_packages()
  project_root <- find_project_root()
  prepared <- prepare_study_level_data(project_root)
  settings <- list(chains = bayesian_constants$full_chains, iter = bayesian_constants$full_iter, warmup = bayesian_constants$full_warmup)
  fits <- list()
  for (response in c("lnRR", "lnCVR")) {
    variance <- if (identical(response, "lnRR")) "var.lnRR" else "var.lnCVR"
    for (extended in c(FALSE, TRUE)) {
      fit <- fit_bayesian_model(prepared$data, response, variance, extended, settings)
      validate_fit(fit, "full")
      fits[[model_label(response, extended)]] <- fit
    }
    # One sensitivity per outcome, conditional on the smoke projection clearing 90 min.
    sensitivity_fit <- fit_bayesian_model(prepared$data, response, variance, TRUE, settings, sensitivity = TRUE)
    validate_fit(sensitivity_fit, "full")
    fits[[model_label(response, TRUE, sensitivity = TRUE)]] <- sensitivity_fit
  }
  diagnostics <- do.call(rbind, lapply(names(fits), function(name) {
    result <- fit_diagnostics(fits[[name]])
    result$model <- name
    result
  }))
  summaries <- do.call(rbind, lapply(fits, extract_posterior_natural_scales))
  predictions <- do.call(rbind, lapply(fits, new_study_predictions, raw_data = prepared$data))
  save_bayesian_artifacts(project_root, fits, summaries, predictions, diagnostics, "full", smoke$projected_full_minutes)
  invisible(list(fits = fits, diagnostics = diagnostics, projected_full_minutes = smoke$projected_full_minutes))
}

run_bayesian_self_test <- function() {
  baseline <- paste(deparse(bayesian_formula(FALSE)), collapse = " ")
  extended <- paste(deparse(bayesian_formula(TRUE)), collapse = " ")
  lnrr_ratio_component <- natural_scale_ratio_component("lnRR")
  lncvr_ratio_component <- natural_scale_ratio_component("lnCVR")
  prediction_contract <- canonical_prediction_draws(
    beta = rep(0, 8L), study_sd = rep(0.2, 8L), residual_sd = rep(0.3, 8L),
    vi_grid = c(median = 0.01, iqr_low = 0.005, iqr_high = 0.02)
  )
  stopifnot(grepl("se\\(sqrt\\(vi\\), sigma = TRUE\\)", baseline), grepl("0 \\+ fertilizer", baseline),
    grepl("\\(1 \\| study_ID\\)", baseline), grepl("\\(0 \\+ fertilizer \\|\\| study_ID\\)", extended),
    grepl("sigma ~ 0 \\+ fertilizer", extended), nrow(bayesian_priors()) == 3L,
    identical(prediction_contract$observed_new_effect$median, prediction_contract$mean_of_1_new_effects$median),
    identical(prediction_contract$latent_new_study_mean$median, prediction_contract$latent_new_study_mean$iqr_low),
    identical(prediction_contract$latent_new_effect$median, prediction_contract$latent_new_effect$iqr_high),
    !identical(bayesian_cache_paths(".", "smoke")$results, bayesian_cache_paths(".", "full")$results),
    identical(lnrr_ratio_component, "mean_response_ratio"),
    identical(lncvr_ratio_component, "coefficient_of_variation_ratio"))
  invisible(TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) {
  if ("--self-test" %in% arguments) run_bayesian_self_test()
  else if ("--regenerate-predictions" %in% arguments) {
    mode <- intersect(arguments, c("smoke", "full"))
    if (length(mode) != 1L) stop("Use --regenerate-predictions with exactly one of smoke or full.")
    regenerate_bayesian_predictions(mode)
  }
  else if ("--smoke" %in% arguments) run_bayesian_smoke()
  else if ("--full" %in% arguments) run_bayesian_full()
  else stop("Use one of --self-test, --smoke, or --full.")
}
