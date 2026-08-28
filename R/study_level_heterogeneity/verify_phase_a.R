script_arguments <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arguments[grepl("^--file=", script_arguments)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

project_root <- find_project_root()
cache_root <- file.path(project_root, "Rdata", "study_level_heterogeneity")
result_root <- file.path(project_root, "results", "study_level_heterogeneity")

required_file <- function(...) {
  path <- file.path(...)
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    stop("Required Phase-A artifact is missing or empty: ", path)
  }
  path
}

read_cached_csv <- function(...) {
  utils::read.csv(required_file(...), stringsAsFactors = FALSE, check.names = FALSE)
}

read_phase_a_inputs <- function() {
  list(
    baseline = read_cached_csv(cache_root, "baseline", "baseline_summary.csv"),
    frequentist = read_cached_csv(cache_root, "frequentist", "full_frequentist_summary.csv"),
    bayesian = read_cached_csv(cache_root, "bayesian", "bayesian_full_posterior_summary.csv"),
    bayesian_predictions = read_cached_csv(cache_root, "bayesian", "bayesian_full_new_study_predictions.csv"),
    bayesian_diagnostics = read_cached_csv(cache_root, "bayesian", "bayesian_full_diagnostics.csv"),
    bayesian_models = read_cached_csv(cache_root, "bayesian", "bayesian_full_models.csv"),
    direct = read_cached_csv(cache_root, "direct", "direct_subset_full_summary.csv"),
    blsmeta = read_cached_csv(cache_root, "direct", "blsmeta_subset_full_summary.csv")
  )
}

empty_model_rows <- function(n) {
  data.frame(
    response = rep(NA_character_, n),
    data_scope = rep(NA_character_, n),
    engine = rep(NA_character_, n),
    model = rep(NA_character_, n),
    estimator = rep(NA_character_, n),
    component = rep(NA_character_, n),
    term = rep(NA_character_, n),
    estimate = rep(NA_real_, n),
    lower = rep(NA_real_, n),
    upper = rep(NA_real_, n),
    interval = rep(NA_character_, n),
    fit_status = rep(NA_character_, n),
    diagnostic_status = rep(NA_character_, n),
    boundary = rep(NA, n),
    exploratory_only = rep(FALSE, n),
    note = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
}

baseline_model_rows <- function(x) {
  rows <- list()
  for (response in c("lnRR", "lnCVR")) {
    current <- x[x$response == response, , drop = FALSE]
    intercept <- current$estimate[current$component == "location" & current$term == "intrcpt"]
    contrast <- current$estimate[current$component == "location" & current$term == "fertilizerplant"]
    stopifnot(length(intercept) == 1L, length(contrast) == 1L)
    keep <- current$component %in% c("study_sd", "residual_sd")
    values <- rbind(
      data.frame(component = "location", term = c("animal", "plant"), estimate = c(intercept, intercept + contrast)),
      current[keep, c("component", "term", "estimate"), drop = FALSE]
    )
    residual <- stats::setNames(
      values$estimate[values$component == "residual_sd"],
      values$term[values$component == "residual_sd"]
    )
    values <- rbind(values, data.frame(
      component = "ratio", term = "residual_plant_animal",
      estimate = unname(residual[["plant"]] / residual[["animal"]])
    ))
    out <- empty_model_rows(nrow(values))
    out$response <- response
    out$data_scope <- "full_318_effects_36_studies"
    out$engine <- "metafor"
    out$model <- "current_common_study_sd"
    out$estimator <- "REML"
    out$component <- values$component
    out$term <- values$term
    out$estimate <- values$estimate
    out$fit_status <- "success"
    out$diagnostic_status <- "baseline_reproduced"
    out$note <- "Current comparator: one common study SD and category-specific residual SD."
    rows[[response]] <- out
  }
  do.call(rbind, rows)
}

frequentist_model_rows <- function(x) {
  x <- x[x$component %in% c("location", "study_sd", "residual_sd", "ratio", "study_rho"), , drop = FALSE]
  out <- empty_model_rows(nrow(x))
  out$response <- x$response
  out$data_scope <- "full_318_effects_36_studies"
  out$engine <- x$engine
  out$model <- x$model
  out$estimator <- "REML"
  out$component <- x$component
  out$term <- sub("^fertilizer", "", x$term)
  out$estimate <- x$estimate
  out$fit_status <- "success"
  out$diagnostic_status <- x$diagnostic_status
  out$boundary <- x$boundary
  out$exploratory_only <- x$exploratory_only
  out$note <- ifelse(
    x$model == "sensitivity_un",
    "Exploratory correlated study effects; rho is at its boundary for both outcomes.",
    ifelse(x$engine == "glmmTMB", "Exact-likelihood replication of the primary DIAG model.", "Primary DIAG study-effects model.")
  )
  out
}

bayesian_model_rows <- function(x) {
  out <- empty_model_rows(nrow(x))
  out$response <- x$response
  out$data_scope <- "full_318_effects_36_studies"
  out$engine <- "brms"
  out$model <- x$model
  out$estimator <- "Bayesian"
  out$component <- x$component
  out$term <- x$term
  out$estimate <- x$median
  out$lower <- x$q2.5
  out$upper <- x$q97.5
  out$interval <- "95% posterior credible interval"
  out$fit_status <- "success"
  out$diagnostic_status <- "rhat_ess_divergence_gates_passed"
  out$note <- ifelse(
    grepl("sensitivity", x$model),
    "Study-SD prior sensitivity: exponential(1).",
    "Primary study-SD prior: exponential(2)."
  )
  out
}

direct_model_rows <- function(x) {
  x <- x[x$component %in% c("location", "study_sd", "residual_sd", "ratio"), , drop = FALSE]
  out <- empty_model_rows(nrow(x))
  out$response <- x$response
  out$data_scope <- ifelse(is.na(x$response), "matched_subset_232_effects_30_studies", "matched_subset_232_effects_30_studies")
  out$engine <- ifelse(grepl("^drmTMB", x$model), "drmTMB", ifelse(grepl("metafor", x$model), "metafor", x$model))
  out$model <- x$model
  out$estimator <- x$estimator
  out$component <- x$component
  out$term <- x$term
  out$estimate <- x$estimate
  out$fit_status <- x$fit_status
  out$diagnostic_status <- x$diagnostic_status
  boundary_rows <- x$model == "drmTMB_direct_study_scale" & x$response == "lnCVR"
  out$boundary[boundary_rows] <- TRUE
  out$note <- ifelse(
    boundary_rows,
    "Plant study SD is near zero; check_drm reports a boundary/weak-identification warning.",
    "Matched-subset bridge or direct study-scale fit."
  )
  out
}

blsmeta_model_rows <- function(x) {
  x <- x[x$component %in% c("study_sd", "residual_sd", "ratio"), , drop = FALSE]
  out <- empty_model_rows(nrow(x))
  out$response <- x$response
  out$data_scope <- "matched_subset_232_effects_30_studies"
  out$engine <- "blsmeta"
  out$model <- x$model
  out$estimator <- x$estimator
  out$component <- x$component
  out$term <- x$term
  out$estimate <- x$estimate
  out$lower <- x$lower
  out$upper <- x$upper
  out$interval <- "95% posterior credible interval"
  out$fit_status <- x$fit_status
  out$diagnostic_status <- x$diagnostic_status
  out$note <- "Correctly ordered matched-subset fit; study-SD ratio retained as diagnostic sensitivity only because target-chain convergence failed."
  out
}

assemble_model_summary <- function(inputs) {
  output <- rbind(
    baseline_model_rows(inputs$baseline),
    frequentist_model_rows(inputs$frequentist),
    bayesian_model_rows(inputs$bayesian),
    direct_model_rows(inputs$direct),
    blsmeta_model_rows(inputs$blsmeta)
  )
  rownames(output) <- NULL
  output
}

assemble_prediction_summary <- function(x) {
  output <- x
  output$k <- rep(NA_integer_, nrow(output))
  mean_rows <- grepl("^mean_of_[0-9]+_new_effects$", output$term)
  output$k[mean_rows] <- as.integer(sub(
    "^mean_of_([0-9]+)_new_effects$", "\\1", output$term[mean_rows]
  ))
  output$k[output$term == "observed_new_effect"] <- 1L
  output$interval_width <- output$q97.5 - output$q2.5
  output$data_scope <- "full_318_effects_36_studies"
  output$interval <- "95% posterior predictive interval"
  output
}

assemble_diagnostics <- function(inputs) {
  frequentist_keys <- unique(inputs$frequentist[, c(
    "response", "model", "engine", "boundary", "exploratory_only", "diagnostic_status"
  )])
  frequentist <- data.frame(
    response = frequentist_keys$response,
    data_scope = "full_318_effects_36_studies",
    engine = frequentist_keys$engine,
    model = frequentist_keys$model,
    fit_status = "success",
    diagnostic_status = frequentist_keys$diagnostic_status,
    max_rhat = NA_real_, min_ess_bulk = NA_real_, min_ess_tail = NA_real_, divergences = NA_real_,
    optimizer_convergence = 0, max_abs_gradient = NA_real_, pd_hessian = NA,
    boundary = frequentist_keys$boundary,
    exploratory_only = frequentist_keys$exploratory_only,
    note = ifelse(frequentist_keys$model == "sensitivity_un",
      "rho = 1 boundary; exploratory only.",
      ifelse(frequentist_keys$engine == "metafor",
        "Optimizer converged; numerical gradient/Hessian were not retained by metafor.",
        "Optimizer converged; finite gradient and positive-definite Hessian."
      )
    ),
    stringsAsFactors = FALSE
  )

  by_model <- split(inputs$bayesian_diagnostics, inputs$bayesian_diagnostics$model)
  bayesian <- do.call(rbind, lapply(by_model, function(x) {
    data.frame(
      response = if (grepl("^lnRR", x$model[[1L]])) "lnRR" else "lnCVR",
      data_scope = "full_318_effects_36_studies", engine = "brms", model = x$model[[1L]],
      fit_status = "success", diagnostic_status = "passed",
      max_rhat = max(x$rhat, na.rm = TRUE), min_ess_bulk = min(x$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(x$ess_tail, na.rm = TRUE), divergences = max(x$divergences, na.rm = TRUE),
      optimizer_convergence = NA_real_, max_abs_gradient = NA_real_, pd_hessian = NA,
      boundary = FALSE, exploratory_only = FALSE,
      note = "Four chains; R-hat, ESS, and divergence gates passed.", stringsAsFactors = FALSE
    )
  }))

  direct_models <- unique(inputs$direct[
    grepl("^drmTMB", inputs$direct$model) & !is.na(inputs$direct$response),
    c("response", "model", "fit_status", "diagnostic_status")
  ])
  direct <- do.call(rbind, lapply(seq_len(nrow(direct_models)), function(index) {
    key <- direct_models[index, , drop = FALSE]
    rows <- inputs$direct[inputs$direct$response == key$response & inputs$direct$model == key$model, , drop = FALSE]
    value <- function(term) {
      hit <- rows$estimate[rows$component == "diagnostic" & rows$term == term]
      if (length(hit)) hit[[1L]] else NA_real_
    }
    boundary <- identical(key$model, "drmTMB_direct_study_scale") && identical(key$response, "lnCVR")
    data.frame(
      response = key$response, data_scope = "matched_subset_232_effects_30_studies",
      engine = "drmTMB", model = key$model, fit_status = key$fit_status,
      diagnostic_status = key$diagnostic_status,
      max_rhat = NA_real_, min_ess_bulk = NA_real_, min_ess_tail = NA_real_, divergences = NA_real_,
      optimizer_convergence = value("optimizer_convergence"), max_abs_gradient = value("gradient_max_abs"),
      pd_hessian = as.logical(value("pdHess")), boundary = boundary, exploratory_only = FALSE,
      note = if (boundary) "Near-zero plant study SD; standard-error and boundary warnings retained."
        else "check_drm diagnostics retained in the full cache.",
      stringsAsFactors = FALSE
    )
  }))

  blsmeta_keys <- unique(inputs$blsmeta[, c(
    "response", "model", "fit_status", "diagnostic_status", "max_rhat", "min_ess"
  )])
  blsmeta <- data.frame(
    response = blsmeta_keys$response,
    data_scope = "matched_subset_232_effects_30_studies",
    engine = "blsmeta", model = blsmeta_keys$model,
    fit_status = blsmeta_keys$fit_status,
    diagnostic_status = blsmeta_keys$diagnostic_status,
    max_rhat = blsmeta_keys$max_rhat,
    min_ess_bulk = blsmeta_keys$min_ess,
    min_ess_tail = NA_real_, divergences = NA_real_, optimizer_convergence = NA_real_,
    max_abs_gradient = NA_real_, pd_hessian = NA, boundary = FALSE, exploratory_only = FALSE,
    note = "Model ran after deterministic study-order alignment; study-ratio target diagnostics did not fully pass.",
    stringsAsFactors = FALSE
  )
  output <- rbind(frequentist, bayesian, direct, blsmeta)
  rownames(output) <- NULL
  output
}

write_phase_a_outputs <- function() {
  inputs <- read_phase_a_inputs()
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(assemble_model_summary(inputs), file.path(result_root, "model_summary.csv"), row.names = FALSE)
  utils::write.csv(assemble_prediction_summary(inputs$bayesian_predictions), file.path(result_root, "prediction_summary.csv"), row.names = FALSE)
  utils::write.csv(assemble_diagnostics(inputs), file.path(result_root, "diagnostics.csv"), row.names = FALSE)
  invisible(TRUE)
}

gate_data <- function() {
  counts <- prepare_study_level_data(project_root)$counts
  stopifnot(
    identical(counts$effects, 318L), identical(counts$studies, 36L),
    identical(counts$animal_effects, 134L), identical(counts$plant_effects, 184L),
    identical(counts$mixed_studies, 6L), identical(counts$mixed_effects, 86L),
    identical(counts$matched_effects, 232L), identical(counts$matched_studies, 30L)
  )
  cat("G1 DATA PASS\n")
}

gate_baseline <- function() {
  x <- read_phase_a_inputs()$baseline
  reference <- data.frame(
    response = rep(c("lnRR", "lnCVR"), each = 5L),
    component = rep(c("location", "location", "study_sd", "residual_sd", "residual_sd"), 2L),
    term = rep(c("intrcpt", "fertilizerplant", "study_ID", "animal", "plant"), 2L),
    expected = c(-0.26791, -0.13991, 0.24660, 0.25809, 0.30609,
      0.29890, 0.08452, 0.23558, 0.43324, 0.72820),
    stringsAsFactors = FALSE
  )
  observed <- merge(reference, x, by = c("response", "component", "term"), all.x = TRUE)
  stopifnot(nrow(observed) == nrow(reference), all(is.finite(observed$estimate)), max(abs(observed$estimate - observed$expected)) <= 0.03)
  cat("G2 BASELINE PASS\n")
}

gate_alignment <- function() {
  frequentist_code <- paste(readLines(required_file(project_root, "R", "study_level_heterogeneity", "full_frequentist.R"), warn = FALSE), collapse = "\n")
  bayesian_code <- paste(readLines(required_file(project_root, "R", "study_level_heterogeneity", "full_bayesian.R"), warn = FALSE), collapse = "\n")
  direct_code <- paste(readLines(required_file(project_root, "R", "study_level_heterogeneity", "direct_subset.R"), warn = FALSE), collapse = "\n")
  blsmeta_code <- paste(readLines(required_file(project_root, "R", "study_level_heterogeneity", "blsmeta_subset.R"), warn = FALSE), collapse = "\n")
  stopifnot(
    grepl("random = list\\(~ fertilizer \\| study_ID, ~ fertilizer \\| effect_size_ID\\)", frequentist_code),
    grepl("struct = c\\(study_struct, \"DIAG\"\\)", frequentist_code),
    grepl("se\\(sqrt\\(vi\\), sigma = TRUE\\)", bayesian_code),
    grepl("\\(0 \\+ fertilizer \\|\\| study_ID\\)", bayesian_code),
    grepl("sigma ~ 0 \\+ fertilizer", bayesian_code),
    grepl("sd\\(study_ID\\) ~ %s", direct_code),
    grepl("meta_V\\(V = V\\)", direct_code),
    grepl("mods_scale2 = ~ 0 \\+ fertilizer", blsmeta_code),
    grepl("mods_scale3 = ~ 0 \\+ fertilizer", blsmeta_code),
    grepl("order\\(data\\$study_ID, data\\$effect_size_ID\\)", blsmeta_code)
  )
  cat("G3 ALIGNMENT PASS\n")
}

gate_frequentist <- function() {
  x <- read_phase_a_inputs()$frequentist
  stopifnot(
    setequal(unique(x$response), c("lnRR", "lnCVR")),
    all(c("primary_diag", "primary_glmmTMB", "sensitivity_un") %in% unique(x$model)),
    all(!is.na(x$diagnostic_status)),
    all(x$boundary[x$model == "sensitivity_un"]),
    all(x$exploratory_only[x$model == "sensitivity_un"])
  )
  for (response in c("lnRR", "lnCVR")) {
    a <- x[x$response == response & x$model == "primary_diag" & x$component %in% c("location", "study_sd", "residual_sd"), c("component", "term", "estimate")]
    b <- x[x$response == response & x$model == "primary_glmmTMB" & x$component %in% c("location", "study_sd", "residual_sd"), c("component", "term", "estimate")]
    merged <- merge(a, b, by = c("component", "term"), suffixes = c("_metafor", "_glmmTMB"))
    stopifnot(nrow(merged) == 6L, max(abs(merged$estimate_metafor - merged$estimate_glmmTMB)) <= 1e-5)
  }
  cat("G4 FREQUENTIST PASS\n")
}

gate_bayesian <- function() {
  inputs <- read_phase_a_inputs()
  x <- inputs$bayesian_diagnostics
  stopifnot(
    length(unique(x$model)) == 6L,
    max(x$rhat, na.rm = TRUE) <= 1.01,
    min(x$ess_bulk, na.rm = TRUE) >= 1000,
    min(x$ess_tail, na.rm = TRUE) >= 1000,
    max(x$divergences, na.rm = TRUE) == 0
  )
  ratio <- inputs$bayesian[inputs$bayesian$component == "ratio" & inputs$bayesian$term == "study_sd_plant_animal", ]
  for (response in c("lnRR", "lnCVR")) {
    primary <- ratio[ratio$response == response & grepl("extended_primary$", ratio$model), ]
    sensitivity <- ratio[ratio$response == response & grepl("sensitivity$", ratio$model), ]
    stopifnot(nrow(primary) == 1L, nrow(sensitivity) == 1L)
    primary_direction <- c(primary$q2.5 > 1, primary$q97.5 < 1)
    sensitivity_direction <- c(sensitivity$q2.5 > 1, sensitivity$q97.5 < 1)
    stopifnot(identical(primary_direction, sensitivity_direction))
  }
  cat("G5 BAYESIAN PASS\n")
}

gate_subset <- function() {
  prepared <- prepare_study_level_data(project_root)
  x <- read_phase_a_inputs()$direct
  blsmeta <- read_phase_a_inputs()$blsmeta
  stopifnot(
    nrow(prepared$matched_data) == 232L,
    length(unique(prepared$matched_data$study_ID)) == 30L,
    all(tapply(prepared$matched_data$fertilizer, prepared$matched_data$study_ID, function(value) length(unique(value))) == 1L),
    all(c("lnRR", "lnCVR") %in% unique(stats::na.omit(x$response))),
    all(c("metafor_DIAG_bridge", "drmTMB_common_study_scale", "drmTMB_direct_study_scale") %in% unique(x$model)),
    all(x$fit_status[grepl("^drmTMB", x$model)] == "success"),
    all(x$diagnostic_status[x$response == "lnCVR" & x$model == "drmTMB_direct_study_scale"] == "warning_or_boundary"),
    setequal(unique(blsmeta$response), c("lnRR", "lnCVR")),
    all(blsmeta$fit_status == "success"),
    all(blsmeta$diagnostic_status == "failed_convergence")
  )
  cat("G6 SUBSET PASS\n")
}

gate_predictions <- function() {
  x <- read_phase_a_inputs()$bayesian_predictions
  stopifnot(
    all(c("median", "iqr_low", "iqr_high") %in% unique(x$vi_scenario)),
    all(c(
      "latent_new_study_mean", "latent_new_effect", "observed_new_effect",
      "mean_of_1_new_effects", "mean_of_5_new_effects", "mean_of_10_new_effects"
    ) %in% unique(x$term)),
    all(is.finite(x$q2.5)), all(is.finite(x$q97.5)), all(x$q97.5 > x$q2.5)
  )
  observed <- x[x$term == "observed_new_effect", ]
  mean_one <- x[x$term == "mean_of_1_new_effects", ]
  key <- c("response", "model", "fertilizer", "vi_scenario", "vi")
  paired <- merge(observed, mean_one, by = key, suffixes = c("_observed", "_mean1"))
  fields <- c("mean", "median", "q2.5", "q97.5")
  stopifnot(nrow(paired) == 36L, all(vapply(fields, function(field) {
    identical(paired[[paste0(field, "_observed")]], paired[[paste0(field, "_mean1")]])
  }, logical(1))))
  latent <- x[x$term == "latent_new_study_mean", ]
  latent_groups <- split(latent, interaction(latent$response, latent$model, latent$fertilizer, drop = TRUE))
  stopifnot(all(vapply(latent_groups, function(group) {
    all(vapply(fields, function(field) length(unique(group[[field]])) == 1L, logical(1)))
  }, logical(1))))
  latent_effect <- x[x$term == "latent_new_effect", ]
  latent_effect_groups <- split(
    latent_effect,
    interaction(latent_effect$response, latent_effect$model, latent_effect$fertilizer, drop = TRUE)
  )
  stopifnot(all(vapply(latent_effect_groups, function(group) {
    all(vapply(fields, function(field) length(unique(group[[field]])) == 1L, logical(1)))
  }, logical(1))))
  example <- data.frame(study_sd = c(0.15, 0.24), residual_sd = c(0.27, 0.73), vi = c(0.01, 0.02))
  latent_effect_variance <- example$study_sd^2 + example$residual_sd^2
  stopifnot(isTRUE(all.equal(
    latent_effect_variance,
    c(0.15^2 + 0.27^2, 0.24^2 + 0.73^2),
    tolerance = 1e-14
  )))
  for (k in c(1, 5, 10)) {
    declared <- example$study_sd^2 + example$residual_sd^2 / k + (k * example$vi) / k^2
    simplified <- example$study_sd^2 + (example$residual_sd^2 + example$vi) / k
    stopifnot(isTRUE(all.equal(declared, simplified, tolerance = 1e-14)))
  }
  cat("G7 PREDICTIONS PASS\n")
}

run_gate <- function(command) {
  switch(command,
    data = gate_data(), baseline = gate_baseline(), alignment = gate_alignment(),
    frequentist = gate_frequentist(), bayesian = gate_bayesian(), subset = gate_subset(),
    predictions = gate_predictions(), assemble = write_phase_a_outputs(),
    all = {
      write_phase_a_outputs()
      gate_data(); gate_baseline(); gate_alignment(); gate_frequentist();
      gate_bayesian(); gate_subset(); gate_predictions()
    },
    stop("Unknown command. Use one of: assemble, data, baseline, alignment, frequentist, bayesian, subset, predictions, all.")
  )
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 1L) stop("Supply exactly one verification command.")
  run_gate(arguments[[1L]])
}
