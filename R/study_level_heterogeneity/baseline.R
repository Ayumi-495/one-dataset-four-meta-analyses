script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

fit_heteroscedastic_baseline <- function(data, response, variance) {
  metafor::rma.mv(
    yi = data[[response]],
    V = data[[variance]],
    mods = ~ fertilizer,
    random = list(~ 1 | study_ID, ~ fertilizer | effect_size_ID),
    struct = "DIAG",
    test = "t",
    method = "REML",
    data = data
  )
}

extract_baseline <- function(model, response) {
  residual_sd <- sqrt(model$tau2)
  names(residual_sd) <- levels(model$data$fertilizer)
  fit_statistics <- model$fit.stats[, "REML"]
  names(fit_statistics) <- rownames(model$fit.stats)
  convergence <- if (is.null(model$opt.res$convergence)) NA_integer_ else {
    as.integer(model$opt.res$convergence)
  }

  list(
    response = response,
    location = stats::setNames(as.numeric(model$beta), rownames(model$beta)),
    study_sd = sqrt(model$sigma2[1]),
    residual_sd = residual_sd,
    convergence = convergence,
    fit = fit_statistics,
    model = model
  )
}

fit_baselines <- function(data) {
  lnrr_model <- fit_heteroscedastic_baseline(data, "lnRR", "var.lnRR")
  lncvr_model <- fit_heteroscedastic_baseline(data, "lnCVR", "var.lnCVR")

  list(
    lnRR = extract_baseline(lnrr_model, "lnRR"),
    lnCVR = extract_baseline(lncvr_model, "lnCVR")
  )
}

baseline_summary_table <- function(baselines) {
  do.call(rbind, lapply(baselines, function(result) {
    rbind(
      data.frame(
        response = result$response,
        component = "location",
        term = names(result$location),
        estimate = unname(result$location),
        row.names = NULL
      ),
      data.frame(
        response = result$response,
        component = "study_sd",
        term = "study_ID",
        estimate = result$study_sd,
        row.names = NULL
      ),
      data.frame(
        response = result$response,
        component = "residual_sd",
        term = names(result$residual_sd),
        estimate = unname(result$residual_sd),
        row.names = NULL
      ),
      data.frame(
        response = result$response,
        component = "convergence",
        term = "optimizer",
        estimate = result$convergence,
        row.names = NULL
      ),
      data.frame(
        response = result$response,
        component = "fit",
        term = names(result$fit),
        estimate = unname(result$fit),
        row.names = NULL
      )
    )
  }))
}

save_baseline_artifacts <- function(prepared, baselines, project_root) {
  cache_dir <- file.path(
    project_root,
    "Rdata",
    "study_level_heterogeneity",
    "baseline"
  )
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list(
    counts = prepared$counts,
    mixed_studies = prepared$mixed_studies,
    baselines = baselines
  )
  saveRDS(results, file.path(cache_dir, "baseline_results.rds"))
  utils::write.csv(
    baseline_summary_table(baselines),
    file.path(cache_dir, "baseline_summary.csv"),
    row.names = FALSE
  )
}

assert_near <- function(actual, expected, tolerance = 0.003) {
  if (length(actual) != length(expected) ||
      any(!is.finite(actual)) ||
      any(abs(actual - expected) > tolerance)) {
    stop("Baseline estimate is outside its expected tolerance.")
  }
}

validate_baseline <- function(result) {
  expected_location <- c("intrcpt", "fertilizerplant")
  expected_residual <- c("animal", "plant")
  valid <- identical(result$convergence, 0L) &&
    identical(names(result$location), expected_location) &&
    identical(names(result$residual_sd), expected_residual) &&
    length(result$study_sd) == 1L &&
    length(result$fit) > 0L &&
    all(is.finite(result$location)) &&
    all(is.finite(result$study_sd)) &&
    all(is.finite(result$residual_sd)) &&
    all(is.finite(result$fit))

  if (!valid) {
    stop("Baseline validation failed.")
  }
}

validate_baselines <- function(baselines) {
  if (!identical(names(baselines), c("lnRR", "lnCVR"))) {
    stop("Baseline validation failed.")
  }
  lapply(baselines, validate_baseline)
  invisible(baselines)
}

validate_smoke <- function(baselines) {
  lnrr <- baselines$lnRR
  lncvr <- baselines$lnCVR
  assert_near(lnrr$location, c(-0.2679, -0.1399))
  assert_near(lnrr$study_sd, 0.2466)
  assert_near(lnrr$residual_sd, c(0.2581, 0.3061))
  assert_near(lncvr$location, c(0.2989, 0.0845))
  assert_near(lncvr$study_sd, 0.2356)
  assert_near(lncvr$residual_sd, c(0.4332, 0.7282))
}

run_baseline <- function(smoke = FALSE) {
  project_root <- find_project_root()
  prepared <- prepare_study_level_data(project_root)
  baselines <- fit_baselines(prepared$data)
  validate_baselines(baselines)

  if (smoke) {
    validate_smoke(baselines)
  }
  save_baseline_artifacts(prepared, baselines, project_root)

  if (smoke) {
    cat("BASELINE SMOKE PASS\n")
  }

  invisible(list(prepared = prepared, baselines = baselines))
}

if (sys.nframe() == 0L) {
  run_baseline(smoke = "--smoke" %in% commandArgs(trailingOnly = TRUE))
}
