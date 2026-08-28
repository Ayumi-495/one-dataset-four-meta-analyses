script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

# This script isolates the matched subset where fertilizer is a study-level
# predictor.  That condition is necessary for `sd(study_ID) ~ fertilizer`.
assert_matched_subset <- function(data) {
  study_fertilizer_count <- tapply(
    data$fertilizer,
    data$study_ID,
    function(x) length(unique(x))
  )
  stopifnot(
    nrow(data) == 232L,
    length(unique(data$study_ID)) == 30L,
    all(study_fertilizer_count == 1L)
  )
  invisible(data)
}

package_status <- function() {
  available <- requireNamespace("blsmeta", quietly = TRUE)
  data.frame(
    package = "blsmeta",
    status = if (available) {
      "available"
    } else {
      "not_installed"
    },
    action = if (available) "available_not_used" else "not_installed_not_used",
    stringsAsFactors = FALSE
  )
}

fit_metafor_bridge <- function(data, response, variance) {
  metafor::rma.mv(
    yi = data[[response]],
    V = data[[variance]],
    mods = ~ fertilizer,
    random = list(~ fertilizer | study_ID, ~ fertilizer | effect_size_ID),
    struct = "DIAG",
    test = "t",
    method = "REML",
    data = data
  )
}

extract_metafor_bridge <- function(model, response) {
  location <- as.numeric(model$beta)
  names(location) <- rownames(model$beta)
  residual_sd <- sqrt(model$gamma2)
  names(residual_sd) <- levels(model$data$fertilizer)
  study_sd <- sqrt(model$tau2)
  names(study_sd) <- levels(model$data$fertilizer)
  list(
    response = response,
    estimator = "REML",
    location = c(
      animal = location[["intrcpt"]],
      plant = location[["intrcpt"]] + location[["fertilizerplant"]]
    ),
    residual_sd = residual_sd,
    study_sd = study_sd,
    residual_sd_ratio_plant_animal = unname(residual_sd[["plant"]] / residual_sd[["animal"]]),
    study_sd_ratio_plant_animal = unname(study_sd[["plant"]] / study_sd[["animal"]]),
    convergence = if (is.null(model$opt.res$convergence)) NA_integer_ else as.integer(model$opt.res$convergence),
    fit = model$fit.stats[, "REML"],
    model = model
  )
}

make_drm_formula <- function(response, study_scale) {
  response_formula <- stats::as.formula(sprintf(
    "%s ~ fertilizer + (1 | study_ID) + drmTMB::meta_V(V = V)",
    response
  ))
  do.call(
    drmTMB::bf,
    list(
      response_formula,
      stats::as.formula("sigma ~ 0 + fertilizer"),
      stats::as.formula(sprintf("sd(study_ID) ~ %s", study_scale))
    )
  )
}

fit_drm_model <- function(data, response, variance, study_scale) {
  drm_data <- data
  drm_data$V <- drm_data[[variance]]
  drmTMB::drmTMB(
    formula = make_drm_formula(response, study_scale),
    family = stats::gaussian(),
    data = drm_data,
    REML = FALSE
  )
}

capture_drm_model <- function(data, response, variance, study_scale, model_label) {
  started <- Sys.time()
  result <- tryCatch(
    fit_drm_model(data, response, variance, study_scale),
    error = function(error) error
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  if (inherits(result, "error")) {
    return(list(
      response = response,
      model = model_label,
      study_scale_formula = study_scale,
      estimator = "ML",
      fit_status = "capability_failure",
      diagnostic_status = "not_assessed",
      error = conditionMessage(result),
      elapsed_seconds = elapsed_seconds,
      fit = NULL
    ))
  }

  list(
    response = response,
    model = model_label,
    study_scale_formula = study_scale,
    estimator = "ML",
    fit_status = "success",
    diagnostic_status = "pending",
    error = NA_character_,
    elapsed_seconds = elapsed_seconds,
    fit = result
  )
}

fixed_scale_by_fertilizer <- function(coefficients, component, direct) {
  coefficient_names <- names(coefficients)
  coefficients <- stats::setNames(as.numeric(coefficients), coefficient_names)
  if (direct) {
    values <- coefficients[grepl("fertilizer", coefficient_names)]
    names(values) <- sub("^.*fertilizer", "", names(values))
  } else {
    intercept <- coefficients[grepl("Intercept", coefficient_names)]
    values <- c(animal = intercept[[1L]], plant = intercept[[1L]])
  }
  exp(values)
}

extract_drm_model <- function(captured) {
  if (!identical(captured$fit_status, "success")) {
    return(captured)
  }

  fit <- captured$fit
  mu <- stats::coef(fit, "mu")
  mu <- as.numeric(mu)
  names(mu) <- names(stats::coef(fit, "mu"))
  location <- c(
    animal = mu[[grep("Intercept", names(mu), value = TRUE)[[1L]]]],
    plant = mu[[grep("Intercept", names(mu), value = TRUE)[[1L]]]] +
      mu[[grep("fertilizerplant", names(mu), value = TRUE)[[1L]]]]
  )
  residual_sd <- fixed_scale_by_fertilizer(
    stats::coef(fit, "sigma"),
    component = "sigma",
    direct = TRUE
  )
  study_sd <- fixed_scale_by_fertilizer(
    stats::coef(fit, "sd(study_ID)"),
    component = "sd(study_ID)",
    direct = identical(captured$study_scale_formula, "0 + fertilizer")
  )
  check <- drmTMB::check_drm(fit)
  gradient <- fit$obj$gr(fit$opt$par)
  check_table <- as.data.frame(check, stringsAsFactors = FALSE)
  check_issues <- check_table[check_table$status != "ok", , drop = FALSE]
  boundary_row <- check_table[
    check_table$check == "random_effect_sd_boundary",
    ,
    drop = FALSE
  ]

  extracted <- c(
    captured,
    list(
      location = location,
      residual_sd = residual_sd,
      study_sd = study_sd,
      residual_sd_ratio_plant_animal = unname(residual_sd[["plant"]] / residual_sd[["animal"]]),
      study_sd_ratio_plant_animal = unname(study_sd[["plant"]] / study_sd[["animal"]]),
      diagnostics = list(
        check_drm_status = if (isTRUE(attr(check, "ok"))) "ok" else "issues_detected",
        check_drm = check_table,
        check_drm_issues = check_issues,
        check_drm_messages = stats::setNames(check_table$message, check_table$check),
        check_drm_ok = isTRUE(attr(check, "ok")),
        random_effect_sd_boundary = list(
          flag = nrow(boundary_row) == 1L && !identical(boundary_row$status, "ok"),
          status = if (nrow(boundary_row) == 1L) boundary_row$status[[1L]] else NA_character_,
          value = if (nrow(boundary_row) == 1L) boundary_row$value[[1L]] else NA_character_,
          message = if (nrow(boundary_row) == 1L) boundary_row$message[[1L]] else NA_character_
        ),
        pdHess = isTRUE(fit$sdr$pdHess),
        optimizer_convergence = fit$opt$convergence,
        gradient_max_abs = max(abs(gradient))
      )
    )
  )
  extracted$diagnostic_status <- if (isTRUE(attr(check, "ok"))) {
    "ok"
  } else if (nrow(boundary_row) == 1L && !identical(boundary_row$status, "ok")) {
    "warning_or_boundary"
  } else {
    "check_drm_issues"
  }
  extracted
}

direct_summary_table <- function(metafor_results, drm_results, packages) {
  rows <- list()
  add_row <- function(response, model, estimator, status, component, term, estimate, detail = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      response = response,
      model = model,
      estimator = estimator,
      status = status,
      component = component,
      term = term,
      estimate = estimate,
      detail = detail,
      stringsAsFactors = FALSE
    )
  }

  for (response in names(metafor_results)) {
    result <- metafor_results[[response]]
    for (term in names(result$location)) add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "location", term, result$location[[term]])
    for (term in names(result$residual_sd)) add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "residual_sd", term, result$residual_sd[[term]])
    for (term in names(result$study_sd)) add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "study_sd", term, result$study_sd[[term]])
    add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "ratio", "residual_plant_animal", result$residual_sd_ratio_plant_animal)
    add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "ratio", "study_plant_animal", result$study_sd_ratio_plant_animal)
    add_row(response, "metafor_DIAG_bridge", result$estimator, "success", "diagnostic", "optimizer_convergence", result$convergence)
  }

  for (result in drm_results) {
    if (!identical(result$fit_status, "success")) {
      add_row(result$response, result$model, result$estimator, result$fit_status, "capability", "error", NA_real_, result$error)
      add_row(result$response, result$model, result$estimator, result$fit_status, "runtime", "elapsed_seconds", result$elapsed_seconds)
      next
    }
    add_row(result$response, result$model, result$estimator, result$fit_status, "formula", "sd(study_ID)", NA_real_, result$study_scale_formula)
    add_row(result$response, result$model, result$estimator, result$fit_status, "runtime", "elapsed_seconds", result$elapsed_seconds)
    for (term in names(result$location)) add_row(result$response, result$model, result$estimator, result$fit_status, "location", term, result$location[[term]])
    for (term in names(result$residual_sd)) add_row(result$response, result$model, result$estimator, result$fit_status, "residual_sd", term, result$residual_sd[[term]])
    for (term in names(result$study_sd)) add_row(result$response, result$model, result$estimator, result$fit_status, "study_sd", term, result$study_sd[[term]])
    add_row(result$response, result$model, result$estimator, result$fit_status, "ratio", "residual_plant_animal", result$residual_sd_ratio_plant_animal)
    add_row(result$response, result$model, result$estimator, result$fit_status, "ratio", "study_plant_animal", result$study_sd_ratio_plant_animal)
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "check_drm_ok", as.numeric(result$diagnostics$check_drm_ok))
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "check_drm_status", NA_real_, result$diagnostics$check_drm_status)
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "pdHess", as.numeric(result$diagnostics$pdHess))
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "optimizer_convergence", result$diagnostics$optimizer_convergence)
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "gradient_max_abs", result$diagnostics$gradient_max_abs)
    add_row(result$response, result$model, result$estimator, result$fit_status, "diagnostic", "random_effect_sd_boundary", as.numeric(result$diagnostics$random_effect_sd_boundary$flag), result$diagnostics$random_effect_sd_boundary$message)
    for (check_row in seq_len(nrow(result$diagnostics$check_drm))) {
      check <- result$diagnostics$check_drm[check_row, , drop = FALSE]
      add_row(
        result$response,
        result$model,
        result$estimator,
        result$fit_status,
        "check_drm",
        check$check[[1L]],
        NA_real_,
        paste0("status=", check$status[[1L]], "; value=", check$value[[1L]], "; ", check$message[[1L]])
      )
    }
    for (issue_row in seq_len(nrow(result$diagnostics$check_drm_issues))) {
      issue <- result$diagnostics$check_drm_issues[issue_row, , drop = FALSE]
      add_row(
        result$response,
        result$model,
        result$estimator,
        result$fit_status,
        "check_drm_issue",
        issue$check[[1L]],
        NA_real_,
        paste0("status=", issue$status[[1L]], "; value=", issue$value[[1L]], "; ", issue$message[[1L]])
      )
    }
  }

  add_row(NA_character_, "blsmeta", NA_character_, packages$status, "package_status", packages$package, NA_real_, packages$action)
  summary <- do.call(rbind, rows)
  summary$fit_status <- summary$status
  summary$diagnostic_status <- "not_assessed"
  for (result in drm_results) {
    rows_for_model <- summary$response == result$response & summary$model == result$model
    summary$diagnostic_status[rows_for_model] <- result$diagnostic_status
  }
  summary
}

validate_direct_result <- function(result) {
  if (!identical(result$fit_status, "success")) return(invisible(result))
  stopifnot(
    identical(names(result$location), c("animal", "plant")),
    identical(names(result$residual_sd), c("animal", "plant")),
    identical(names(result$study_sd), c("animal", "plant")),
    all(is.finite(result$location)),
    all(is.finite(result$residual_sd)),
    all(is.finite(result$study_sd)),
    is.finite(result$diagnostics$gradient_max_abs),
    result$diagnostic_status %in% c("ok", "warning_or_boundary", "check_drm_issues")
  )
  invisible(result)
}

validate_artifact_payload <- function(results, summary, run_label) {
  expected_responses <- if (identical(run_label, "full")) c("lnRR", "lnCVR") else "lnRR"
  expected_drm <- if (identical(run_label, "full")) {
    c(
      "lnRR_drmTMB_common_study_scale",
      "lnRR_drmTMB_direct_study_scale",
      "lnCVR_drmTMB_common_study_scale",
      "lnCVR_drmTMB_direct_study_scale"
    )
  } else {
    "lnRR_drmTMB_direct_study_scale"
  }
  expected_direct_models <- sub("^[^_]+_", "", expected_drm)
  expected_summary_models <- c("metafor_DIAG_bridge", unique(expected_direct_models))
  direct_rows <- summary[
    !is.na(summary$response) & summary$model %in% expected_direct_models,
    ,
    drop = FALSE
  ]
  summary_direct_keys <- unique(paste(direct_rows$response, direct_rows$model, sep = "_"))
  summary_models <- unique(summary$model[!is.na(summary$response)])
  required_components <- c(
    "formula", "runtime", "location", "residual_sd", "study_sd",
    "ratio", "diagnostic", "check_drm"
  )
  stopifnot(
    identical(names(results$metafor_bridge), expected_responses),
    identical(names(results$drmTMB), expected_drm),
    all(c("fit_status", "diagnostic_status") %in% names(summary)),
    setequal(unique(summary$response[!is.na(summary$response)]), expected_responses),
    setequal(summary_models, expected_summary_models),
    setequal(summary_direct_keys, expected_drm),
    nrow(direct_rows) > 0L
  )

  for (result_name in expected_drm) {
    result <- results$drmTMB[[result_name]]
    model_rows <- direct_rows[
      direct_rows$response == result$response & direct_rows$model == result$model,
      ,
      drop = FALSE
    ]
    stopifnot(
      nrow(model_rows) > 0L,
      all(model_rows$fit_status == result$fit_status),
      all(model_rows$diagnostic_status == result$diagnostic_status)
    )
    if (identical(result$fit_status, "success")) {
      stopifnot(all(required_components %in% model_rows$component))
    } else {
      stopifnot(all(c("capability", "runtime") %in% model_rows$component))
    }
  }
  invisible(TRUE)
}

validate_artifact_negative_self_test <- function(results, summary, run_label) {
  malformed <- summary[summary$model == "metafor_DIAG_bridge", , drop = FALSE]
  rejected <- inherits(
    tryCatch(
      validate_artifact_payload(results, malformed, run_label),
      error = function(error) error
    ),
    "error"
  )
  if (!rejected) stop("Artifact validator accepted a malformed CSV without drmTMB rows.")
  invisible(TRUE)
}

validate_persisted_artifacts <- function(results_path, summary_path, run_label) {
  persisted_results <- readRDS(results_path)
  persisted_summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  validate_artifact_payload(persisted_results, persisted_summary, run_label)
}

save_direct_artifacts <- function(prepared, metafor_results, drm_results, packages, project_root, run_label) {
  cache_dir <- file.path(project_root, "Rdata", "study_level_heterogeneity", "direct")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  results <- list(
    subset_counts = list(effects = nrow(prepared$matched_data), studies = length(unique(prepared$matched_data$study_ID))),
    package_status = packages,
    metafor_bridge = metafor_results,
    drmTMB = drm_results
  )
  summary <- direct_summary_table(metafor_results, drm_results, packages)
  validate_artifact_payload(results, summary, run_label)
  validate_artifact_negative_self_test(results, summary, run_label)

  results_path <- file.path(cache_dir, sprintf("direct_subset_%s_results.rds", run_label))
  summary_path <- file.path(cache_dir, sprintf("direct_subset_%s_summary.csv", run_label))
  results_temp <- tempfile(pattern = "direct_subset_results_", tmpdir = cache_dir, fileext = ".rds")
  summary_temp <- tempfile(pattern = "direct_subset_summary_", tmpdir = cache_dir, fileext = ".csv")
  on.exit(unlink(c(results_temp, summary_temp)), add = TRUE)

  saveRDS(results, results_temp)
  utils::write.csv(
    summary,
    summary_temp,
    row.names = FALSE
  )
  stopifnot(file.rename(results_temp, results_path), file.rename(summary_temp, summary_path))
  validate_persisted_artifacts(results_path, summary_path, run_label)
}

run_direct_subset <- function(smoke = FALSE) {
  project_root <- find_project_root()
  prepared <- prepare_study_level_data(project_root)
  data <- prepared$matched_data
  assert_matched_subset(data)
  packages <- package_status()

  responses <- if (smoke) "lnRR" else c("lnRR", "lnCVR")
  variances <- if (smoke) "var.lnRR" else c("var.lnRR", "var.lnCVR")
  names(variances) <- responses

  metafor_results <- lapply(responses, function(response) {
    extract_metafor_bridge(fit_metafor_bridge(data, response, variances[[response]]), response)
  })
  names(metafor_results) <- responses

  drm_specs <- if (smoke) {
    list(list(response = "lnRR", variance = "var.lnRR", study_scale = "0 + fertilizer", model = "drmTMB_direct_study_scale"))
  } else {
    unlist(lapply(responses, function(response) list(
      list(response = response, variance = variances[[response]], study_scale = "1", model = "drmTMB_common_study_scale"),
      list(response = response, variance = variances[[response]], study_scale = "0 + fertilizer", model = "drmTMB_direct_study_scale")
    )), recursive = FALSE)
  }
  drm_results <- lapply(drm_specs, function(specification) {
    extract_drm_model(capture_drm_model(
      data = data,
      response = specification$response,
      variance = specification$variance,
      study_scale = specification$study_scale,
      model_label = specification$model
    ))
  })
  names(drm_results) <- vapply(drm_specs, function(x) paste(x$response, x$model, sep = "_"), character(1))
  lapply(drm_results, validate_direct_result)

  save_direct_artifacts(
    prepared = prepared,
    metafor_results = metafor_results,
    drm_results = drm_results,
    packages = packages,
    project_root = project_root,
    run_label = if (smoke) "smoke" else "full"
  )

  if (smoke) {
    # An unsupported direct model is retained as a capability finding, not
    # silently replaced by another package or likelihood.
    stopifnot(
      identical(packages$status, "not_installed"),
      identical(names(metafor_results), "lnRR"),
      drm_results[[1L]]$fit_status %in% c("success", "capability_failure")
    )
    cat("DIRECT SMOKE PASS\n")
  }

  invisible(list(
    prepared = prepared,
    package_status = packages,
    metafor_bridge = metafor_results,
    drmTMB = drm_results
  ))
}

if (sys.nframe() == 0L) {
  run_direct_subset(smoke = "--smoke" %in% commandArgs(trailingOnly = TRUE))
}
