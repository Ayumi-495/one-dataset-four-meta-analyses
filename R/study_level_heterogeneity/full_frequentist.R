script_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_argument[grepl("^--file=", script_argument)])
if (length(script_path) == 1L) {
  source(file.path(dirname(normalizePath(script_path)), "common.R"))
} else {
  source(file.path("R", "study_level_heterogeneity", "common.R"))
}

fertilizer_levels <- function(data) {
  levels(data$fertilizer)
}

fit_metafor_location_scale <- function(data, response, variance, study_struct = "DIAG") {
  metafor::rma.mv(
    yi = data[[response]],
    V = data[[variance]],
    mods = ~ 0 + fertilizer,
    random = list(~ fertilizer | study_ID, ~ fertilizer | effect_size_ID),
    struct = c(study_struct, "DIAG"),
    test = "t",
    method = "REML",
    data = data
  )
}

extract_metafor_hessian <- function(model) {
  hessian <- model$hessian
  if (!is.matrix(hessian) || !all(is.finite(hessian))) {
    return(list(
      available = FALSE,
      status = "unavailable",
      reason = paste(
        "metafor did not retain a finite Hessian in this fit;",
        "no Hessian-scale diagnostic is available"
      )
    ))
  }

  eigenvalues <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  list(
    available = TRUE,
    status = if (all(eigenvalues > 0)) "positive_definite" else "not_positive_definite",
    reason = NA_character_,
    source_scale = "optimizer_parameter_scale",
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    positive_definite = all(eigenvalues > 0)
  )
}

extract_metafor_fit <- function(model, response, study_struct) {
  categories <- fertilizer_levels(model$data)
  location <- as.numeric(model$beta)
  names(location) <- rownames(model$beta)
  study_sd <- sqrt(model$tau2)
  residual_sd <- sqrt(model$gamma2)
  names(study_sd) <- categories
  names(residual_sd) <- categories
  fit_statistics <- model$fit.stats[, "REML"]
  names(fit_statistics) <- rownames(model$fit.stats)
  optimizer <- model$opt.res
  optimizer_diagnostics <- list(
    available = TRUE,
    status = if (identical(optimizer$convergence, 0L)) "converged" else "not_converged",
    reason = NA_character_,
    convergence = if (is.null(optimizer$convergence)) NA_integer_ else {
      as.integer(optimizer$convergence)
    },
    message = if (is.null(optimizer$message)) NA_character_ else {
      as.character(optimizer$message)
    },
    iterations = if (is.null(optimizer$iterations)) NA_integer_ else {
      as.integer(optimizer$iterations)
    },
    objective = if (is.null(optimizer$objective)) NA_real_ else {
      as.numeric(optimizer$objective)
    }
  )
  gradient_diagnostics <- list(
    available = FALSE,
    status = "unavailable",
    reason = "metafor::rma.mv() did not retain a numerical optimizer gradient"
  )
  natural_sd <- c(study_sd, residual_sd)
  boundary_tolerance <- 1e-6
  rho <- if (identical(study_struct, "UN")) as.numeric(model$rho) else NA_real_

  list(
    engine = "metafor",
    response = response,
    study_structure = study_struct,
    location = location,
    study_sd = study_sd,
    residual_sd = residual_sd,
    study_sd_ratio_plant_animal = unname(study_sd["plant"] / study_sd["animal"]),
    residual_sd_ratio_plant_animal = unname(residual_sd["plant"] / residual_sd["animal"]),
    rho = rho,
    fit = fit_statistics,
    diagnostics = list(
      natural_sd_scale = TRUE,
      optimizer = optimizer_diagnostics,
      gradient = gradient_diagnostics,
      hessian = extract_metafor_hessian(model),
      boundary = list(
        available = TRUE,
        status = if (any(natural_sd <= boundary_tolerance) ||
          (!is.na(rho) && abs(rho) >= 0.999)) "at_or_near_boundary" else "clear",
        reason = NA_character_,
        tolerance = boundary_tolerance,
        minimum_sd = min(natural_sd),
        sd_at_or_below_tolerance = any(natural_sd <= boundary_tolerance),
        correlation_at_or_near_boundary = !is.na(rho) && abs(rho) >= 0.999,
        at_or_near_boundary = any(natural_sd <= boundary_tolerance) ||
          (!is.na(rho) && abs(rho) >= 0.999)
      )
    ),
    model = model
  )
}

prepare_glmmtmb_data <- function(data, variance) {
  glmm_data <- data
  glmm_data$study_ID <- factor(glmm_data$study_ID)
  glmm_data$effect_size_ID <- factor(glmm_data$effect_size_ID)
  glmm_data$g <- factor(rep("known_sampling", nrow(glmm_data)))
  V <- diag(glmm_data[[variance]])
  effect_size_levels <- levels(glmm_data$effect_size_ID)
  dimnames(V) <- list(effect_size_levels, effect_size_levels)

  if (!identical(rownames(V), effect_size_levels) ||
      !identical(colnames(V), effect_size_levels) ||
      ncol(V) != nlevels(glmm_data$effect_size_ID)) {
    stop("Known-variance matrix is not aligned to effect_size_ID levels.")
  }

  list(
    data = glmm_data,
    V = V,
    variance_alignment = list(
      verified = TRUE,
      V_dimension = dim(V),
      effect_size_levels = effect_size_levels,
      formula_term_order = "equalto(0 + effect_size_ID | g, V)"
    )
  )
}

fit_glmmtmb_location_scale <- function(data, response, variance) {
  prepared <- prepare_glmmtmb_data(data, variance)
  V <- prepared$V
  model_formula <- stats::as.formula(paste(
    response,
    "~ 0 + fertilizer + diag(0 + fertilizer | study_ID) +",
    "equalto(0 + effect_size_ID | g, V)"
  ), env = environment())
  model <- glmmTMB::glmmTMB(
    model_formula,
    dispformula = ~ 0 + fertilizer,
    REML = TRUE,
    data = prepared$data
  )
  list(model = model, variance_alignment = prepared$variance_alignment)
}

extract_glmmtmb_fit <- function(fitted, response) {
  model <- fitted$model
  study_matrix <- glmmTMB::VarCorr(model)$cond$study_ID
  study_sd <- sqrt(diag(study_matrix))
  names(study_sd) <- sub("^fertilizer", "", names(study_sd))
  residual_sd <- exp(glmmTMB::fixef(model)$disp)
  names(residual_sd) <- sub("^fertilizer", "", names(residual_sd))
  formula_text <- paste(deparse(stats::formula(model)), collapse = " ")
  natural_sd <- c(study_sd, residual_sd)
  boundary_tolerance <- 1e-6

  list(
    engine = "glmmTMB",
    response = response,
    status = "fit",
    formula = formula_text,
    location = glmmTMB::fixef(model)$cond,
    study_sd = study_sd,
    residual_sd = residual_sd,
    study_sd_ratio_plant_animal = unname(study_sd["plant"] / study_sd["animal"]),
    residual_sd_ratio_plant_animal = unname(residual_sd["plant"] / residual_sd["animal"]),
    known_variance_alignment = fitted$variance_alignment,
    diagnostics = list(
      natural_sd_scale = TRUE,
      optimizer = list(
        available = TRUE,
        status = if (identical(model$fit$convergence, 0L)) "converged" else "not_converged",
        reason = NA_character_,
        convergence = model$fit$convergence,
        message = model$fit$message,
        iterations = model$fit$iterations,
        evaluations = model$fit$evaluations
      ),
      gradient = list(
        available = !is.null(model$sdr$gradient.fixed) &&
          all(is.finite(model$sdr$gradient.fixed)),
        status = if (!is.null(model$sdr$gradient.fixed) &&
          all(is.finite(model$sdr$gradient.fixed))) "available" else "unavailable",
        reason = if (!is.null(model$sdr$gradient.fixed) &&
          all(is.finite(model$sdr$gradient.fixed))) NA_character_ else {
          "glmmTMB did not retain a finite fixed-parameter gradient"
        },
        maximum_absolute_gradient = if (!is.null(model$sdr$gradient.fixed) &&
          all(is.finite(model$sdr$gradient.fixed))) max(abs(model$sdr$gradient.fixed)) else NA_real_
      ),
      hessian = list(
        available = !is.null(model$sdr$pdHess),
        status = if (isTRUE(model$sdr$pdHess)) "positive_definite" else "not_positive_definite",
        reason = if (is.null(model$sdr$pdHess)) {
          "glmmTMB did not retain a positive-definiteness Hessian diagnostic"
        } else NA_character_,
        pdHess = model$sdr$pdHess
      ),
      boundary = list(
        available = TRUE,
        status = if (any(natural_sd <= boundary_tolerance)) "at_or_near_boundary" else "clear",
        reason = NA_character_,
        tolerance = boundary_tolerance,
        minimum_sd = min(natural_sd),
        at_or_below_tolerance = any(natural_sd <= boundary_tolerance),
        at_or_near_boundary = any(natural_sd <= boundary_tolerance)
      )
    ),
    model = model
  )
}

fit_glmmtmb_primary <- function(data, response, variance) {
  result <- tryCatch(
    {
      fitted <- fit_glmmtmb_location_scale(data, response, variance)
      extract_glmmtmb_fit(fitted, response)
    },
    error = function(error) {
      list(
        engine = "glmmTMB",
        response = response,
        status = "unavailable_or_failed",
        reason = conditionMessage(error)
      )
    }
  )
  result
}

validate_metafor_result <- function(result) {
  expected_location <- c("fertilizeranimal", "fertilizerplant")
  expected_categories <- c("animal", "plant")
  valid <- identical(names(result$location), expected_location) &&
    identical(names(result$study_sd), expected_categories) &&
    identical(names(result$residual_sd), expected_categories) &&
    all(is.finite(result$location)) &&
    all(is.finite(result$study_sd)) &&
    all(is.finite(result$residual_sd)) &&
    all(is.finite(result$fit)) &&
    identical(result$diagnostics$optimizer$convergence, 0L)
  if (!valid) {
    stop("Metafor location-scale validation failed.")
  }
}

validate_glmmtmb_result <- function(result) {
  valid <- identical(result$status, "fit") &&
    all(is.finite(result$study_sd)) &&
    all(is.finite(result$residual_sd)) &&
    isTRUE(result$known_variance_alignment$verified) &&
    isTRUE(result$diagnostics$hessian$pdHess) &&
    identical(as.integer(result$diagnostics$optimizer$convergence), 0L)
  if (!valid) {
    stop("Exact glmmTMB location-scale validation failed.")
  }
}

agreement_tolerance <- 1e-5

validate_primary_agreement <- function(metafor_result, glmmtmb_result, tolerance = agreement_tolerance) {
  validate_glmmtmb_result(glmmtmb_result)
  components <- c("location", "study_sd", "residual_sd")
  maximum_absolute_differences <- vapply(components, function(component) {
    metafor_estimates <- metafor_result[[component]]
    glmmtmb_estimates <- glmmtmb_result[[component]]
    if (!identical(names(metafor_estimates), names(glmmtmb_estimates))) {
      stop("metafor and glmmTMB component names do not align for agreement checking.")
    }
    max(abs(metafor_estimates - glmmtmb_estimates))
  }, numeric(1))
  if (any(maximum_absolute_differences > tolerance)) {
    stop("metafor and glmmTMB primary estimates exceed the agreement tolerance.")
  }
  list(
    available = TRUE,
    status = "within_tolerance",
    tolerance = tolerance,
    maximum_absolute_differences = maximum_absolute_differences
  )
}

diagnostic_status <- function(result) {
  diagnostics <- result$diagnostics
  paste(
    diagnostics$optimizer$status,
    diagnostics$gradient$status,
    diagnostics$hessian$status,
    diagnostics$boundary$status,
    sep = ";"
  )
}

result_summary_rows <- function(result, model_label) {
  exploratory_only <- identical(model_label, "sensitivity_un")
  if (!identical(result$status %||% "fit", "fit")) {
    return(data.frame(
      model = model_label,
      engine = result$engine,
      response = result$response,
      component = "status",
      term = "unavailable_or_failed",
      estimate = NA_real_,
      value = result$reason,
      boundary = NA,
      exploratory_only = exploratory_only,
      diagnostic_status = "unavailable_or_failed",
      stringsAsFactors = FALSE
    ))
  }

  row_metadata <- data.frame(
    model = model_label,
    engine = result$engine,
    response = result$response,
    boundary = result$diagnostics$boundary$at_or_near_boundary,
    exploratory_only = exploratory_only,
    diagnostic_status = diagnostic_status(result),
    stringsAsFactors = FALSE
  )
  make_rows <- function(component, term, estimate) {
    data.frame(
      model = row_metadata$model,
      engine = row_metadata$engine,
      response = row_metadata$response,
      component = component,
      term = term,
      estimate = estimate,
      value = NA_character_,
      boundary = row_metadata$boundary,
      exploratory_only = row_metadata$exploratory_only,
      diagnostic_status = row_metadata$diagnostic_status,
      stringsAsFactors = FALSE
    )
  }

  rows <- rbind(
    make_rows("location", names(result$location), unname(result$location)),
    make_rows("study_sd", names(result$study_sd), unname(result$study_sd)),
    make_rows("residual_sd", names(result$residual_sd), unname(result$residual_sd)),
    make_rows(
      "ratio",
      c("study_plant_animal", "residual_plant_animal"),
      c(result$study_sd_ratio_plant_animal, result$residual_sd_ratio_plant_animal)
    )
  )
  rho <- result$rho %||% NA_real_
  if (!is.na(rho)) {
    rows <- rbind(rows, make_rows("study_rho", "animal_plant", rho))
  }
  rows
}

`%||%` <- function(x, y) if (is.null(x)) y else x

full_frequentist_summary_table <- function(results) {
  rows <- list()
  for (response in names(results$primary$metafor)) {
    rows[[length(rows) + 1L]] <- result_summary_rows(
      results$primary$metafor[[response]], "primary_diag"
    )
    rows[[length(rows) + 1L]] <- result_summary_rows(
      results$primary$glmmTMB[[response]], "primary_glmmTMB"
    )
    if (!is.null(results$sensitivity$metafor[[response]])) {
      rows[[length(rows) + 1L]] <- result_summary_rows(
        results$sensitivity$metafor[[response]], "sensitivity_un"
      )
    }
  }
  do.call(rbind, rows)
}

frequentist_artifact_paths <- function(project_root, smoke) {
  cache_dir <- file.path(
    project_root,
    "Rdata",
    "study_level_heterogeneity",
    "frequentist"
  )
  filename_prefix <- if (smoke) "smoke_frequentist" else "full_frequentist"
  list(
    cache_dir = cache_dir,
    results = file.path(cache_dir, paste0(filename_prefix, "_results.rds")),
    summary = file.path(cache_dir, paste0(filename_prefix, "_summary.csv"))
  )
}

save_full_frequentist_artifacts <- function(prepared, results, project_root, smoke) {
  paths <- frequentist_artifact_paths(project_root, smoke)
  cache_dir <- paths$cache_dir
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(counts = prepared$counts, mixed_studies = prepared$mixed_studies, results = results),
    paths$results
  )
  utils::write.csv(
    full_frequentist_summary_table(results),
    paths$summary,
    row.names = FALSE
  )
}

run_full_frequentist <- function(smoke = FALSE) {
  project_root <- find_project_root()
  prepared <- prepare_study_level_data(project_root)
  response_specs <- if (smoke) {
    list(lnRR = "var.lnRR")
  } else {
    list(lnRR = "var.lnRR", lnCVR = "var.lnCVR")
  }
  primary_metafor <- lapply(names(response_specs), function(response) {
    extract_metafor_fit(
      fit_metafor_location_scale(prepared$data, response, response_specs[[response]]),
      response,
      "DIAG"
    )
  })
  names(primary_metafor) <- names(response_specs)
  lapply(primary_metafor, validate_metafor_result)

  primary_glmmtmb <- lapply(names(response_specs), function(response) {
    fit_glmmtmb_primary(prepared$data, response, response_specs[[response]])
  })
  names(primary_glmmtmb) <- names(response_specs)

  sensitivity_metafor <- if (smoke) {
    list()
  } else {
    fitted <- lapply(names(response_specs), function(response) {
      extract_metafor_fit(
        fit_metafor_location_scale(prepared$data, response, response_specs[[response]], "UN"),
        response,
        "UN"
      )
    })
    names(fitted) <- names(response_specs)
    lapply(fitted, validate_metafor_result)
    fitted
  }

  results <- list(
    primary = list(metafor = primary_metafor, glmmTMB = primary_glmmtmb),
    sensitivity = list(metafor = sensitivity_metafor)
  )
  primary_agreement <- lapply(names(response_specs), function(response) {
    validate_primary_agreement(
      results$primary$metafor[[response]],
      results$primary$glmmTMB[[response]]
    )
  })
  names(primary_agreement) <- names(response_specs)
  results$validation <- list(
    glmmTMB_primary = lapply(results$primary$glmmTMB, validate_glmmtmb_result),
    metafor_glmmTMB_agreement = primary_agreement
  )
  save_full_frequentist_artifacts(prepared, results, project_root, smoke)

  if (smoke) {
    cat("FREQUENTIST SMOKE PASS\n")
  }
  invisible(list(prepared = prepared, results = results))
}

if (sys.nframe() == 0L) {
  run_full_frequentist(smoke = "--smoke" %in% commandArgs(trailingOnly = TRUE))
}
