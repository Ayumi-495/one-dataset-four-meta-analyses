# Reproducible study-level heterogeneity extension for the revision tutorial.
#
# Run from revision/ with:
#   Rscript R/study_level_heterogeneity.R --fit-bayesian
#
# The script writes Rdata/study_level_artifacts.rds. Rendering only reads that
# artifact; it never substitutes a hand-entered table for a fitted result.

# This is the primary model specification recovered from commit 1b759ac.
# It is declared here rather than tuned to match a reported ratio: the same
# prior and sampler settings were used by the historical script that generated
# the manuscript's lnCVR study-SD ratio.
historical_primary_spec <- list(
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  seed_base = 20260827L,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  backend = "cmdstanr",
  priors = c(
    brms::set_prior("normal(0, 1)", class = "b"),
    brms::set_prior("normal(-1, 1)", class = "b", dpar = "sigma"),
    brms::set_prior("exponential(2)", class = "sd", group = "study_ID")
  )
)

# This prespecified sensitivity changes only the prior on the category-specific
# study SDs.  It deliberately is not chosen to reproduce a reported interval.
study_sd_prior_sensitivity_spec <- historical_primary_spec
study_sd_prior_sensitivity_spec$priors <- c(
  brms::set_prior("normal(0, 1)", class = "b"),
  brms::set_prior("normal(-1, 1)", class = "b", dpar = "sigma"),
  brms::set_prior("exponential(1)", class = "sd", group = "study_ID")
)

historical_fit_seed <- function(response) {
  historical_primary_spec$seed_base + 10L +
    if (identical(response, "lnCVR")) 100L else 0L
}

study_level_formula <- function() {
  brms::bf(
    y | se(sqrt(vi), sigma = TRUE) ~ 0 + fertilizer + (0 + fertilizer || study_ID),
    sigma ~ 0 + fertilizer
  )
}

package_versions <- function(packages) {
  stats::setNames(vapply(packages, function(package) {
    as.character(utils::packageVersion(package))
  }, character(1)), packages)
}

required_packages <- c("dplyr", "metafor", "brms", "posterior", "cmdstanr", "glmmTMB", "drmTMB", "here")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required packages before running this script: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

prepare_study_level_data <- function() {
  datfull <- read.csv(here::here("data", "ponisio2014dataset.csv"), header = TRUE)
  datfull$study_ID <- as.numeric(as.factor(paste(
    datfull$Author, datfull$Year, datfull$Journal, sep = "_"
  )))
  datfull$effect_size_ID <- seq_len(nrow(datfull))
  dat <- datfull |>
    dplyr::filter(Org.fertilizer.type %in% c("plant", "animal"),
                  Crop.type == "cereals") |>
    droplevels()
  dat <- metafor::escalc(
    measure = "VR", m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
    m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv, data = dat,
    var.names = c("lnVR", "var.lnVR")
  )
  dat <- metafor::escalc(
    measure = "CVR", m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
    m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv, data = dat,
    var.names = c("lnCVR", "var.lnCVR"), correct = TRUE
  )
  dat$fertilizer <- relevel(factor(dat$Org.fertilizer.type), ref = "animal")
  study_categories <- dat |>
    dplyr::distinct(study_ID, fertilizer) |>
    dplyr::count(study_ID, name = "n_categories")
  matched <- dat |>
    dplyr::inner_join(study_categories, by = "study_ID") |>
    dplyr::filter(n_categories == 1L) |>
    dplyr::select(-n_categories)
  counts <- list(
    full_effects = nrow(dat), full_studies = dplyr::n_distinct(dat$study_ID),
    mixed_studies = sum(study_categories$n_categories == 2L),
    mixed_effects = sum(dat$study_ID %in% study_categories$study_ID[study_categories$n_categories == 2L]),
    matched_effects = nrow(matched), matched_studies = dplyr::n_distinct(matched$study_ID)
  )
  stopifnot(identical(unname(unlist(counts)), c(318L, 36L, 6L, 86L, 232L, 30L)))
  list(full = dat, matched = matched, counts = counts)
}

fit_metafor_category_model <- function(data, response, variance) {
  metafor::rma.mv(
    yi = data[[response]], V = data[[variance]], mods = ~ 0 + fertilizer,
    random = list(~ fertilizer | study_ID, ~ fertilizer | effect_size_ID),
    struct = c("DIAG", "DIAG"), test = "t", method = "REML", data = data
  )
}

metafor_ratio_table <- function(fit, response, scope, model) {
  # `rma.mv` stores the category-specific study variances in tau2 and
  # effect-size-level residual variances in gamma2. The fixed order is asserted
  # so a changed metafor representation cannot silently create a misleading table.
  stopifnot(length(fit$tau2) == 2L, length(fit$gamma2) == 2L)
  data.frame(
    response = response, data_scope = scope, engine = "metafor", model = model,
    component = c("study_sd_ratio", "residual_sd_ratio"),
    estimate = c(sqrt(fit$tau2[2] / fit$tau2[1]), sqrt(fit$gamma2[2] / fit$gamma2[1])),
    diagnostic_status = "not_assessed", boundary = FALSE,
    stringsAsFactors = FALSE
  )
}

fit_glmmtmb_category_model <- function(data, response, variance) {
  model_data <- data
  model_data$study_ID <- factor(model_data$study_ID)
  model_data$effect_size_ID <- factor(model_data$effect_size_ID)
  model_data$g <- factor(rep("known_sampling", nrow(model_data)))
  V <- diag(model_data[[variance]])
  effect_levels <- levels(model_data$effect_size_ID)
  dimnames(V) <- list(effect_levels, effect_levels)
  stopifnot(identical(rownames(V), effect_levels), identical(colnames(V), effect_levels))
  model_formula <- stats::as.formula(paste(
    response, "~ 0 + fertilizer + diag(0 + fertilizer | study_ID) +",
    "equalto(0 + effect_size_ID | g, V)"
  ), env = environment())
  glmmTMB::glmmTMB(
    model_formula, dispformula = ~ 0 + fertilizer, REML = TRUE, data = model_data
  )
}

glmmtmb_ratio_table <- function(fit, response) {
  study_sd <- sqrt(diag(glmmTMB::VarCorr(fit)$cond$study_ID))
  names(study_sd) <- sub("^fertilizer", "", names(study_sd))
  residual_sd <- exp(glmmTMB::fixef(fit)$disp)
  names(residual_sd) <- sub("^fertilizer", "", names(residual_sd))
  stopifnot(all(c("animal", "plant") %in% names(study_sd)),
            all(c("animal", "plant") %in% names(residual_sd)),
            identical(fit$fit$convergence, 0L), isTRUE(fit$sdr$pdHess))
  data.frame(
    response = response, data_scope = "full_318_effects_36_studies",
    engine = "glmmTMB", model = "primary_glmmTMB",
    component = c("study_sd_ratio", "residual_sd_ratio"),
    estimate = c(study_sd["plant"] / study_sd["animal"],
                 residual_sd["plant"] / residual_sd["animal"]),
    diagnostic_status = "converged", boundary = FALSE,
    stringsAsFactors = FALSE
  )
}

fit_drmtmb_direct_model <- function(data, response, variance) {
  model_data <- data.frame(
    response = data[[response]], V = data[[variance]],
    fertilizer = data$fertilizer, study_ID = data$study_ID
  )
  drmTMB::drmTMB(
    drmTMB::bf(
      response ~ fertilizer + (1 | study_ID) + drmTMB::meta_V(V = V),
      sigma ~ 0 + fertilizer,
      sd(study_ID) ~ 0 + fertilizer
    ),
    family = gaussian(), data = model_data, REML = FALSE
  )
}

drmtmb_ratio_table <- function(fit, response) {
  coefficients <- summary(fit)$coefficients
  term <- function(name) coefficients[name, "estimate"]
  plant_study <- term("sd(study_ID):fertilizerplant")
  boundary <- plant_study < log(1e-3)
  data.frame(
    response = response, data_scope = "matched_subset_232_effects_30_studies",
    engine = "drmTMB", model = "drmTMB_direct_study_scale",
    component = c("study_sd_ratio", "residual_sd_ratio"),
    estimate = c(
      exp(plant_study - term("sd(study_ID):fertilizeranimal")),
      exp(term("sigma:fertilizerplant") - term("sigma:fertilizeranimal"))
    ),
    diagnostic_status = if (summary(fit)$convergence == 0L) "ok" else "failed",
    # The weakly identified boundary is a study-scale issue.  Do not attach it
    # to the separate residual-SD contrast merely because both appear in a
    # two-row display table.
    boundary = c(boundary, FALSE),
    stringsAsFactors = FALSE
  )
}

fit_brms_category_model <- function(data, response, variance, seed = historical_fit_seed(response),
                                    specification = historical_primary_spec) {
  model_data <- data.frame(
    y = data[[response]], vi = data[[variance]], fertilizer = data$fertilizer,
    study_ID = data$study_ID
  )
  brms::brm(
    study_level_formula(), prior = specification$priors,
    family = gaussian(), data = model_data,
    chains = specification$chains, cores = specification$chains,
    iter = specification$iter, warmup = specification$warmup, seed = seed,
    backend = specification$backend,
    control = list(adapt_delta = specification$adapt_delta,
                   max_treedepth = specification$max_treedepth), refresh = 1000
  )
}

draw_column <- function(draws, pattern) {
  hits <- grep(pattern, names(draws), value = TRUE)
  if (length(hits) != 1L) {
    stop("Expected exactly one posterior column matching ", pattern,
         "; found: ", paste(hits, collapse = ", "), call. = FALSE)
  }
  hits
}

ratio_summary <- function(numerator, denominator, response, component) {
  ratio <- numerator / denominator
  limits <- stats::quantile(ratio, c(0.5, 0.025, 0.975), names = FALSE)
  data.frame(response = response, component = component, estimate = limits[1],
             lower = limits[2], upper = limits[3], stringsAsFactors = FALSE)
}

brms_artifacts <- function(fit, response, data, variance, model = "extended_primary") {
  draws <- posterior::as_draws_df(fit)
  study_animal <- draw_column(draws, "^sd_study_ID__fertilizeranimal$")
  study_plant <- draw_column(draws, "^sd_study_ID__fertilizerplant$")
  # `b_sigma_*` is on the log-SD scale. Convert it once before using these
  # draws in either ratios or prediction intervals.
  residual_animal <- exp(draws[[draw_column(draws, "^b_sigma_fertilizeranimal$")]])
  residual_plant <- exp(draws[[draw_column(draws, "^b_sigma_fertilizerplant$")]])
  ratio_table <- rbind(
    ratio_summary(draws[[study_plant]], draws[[study_animal]], response, "study_sd_ratio"),
    ratio_summary(residual_plant, residual_animal, response, "residual_sd_ratio")
  )
  ratio_table$model <- model
  ratio_table$engine <- "brms"
  ratio_table$data_scope <- "full_318_effects_36_studies"
  ratio_table <- ratio_table[, c("response", "data_scope", "engine", "model", "component", "estimate", "lower", "upper")]

  intercept <- draw_column(draws, "^b_fertilizeranimal$")
  plant_mean <- draw_column(draws, "^b_fertilizerplant$")
  set.seed(if (response == "lnRR") 20260828L else 20260829L)
  prediction_rows <- lapply(c("animal", "plant"), function(category) {
    mean_draw <- draws[[if (category == "animal") intercept else plant_mean]]
    study_sd <- draws[[if (category == "animal") study_animal else study_plant]]
    residual_sd <- if (category == "animal") residual_animal else residual_plant
    # Sampling variance is a response-by-fertiliser scenario, not one median
    # pooled across categories within a response.
    median_vi <- stats::median(data[[variance]][data$fertilizer == category])
    simulate_target <- function(k = NULL) {
      target_variance <- if (is.null(k)) {
        study_sd^2 + residual_sd^2
      } else {
        study_sd^2 + residual_sd^2 / k + median_vi / k
      }
      mean_draw + stats::rnorm(length(mean_draw), 0, sqrt(target_variance))
    }
    make_row <- function(term, values, k = NA_integer_) {
      q <- stats::quantile(values, c(0.5, 0.025, 0.975), names = FALSE)
      data.frame(response = response, fertilizer = category, term = term,
                 vi_scenario = "median_response_by_fertilizer", vi = median_vi, k = k,
                 median = q[1], q2.5 = q[2], q97.5 = q[3], stringsAsFactors = FALSE)
    }
    rbind(
      make_row("latent_new_effect", simulate_target()),
      make_row("observed_new_effect", simulate_target(1L), 1L),
      make_row("mean_of_1_new_effects", simulate_target(1L), 1L),
      make_row("mean_of_5_new_effects", simulate_target(5L), 5L),
      make_row("mean_of_10_new_effects", simulate_target(10L), 10L)
    )
  })
  diagnostics <- posterior::summarise_draws(
    posterior::as_draws_array(fit), posterior::rhat, posterior::ess_bulk
  )
  rhat_column <- grep("rhat$", names(diagnostics), value = TRUE)
  ess_column <- grep("ess_bulk$", names(diagnostics), value = TRUE)
  stopifnot(length(rhat_column) == 1L, length(ess_column) == 1L)
  sampler <- posterior::as_draws_df(brms::nuts_params(fit))
  list(
    ratios = ratio_table,
    predictions = do.call(rbind, prediction_rows),
    diagnostics = data.frame(
      response = response,
      max_rhat = max(diagnostics[[rhat_column]], na.rm = TRUE),
      min_bulk_ess = min(diagnostics[[ess_column]], na.rm = TRUE),
      divergences = sum(sampler$Parameter == "divergent__" & sampler$Value == 1),
      stringsAsFactors = FALSE
    ),
    draws = draws
  )
}

# Optional direct Bayesian study-scale sensitivity route. It is separate from
# the primary artifact: blsmeta uses JAGS and must clear its own post-warmup
# diagnostics before its estimates are interpreted.
prepare_blsmeta_subset <- function() {
  data <- prepare_study_level_data()$matched
  data <- data[order(data$study_ID, data$effect_size_ID), , drop = FALSE]
  data$es_blsmeta <- seq_len(nrow(data))
  data$study_blsmeta <- blsmeta::make_study_id(data$study_ID)
  first_rows <- !duplicated(data$study_blsmeta)
  stopifnot(
    nrow(data) == 232L,
    length(unique(data$study_blsmeta)) == 30L,
    identical(as.integer(data$study_blsmeta[first_rows]), seq_len(30L)),
    all(tapply(data$fertilizer, data$study_blsmeta,
               function(x) length(unique(x))) == 1L)
  )
  data
}

fit_blsmeta_direct_model <- function(data, response, variance, iter, warmup, seed) {
  fit_data <- data
  fit_data$y_blsmeta <- fit_data[[response]]
  fit_data$v_blsmeta <- fit_data[[variance]]
  set.seed(seed)
  fit <- blsmeta::blsmeta(
    yi = y_blsmeta, vi = v_blsmeta,
    es_id = es_blsmeta, study_id = study_blsmeta,
    mods = ~ fertilizer,
    mods_scale2 = ~ 0 + fertilizer,
    mods_scale3 = ~ 0 + fertilizer,
    iter = iter, warmup = warmup, chains = 4, data = fit_data
  )
  attr(fit, "blsmeta_metadata") <- list(
    response = response, iter = iter, warmup = warmup, seed = seed,
    package_version = as.character(utils::packageVersion("blsmeta"))
  )
  fit
}

blsmeta_post_warmup <- function(fit) {
  first_iteration <- stats::start(fit$posterior_samples)[[1L]]
  stats::window(fit$posterior_samples,
                start = first_iteration + attr(fit, "blsmeta_metadata")$warmup)
}

summarise_blsmeta_sensitivity <- function(fit) {
  metadata <- attr(fit, "blsmeta_metadata")
  samples <- blsmeta_post_warmup(fit)
  draws <- do.call(rbind, lapply(samples, as.matrix))
  residual_ratio <- exp(draws[, "gamma[2]"] - draws[, "gamma[1]"])
  study_ratio <- exp(draws[, "eta[2]"] - draws[, "eta[1]"])
  ratio_chains <- coda::mcmc.list(lapply(samples, function(chain) {
    chain <- as.matrix(chain)
    coda::mcmc(cbind(
      residual_sd_ratio = exp(chain[, "gamma[2]"] - chain[, "gamma[1]"]),
      study_sd_ratio = exp(chain[, "eta[2]"] - chain[, "eta[1]"])
    ))
  }))
  rhat <- coda::gelman.diag(ratio_chains, multivariate = FALSE)$psrf[, 1L]
  ess <- coda::effectiveSize(ratio_chains)
  status <- if (max(rhat) <= 1.01 && min(ess) >= 400) "passed" else "failed_convergence"
  make_row <- function(component, values, index) {
    limits <- stats::quantile(values, c(0.5, 0.025, 0.975), names = FALSE)
    data.frame(
      response = metadata$response, engine = "blsmeta",
      model = "blsmeta_direct_study_scale", component = component,
      estimate = limits[1], lower = limits[2], upper = limits[3],
      ratio_rhat = rhat[index], ratio_ess = ess[index],
      diagnostic_status = status, package_version = metadata$package_version,
      stringsAsFactors = FALSE
    )
  }
  rbind(
    make_row("residual_sd_ratio", residual_ratio, "residual_sd_ratio"),
    make_row("study_sd_ratio", study_ratio, "study_sd_ratio")
  )
}

run_blsmeta_sensitivity <- function() {
  if (!requireNamespace("blsmeta", quietly = TRUE) || !requireNamespace("coda", quietly = TRUE)) {
    stop("--fit-blsmeta requires blsmeta and coda. See revision/README.md.", call. = FALSE)
  }
  data <- prepare_blsmeta_subset()
  settings <- list(
    lnRR = list(variance = "var.lnRR", iter = 20000L, warmup = 5000L, seed = 20260830L),
    lnCVR = list(variance = "var.lnCVR", iter = 5000L, warmup = 1000L, seed = 20260828L)
  )
  fits <- lapply(names(settings), function(response) {
    configuration <- settings[[response]]
    fit_blsmeta_direct_model(
      data, response, configuration$variance, configuration$iter,
      configuration$warmup, configuration$seed
    )
  })
  names(fits) <- names(settings)
  list(
    fits = fits,
    summary = do.call(rbind, lapply(fits, summarise_blsmeta_sensitivity)),
    settings = settings
  )
}

run_study_level_pipeline <- function(fit_bayesian = FALSE) {
  prepared <- prepare_study_level_data()
  metafor_fits <- list(
    lnRR_full = fit_metafor_category_model(prepared$full, "lnRR", "var.lnRR"),
    lnCVR_full = fit_metafor_category_model(prepared$full, "lnCVR", "var.lnCVR"),
    lnRR_bridge = fit_metafor_category_model(prepared$matched, "lnRR", "var.lnRR"),
    lnCVR_bridge = fit_metafor_category_model(prepared$matched, "lnCVR", "var.lnCVR")
  )
  glmmtmb_fits <- list(
    lnRR = fit_glmmtmb_category_model(prepared$full, "lnRR", "var.lnRR"),
    lnCVR = fit_glmmtmb_category_model(prepared$full, "lnCVR", "var.lnCVR")
  )
  direct_fits <- list(
    lnRR = fit_drmtmb_direct_model(prepared$matched, "lnRR", "var.lnRR"),
    lnCVR = fit_drmtmb_direct_model(prepared$matched, "lnCVR", "var.lnCVR")
  )
  frequentist <- rbind(
    metafor_ratio_table(metafor_fits$lnRR_full, "lnRR", "full_318_effects_36_studies", "primary_diag"),
    metafor_ratio_table(metafor_fits$lnCVR_full, "lnCVR", "full_318_effects_36_studies", "primary_diag"),
    metafor_ratio_table(metafor_fits$lnRR_bridge, "lnRR", "matched_subset_232_effects_30_studies", "metafor_DIAG_bridge"),
    metafor_ratio_table(metafor_fits$lnCVR_bridge, "lnCVR", "matched_subset_232_effects_30_studies", "metafor_DIAG_bridge"),
    glmmtmb_ratio_table(glmmtmb_fits$lnRR, "lnRR"),
    glmmtmb_ratio_table(glmmtmb_fits$lnCVR, "lnCVR"),
    drmtmb_ratio_table(direct_fits$lnRR, "lnRR"),
    drmtmb_ratio_table(direct_fits$lnCVR, "lnCVR")
  )
  output <- list(counts = prepared$counts, full_data = prepared$full,
                 matched_data = prepared$matched, metafor_fits = metafor_fits,
                 glmmtmb_fits = glmmtmb_fits,
                 direct_fits = direct_fits,
                 frequentist = frequentist, generated_at = Sys.time(),
                 software = utils::sessionInfo(),
                 model_provenance = list(
                   primary_specification = "recovered_commit_1b759ac",
                   formula = paste(deparse(study_level_formula()), collapse = "\n"),
                   priors = as.data.frame(historical_primary_spec$priors),
                   mcmc = historical_primary_spec[names(historical_primary_spec) != "priors"],
                   package_versions = c(
                     package_versions(c("metafor", "brms", "posterior", "cmdstanr")),
                     CmdStan = as.character(cmdstanr::cmdstan_version())
                   )
                 ))
  if (fit_bayesian) {
    fit_plan <- list(
      lnRR_primary = list(response = "lnRR", variance = "var.lnRR",
                          specification = historical_primary_spec,
                          model = "extended_primary", seed_offset = 0L),
      lnCVR_primary = list(response = "lnCVR", variance = "var.lnCVR",
                           specification = historical_primary_spec,
                           model = "extended_primary", seed_offset = 0L),
      lnRR_sensitivity = list(response = "lnRR", variance = "var.lnRR",
                              specification = study_sd_prior_sensitivity_spec,
                              model = "extended_sd_prior_sensitivity", seed_offset = 1L),
      lnCVR_sensitivity = list(response = "lnCVR", variance = "var.lnCVR",
                               specification = study_sd_prior_sensitivity_spec,
                               model = "extended_sd_prior_sensitivity", seed_offset = 1L)
    )
    fits <- lapply(fit_plan, function(plan) {
      fit_brms_category_model(
        prepared$full, plan$response, plan$variance,
        seed = historical_fit_seed(plan$response) + plan$seed_offset,
        specification = plan$specification
      )
    })
    derived <- Map(function(fit, plan) {
      brms_artifacts(fit, plan$response, prepared$full, plan$variance, plan$model)
    }, fits, fit_plan)
    output$brms_fits <- fits
    output$ratios <- do.call(rbind, lapply(derived, `[[`, "ratios"))
    output$predictions <- do.call(rbind, lapply(derived[names(fit_plan) %in% c("lnRR_primary", "lnCVR_primary")], `[[`, "predictions"))
    output$diagnostics <- do.call(rbind, lapply(derived, `[[`, "diagnostics"))
    output$diagnostics$model <- vapply(fit_plan, `[[`, character(1), "model")
  }
  output
}

run_lncvr_monte_carlo_audit <- function(n_replicates = 5L) {
  prepared <- prepare_study_level_data()
  # These consecutive seeds were fixed before fitting and are solely used to
  # quantify Monte Carlo variation under the recovered historical model.
  seeds <- historical_fit_seed("lnCVR") + seq_len(n_replicates)
  runs <- lapply(seeds, function(seed) {
    fit <- fit_brms_category_model(prepared$full, "lnCVR", "var.lnCVR", seed = seed)
    ratio <- subset(brms_artifacts(fit, "lnCVR", prepared$full, "var.lnCVR")$ratios,
                    component == "study_sd_ratio")
    cbind(seed = seed, ratio)
  })
  do.call(rbind, runs)
}

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L) {
  if ("--mc-audit" %in% args) {
    dir.create(here::here("Rdata"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(run_lncvr_monte_carlo_audit(),
            here::here("Rdata", "lncvr_monte_carlo_audit.rds"))
    message("LNCVR MONTE CARLO AUDIT PASSED")
  } else if ("--fit-blsmeta" %in% args) {
    dir.create(here::here("Rdata"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(run_blsmeta_sensitivity(),
            here::here("Rdata", "blsmeta_sensitivity.rds"))
    message("BLSMETA SENSITIVITY COMPLETED")
  } else {
    artifacts <- run_study_level_pipeline(fit_bayesian = "--fit-bayesian" %in% args)
    dir.create(here::here("Rdata"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(artifacts, here::here("Rdata", "study_level_artifacts.rds"))
    message("STUDY LEVEL PIPELINE PASSED")
  }
}
