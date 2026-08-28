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

required_packages <- c("dplyr", "metafor", "brms", "posterior", "cmdstanr", "drmTMB", "here")
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
    boundary = boundary,
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

brms_artifacts <- function(fit, response, vi) {
  draws <- posterior::as_draws_df(fit)
  study_animal <- draw_column(draws, "^sd_study_ID__fertilizeranimal$")
  study_plant <- draw_column(draws, "^sd_study_ID__fertilizerplant$")
  residual_animal <- draw_column(draws, "^b_sigma_fertilizeranimal$")
  residual_plant <- draw_column(draws, "^b_sigma_fertilizerplant$")
  ratio_table <- rbind(
    ratio_summary(draws[[study_plant]], draws[[study_animal]], response, "study_sd_ratio"),
    # brms models sigma on the log scale; exponentiate before forming an SD ratio.
    ratio_summary(exp(draws[[residual_plant]]), exp(draws[[residual_animal]]), response, "residual_sd_ratio")
  )
  ratio_table$model <- "extended_primary"
  ratio_table$engine <- "brms"
  ratio_table$data_scope <- "full_318_effects_36_studies"
  ratio_table <- ratio_table[, c("response", "data_scope", "engine", "model", "component", "estimate", "lower", "upper")]

  intercept <- draw_column(draws, "^b_fertilizeranimal$")
  plant_mean <- draw_column(draws, "^b_fertilizerplant$")
  set.seed(if (response == "lnRR") 20260828L else 20260829L)
  prediction_rows <- lapply(c("animal", "plant"), function(category) {
    mean_draw <- draws[[if (category == "animal") intercept else plant_mean]]
    study_sd <- draws[[if (category == "animal") study_animal else study_plant]]
    residual_sd <- draws[[if (category == "animal") residual_animal else residual_plant]]
    one <- mean_draw + stats::rnorm(length(mean_draw), 0, sqrt(study_sd^2 + residual_sd^2 + stats::median(vi)))
    ten <- mean_draw + stats::rnorm(length(mean_draw), 0, sqrt(study_sd^2 + residual_sd^2 / 10 + stats::median(vi) / 10))
    make_row <- function(term, values, k) {
      q <- stats::quantile(values, c(0.5, 0.025, 0.975), names = FALSE)
      data.frame(response = response, fertilizer = category, term = term, k = k,
                 median = q[1], q2.5 = q[2], q97.5 = q[3], stringsAsFactors = FALSE)
    }
    rbind(make_row("observed_new_effect", one, 1L), make_row("mean_of_10_new_effects", ten, 10L))
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

run_study_level_pipeline <- function(fit_bayesian = FALSE) {
  prepared <- prepare_study_level_data()
  metafor_fits <- list(
    lnRR_full = fit_metafor_category_model(prepared$full, "lnRR", "var.lnRR"),
    lnCVR_full = fit_metafor_category_model(prepared$full, "lnCVR", "var.lnCVR"),
    lnRR_bridge = fit_metafor_category_model(prepared$matched, "lnRR", "var.lnRR"),
    lnCVR_bridge = fit_metafor_category_model(prepared$matched, "lnCVR", "var.lnCVR")
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
    drmtmb_ratio_table(direct_fits$lnRR, "lnRR"),
    drmtmb_ratio_table(direct_fits$lnCVR, "lnCVR")
  )
  output <- list(counts = prepared$counts, full_data = prepared$full,
                 matched_data = prepared$matched, metafor_fits = metafor_fits,
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
    fits <- list(
      lnRR = fit_brms_category_model(prepared$full, "lnRR", "var.lnRR"),
      lnCVR = fit_brms_category_model(prepared$full, "lnCVR", "var.lnCVR")
    )
    derived <- Map(function(fit, response, variance) brms_artifacts(fit, response, prepared$full[[variance]]),
                   fits, c("lnRR", "lnCVR"), c("var.lnRR", "var.lnCVR"))
    output$brms_fits <- fits
    output$ratios <- do.call(rbind, lapply(derived, `[[`, "ratios"))
    output$predictions <- do.call(rbind, lapply(derived, `[[`, "predictions"))
    output$diagnostics <- do.call(rbind, lapply(derived, `[[`, "diagnostics"))
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
    ratio <- subset(brms_artifacts(fit, "lnCVR", prepared$full$var.lnCVR)$ratios,
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
  } else {
    artifacts <- run_study_level_pipeline(fit_bayesian = "--fit-bayesian" %in% args)
    dir.create(here::here("Rdata"), showWarnings = FALSE, recursive = TRUE)
    saveRDS(artifacts, here::here("Rdata", "study_level_artifacts.rds"))
    message("STUDY LEVEL PIPELINE PASSED")
  }
}
