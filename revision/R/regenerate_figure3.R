library(tidyverse)
library(metafor)
library(orchaRd)
library(brms)
library(tidybayes)
library(patchwork)
library(here)

dir.create(here("Rdata"), showWarnings = FALSE)
dir.create(here("figures"), showWarnings = FALSE)

datfull <- read.csv(here("data", "ponisio2014dataset.csv"), header = TRUE)
datfull$study_ID <- as.numeric(as.factor(
  paste(datfull$Author, datfull$Year, datfull$Journal, sep = "_")
))
datfull$effect_size_ID <- seq_len(nrow(datfull))

dat <- datfull %>%
  filter(Org.fertilizer.type %in% c("plant", "animal"),
         Crop.type == "cereals") %>%
  droplevels()

stopifnot(identical(colnames(model.matrix(~ Org.fertilizer.type, dat))[2],
                    "Org.fertilizer.typeplant"))

dat <- escalc(
  measure = "CVR",
  m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
  m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv,
  data = dat,
  var.names = c("lnCVR", "var.lnCVR"),
  correct = TRUE
)

fit_metafor <- function(response, variance) {
  rma.mv(
    yi = dat[[response]],
    V = dat[[variance]],
    mods = ~ Org.fertilizer.type,
    random = list(~ 1 | study_ID,
                  ~ Org.fertilizer.type | effect_size_ID),
    data = dat,
    test = "t",
    method = "REML",
    struct = "DIAG"
  )
}

make_orchard <- function(model, xlab, colour) {
  orchard_plot(
    model,
    mod = "Org.fertilizer.type",
    xlab = xlab,
    group = "study_ID",
    trunk.size = 0.6,
    branch.size = 5,
    angle = 0,
    k.pos = "left",
    # mod.order is bottom-to-top when flip = TRUE: Plant below, Animal above.
    # Reordering inside orchard_plot keeps estimates, raw data, intervals,
    # category labels, and k/study annotations aligned.
    mod.order = c("Plant", "Animal")
  ) +
    scale_colour_manual(values = rep("grey20", 8)) +
    scale_fill_manual(values = rep(colour, 2))
}

mod_lnRR_mr2 <- fit_metafor("lnRR", "var.lnRR")
mod_lnCVR_mr2 <- fit_metafor("lnCVR", "var.lnCVR")
orchard_mr_lnRR <- make_orchard(mod_lnRR_mr2, "lnRR", "#CC6677")
orchard_mr_lnCVR <- make_orchard(mod_lnCVR_mr2, "lnCVR", "#117733")

group_counts <- dat %>%
  group_by(Org.fertilizer.type) %>%
  summarise(k = n(), studies = n_distinct(study_ID), .groups = "drop")
stopifnot(
  identical(group_counts$k[group_counts$Org.fertilizer.type == "animal"], 134L),
  identical(group_counts$studies[group_counts$Org.fertilizer.type == "animal"], 26L),
  identical(group_counts$k[group_counts$Org.fertilizer.type == "plant"], 184L),
  identical(group_counts$studies[group_counts$Org.fertilizer.type == "plant"], 16L)
)

validate_orchard_annotations <- function(plot) {
  built_layers <- ggplot_build(plot)$data
  annotation_layer <- built_layers[
    vapply(built_layers, function(layer) "label" %in% names(layer), logical(1))
  ][[1]]
  annotation_layer <- annotation_layer[order(annotation_layer$x), ]
  stopifnot(
    grepl("184.*16", annotation_layer$label[1]), # Plant, bottom
    grepl("134.*26", annotation_layer$label[2])  # Animal, top
  )
}
validate_orchard_annotations(orchard_mr_lnRR)
validate_orchard_annotations(orchard_mr_lnCVR)

saveRDS(orchard_mr_lnRR, here("Rdata", "orchard_mr_lnRR.rds"))
saveRDS(orchard_mr_lnCVR, here("Rdata", "orchard_mr_lnCVR.rds"))

fit_brms_location_scale <- function(response, variance, seed) {
  vcv <- diag(dat[[variance]])
  rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID
  formula <- bf(
    as.formula(paste0(
      response,
      " ~ 1 + Org.fertilizer.type + (1 | study_ID) + ",
      "(1 | gr(effect_size_ID, cov = vcv))"
    )),
    sigma ~ 1 + Org.fertilizer.type
  )
  prior <- default_prior(
    formula,
    data = dat,
    data2 = list(vcv = vcv),
    family = gaussian()
  )
  prior$prior[5] <- "constant(1)"
  brm(
    formula = formula,
    data = dat,
    data2 = list(vcv = vcv),
    chains = 2,
    cores = 2,
    iter = 6000,
    warmup = 3000,
    prior = prior,
    seed = seed,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )
}

fit_ls_ma1 <- fit_brms_location_scale("lnRR", "var.lnRR", 20260820)
fit_ls_ma2 <- fit_brms_location_scale("lnCVR", "var.lnCVR", 20260821)
saveRDS(fit_ls_ma1, here("Rdata", "fit_ls_ma1.rds"))
saveRDS(fit_ls_ma2, here("Rdata", "fit_ls_ma2.rds"))

rename_vars <- function(variable) {
  variable <- gsub("b_Intercept", "b_l_intercept", variable)
  variable <- gsub("b_sigma_Intercept", "b_s_intercept", variable)
  variable <- gsub("b_Org.fertilizer.typeplant", "b_l_contrast", variable)
  variable <- gsub("b_sigma_Org.fertilizer.typeplant", "b_s_contrast", variable)
  variable <- gsub("sd_study_ID__Intercept", "sd_study_ID", variable)
  gsub("sigma", "sd_effect_ID", variable)
}

visualize_fixed_effects <- function(model) {
  variables <- get_variables(model)
  fixed_effect_vars <- variables[grep("^b_", variables)]
  model %>%
    spread_draws(!!!syms(fixed_effect_vars)) %>%
    pivot_longer(
      cols = all_of(fixed_effect_vars),
      names_to = ".variable",
      values_to = ".value"
    ) %>%
    mutate(
      .variable = factor(
        rename_vars(.variable),
        levels = rev(c("b_l_intercept", "b_l_contrast",
                       "b_s_intercept", "b_s_contrast"))
      )
    ) %>%
    ggplot(aes(x = .value, y = .variable)) +
    stat_halfeye(
      normalize = "xy",
      point_interval = "mean_qi",
      fill = "lightcyan3",
      color = "lightcyan4"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
    labs(y = "Fixed effects", x = "Posterior values") +
    theme_classic()
}

visualize_random_effects <- function(model) {
  variables <- get_variables(model)
  random_effect_vars <- variables[grep("^sd_study", variables)]
  model %>%
    spread_draws(!!!syms(random_effect_vars)) %>%
    pivot_longer(
      cols = all_of(random_effect_vars),
      names_to = ".variable",
      values_to = ".value"
    ) %>%
    mutate(.variable = rename_vars(.variable)) %>%
    ggplot(aes(x = .value, y = .variable)) +
    stat_halfeye(
      normalize = "xy",
      point_interval = "mean_qi",
      fill = "olivedrab3",
      color = "olivedrab4"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
    labs(y = "Random effects (SD)", x = "Posterior values") +
    theme_classic()
}

figure3 <- (
  (orchard_mr_lnRR + orchard_mr_lnCVR) /
    (visualize_fixed_effects(fit_ls_ma1) + visualize_fixed_effects(fit_ls_ma2)) /
    (visualize_random_effects(fit_ls_ma1) + visualize_random_effects(fit_ls_ma2)) +
    plot_layout(heights = c(1, 1.5, 0.5)) +
    plot_annotation(tag_levels = "A")
) & theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(here("figures", "Figure3.png"), figure3,
       width = 12, height = 11, units = "in", dpi = 300)
ggsave(here("figures", "Figure3.pdf"), figure3,
       width = 12, height = 11, units = "in")

print(figure3)
print(summary(fit_ls_ma1))
print(summary(fit_ls_ma2))
sessionInfo()
