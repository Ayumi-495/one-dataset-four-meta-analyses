library(tidyverse)
library(metafor)
library(orchaRd)
library(brms)
library(tidybayes)
library(patchwork)
library(here)

dir.create(here("figures"), showWarnings = FALSE)

orchard_ma_lnRR <- readRDS(here("Rdata", "orchard_ma_lnRR.rds"))
orchard_ma_lnCVR <- readRDS(here("Rdata", "orchard_ma_lnCVR.rds"))
fit_ma1 <- readRDS(here("Rdata", "fit_ma1.rds"))
fit_ma2 <- readRDS(here("Rdata", "fit_ma2.rds"))

# Recalculate lnCVR I2 from the explicitly bias-corrected effect sizes so that
# the plotted annotation cannot silently retain the pre-metafor-5.0 default.
datfull <- read.csv(here("data", "ponisio2014dataset.csv"), header = TRUE)
datfull$study_ID <- as.numeric(as.factor(
  paste(datfull$Author, datfull$Year, datfull$Journal, sep = "_")
))
datfull$effect_size_ID <- seq_len(nrow(datfull))
dat <- datfull %>%
  filter(Org.fertilizer.type %in% c("plant", "animal"),
         Crop.type == "cereals") %>%
  droplevels()
dat <- escalc(
  measure = "CVR",
  m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
  m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv,
  data = dat,
  var.names = c("lnCVR", "var.lnCVR"),
  correct = TRUE
)
mod_lnCVR <- rma.mv(
  lnCVR,
  var.lnCVR,
  random = list(~ 1 | study_ID, ~ 1 | effect_size_ID),
  data = dat,
  test = "t",
  method = "REML"
)
i2_lnCVR <- i2_ml(mod_lnCVR)
stopifnot(
  isTRUE(all.equal(unname(i2_lnCVR["I2_Total"]), 55.863504, tolerance = 1e-6)),
  isTRUE(all.equal(unname(i2_lnCVR["I2_study_ID"]), 9.357345, tolerance = 1e-6)),
  isTRUE(all.equal(unname(i2_lnCVR["I2_effect_size_ID"]), 46.506159,
                   tolerance = 1e-6))
)

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
    mutate(.variable = rename_vars(.variable)) %>%
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
  random_effect_vars <- variables[grep("^sd_study|^sigma", variables)]
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

p1 <- orchard_ma_lnRR +
  annotate(
    geom = "text", x = 1.5, y = -2,
    label = paste0("italic(I)^{2} == ", 98.11, "*\"%\""),
    color = "black", parse = TRUE, size = 4
  )
p2 <- orchard_ma_lnCVR +
  annotate(
    geom = "text", x = 1.5, y = -3,
    label = paste0("italic(I)^{2} == ", round(i2_lnCVR["I2_Total"], 2),
                   "*\"%\""),
    color = "black", parse = TRUE, size = 4
  )

figure2 <- (
  (p1 + p2) /
    (visualize_fixed_effects(fit_ma1) + visualize_fixed_effects(fit_ma2)) /
    (visualize_random_effects(fit_ma1) + visualize_random_effects(fit_ma2)) +
    plot_annotation(tag_levels = "A")
) & theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(here("figures", "Figure2.png"), figure2,
       width = 12, height = 11, units = "in", dpi = 300)
ggsave(here("figures", "Figure2.pdf"), figure2,
       width = 12, height = 11, units = "in")

print(i2_lnCVR)
print(figure2)
sessionInfo()
