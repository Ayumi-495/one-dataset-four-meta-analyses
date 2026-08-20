library(tidyverse)
library(metafor)
library(orchaRd)
library(brms)
library(here)

dir.create(here("Rdata"), showWarnings = FALSE)

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
  measure = "VR",
  m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
  m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv,
  data = dat,
  var.names = c("lnVR", "var.lnVR")
)
dat <- escalc(
  measure = "CVR",
  m1i = Mean.Org, sd1i = SD.Org, n1i = N.Org,
  m2i = Mean.conv, sd2i = SD.conv, n2i = N.Conv,
  data = dat,
  var.names = c("lnCVR", "var.lnCVR"),
  correct = TRUE
)

make_meta_model <- function(response, variance, moderator = FALSE) {
  rma.mv(
    yi = dat[[response]],
    V = dat[[variance]],
    mods = if (moderator) ~ Org.fertilizer.type else ~ 1,
    random = list(~ 1 | study_ID, ~ 1 | effect_size_ID),
    data = dat,
    test = "t",
    method = "REML"
  )
}

make_overall_orchard <- function(model, xlab, colour) {
  result <- unclass(mod_results(model, group = "study_ID"))
  result$mod_table$name <- "Overall"
  result$data$moderator <- "Overall"
  class(result) <- "orchard"
  orchard_plot(
    result,
    xlab = xlab,
    group = "study_ID",
    trunk.size = 0.6,
    branch.size = 5,
    angle = 0,
    k.pos = "left"
  ) +
    scale_colour_manual(values = rep("grey20", 8)) +
    scale_fill_manual(values = colour)
}

specs <- list(
  lnRR = c("lnRR", "var.lnRR", "#CC6677"),
  lnCVR = c("lnCVR", "var.lnCVR", "#117733"),
  lnVR = c("lnVR", "var.lnVR", "#88CCEE")
)

for (name in names(specs)) {
  spec <- specs[[name]]
  overall <- make_meta_model(spec[1], spec[2], FALSE)
  meta_regression <- make_meta_model(spec[1], spec[2], TRUE)
  saveRDS(
    make_overall_orchard(overall, spec[1], spec[3]),
    here("Rdata", paste0("orchard_ma_", name, ".rds"))
  )
  saveRDS(
    meta_regression,
    here("Rdata", paste0("mod_", name, "_mr_metafor.rds"))
  )
}

mod_lnVR_mr2 <- rma.mv(
  lnVR,
  var.lnVR,
  mods = ~ Org.fertilizer.type,
  random = list(~ 1 | study_ID,
                ~ Org.fertilizer.type | effect_size_ID),
  data = dat,
  test = "t",
  method = "REML",
  struct = "DIAG"
)
orchard_mr_lnVR <- orchard_plot(
  mod_lnVR_mr2,
  mod = "Org.fertilizer.type",
  xlab = "lnVR",
  group = "study_ID",
  trunk.size = 0.6,
  branch.size = 5,
  angle = 0
) +
  scale_colour_manual(values = rep("grey20", 8)) +
  scale_fill_manual(values = rep("#88CCEE", 2))
saveRDS(orchard_mr_lnVR, here("Rdata", "orchard_mr_lnVR.rds"))

fit_brms_meta <- function(response, variance, moderator, seed) {
  vcv <- diag(dat[[variance]])
  rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID
  rhs <- if (moderator) "1 + Org.fertilizer.type" else "1"
  formula <- bf(as.formula(paste0(
    response, " ~ ", rhs,
    " + (1 | study_ID) + (1 | gr(effect_size_ID, cov = vcv))"
  )))
  prior <- default_prior(
    formula,
    data = dat,
    data2 = list(vcv = vcv),
    family = gaussian()
  )
  prior$prior[if (moderator) 5 else 3] <- "constant(1)"
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

fits <- list(
  fit_ma1 = fit_brms_meta("lnRR", "var.lnRR", FALSE, 20260822),
  fit_mr1 = fit_brms_meta("lnRR", "var.lnRR", TRUE, 20260823),
  fit_ma2 = fit_brms_meta("lnCVR", "var.lnCVR", FALSE, 20260824),
  fit_ma3 = fit_brms_meta("lnVR", "var.lnVR", FALSE, 20260825)
)
for (name in names(fits)) {
  saveRDS(fits[[name]], here("Rdata", paste0(name, ".rds")))
  print(summary(fits[[name]]))
}

sessionInfo()
