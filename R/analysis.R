# relevant webpages for meta-analysis code:

#  https://itchyshin.github.io/Meta-analysis_tutorial/
# https://itchyshin.github.io/location-scale_meta-analysis/#figure-code


# packages

library(tidybayes)
library(tidyverse)
library(orchaRd)
library(metafor)
#library(glmmTMB)
library(here)
library(brms)
# build package from branch: 
# required to do meta-analysis with glmmTMB
#remove.packages("glmmTMB")
#remotes::install_github("coraliewilliams/glmmTMB", ref="equalto_covstruc", subdir="glmmTMB", force = TRUE) 
library(glmmTMB)
glmmTMB:::.valid_covstruct


datfull <- read.csv(here::here("data", "ponisio2014dataset.csv"), header = TRUE)

dim(datfull)
head(datfull)

# create ID using Author + Year + Journal - gives unique identifier for each study

datfull$study_ID <- as.numeric(as.factor(paste(datfull$Author, datfull$Year, datfull$Journal, sep = "_")))
datfull$effect_size_ID <- 1:nrow(datfull)

# reduce data - Org.fertilizer.type = either plant or animal 

dat <- datfull %>%
  filter(Org.fertilizer.type %in% c("plant", "animal")) %>%
  droplevels()

# reduce data - and also Crop.type = cereals

dat <- dat %>%
  filter(Crop.type == "cereals") %>%
  droplevels()

dim(dat) # 611 effect sizes
head(datfull)


# getting effect size - lnVR and lnCVR using escalc from metafor

dat <- escalc(measure = "VR",
                      m1i = Mean.Org,
                      sd1i = SD.Org,
                      n1i = N.Org,
                      m2i = Mean.conv,
                      sd2i = SD.conv,
                      n2i = N.Conv,
                      data = dat,
                      var.names = c("lnVR", "var.lnVR"))

dat <- escalc(measure = "CVR",
              m1i = Mean.Org,
              sd1i = SD.Org,
              n1i = N.Org,
              m2i = Mean.conv,
              sd2i = SD.conv,
              n2i = N.Conv,
              data = dat,
              var.names = c("lnCVR", "var.lnCVR"))

############################
# meta-analysis with metafor
############################

# lnRR
mod_lnRR <- rma.mv(lnRR,
                    var.lnRR,
                    random = list(~ 1 | study_ID,
                                  ~ 1 | effect_size_ID),
                    data = dat,
                    test = "t",
                    method = "REML")

summary(mod_lnRR)
i2_ml(mod_lnRR)

orchard_plot(mod_lnRR,  xlab = "lnRR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = "#CC6677") 

# adding overall to orchard plot
main <- mod_results(mod_lnRR, group = "study_ID")
main <- unclass(main)
main$mod_table$name <- "Overall" 
main$data$moderator <- "Overall"
class(main) <- "orchard"

orchard_ma_lnRR <-orchard_plot(main,  xlab = "lnRR", 
                               group = "study_ID",
                               trunk.size = 0.6,
                               branch.size = 5,
                               angle = 0,
                               k.pos = "left") +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = "#CC6677") 

saveRDS(orchard_ma_lnRR, here::here("Rdata", "orchard_ma_lnRR.rds"))


#lnCVR
mod_lnCVR <- rma.mv(lnCVR,
                   var.lnCVR,
                   random = list(~ 1 | study_ID,
                                 ~ 1 | effect_size_ID),
                   data = dat,
                   test = "t",
                   method = "REML")
summary(mod_lnCVR)
i2_ml(mod_lnCVR)

orchard_plot(mod_lnCVR,  xlab = "lnCVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = "#117733")


# adding overall to orchard plot
main <- mod_results(mod_lnCVR, group = "study_ID")
main <- unclass(main)
main$mod_table$name <- "Overall" 
main$data$moderator <- "Overall"
class(main) <- "orchard"

orchard_ma_lnCVR <-orchard_plot(main,  xlab = "lnCVR", 
                               group = "study_ID",
                               trunk.size = 0.6,
                               branch.size = 5,
                               angle = 0) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = "#117733")

saveRDS(orchard_ma_lnCVR, here::here("Rdata", "orchard_ma_lnCVR.rds"))

#lnVR
mod_lnVR <- rma.mv(lnVR,
                   var.lnVR,
                   random = list(~ 1 | study_ID,
                                 ~ 1 | effect_size_ID),
                   data = dat,
                   test = "t",
                   method = "REML")
summary(mod_lnVR)
i2_ml(mod_lnVR)

orchard_plot(mod_lnVR,  xlab = "lnVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = "#88CCEE")


# meta-regression using Org.fertilizer.type
# lnRR

mod_lnRR_mr <- rma.mv(lnRR,
                    var.lnRR,
                    mods = ~ Org.fertilizer.type,
                    random = list(~ 1 | study_ID,
                                  ~ 1 | effect_size_ID),
                    data = dat,
                    test = "t",
                    method = "REML")
summary(mod_lnRR_mr)
r2_ml(mod_lnRR_mr)

orchard_plot(mod_lnRR_mr, 
             mod = "Org.fertilizer.type",
             xlab = "lnRR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#CC6677", "#CC6677"))

# lnCVR
mod_lnCVR_mr <- rma.mv(lnCVR,
                    var.lnCVR,
                    mods = ~ Org.fertilizer.type,
                    random = list(~ 1 | study_ID,
                                  ~ 1 | effect_size_ID),
                    data = dat,
                    test = "t",
                    method = "REML")
summary(mod_lnCVR_mr)
r2_ml(mod_lnCVR_mr)
orchard_plot(mod_lnCVR_mr, 
             mod = "Org.fertilizer.type",
             xlab = "lnCVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#117733", "#117733"))

# lnVR
mod_lnVR_mr <- rma.mv(lnVR,
                      var.lnVR,
                      mods = ~ Org.fertilizer.type,
                      random = list(~ 1 | study_ID,
                                    ~ 1 | effect_size_ID),
                      data = dat,
                      test = "t",
                      method = "REML")
summary(mod_lnVR_mr)
r2_ml(mod_lnVR_mr)

orchard_plot(mod_lnVR_mr, 
             mod = "Org.fertilizer.type",
             xlab = "lnVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#88CCEE", "#88CCEE"))


# heteroscad model
# lnRR

mod_lnRR_mr2 <- rma.mv(lnRR,
                       var.lnRR,
                       mods = ~ Org.fertilizer.type,
                       random = list(~ 1 | study_ID,
                                     ~ Org.fertilizer.type | effect_size_ID),
                       data = dat,
                       test = "t",
                       method = "REML",
                       struct =  'DIAG')

summary(mod_lnRR_mr2)

orchard_plot(mod_lnRR_mr2, 
             mod = "Org.fertilizer.type",
             xlab = "lnRR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#CC6677", "#CC6677"))


orchard_mr_lnRR <- orchard_plot(mod_lnRR_mr2, 
                                mod = "Org.fertilizer.type",
                                xlab = "lnRR", 
                                group = "study_ID",
                                trunk.size = 0.6,
                                branch.size = 5,
                                angle = 0,
                                k.pos = "left") +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#CC6677", "#CC6677"))

saveRDS(orchard_mr_lnRR, here::here("Rdata", "orchard_mr_lnRR.rds"))


# lnCVR
mod_lnCVR_mr2 <- rma.mv(lnCVR,
                       var.lnCVR,
                       mods = ~ Org.fertilizer.type,
                       random = list(~ 1 | study_ID,
                                     ~ Org.fertilizer.type | effect_size_ID),
                       data = dat,
                       test = "t",
                       method = "REML",
                       struct =  'DIAG')
summary(mod_lnCVR_mr2)
orchard_plot(mod_lnCVR_mr2, 
             mod = "Org.fertilizer.type",
             xlab = "lnCVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#117733", "#117733"))

orchard_mr_lnCVR <- orchard_plot(mod_lnCVR_mr2, 
                                mod = "Org.fertilizer.type",
                                xlab = "lnCVR", 
                                group = "study_ID",
                                trunk.size = 0.6,
                                branch.size = 5,
                                angle = 0) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#117733", "#117733"))

saveRDS(orchard_mr_lnCVR, here::here("Rdata", "orchard_mr_lnCVR.rds"))

# lnVR
mod_lnVR_mr2 <- rma.mv(lnVR,
                       var.lnVR,
                       mods = ~ Org.fertilizer.type,
                       random = list(~ 1 | study_ID,
                                     ~ Org.fertilizer.type | effect_size_ID),
                       data = dat,
                       test = "t",
                       method = "REML",
                       struct =  'DIAG')
summary(mod_lnVR_mr2)

orchard_plot(mod_lnVR_mr2, 
             mod = "Org.fertilizer.type",
             xlab = "lnVR", 
             group = "study_ID",
             trunk.size = 0.6,
             branch.size = 5) +
  scale_colour_manual(values = rep("grey20",8)) +
  scale_fill_manual(values = c("#88CCEE", "#88CCEE"))


########
# location-scale model with metafor - although it cannot take no random effects for variance part

# lnRR
mod_lnRR_ls <- rma(yi = lnRR,
                   vi = var.lnRR,
                       mods = ~ Org.fertilizer.type,
                       scale = ~ Org.fertilizer.type,
                       data = dat,
                       test = "t",
                       method = "REML")

# note scale is variance not SD
summary(mod_lnRR_ls)

# lnCVR
mod_lnCVR_ls <- rma(yi = lnCVR,
                   vi = var.lnCVR,
                   mods = ~ Org.fertilizer.type,
                   scale = ~ Org.fertilizer.type,
                   data = dat,
                   test = "t",
                   method = "REML")
summary(mod_lnCVR_ls)

# lnVR
mod_lnVR_ls <- rma(yi = lnVR,
                   vi = var.lnVR,
                   mods = ~ Org.fertilizer.type,
                   scale = ~ Org.fertilizer.type,
                   data = dat,
                   test = "t",
                   method = "REML")

summary(mod_lnVR_ls)

#############################
# meta-analysis with brms 
#############################

#' @Shinichi - the following code block is duplicated later in the script.
#' Please remove this block if not needed.

# lnRR
######

# vcv <- diag(dat$var.lnRR)
# 
# rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID
# 
# form_ma1 <- bf(lnRR ~ 1 +
#                  (1|study_ID) + # this is u_l (the between-study effect)
#                  (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
# )
# 
# # generate default priors
# prior_ma1 <- default_prior(form_ma1, 
#                            data=dat, 
#                            data2=list(vcv=vcv),
#                            family=gaussian())
# prior_ma1$prior[3] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
# prior_ma1
# 
# # fitting model
# fit_ma1 <- brm(
#   formula = form_ma1,
#   data = dat,
#   data2 = list(vcv=vcv),
#   chains = 2,
#   cores = 2,
#   iter = 6000,
#   warmup = 3000,
#   prior = prior_ma1,
#   control = list(adapt_delta=0.95, max_treedepth=15)
# )
# 
# summary(fit_ma1)
# 
# # save this as rds
# saveRDS(fit_ma1, here::here("Rdata", "fit_ma1.rds"))
# 
# # meta-regression
# 
# form_mr1 <- bf(lnRR ~ 1 + Org.fertilizer.type +
#                   (1|study_ID) + # this is u_l (the between-study effect)
#                   (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
# )
# 
# # generate default priors
# prior_mr1 <- default_prior(form_mr1,
#                            data=dat, 
#                            data2=list(vcv=vcv),
#                            family=gaussian())
# prior_mr1$prior[5] = "constant(1)" # meta-analysis assumes
# 
# 
# # fitting model
# fit_mr1 <- brm(
#   formula = form_mr1,
#   data = dat,
#   data2 = list(vcv=vcv),
#   chains = 2,
#   cores = 2,
#   iter = 6000,
#   warmup = 3000,
#   prior = prior_mr1,
#   control = list(adapt_delta=0.95, max_treedepth=15)
# )
# summary(fit_mr1)
# 
# # save this as rds
# saveRDS(fit_mr1, here::here("Rdata", "fit_mr1.rds"))
# 
# # location-scale meta-regression: Model 1
# 
# form_ls_ma1 <- bf(lnRR ~ 1 + Org.fertilizer.type +
#                     (1|study_ID) + # this is u_l (the between-study effect)
#                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
#                   sigma ~ 1 +  Org.fertilizer.type
# )
# 
# 
# # generate default priors
# prior_ls_ma1 <- default_prior(form_ls_ma1, 
#                               data=dat, 
#                               data2=list(vcv=vcv),
#                               family=gaussian())
# 
# prior_ls_ma1$prior[5] = "constant(1)" # meta-analysis assumes
# 
# # fitting model
# fit_ls_ma1 <- brm(
#   formula = form_ls_ma1,
#   data = dat,
#   data2 = list(vcv=vcv),
#   chains = 2,
#   cores = 2,
#   iter = 6000,
#   warmup = 3000,
#   prior = prior_ls_ma1,
#   control = list(adapt_delta=0.95, max_treedepth=15)
# )
# summary(fit_ls_ma1)
# 
# # save this as rds
# saveRDS(fit_ls_ma1, here::here("Rdata", "fit_ls_ma1.rds"))
# 
# # location-scale meta-regression: Model 2
# 
# form_ls_ma1b <- bf(lnRR ~ 1 + Org.fertilizer.type +
#                     (1|study_ID) + # this is u_l (the between-study effect)
#                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
#                   sigma ~ 1 +  Org.fertilizer.type + 
#                     (1|study_ID) # this is u_s (between-study effect for variance)
# )
# 
# 
# # generate default priors
# prior_ls_ma1b <- default_prior(form_ls_ma1b, 
#                               data=dat, 
#                               data2=list(vcv=vcv),
#                               family=gaussian())
# 
# prior_ls_ma1b$prior[5] = "constant(1)" # meta-analysis assumes
# 
# # fitting model
# fit_ls_ma1b <- brm(
#   formula = form_ls_ma1b,
#   data = dat,
#   data2 = list(vcv=vcv),
#   chains = 2,
#   cores = 2,
#   iter = 9000,
#   warmup = 6000,
#   prior = prior_ls_ma1b,
#   control = list(adapt_delta=0.99, max_treedepth=15)
# )
# summary(fit_ls_ma1b)
# 
# # save this as rds
# saveRDS(fit_ls_ma1b, here("Rdata", "fit_ls_ma1b.rds"))
# 
# # location-scale meta-regression: Model 3
# 
# form_ls_ma1c <- bf(lnRR ~ 1 + Org.fertilizer.type +
#                      (1|p|study_ID) + # this is u_l (the between-study effect)
#                      (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
#                    sigma ~ 1 +  Org.fertilizer.type + 
#                      (1|p|study_ID) # this is u_s (between-study effect for variance)
# )
# 
# 
# # generate default priors
# prior_ls_ma1c <- default_prior(form_ls_ma1c, 
#                                data=dat, 
#                                data2=list(vcv=vcv),
#                                family=gaussian())
#
#' @Shinchi - [below line is incorrect - should be 7th prior not 5th]
# prior_ls_ma1c$prior[5] = "constant(1)" # meta-analysis assumes 
# 
# # fitting model
# fit_ls_ma1c <- brm(
#   formula = form_ls_ma1c,
#   data = dat,
#   data2 = list(vcv=vcv),
#   chains = 2,
#   cores = 2,
#   iter = 9000,
#   warmup = 6000,
#   prior = prior_ls_ma1c,
#   control = list(adapt_delta=0.99, max_treedepth=15)
# )
# summary(fit_ls_ma1c)
# 
# # save this as rds
# saveRDS(fit_ls_ma1c, here("Rdata", "fit_ls_ma1c.rds"))


# lnRR
vcv <- diag(dat$var.lnRR)

rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID

form_ma1 <- bf(lnRR ~ 1 +
                 (1|study_ID) + # this is u_l (the between-study effect)
                 (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_ma1 <- default_prior(form_ma1, 
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_ma1$prior[3] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior_ma1

# fitting model
fit_ma1 <- brm(
  formula = form_ma1,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ma1,
  control = list(adapt_delta=0.95, max_treedepth=15)
)

summary(fit_ma1)

# save this as rds
saveRDS(fit_ma1, here::here("Rdata", "fit_ma1.rds"))

# meta-regression

form_mr1 <- bf(lnRR ~ 1 + Org.fertilizer.type +
                  (1|study_ID) + # this is u_l (the between-study effect)
                  (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_mr1 <- default_prior(form_mr1,
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_mr1$prior[5] = "constant(1)" # meta-analysis assumes


# fitting model
fit_mr1 <- brm(
  formula = form_mr1,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_mr1,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_mr1)

# save this as rds
saveRDS(fit_mr1, here::here("Rdata", "fit_mr1.rds"))

# location-scale meta-regression: Model 1

form_ls_ma1 <- bf(lnRR ~ 1 + Org.fertilizer.type +
                    (1|study_ID) + # this is u_l (the between-study effect)
                    (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                  sigma ~ 1 +  Org.fertilizer.type
)


# generate default priors
prior_ls_ma1 <- default_prior(form_ls_ma1, 
                              data=dat, 
                              data2=list(vcv=vcv),
                              family=gaussian())

prior_ls_ma1$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma1 <- brm(
  formula = form_ls_ma1,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ls_ma1,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_ls_ma1)

# save this as rds
saveRDS(fit_ls_ma1, here("Rdata", "fit_ls_ma1.rds"))

# location-scale meta-regression: Model 2

form_ls_ma1b <- bf(lnRR ~ 1 + Org.fertilizer.type +
                    (1|study_ID) + # this is u_l (the between-study effect)
                    (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                  sigma ~ 1 +  Org.fertilizer.type + 
                    (1|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma1b <- default_prior(form_ls_ma1b, 
                              data=dat, 
                              data2=list(vcv=vcv),
                              family=gaussian())

prior_ls_ma1b$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma1b <- brm(
  formula = form_ls_ma1b,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 9000,
  warmup = 6000,
  prior = prior_ls_ma1b,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma1b)

# save this as rds
saveRDS(fit_ls_ma1b, here::here("Rdata", "fit_ls_ma1b.rds"))

# location-scale meta-regression: Model 3

form_ls_ma1c <- bf(lnRR ~ 1 + Org.fertilizer.type +
                     (1|p|study_ID) + # this is u_l (the between-study effect)
                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                   sigma ~ 1 +  Org.fertilizer.type + 
                     (1|p|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma1c <- default_prior(form_ls_ma1c, 
                               data=dat, 
                               data2=list(vcv=vcv),
                               family=gaussian())

prior_ls_ma1c$prior[7] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma1c <- brm(
  formula = form_ls_ma1c,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 12000,
  warmup = 9000,
  prior = prior_ls_ma1c,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma1c)

# save this as rds
saveRDS(fit_ls_ma1c, here::here("Rdata", "fit_ls_ma1c.rds"))

# lnCVR
#######

vcv <- diag(dat$var.lnCVR)

rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID

form_ma2 <- bf(lnCVR ~ 1 +
                 (1|study_ID) + # this is u_l (the between-study effect)
                 (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_ma2 <- default_prior(form_ma2, 
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_ma2$prior[3] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior_ma2

# fitting model
fit_ma2 <- brm(
  formula = form_ma2,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ma2,
  control = list(adapt_delta=0.95, max_treedepth=15)
)

summary(fit_ma2)

# save this as rds
saveRDS(fit_ma2, here::here("Rdata", "fit_ma2.rds"))

# meta-regression

form_mr2 <- bf(lnCVR ~ 1 + Org.fertilizer.type +
                 (1|study_ID) + # this is u_l (the between-study effect)
                 (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_mr2 <- default_prior(form_mr2,
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_mr2$prior[5] = "constant(1)" # meta-analysis assumes


# fitting model
fit_mr2 <- brm(
  formula = form_mr2,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_mr2,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_mr2)

# save this as rds
saveRDS(fit_mr2, here::here("Rdata", "fit_mr2.rds"))

# location-scale meta-regression: Model 1

form_ls_ma2 <- bf(lnCVR ~ 1 + Org.fertilizer.type +
                    (1|study_ID) + # this is u_l (the between-study effect)
                    (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                  sigma ~ 1 +  Org.fertilizer.type
)


# generate default priors
prior_ls_ma2 <- default_prior(form_ls_ma2, 
                              data=dat, 
                              data2=list(vcv=vcv),
                              family=gaussian())

prior_ls_ma2$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma2 <- brm(
  formula = form_ls_ma2,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ls_ma2,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_ls_ma2)

# save this as rds
saveRDS(fit_ls_ma2, here::here("Rdata", "fit_ls_ma2.rds"))

# location-scale meta-regression: Model 2

form_ls_ma2b <- bf(lnCVR ~ 1 + Org.fertilizer.type +
                     (1|study_ID) + # this is u_l (the between-study effect)
                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                   sigma ~ 1 +  Org.fertilizer.type + 
                     (1|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma2b <- default_prior(form_ls_ma2b, 
                               data=dat, 
                               data2=list(vcv=vcv),
                               family=gaussian())

prior_ls_ma2b$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma2b <- brm(
  formula = form_ls_ma2b,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 9000,
  warmup = 6000,
  prior = prior_ls_ma2b,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma2b)

# save this as rds
saveRDS(fit_ls_ma2b, here::here("Rdata", "fit_ls_ma2b.rds"))

# location-scale meta-regression: Model 3

form_ls_ma2c <- bf(lnCVR ~ 1 + Org.fertilizer.type +
                     (1|p|study_ID) + # this is u_l (the between-study effect)
                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                   sigma ~ 1 +  Org.fertilizer.type + 
                     (1|p|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma2c <- default_prior(form_ls_ma2c, 
                               data=dat, 
                               data2=list(vcv=vcv),
                               family=gaussian())

prior_ls_ma2c$prior[7] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma2c <- brm(
  formula = form_ls_ma2c,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 9000,
  warmup = 6000,
  prior = prior_ls_ma2c,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma2c)

# save this as rds
saveRDS(fit_ls_ma2c, here::here("Rdata", "fit_ls_ma2c.rds"))


# lnVR
######

vcv <- diag(dat$var.lnVR)

rownames(vcv) <- colnames(vcv) <- dat$effect_size_ID

form_ma3 <- bf(lnVR ~ 1 +
                 (1|study_ID) + # this is u_l (the between-study effect)
                 (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_ma3 <- default_prior(form_ma3, 
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_ma3$prior[3] = "constant(1)" # meta-analysis assumes sampling variance is known so fixing this to 1
prior_ma3

# fitting model
fit_ma3 <- brm(
  formula = form_ma3,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ma3,
  control = list(adapt_delta=0.95, max_treedepth=15)
)

summary(fit_ma3)

# save this as rds
saveRDS(fit_ma3, here::here("Rdata", "fit_ma3.rds"))

# meta-regression

form_mr3 <- bf(lnRR ~ 1 + Org.fertilizer.type +
                 (1|study_ID) + # this is u_l (the between-study effect)
                 (1|gr(effect_size_ID, cov = vcv)) # this is m (sampling error)
)

# generate default priors
prior_mr3 <- default_prior(form_mr3,
                           data=dat, 
                           data2=list(vcv=vcv),
                           family=gaussian())
prior_mr3$prior[5] = "constant(1)" # meta-analysis assumes


# fitting model
fit_mr3 <- brm(
  formula = form_mr1,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_mr3,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_mr3)

# save this as rds
saveRDS(fit_mr3, here::here("Rdata", "fit_mr3.rds"))

# location-scale meta-regression: Model 1

form_ls_ma3 <- bf(lnRR ~ 1 + Org.fertilizer.type +
                    (1|study_ID) + # this is u_l (the between-study effect)
                    (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                  sigma ~ 1 +  Org.fertilizer.type
)


# generate default priors
prior_ls_ma3 <- default_prior(form_ls_ma3, 
                              data=dat, 
                              data2=list(vcv=vcv),
                              family=gaussian())

prior_ls_ma3$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma3 <- brm(
  formula = form_ls_ma3,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 6000,
  warmup = 3000,
  prior = prior_ls_ma3,
  control = list(adapt_delta=0.95, max_treedepth=15)
)
summary(fit_ls_ma3)

# save this as rds
saveRDS(fit_ls_ma3, here::here("Rdata", "fit_ls_ma3.rds"))

# location-scale meta-regression: Model 2

form_ls_ma3b <- bf(lnRR ~ 1 + Org.fertilizer.type +
                     (1|study_ID) + # this is u_l (the between-study effect)
                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                   sigma ~ 1 +  Org.fertilizer.type + 
                     (1|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma3b <- default_prior(form_ls_ma3b, 
                               data=dat, 
                               data2=list(vcv=vcv),
                               family=gaussian())

prior_ls_ma3b$prior[5] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma3b <- brm(
  formula = form_ls_ma1b,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 9000,
  warmup = 6000,
  prior = prior_ls_ma3b,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma3b)

# save this as rds
saveRDS(fit_ls_ma3b, here::here("Rdata", "fit_ls_ma1b.rds"))

# location-scale meta-regression: Model 3

form_ls_ma3c <- bf(lnRR ~ 1 + Org.fertilizer.type +
                     (1|p|study_ID) + # this is u_l (the between-study effect)
                     (1|gr(effect_size_ID, cov = vcv)), # this is m (sampling error)
                   sigma ~ 1 +  Org.fertilizer.type + 
                     (1|p|study_ID) # this is u_s (between-study effect for variance)
)


# generate default priors
prior_ls_ma3c <- default_prior(form_ls_ma3c, 
                               data=dat, 
                               data2=list(vcv=vcv),
                               family=gaussian())

prior_ls_ma3c$prior[7] = "constant(1)" # meta-analysis assumes

# fitting model
fit_ls_ma3c <- brm(
  formula = form_ls_ma3c,
  data = dat,
  data2 = list(vcv=vcv),
  chains = 2,
  cores = 2,
  iter = 9000,
  warmup = 6000,
  prior = prior_ls_ma1c,
  control = list(adapt_delta=0.99, max_treedepth=15)
)
summary(fit_ls_ma3c)

# save this as rds
saveRDS(fit_ls_ma3c, here("Rdata", "fit_ls_ma1c.rds"))


#############################
# meta-analysis using glmmTMB
#############################

# lnRR
vcv <- diag(dat$var.lnRR)

row.names(vcv) <- colnames(vcv) <- dat$effect_size_ID

# this is necessary step (but can it make it automatic)
dat$effect_size_ID <- factor(dat$effect_size_ID)

g <- rep(1, nrow(dat))

# Model 1
glmm_lnRR1 <- glmmTMB(lnRR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                      equalto(0 + effect_size_ID|g, vcv),
                  dispformula = ~ 1 + Org.fertilizer.type,
                  data = dat,
                  REML = TRUE)

summary(glmm_lnRR1)

# Model 2
glmm_lnRR2 <- glmmTMB(lnRR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                        equalto(0 + effect_size_ID|g, vcv),
                      dispformula = ~ 1 + Org.fertilizer.type + (1 | study_ID),
                      data = dat,
                      REML = TRUE)

summary(glmm_lnRR2)

# lnCVR
vcv <- diag(dat$var.lnCVR)

row.names(vcv) <- colnames(vcv) <- dat$effect_size_ID

# this is necessary step (but can it make it automatic)
dat$effect_size_ID <- factor(dat$effect_size_ID)
g <- rep(1, nrow(dat))

# Model 1
glmm_lnCVR1 <- glmmTMB(lnCVR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                        equalto(0 + effect_size_ID|g, vcv),
                      dispformula = ~ 1 + Org.fertilizer.type,
                      data = dat,
                      REML = TRUE)
summary(glmm_lnCVR1)

# Model 2
glmm_lnCVR2 <- glmmTMB(lnCVR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                        equalto(0 + effect_size_ID|g, vcv),
                      dispformula = ~ 1 + Org.fertilizer.type + (1 | study_ID),
                      data = dat,
                      REML = TRUE)
summary(glmm_lnCVR2)

# lnVR
vcv <- diag(dat$var.lnVR)
row.names(vcv) <- colnames(vcv) <- dat$effect_size_ID

# this is necessary step (but can it make it automatic)
dat$effect_size_ID <- factor(dat$effect_size_ID)
g <- rep(1, nrow(dat))

# Model 1
glmm_lnVR1 <- glmmTMB(lnVR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                        equalto(0 + effect_size_ID|g, vcv),
                      dispformula = ~ 1 + Org.fertilizer.type,
                      data = dat,
                      REML = TRUE)
summary(glmm_lnVR1)

# Model 2
glmm_lnVR2 <- glmmTMB(lnVR ~ 1 + Org.fertilizer.type + (1 | study_ID) + 
                        equalto(0 + effect_size_ID|g, vcv),
                      dispformula = ~ 1 + Org.fertilizer.type + (1 | study_ID),
                      data = dat,
                      REML = TRUE)
summary(glmm_lnVR2)



