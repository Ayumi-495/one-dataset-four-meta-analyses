library(tidyverse)
library(patchwork)
library(here)
library(brms)
library(tidybayes)
#library(bayesplot)
library(stringr)


# reading in data

orchard_ma_lnRR <- readRDS(here::here("Rdata", "orchard_ma_lnRR.rds"))
orchard_ma_lnCVR <- readRDS(here::here("Rdata", "orchard_ma_lnCVR.rds"))
orchard_mr_lnRR <- readRDS(here::here("Rdata", "orchard_mr_lnRR.rds"))
orchard_mr_lnCVR <- readRDS(here::here("Rdata", "orchard_mr_lnCVR.rds"))

fit_ma1 <- readRDS(here::here("Rdata", "fit_ma1.rds"))
fit_ma2 <- readRDS(here::here("Rdata", "fit_ma2.rds"))
fit_ls_ma1 <- readRDS(here::here("Rdata", "fit_ls_ma1.rds"))
fit_ls_ma2 <- readRDS(here::here("Rdata", "fit_ls_ma2.rds"))
                      
#fit_ls_ma2b <- readRDS(here::here("Rdata", "fit_ls_ma2b.rds"))                      

# Function to get variable names dynamically
get_variables_dynamic <- function(model, patterns) {
  variables <- get_variables(model)
  idx <- unique(unlist(lapply(patterns, function(p) grep(p, variables))))
  if (length(idx) == 0) return(character(0))
  variables[idx]
}

# preparation
rename_vars <- function(variable) {
  # fixed effects
  variable <- gsub("b_Intercept", "b_l_intercept", variable)
  variable <- gsub("b_sigma_Intercept", "b_s_intercept", variable)  
  variable <- gsub("b_Org.fertilizer.typeplant", "b_l_contrast", variable)
  variable <- gsub("b_sigma_Org.fertilizer.typeplant", "b_s_contrast", variable)            
  # random effects
  variable <- gsub("sd_study_ID__Intercept", "sd_study_ID", variable)  
  variable <- gsub("sigma", "sd_effect_ID", variable)  
  return(variable)
}

# Function to visualize fixed effects
visualize_fixed_effects <- function(model) {
  fixed_effect_vars <- get_variables_dynamic(model, "^b_")
  if (length(fixed_effect_vars) == 0) {
    message("No fixed effects found")
    return(NULL)
  }
  
  tryCatch({
    fixed_effects_samples <- model %>%
      spread_draws(!!!syms(fixed_effect_vars)) %>%
      pivot_longer(cols = all_of(fixed_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable))
    
    ggplot(fixed_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "lightcyan3", 
        color = "lightcyan4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Fixed effects", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_fixed_effects: ", e$message)
    return(NULL)
  })
}

# Function to visualize random effects
visualize_random_effects <- function(model) {
  random_effect_vars <- get_variables_dynamic(fit_ma2, c("^sd_study", "^sigma"))
  #random_effect_vars <- random_effect_vars[random_effect_vars != "sd_es_ID__Intercept"]
  if (length(random_effect_vars) == 0) {
    message("No random effects found")
    return(NULL)
  }
  
  tryCatch({
    random_effects_samples <- model %>%
      spread_draws(!!!syms(random_effect_vars)) %>%
      pivot_longer(cols = all_of(random_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable)) #%>%
    #mutate(.value = .value)  # leave SD as it is
    
    ggplot(random_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "olivedrab3", 
        color = "olivedrab4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Random effects (SD)", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_random_effects: ", e$message)
    return(NULL)
  })
}

# Function to visualize random effects
visualize_random_effects2 <- function(model) {
  random_effect_vars <- get_variables_dynamic(fit_ma2, c("^sd_study"))
  #random_effect_vars <- random_effect_vars[random_effect_vars != "sd_es_ID__Intercept"]
  if (length(random_effect_vars) == 0) {
    message("No random effects found")
    return(NULL)
  }
  
  tryCatch({
    random_effects_samples <- model %>%
      spread_draws(!!!syms(random_effect_vars)) %>%
      pivot_longer(cols = all_of(random_effect_vars), names_to = ".variable", values_to = ".value") %>%
      mutate(.variable = rename_vars(.variable)) #%>%
    #mutate(.value = .value)  # leave SD as it is
    
    ggplot(random_effects_samples, aes(x = .value, y = .variable)) +
      stat_halfeye(
        normalize = "xy", 
        point_interval = "mean_qi", 
        fill = "olivedrab3", 
        color = "olivedrab4"
      ) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#005") +
      labs(y = "Random effects (SD)", x = "Posterior values") +
      theme_classic()
  }, error = function(e) {
    message("Error in visualize_random_effects: ", e$message)
    return(NULL)
  })
}

#############
# Figure 1
#############

p1 <- orchard_ma_lnRR +  
  annotate(geom="text", x= 1.5, y= -2, label= paste0("italic(I)^{2} == ", 98.11, "*\"%\""), 
           color="black", parse = TRUE, size = 4)

p2 <- orchard_ma_lnCVR +  
  annotate(geom="text", x= 1.5, y= -3, label= paste0("italic(I)^{2} == ", 56.53, "*\"%\""), 
           color="black", parse = TRUE, size = 4)

p3 <- visualize_fixed_effects(fit_ma1)

p4 <- visualize_fixed_effects(fit_ma2) 

p5 <- visualize_random_effects(fit_ma1)

p6 <- visualize_random_effects(fit_ma2)


(p1 + p2) / (p3 + p4) / (p5 + p6) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 16, face = "bold"))

#############
# Figure 2
#############

p7 <- orchard_mr_lnRR #+  
  #annotate(geom="text", x= 1.5, y= -2, label= paste0("italic(R)^{2} == ", 3.24, "*\"%\""), 
  #         color="black", parse = TRUE, size = 4)

p8 <- orchard_mr_lnCVR #+  
  #annotate(geom="text", x= 1.5, y= -3, label= paste0("italic(R)^{2} == ", 2.62, "*\"%\""), 
  #         color="black", parse = TRUE, size = 4)

p9 <- visualize_fixed_effects(fit_ls_ma1)

p10 <- visualize_fixed_effects(fit_ls_ma2)

p11 <- visualize_random_effects2(fit_ls_ma1)

p12 <- visualize_random_effects2(fit_ls_ma2)

((p7 + p8) / (p9 + p10) / (p11 + p12) +
    plot_layout(heights = c(1, 1.5, 0.5)) +
    plot_annotation(tag_levels = 'A')) &
  theme(plot.tag = element_text(size = 16, face = "bold"))
