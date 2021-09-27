library(tidyverse)
library(ggforestplot)
library(ggforce)
library(forcats)
source("forestplot2.R")

combined_results <- read.csv("all_mr_results.csv", stringsAsFactors = FALSE)
#Rename values
combined_results <- combined_results[combined_results$method %in% c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode"),] 
combined_results[combined_results == 'selenium meta'] <- 'Se meta-analysis'
combined_results[combined_results == 'selenium_alspac_subsignificant'] <- 'Se ALSPAC subsignificant'
combined_results[combined_results == 'selenium_qimr_subsignificant'] <- 'Se QIMR subsignificant'
combined_results[combined_results == 'Hospitalized_(ver_non-hospitalised)'] <- "Hospitalized (ver. non-hospitalized)"
combined_results[combined_results == 'Hospitalized_(ver_population)'] <- "Hospitalized (ver. population)"
combined_results[combined_results == 'SARS-CoV-2_infection'] <- "SARS-CoV-2 infection"
combined_results[combined_results == 'v.severe_COVID-19'] <- "Very severe COVID-19"
combined_results[combined_results == 'zinc_subsignificant'] <- "Zn subsignificant"
combined_results[combined_results == 'zinc'] <- "Zinc"
combined_results[combined_results == 'copper_subsignificant'] <- "Cu subsignificant"
combined_results[combined_results == 'copper'] <- "Cu"
combined_results[combined_results == 'phylloquinone'] <- "vit. K1 subsignificant"
combined_results$outcome <- factor(combined_results$outcome)
combined_results$outcome <- factor(combined_results$outcome, levels=c("SARS-CoV-2 infection", "Hospitalized (ver. non-hospitalized)", "Hospitalized (ver. population)", "Very severe COVID-19"))
#combined_results$method <- factor(combined_results$method, levels=rev(c("Inverse variance weighted", "MR Egger", "Weighted median", "Weighted mode")))
#Subset to selenium only analyses
selens <- c("Se meta-analysis", "Se ALSPAC subsignificant", "Se QIMR subsignificant")
selenium <- combined_results[combined_results$exposure %in% selens,]
forestplot3(
  df = selenium,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in selenium concentration",
  colour = method,
  shape = method,
  xtickbreaks = c(0.6, 0.7, 0.8, 0.9, 1.1, 1.3, 1.8),
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )

#Selenium - outlier results
selenium <- combined_results[combined_results$exposure== "selenium toe-nail",]
forestplot(
  df = selenium,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in selenium concentration",
  colour = method,
  shape = method,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )


zincs <- c("Zinc", "Zn subsignificant")
zinc <- combined_results[combined_results$exposure %in% zincs,]

forestplot2(
  df = zinc,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in zinc concentration",
  colour = method,
  shape = method,
  xtickbreaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8),
  logodds = TRUE,
  my_levels=c("cu", "Zn subsignificant")
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )

coppers <- c("Cu", "Cu subsignificant")
copper <- combined_results[combined_results$exposure %in% coppers,]
copper$exposure <- factor(copper$exposure, levels=c("Cu", "Cu subsignificant"))

forestplot2(
  df = copper,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in copper concentration",
  xtickbreaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6),
  colour = method,
  shape = method,
  logodds = TRUE,
  my_levels=c("cu", "Cu subsignificant")
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )

vitamink <- combined_results[combined_results$exposure == "vit. K1 subsignificant",]
forestplot3(
  df = vitamink,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in vitamin K1 concentration",
  colour = method,
  shape = method,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )