library(tidyverse)
library(ggforestplot)
library(ggforce)

combined_results <- read.csv("all_mr_results.csv")
combined_results$outcome <- factor(combined_results$outcome)
#Subset to selenium only analyses
selens <- c("selenium meta", "selenium_alspac_subsignificant", "selenium_qimr_subsignificant")
selenium <- combined_results[combined_results$exposure %in% selens,]

forestplot(
  df = selenium,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in selenium concentration",
  colour = method,
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
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )


zincs <- c("zinc", "zinc_subsignificant")
zinc <- combined_results[combined_results$exposure %in% zincs,]

forestplot(
  df = zinc,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in zinc concentration",
  colour = method,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )

coppers <- c("copper", "copper_subsignificant")
copper <- combined_results[combined_results$exposure %in% coppers,]

forestplot(
  df = copper,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in copper concentration",
  colour = method,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )

vitamink <- combined_results[combined_results$exposure == "phylloquinone",]

forestplot(
  df = vitamink,
  name = exposure,
  estimate = b,
  se = se, 
  pvalue = pval,
  psignif = 0.05,
  xlab = "odds-ratio for COVID-19 outcome (95% CI) per 1-SD increment in vitamin K concentration",
  colour = method,
  logodds = TRUE
) +
  ggforce::facet_col(
    facets = ~outcome,
    scales = "free_y",
    space = "free"
  )