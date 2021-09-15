library(data.table)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
library(gwasglue)
library(ieugwasr)
library(ggplot2)
library(MendelianRandomization)


simplemr <- function(instrument_file, GWAS_file, GWAS_phenotype, exposure) {

	output_file <- paste(exposure, "vs", GWAS_phenotype, sep="_")

	exp_dat <- read_exposure_data(
	    filename = instrument_file,
	    sep = "\t",
	    snp_col = "SNP",
	    beta_col = "beta",
	    se_col = "se",
	    effect_allele_col = "effect_allele",
	    other_allele_col = "other_allele",
	    eaf_col = "eaf",
	    pval_col = "pval",
	    units_col = "units",
	    gene_col = "Gene",
	    samplesize_col = "samplesize",
	    phenotype_col = "Phenotype"
	)

	outcome_dat <- read_outcome_data(
	    snps = exp_dat$SNP,
	    filename = GWAS_file,
	    sep = "\t",
	    snp_col = "SNP",
	    beta_col = "beta.outcome",
	    se_col = "se.outcome",
	    effect_allele_col = "effect_allele.outcome",
	    other_allele_col = "other_allele.outcome",
	    eaf_col = "eaf.outcome",
	    pval_col = "pval.outcome"
	)

	outcome_dat$outcome <- GWAS_phenotype

	dat <- harmonise_data(
	    exposure_dat = exp_dat, 
	    outcome_dat = outcome_dat
	)

	mr <- mr(dat)
	odds_ratio_mr <- generate_odds_ratios(mr)
	write.csv(odds_ratio_mr, paste0(output_file, "_mr.csv"), quote=T, row.names=F)

	heterogen <- mr_heterogeneity(dat)
	write.csv(heterogen, paste0(output_file, "_heterogeneity.csv"), quote=T, row.names=F)

	pleio <- mr_pleiotropy_test(dat)
	write.csv(pleio, paste0(output_file, "_pleiotropy.csv"), quote=T, row.names=F)

	single_snp <- mr_singlesnp(dat)
	odds_ratio_single_snp <- generate_odds_ratios(single_snp)
	write.csv(odds_ratio_single_snp, paste0(output_file, "_singlesnp.csv"), quote=T, row.names=F)

	leave1out <- mr_leaveoneout(dat)
	odds_ratio_leave1out <- generate_odds_ratios(leave1out)
	write.csv(odds_ratio_leave1out, paste0(output_file, "_leave1out.csv"), quote=T, row.names=F)

	p1 <- mr_scatter_plot(mr, dat)
	ggsave(p1[[1]], file=paste0(output_file, "_scatterplot.pdf"), width=7, height=7)
	#ggsave(p1[[1]], file=paste0(output_file, "_scatterplot.png"), width=7, height=7)

	p2 <- mr_forest_plot(single_snp)
	ggsave(p2[[1]], file=paste0(output_file, "_forestplot.pdf"), width=7, height=7)
	#ggsave(p2[[1]], file=paste0(output_file, "_forestplot.png"), width=7, height=7)

	p3 <- mr_leaveoneout_plot(leave1out)
	ggsave(p3[[1]], file=paste0(output_file, "_leave1out.pdf"), width=7, height=7)
	#ggsave(p3[[1]], file=paste0(output_file, "_leave1out.png"), width=7, height=7)

	p4 <- mr_funnel_plot(single_snp)
	ggsave(p4[[1]], file=paste0(output_file, "_singlesnp.pdf"), width=7, height=7)
	#ggsave(p4[[1]], file=paste0(output_file, "_singlesnp.png"), width=7, height=7)
}

mr_cochraneq <- function(instrument_file, GWAS_file, GWAS_phenotype, exposure) {

	output_file <- paste(exposure, "vs", GWAS_phenotype, sep="_")

	exp_dat <- read_exposure_data(
	    filename = instrument_file,
	    sep = "\t",
	    snp_col = "SNP",
	    beta_col = "beta",
	    se_col = "se",
	    effect_allele_col = "effect_allele",
	    other_allele_col = "other_allele",
	    eaf_col = "eaf",
	    pval_col = "pval",
	    units_col = "units",
	    gene_col = "Gene",
	    samplesize_col = "samplesize",
	    phenotype_col = "Phenotype"
	)

	outcome_dat <- read_outcome_data(
	    snps = exp_dat$SNP,
	    filename = GWAS_file,
	    sep = "\t",
	    snp_col = "SNP",
	    beta_col = "beta.outcome",
	    se_col = "se.outcome",
	    effect_allele_col = "effect_allele.outcome",
	    other_allele_col = "other_allele.outcome",
	    eaf_col = "eaf.outcome",
	    pval_col = "pval.outcome"
	)

	outcome_dat$outcome <- GWAS_phenotype

	dat <- harmonise_data(
	    exposure_dat = exp_dat, 
	    outcome_dat = outcome_dat
	)

	#Convert to MendelianRandomization format.
	dat2 <- dat_to_MRInput(dat)
	all_ivw_results <- MendelianRandomization::mr_ivw(dat2[[1]])
	cochranesq <- all_ivw_results$Heter.Stat
	write.table(cochranesq, paste0(output_file, "_cochraneq.tsv"), quote=F)
	print(all_ivw_results)
}

instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "copper"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "copper_subsignificant.txt"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "zinc"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "zinc_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_meta"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_toenail"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_subsignificant_alspac"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_subsignificant_qimr"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "vitamin_K1"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "copper"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "copper_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "zinc"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "zinc_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_meta"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_toenail"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_subsignificant_qimr"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_subsignificant_alspac"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "vitamin_K1"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)



instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "copper_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "copper"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "zinc"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "zinc_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_meta"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_toenail"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_subsignificant_alspac"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_subsignificant_qimr.txt"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "vitamin_K1"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "copper"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "copper_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "zinc"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "zinc_subsignificant"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_meta"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_toenail"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_subsignificant_qimr"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_subsignificant_alspac"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "vitamin_K1"
mr_cochraneq(instrument_file, GWAS_file, GWAS_phenotype, exposure)
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


loadFile <- function(x) {
  print (x)
  df <- read.csv(x, header=T, stringsAsFactors=F,row.names=NULL)
  return (df)
}

 all <- list.files(pattern="\\mr.csv")
 all_normal <- lapply(all, loadFile)
 all_normal_together <- do.call(rbind,all_normal)

 write.csv(all_normal_together, "all_mr_results.csv", row.names=F)