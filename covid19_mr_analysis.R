library(data.table)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
library(gwasglue)
library(ieugwasr)
library(ggplot2)


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

	outcome_dat$Phenotype <- GWAS_phenotype

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

instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "copper"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "copper_subsignificant.txt"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "zinc"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "zinc_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_meta"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_subsignificant_alspac"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "selenium_subsignificant_qimr"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "SARS-CoV-2_infection"
exposure <- "vitamin_K1"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "copper"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "copper_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "zinc"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "zinc_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_meta"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_toenail"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_subsignificant_qimr"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "selenium_subsignificant_alspac"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_population)"
exposure <- "vitamin_K1"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)



instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "copper_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "copper"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "zinc"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "zinc_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_meta"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_toenail"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_subsignificant_alspac"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "selenium_subsignificant_qimr.txt"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "Hospitalized_(ver_non-hospitalised)"
exposure <- "vitamin_K1"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "copper_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "copper"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "copper_subsignificant.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "copper_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "zinc"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "zinc_subsignificant.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "zinc_subsignificant"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium1_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_meta"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium2_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_toenail"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

instrument_file <- "selenium_subsignificant_qimr.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_subsignificant_qimr"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "selenium_subsignificant_alspac.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "selenium_subsignificant_alspac"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)


instrument_file <- "vitamink_instruments.txt"
GWAS_file <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
GWAS_phenotype <- "v.severe_COVID-19"
exposure <- "vitamin_K1"
simplemr(instrument_file, GWAS_file, GWAS_phenotype, exposure)

