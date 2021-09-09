library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
library(ieugwasr)




convert2outcome <- function(gwas, outcome, output) {

gwas_outcome_format <-
  vroom(gwas,       # vroom is faster than fread!
        #  only read in columns that we need and 
        col_select = c("rsid","all_inv_var_meta_beta","all_inv_var_meta_sebeta",
                       "ALT","REF","all_meta_AF",
                       "all_inv_var_meta_p")) %>% 
  # format data into the 'outcome' format right away
  format_data(., type = "outcome",
              snp_col = "rsid",
              beta_col = "all_inv_var_meta_beta",
              se_col = "all_inv_var_meta_sebeta",
              effect_allele_col = "ALT",
              other_allele_col = "REF",
              eaf_col = "all_meta_AF",
              pval_col  = "all_inv_var_meta_p") %>% 
  # store trait name in the data frame (this will make later analysis easier)
  mutate(outcome = outcome)

# save this tidy and formatted file so that it cam be used directly in MR analysis later
vroom_write(gwas_outcome_format, path = output)

}

# raw file location
raw_gwas_file <- 'COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.txt.gz'
my_output <- 'COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz'

convert2outcome(raw_gwas_file, "A2_severe", my_output)


raw_gwas_file <- 'COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.txt.gz'
my_output <- 'COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz'

convert2outcome(raw_gwas_file, "B1_hospitalized", my_output)

raw_gwas_file <- 'COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.txt.gz'
my_output <- 'COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz'

convert2outcome(raw_gwas_file, "B2_hospitalized", my_output)

raw_gwas_file <- 'COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.txt.gz'
my_output <- 'COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz'

convert2outcome(raw_gwas_file, "C2_hospitalized", my_output)

# raw file location
raw_gwas_file <- 'COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.txt.gz'
my_output <- 'COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.txt.gz'

convert2outcome(raw_gwas_file, "A2_severe", my_output)


raw_gwas_file <- 'COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.txt.gz'
my_output <- 'COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.txt.gz'

convert2outcome(raw_gwas_file, "B1_hospitalized", my_output)

raw_gwas_file <- 'COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.txt.gz'
my_output <- 'COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.txt.gz'

convert2outcome(raw_gwas_file, "B2_hospitalized", my_output)

raw_gwas_file <- 'COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.txt.gz'
my_output <- 'COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.txt.gz'

convert2outcome(raw_gwas_file, "C2_hospitalized", my_output)