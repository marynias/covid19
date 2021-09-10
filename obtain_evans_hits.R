library("gwasvcf")
library("dplyr")
library("gwasglue")

open_gwas <- "/mnt/storage/home/qh18484/opengwas/public/"
gwasvcf::set_bcftools("/mnt/storage/software/apps/BCFTOOLS/bcftools/bcftools")
ldref = "/mnt/storage/home/qh18484/scratch/pleiotropy/data/plink_ref/EUR_all_1k"
plink_bin = genetics.binaRies::get_plink_binary()

ld_clump_local2 <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin, my_prefix){
  # Make textfile
  shell <- "sh"
  fn <- paste0(my_prefix, "infile.txt")
  write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)
  
  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --clump ", shQuote(fn, type=shell), 
    " --clump-p1 ", clump_p, 
    " --clump-p2 ", 1, 
    " --clump-r2 ", clump_r2, 
    " --clump-kb ", clump_kb, 
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2)
  res <- read.table(paste(fn, ".clumped", sep=""), header=T)
  #unlink(paste(fn, "*", sep=""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
  if(nrow(y) > 0)  {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
  }
  return(res)

}

#zinc ieu-a-1079
zinc <- "ieu-a-1079"
#Establish other top zinc hits and clump them.

gwas_input <- paste0(zinc, "/", zinc, ".vcf.gz")
full_path <- paste0(open_gwas, gwas_input)
vcf <- query_gwas(full_path, pval=0.00001)
out <- gwasvcf_to_TwoSampleMR(vcf, type="outcome")
write.table(out, "zinc-ieu-a-1079.txt", quote=F, row.names=F)

zinc_hits <- vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% dplyr::mutate(pval=10^{-LP})

# Perform clumping
for_clumping <- zinc_hits %>%
  dplyr::select(ID, pval) %>%
  dplyr::rename(rsid = ID)

clumped <- ld_clump_local2(for_clumping, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.00001, "zinc_weak")

copper <- "ieu-a-1073"
gwas_input <- paste0(copper, "/", copper, ".vcf.gz")
full_path <- paste0(open_gwas, gwas_input)
vcf <- query_gwas(full_path, pval=0.00001)
out <- gwasvcf_to_TwoSampleMR(vcf, type="outcome")
write.table(out, "copper-ieu-a-1073.txt", quote=F, row.names=F)

copper_hits <- vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% dplyr::mutate(pval=10^{-LP})

# Perform clumping
for_clumping <- copper_hits %>%
  dplyr::select(ID, pval) %>%
  dplyr::rename(rsid = ID)

clumped <- ld_clump_local2(for_clumping, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.00001, "copper_weak")

##ieu-a-1075 - ALSPAC Selenium
selenium <- "ieu-a-1075"
gwas_input <- paste0(selenium, "/", selenium, ".vcf.gz")
full_path <- paste0(open_gwas, gwas_input)
vcf <- query_gwas(full_path, pval=0.00001)
out <- gwasvcf_to_TwoSampleMR(vcf, type="outcome")
write.table(out, "selenium-ieu-a-1075.txt", quote=F, row.names=F)

selenium_hits <- vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% dplyr::mutate(pval=10^{-LP})

# Perform clumping
for_clumping <- selenium_hits %>%
  dplyr::select(ID, pval) %>%
  dplyr::rename(rsid = ID)

clumped <- ld_clump_local2(for_clumping, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.00001, "selenium_weak_alspac")