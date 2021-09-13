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
out <- gwasvcf_to_TwoSampleMR(vcf, type="exposure")
write.table(out, "zinc-ieu-a-1079.txt", quote=F, row.names=F, sep="\t")

mr_loci <- out %>% dplyr::select(chr.outcome, pos.outcome, other_allele.outcome, effect_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome, samplesize.outcome, SNP) %>% dplyr::rename(chromosome=chr.outcome, position=pos.outcome, effect_allele=other_allele.outcome, non_effect_allele=effect_allele.outcome, eaf=eaf.outcome, variant_ID=SNP, beta=beta.outcome, se=se.outcome, pval=pval.outcome, N=samplesize.outcome) %>% dplyr::select(variant_ID, chromosome, position, effect_allele, non_effect_allele, eaf, beta, se, pval, N)

zn_select_hits <- c("rs2120019", "rs11638477", "rs13273360", "rs7148590", "rs11232535", "rs11763353", "rs4333127", "rs10484101", "rs10931753", "rs7569234", "rs8099461", "rs17511001", "rs6545343", "rs17097781")
mr_loci <- mr_loci[mr_loci$variant_ID %in% zn_select_hits,]
write.table(mr_loci, "zinc_evans_weak_hits.txt", quote=F, row.names=F, sep="\t")

final_hits <- c("rs7569234", "rs4333127", "rs17511001", "rs11763353", "rs13273360", "rs11232535", "rs7148590", "rs10484101", "rs17097781", "rs2120019", "rs11638477", "rs8099461")
final_hits_out <- out[out$SNP %in% final_hits,]
write.table(final_hits_out, "zinc_evans_weak_hits.exposure", quote=F, row.names=F, sep="\t")

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
out <- gwasvcf_to_TwoSampleMR(vcf, type="exposure")
write.table(out, "copper-ieu-a-1073.txt", quote=F, row.names=F, sep="\t")

mr_loci <- out %>% dplyr::select(chr.outcome, pos.outcome, other_allele.outcome, effect_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome, samplesize.outcome, SNP) %>% dplyr::rename(chromosome=chr.outcome, position=pos.outcome, effect_allele=other_allele.outcome, non_effect_allele=effect_allele.outcome, eaf=eaf.outcome, variant_ID=SNP, beta=beta.outcome, se=se.outcome, pval=pval.outcome, N=samplesize.outcome) %>% dplyr::select(variant_ID, chromosome, position, effect_allele, non_effect_allele, eaf, beta, se, pval, N)

cu_select_hits <- c("rs2769264", "rs2769270", "rs1175550", "rs10014072", "rs12153606", "rs12582659", "rs3857536", "rs9324493", "rs764560", "rs13074172", "rs572585", "rs7206796")
mr_loci <- mr_loci[mr_loci$variant_ID %in% cu_select_hits,]
write.table(mr_loci, "copper_evans_weak_hits.txt", quote=F, row.names=F, sep="\t")

final_hits <- c("rs1175550", "rs2769264", "rs2769270", "rs13074172", "rs10014072", "rs12153606", "rs9324493")
final_hits_out <- out[out$SNP %in% final_hits,]
write.table(final_hits_out, "copper_evans_weak_hits.exposure", quote=F, row.names=F, sep="\t")

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
out <- gwasvcf_to_TwoSampleMR(vcf, type="exposure")
write.table(out, "selenium-ieu-a-1075.txt", quote=F, row.names=F, sep="\t")

mr_loci <- out %>% dplyr::select(chr.outcome, pos.outcome, other_allele.outcome, effect_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome, samplesize.outcome, SNP) %>% dplyr::rename(chromosome=chr.outcome, position=pos.outcome, effect_allele=other_allele.outcome, non_effect_allele=effect_allele.outcome, eaf=eaf.outcome, variant_ID=SNP, beta=beta.outcome, se=se.outcome, pval=pval.outcome, N=samplesize.outcome) %>% dplyr::select(variant_ID, chromosome, position, effect_allele, non_effect_allele, eaf, beta, se, pval, N)

se_select_hits <- c("rs921943", "rs694290", "rs11948804", "rs12951643", "rs2631524", "rs6823178", "rs3770549", "rs9609603", "rs667707", "rs1521328", "rs17771017", "rs1153734", "rs9348626")
mr_loci <- mr_loci[mr_loci$variant_ID %in% se_select_hits,]
write.table(mr_loci, "selenium_alspac_weak_hits.txt", quote=F, row.names=F, sep="\t")

final_hits <- c("rs17771017", "rs1153734", "rs1521328", "rs6823178", "rs921943", "rs694290", "rs11948804", "rs9348626", "rs667707", "rs2631524", "rs12951643", "rs9609603")
final_hits_out <- out[out$SNP %in% final_hits,]
write.table(final_hits_out, "selenium_alspac_weak_hits.exposure", quote=F, row.names=F, sep="\t")

selenium_hits <- vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% dplyr::mutate(pval=10^{-LP})

# Perform clumping
for_clumping <- selenium_hits %>%
  dplyr::select(ID, pval) %>%
  dplyr::rename(rsid = ID)

clumped <- ld_clump_local2(for_clumping, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.00001, "selenium_weak_alspac")

##ieu-a-1077 - QIMR Selenium
selenium <- "ieu-a-1077"
gwas_input <- paste0(selenium, "/", selenium, ".vcf.gz")
full_path <- paste0(open_gwas, gwas_input)
vcf <- query_gwas(full_path, pval=0.00001)
out <- gwasvcf_to_TwoSampleMR(vcf, type="exposure")
write.table(out, "selenium-ieu-a-1077.txt", quote=F, row.names=F, sep="\t")

mr_loci <- out %>% dplyr::select(chr.outcome, pos.outcome, other_allele.outcome, effect_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome, samplesize.outcome, SNP) %>% dplyr::rename(chromosome=chr.outcome, position=pos.outcome, effect_allele=other_allele.outcome, non_effect_allele=effect_allele.outcome, eaf=eaf.outcome, variant_ID=SNP, beta=beta.outcome, se=se.outcome, pval=pval.outcome, N=samplesize.outcome) %>% dplyr::select(variant_ID, chromosome, position, effect_allele, non_effect_allele, eaf, beta, se, pval, N)

se_select_hits <- c("rs7700970", "rs4950779", "rs10023369", "rs11779526", "rs478651", "rs3785832", "rs7163368", "rs713550", "rs10475981", "rs9309294", "rs11114989", "rs1492483", "rs2461362", "rs4467633", "rs566108", "rs4779561")
mr_loci <- mr_loci[mr_loci$variant_ID %in% se_select_hits,]
write.table(mr_loci, "selenium_qimr_weak_hits.txt", quote=F, row.names=F, sep="\t")

final_hits <- c("rs566108", "rs4950779", "rs9309294", "rs2461362", "rs1492483", "rs10023369", "rs4467633", "rs478651", "rs7700970", "rs10475981", "rs713550", "rs11114989", "rs4779561", "rs7163368", "rs3785832")
final_hits_out <- out[out$SNP %in% final_hits,]
write.table(final_hits_out, "selenium_qimr_weak_hits.exposure", quote=F, row.names=F, sep="\t")

selenium_hits <- vcf_to_granges(vcf) %>% dplyr::as_tibble() %>% dplyr::mutate(pval=10^{-LP})

# Perform clumping
for_clumping <- selenium_hits %>%
  dplyr::select(ID, pval) %>%
  dplyr::rename(rsid = ID)

clumped <- ld_clump_local2(for_clumping, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.00001, "selenium_weak_qimr")