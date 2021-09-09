library(dplyr)

#Clump our exposures of interest.
selenium1 <- read.delim("selenium1.txt", header=T, stringsAsFactors=F, row.names=NULL, sep="\t")
selenium2 <- read.delim("selenium2.txt", header=T, stringsAsFactors=F, row.names=NULL, sep="\t")
vitk1 <- read.delim("vitamink1.txt", header=T, stringsAsFactors=F, row.names=NULL, sep="\t")

ld_clump_local2 <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin, my_prefix){
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- paste0(my_prefix, "infile.txt")
  write.table(data.frame(SNP=dat[["variant.id"]], P=dat[["p.value"]]), file=fn, row.names=F, col.names=T, quote=F)
  
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
ldref = "/mnt/storage/home/qh18484/scratch/pleiotropy/data/plink_ref/EUR_all_1k"
plink_bin = genetics.binaRies::get_plink_binary()
clumped_sel1 <- ld_clump_local2(selenium1, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.001, "selenium1_")
clumped_sel2 <- ld_clump_local2(selenium2, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.001, "selenium2_")
clumped_vitk1 <- ld_clump_local2(vitk1, bfile=ldref, plink_bin=plink_bin, clump_kb = 10000, clump_r2 = 0.05, clump_p = 0.001, "vitk1_")