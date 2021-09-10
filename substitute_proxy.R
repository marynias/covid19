library("vroom")

####Add "fake" exposure SNPs with the betas, se etc. lifted from their proxies.

convert_proxies <- function (my_gwas, output) {
	inst <- my_gwas[my_gwas$SNP=="rs73118350",]
	inst$SNP <- "rs2802728"
	inst$effect_allele.outcome <- "C"
	inst$other_allele.outcome <- "T" 
	inst$eaf.outcome <- 0.11

	my_gwas <- rbind(my_gwas, inst)

	inst <- my_gwas[my_gwas$SNP=="rs10895398",]
	inst$SNP <- "rs313426"
	inst$eaf.outcome <- 0.327

	my_gwas <- rbind(my_gwas, inst)

	inst <- my_gwas[my_gwas$SNP=="rs2453868",]
	inst$SNP <- "rs1532423"
	inst$effect_allele.outcome <- "G"
	inst$other_allele.outcome <- "A" 
	inst$eaf.outcome <- 0.59

	my_gwas <- rbind(my_gwas, inst)
	vroom_write(my_gwas, path = output)
}

my_gwas <- vroom("COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz")
my_output <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz")
my_output <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz")
my_output <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.outcome.txt.gz")
my_output <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.txt.gz")
my_output <- "COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.txt.gz")
my_output <- "COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.txt.gz")
my_output <- "COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)

my_gwas <- vroom("COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.txt.gz")
my_output <- "COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.outcome.prox.txt.gz"
convert_proxies(my_gwas, my_output)
#echo