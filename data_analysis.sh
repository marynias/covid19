  
#!/bin/bash
HOME=/mnt/storage/home/qh18484/scratch
scripts=/mnt/storage/home/qh18484/bin/covid19
analysis=$HOME/covid19

cd $analysis

#Download outcome data from COVID-19 HGI:

#Very severe COVID-19 round 5
#A2_ALL_eur
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.b37.txt.gz.tbi
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_A2_ALL_eur_leave_23andme_20210107.txt.gz.tbi

#Hospitalized COVID-19 round 5
#B1_ALL_eur - hospitalized vs non-hospitalized COVID
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.b37.txt.gz.tbi
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B1_ALL_eur_leave_23andme_20210107.txt.gz.tbi

#Hospitalized COVID-19 round 5
#B2_ALL_eur - hospitalized vs population
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.txt.gz.tbi
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.txt.gz.tbi


#Hospitalized COVID-19 round 5
#C2_ALL_eur - Covid-19 vs population
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.b37.txt.gz.tbi
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.txt.gz
wget https://storage.googleapis.com/covid19-hg-public/20201215/results/20210107/COVID19_HGI_C2_ALL_eur_leave_23andme_20210107.txt.gz.tbi

#Downloading results for round 6 from the GS bucket at gs://covid19-hg-public/20210415/results/20210607. However, turns out that EUR files have not been uploaded yet after all and so contacted the authors about their upload.


#Convert the outcome data to Two Sample MR format.
Rscript --vanilla $scripts/format_covid19_outcome.R

#Clump instruments for selenium and vitamin K1.
Rscript --vanilla $scripts/clump_exposures.



