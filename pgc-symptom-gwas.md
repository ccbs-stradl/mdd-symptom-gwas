---
title: PGC Symptom GWAS 
author: Bradley Jermy, Mark Adams 
output:
  html_document:
    toc: TRUE
    code_folding: hide
    number_sections: TRUE
    df_print: kable
    keep_md: true
  md_document:
    variant: markdown_github
---

Run GWAS of MDD symptoms on the [LISA server](http://geneticcluster.org/). Because these commands need to be run on a specific cluster with specific group permissions, this notebook documents the commands that were run rather than reproducing the whole analysis internally, with `eval=FALSE`.

# Directory setup

Softlink to all MDD files. Replace `PATH` with the actual path to the MDD Wave2 directory (listed in `README.mdd2sum`):


```bash
ln -s /home/PATH/* /home/USER/ 
```

# Install [Ricopili](https://sites.google.com/a/broadinstitute.org/ricopili/)


```bash
#Install and unpack ricopili - follow this link for custom installation https://docs.google.com/document/d/14aa-oeT5hF541I8hHsDAL_42oyvlHRC5FWR7gir4xco/edit#heading=h.y8igfs7neh22
wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.2019_Jun_25.001.tar.gz

tar -xvzf rp_bin.2019_Jun_25.001.tar.gz

#Create custom installation file in ricopili

#TEST!!!  Postimp RICOPILI code - run this home directory after all softlinks have been set-up
postimp_navi --out MDD1_case_con --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --refiex refiex.mddw2v01
```

# Symptom phenotype files

Launch an interactive session


```bash
#Go into interactive mode
srun -n 16 -t 1:00:00 --pty bash -il

#Load R

module load R/3.4.3

R
```

Load libraries and configure input directories


```r
library(yaml)
library(gdata)
library(stringr)
library(plyr)
library(dplyr)
library(tidyr)

# load configuration file (see config-example.yaml')
config <- yaml.load_file('config.yaml')

# input and output directories and files
mdd_v1_dir <- config$data$pgc$mdd$v1
output_dir <- config$data$pgc$mdd$output
ukb_remove <- config$data$pgc$mdd$ukb_remove
```

## Phenotype files 

Extract 9 pheno files - One for each symptom.


```r
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
studies <- c("RADIANT", "STR", "ROTTERDAM", "QIMR", "NESDANTR", "GSK", "GENRED", "GENRED2", "EDINBURGH", "CoLaus", "COFAMS", "BonnMH")
study_noupdate <- c("STARD", "SHIP0", "MARS", "Janssen", "GENPOD")

dir.create(output_dir, recursive=TRUE)


updated <- data.frame()
for (i in studies) {
mdd <- read.xls(file.path(mdd_v1_dir, paste0("secondary_phenotypes/1216_update/",i,"_PGC2updated.xls")), header=T)
updated <- rbind.fill(updated, mdd)
}

noupdate <- data.frame()
for (j in study_noupdate) {
mdd <- read.xls(file.path(mdd_v1_dir, paste0("secondary_phenotypes/1216_update/",j,"_PGC2.xls")), header=T)
noupdate <- rbind.fill(noupdate, mdd)
}

for (i in symptoms){
full <- rbind.fill(updated,noupdate)
full <- full[,c("ID1","ID2","Study","Sub.study", i)]
write.table(full,file.path(output_dir, paste0(i,"_PGCMDD2.tsv", sep="")), sep="\t", row.names=F, quote=F)
}
```

## Cases and controls symptom data 

Subset by studies that have symptom data for cases and controls


```r
for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,'_PGCMDD2.tsv', sep="")), sep="\t", header=T)
keepStudy <- c("COF", "Colaus", "NESDA/NTR", "ROTTERDAM", "SHIP-0", "RS")
#Change pheno file for 1=controls, 2=cases and -9=unknowns as stated on RICOPILI
MDD[,i][!MDD$Study %in% keepStudy] <- -9 #Remove cases and controls not in the correct cohort
MDD[,i][is.na(MDD[,i])] <- -9  #Recode NAs to -9
MDD[,i][!MDD[,i] %in% c(0, 1)] <- -9  #Recode values that are not 0 or 1 to -9
MDD[,i][MDD[,i]==1] <- 2 
MDD[,i][MDD[,i]==0] <- 1
print(table(MDD[[i]]))
write.table(MDD, file.path(output_dir, paste0(i,"_CasesControls_PGCMDD2_recoded.tsv", sep="")), sep="\t", row.names=F, quote=F)
}
```

## Cases-only symptom data

Subset by studies that have symptom data for cases only 


```r
for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,'_PGCMDD2.tsv', sep="")), sep="\t", header=T)
keepStudy <- c("STARD", "RADIANT", "STR", "QIMR_317370K", "QIMR_610660K", "QIMR_Omni", "MARS", "GSK", "genpod", "1", "2", NA)
#Change pheno file for 1=controls, 2=cases and -9=unknowns as stated on RICOPILI
MDD[,i][!MDD$Study %in% keepStudy] <- -9 #Remove cases and controls not in the correct cohort
MDD[,i][is.na(MDD[,i])] <- -9  #Recode NAs to -9
MDD[,i][!MDD[,i] %in% c(0, 1)] <- -9  #Recode values that are not 0 or 1 to -9
MDD[,i][MDD[,i]==1] <- 2 
MDD[,i][MDD[,i]==0] <- 1 
print(table(MDD[[i]]))
write.table(MDD, file.path(output_dir, paste0(i,"_CasesOnly_PGCMDD2_recoded.tsv", sep="")), sep="\t", row.names=F, quote=F)
}
```

## Overlap with UK Biobank 

Remove overlapping samples with UKB - remove according to ID1 and ID2 


```r
UKB <- read.table(ukb_remove, sep="\t", header=F)
UKB$ID1adj <- gsub("\\*.*","",UKB$V2) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,"_CasesControls_PGCMDD2_recoded.tsv", sep="")), sep="\t", header=T)
MDD$ID1adj <- gsub("\\*.*","",MDD$ID1) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
print(dim(MDD))
#filter MDD for IDs not in the UKB sample
MDD2 <- MDD %>%
filter(!MDD$ID2 %in% UKB$V3 & !MDD$ID1adj %in% UKB$V2)
print(dim(MDD2))
# check how many overlapped that had phenotypes
MDD %>%
  filter(MDD$ID2 %in% UKB$V3 & MDD$ID1adj %in% UKB$ID1adj) %>%
  filter_at(i, any_vars(. != -9)) %>%
  group_by(Study) %>%
  tally()
#293 removed due to overlap with UKB - 205 from RADIANT which makes sense
MDD_pheno <- MDD2[,c(1,2,5)]
write.table(MDD_pheno, file.path(output_dir, paste0(i,"_CasesControls_PGCMDD2_final.tsv", sep="")), sep="\t", row.names=F, col.names=F, quote=F)
}

for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,"_CasesOnly_PGCMDD2_recoded.tsv", sep="")), sep="\t", header=T)
MDD$ID1adj <- gsub("\\*.*","",MDD$ID1) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
print(dim(MDD))
#filter MDD for IDs not in the UKB sample
MDD2 <- MDD %>%
filter(!MDD$ID2 %in% UKB$V3 & !MDD$ID1adj %in% UKB$V2)
print(dim(MDD2))
MDD %>%
  filter(MDD$ID2 %in% UKB$V3 & MDD$ID1adj %in% UKB$ID1adj) %>%
  filter_at(i, any_vars(. != -9)) %>%
  group_by(Study) %>%
  tally()
#293 removed due to overlap with UKB - 205 from RADIANT which makes sense
MDD_pheno <- MDD2[,c(1,2,5)]
write.table(MDD_pheno, file.path(output_dir, paste0(i,"_CasesOnly_PGCMDD2_final.tsv", sep="")), sep="\t", row.names=F, col.names=F, quote=F)
}
```

# GWAS Pipeline

Run Ricopoli


```bash

cat > datasets_info << EOF
mdd_cof3_eur_sr-qc.hg19.ch.fl
mdd_col3_eur_sr-qc.hg19.ch.fl
mdd_rot4_eur_sr-qc.hg19.ch.fl
mdd_nes1_eur_sr-qc2.hg19.ch.fl
mdd_shp0_eur_sr-qc.hg19.ch.fl
mdd_boma_eur_sr-qc.hg19.ch.fl
mdd_gens_eur_sr-qc.hg19.ch.fl
mdd_gep3_eur_sr-qc.hg19.ch.fl
mdd_grnd_eur_sr-qc.hg19.ch.fl
mdd_gsk2_eur_sr-qc.hg19.ch.fl
mdd_mmi2_eur_sr-qc.hg19.ch.fl
mdd_mmo4_eur_sr-qc.hg19.ch.fl
mdd_qi3c_eur_sr-qc.hg19.ch.fl
mdd_qi6c_eur_sr-qc.hg19.ch.fl
mdd_qio2_eur_sr-qc.hg19.ch.fl
mdd_rad3_eur_sr-qc.hg19.ch.fl
mdd_rage_eur_sr-qc.hg19.ch.fl
mdd_rai2_eur_sr-qc.hg19.ch.fl
mdd_rau2_eur_sr-qc.hg19.ch.fl
mdd_rde4_eur_sr-qc2.hg19.ch.fl
mdd_stm2_eur_sr-qc.hg19.ch.fl
mdd_twg2_eur_sr-qc.hg19.ch.fl
EOF

#Run analysis
postimp_navi --out MDD1_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD1_PGCMDD2_final.txt #Done
postimp_navi --out MDD2_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD2_PGCMDD2_final.txt #Done
postimp_navi --out MDD3a_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3a_PGCMDD2_final.txt #Done
postimp_navi --out MDD3b_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3b_PGCMDD2_final.txt
postimp_navi --out MDD4a_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4a_PGCMDD2_final.txt
postimp_navi --out MDD4b_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4b_PGCMDD2_final.txt
postimp_navi --out MDD5a_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5a_PGCMDD2_final.txt
postimp_navi --out MDD5b_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5b_PGCMDD2_final.txt
postimp_navi --out MDD6_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD6_PGCMDD2_final.txt
postimp_navi --out MDD7_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD7_PGCMDD2_final.txt
postimp_navi --out MDD8_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD8_PGCMDD2_final.txt
postimp_navi --out MDD9_CasesAllCohorts --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD9_PGCMDD2_final.txt

```


