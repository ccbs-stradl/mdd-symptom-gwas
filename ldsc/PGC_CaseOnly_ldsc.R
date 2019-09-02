# Make LDSC covariance matrix from PGC Case-Only symptom GWAS

library(GenomicSEM)
library(readr)
library(dplyr)

pgc_caseonly_dir <- 'sumstats/PGC/CaseOnly'

# number of cases and controls for each GWAS
pgc_caseonly_meta <- read_csv(file.path(pgc_caseonly_dir, 'CaseOnly.ldsc.csv'))

# munged sumstats filepaths
traits <- file.path(pgc_caseonly_dir, paste0('daner_', pgc_caseonly_meta$Symptom, '_CasesOnly.meta.ldsc.hm3.sumstats.gz'))

trait_names <- c('Mood', 'Anhed', 'WtLoss', 'WtGain', 'Insomn', 'Hyprsomn', 'Agit', 'Ret', 'Fatigue', 'Worth', 'Cnctr', 'Suicid')

# calculate sample and approximate pop prevalences
sample_prev <- with(pgc_caseonly_meta, Ncases / (Ncases + Ncontrols))
pop_prev <- sample_prev * 0.15

ld <- wld <- 'sumstats/reference/eur_w_ld_chr/'

# genomic covariance matrix
pgc_caseonly_cov <- ldsc(traits, sample_prev, pop_prev, ld, wld, trait_names)

# save out as RDS object and deparsed R code
saveRDS(pgc_caseonly_cov, 'ldsc/pgc_mdd_dsm_caseonly_covstruct.rds')
dput(pgc_caseonly_cov, 'ldsc/pgc_mdd_dsm_caseonly_covstruct.deparse.R', control=c('all', 'digits17'))

# check for exact match of deparsed object
identical(dget('ldsc/pgc_mdd_dsm_caseonly_covstruct.deparse.R'), pgc_caseonly_cov)
