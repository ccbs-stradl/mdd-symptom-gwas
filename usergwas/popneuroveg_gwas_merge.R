# Common factor GWAS merge

library(GenomicSEM)
library(dplyr)
library(readr)

covstruct_prefix <- 'agds_pgc.alspac_ukb.neuroveg.covstruct'
sumstats_prefix <- 'agds_pgc.alspac_ukb.neuroveg.sumstats'
gwas_prefix <- 'agds_pgc.alspac_ukb.neuroveg.sumstats'

covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))
gwas_rds_dir <- here::here('meta', 'usergwas', gwas_prefix)

sumstats_txt <- here::here('meta', 'usergwas', paste0(gwas_prefix, '.txt.gz'))

symptoms_covstruct <- dget(covstruct_r)

symptoms_gwas <-
bind_rows(lapply(list.files(gwas_rds_dir, pattern=gwas_prefix, full.names=TRUE), readRDS)) %>%
as_tibble() %>%
arrange(CHR, BP) %>%
filter(fail == 0) %>%
select(-i, -lhs, -op, -rhs, -fail, -warning)

write_tsv(symptoms_gwas, sumstats_txt)