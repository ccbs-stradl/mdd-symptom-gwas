library(dplyr)
library(readr)

args <- commandArgs(TRUE)

prefix <- args[1]

gwas_prefix <- paste(prefix, 'sumstats', sep='.')

gwas_rds_dir <- here::here('meta', 'usergwas', gwas_prefix)

sumstats_txt <- here::here('meta', 'usergwas', paste0(gwas_prefix, '.txt.gz'))

symptoms_gwas <-
bind_rows(lapply(list.files(gwas_rds_dir, pattern=gwas_prefix, full.names=TRUE), readRDS)) %>%
as_tibble() %>%
arrange(CHR, BP)

write_tsv(symptoms_gwas, sumstats_txt)