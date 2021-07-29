# convert ALSPAC sumstats to daner format

library(dplyr)
library(readr)
library(ggplot2)

args <- commandArgs(TRUE)

sumstats_file <- args[1]
out <- args[2]

sumstats <- read_table2(sumstats_file)

# number of cases and controls
n_cases <- max(sumstats$n_cases)
n_controls <- max(sumstats$n_controls)

# frequency column names
frq_a_col <- paste0('FRQ_A_', n_cases)
frq_u_col <- paste0('FRQ_U_', n_controls)

# remove extreme BETA/SEs
se_cutoff <- 2

sumstats_qc <- sumstats %>%
mutate(qc=if_else(SE <= se_cutoff, true='include', false='exclude'))


daner <-
sumstats %>%
filter(SE <= se_cutoff) %>%
transmute(CHR=as.numeric(CHR), SNP, BP,
A1=effect_allele, A2=noneffect_allele,
!!frq_a_col:=AF_coded_cases,
!!frq_u_col:=AF_coded_controls,
INFO=info,
OR=exp(BETA),
SE,
P)

out_daner <- paste(out, 'txt.gz', sep='.')
out_qc_png <- paste(out, 'qc.png', sep='.')
out_or_png <- paste(out, 'or.png', sep='.')

write_tsv(daner, out_daner)

ggplot(sumstats_qc, aes(x=abs(BETA), fill=qc)) +
geom_histogram(bins=1000) +
scale_y_log10()
ggsave(out_qc_png)

ggplot(daner, aes(x=log(OR))) +
geom_histogram(bins=100)
ggsave(out_or_png)
