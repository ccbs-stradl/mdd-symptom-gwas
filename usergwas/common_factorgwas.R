# Common factor GWAS

library(GenomicSEM)

k <- 1

covstruct_prefix <- 'agds_pgc.alspac_ukb.common.covstruct'
sumstats_prefix <- 'agds_pgc.alspac_ukb.common.sumstats'
gwas_prefix <- 'agds_pgc.alspac_ukb.common.sumstats'

covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))
gwas_rds <- here::here('meta', 'usergwas', gwas_prefix, paste(gwas_prefix, k, 'rds', sep='.'))

dir.create(dirname(gwas_rds))

symptoms_covstruct <- dget(covstruct_r)

symptoms_sumstats <- readRDS(sumstats_rds)

symptoms_gwas <- commonfactorGWAS(covstruc=symptoms_covstruct, SNPs=symptoms_sumstats[1:10000,], cores=8)
				 
saveRDS(symptoms_gwas, gwas_rds)