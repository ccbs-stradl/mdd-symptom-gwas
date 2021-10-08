# Common factor GWAS

library(GenomicSEM)

k <- 1

covstruct_prefix <- 'agds_pgc.alspac_ukb.common.covstruct'
sumstats_prefix <- 'agds_pgc.alspac_ukb.common.sumstats'
gwas_prefix <- 'agds_pgc.alspac_ukb.common_bif.sumstats'

covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))
gwas_rds <- here::here('meta', 'usergwas', gwas_prefix, paste(gwas_prefix, k, 'rds', sep='.'))

dir.create(dirname(gwas_rds))

symptoms_covstruct <- dget(covstruct_r)

symptoms_sumstats <- readRDS(sumstats_rds)

clin_pop_bif.model <- "
A2 =~ PopDep + PopAnh + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A =~ ClinAppInc + ClinSui + PopDep + PopAnh + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A ~~ 0*A2
"

clin_pop_bif.model <- "F1 =~ PopDep + PopAnh + PopSui
F1 ~ SNP"

clin_pop_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop_bif.model)


SNPs <- symptoms_sumstats[1:10,] %>%
	select(SNP, CHR, BP, MAF, A1, A2,
		   ends_with('PopDep'), ends_with('PopAnh'), ends_with('PopSui'))

symptoms_gwas <- userGWAS(covstruc=symptoms_covstruct,
						 SNPs=SNPs,
						 model=clin_pop_bif.model,
						 sub=c("F1 ~ SNP"),
						 cores=2)
				 
saveRDS(symptoms_gwas, gwas_rds)