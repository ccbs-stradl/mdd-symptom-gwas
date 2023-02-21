library(readr)
library(dplyr)

sumstats <- read_tsv('all.common.sumstats.txt.gz')

effect_loci <- read_tsv('FUMA_job233308/GenomicRiskLoci.txt')
effect_genes <- read_tsv('FUMA_job233308/genes.txt')

het_loci <- read_tsv('FUMA_job233309/GenomicRiskLoci.txt')
het_genes <- read_tsv('FUMA_job233309/genes.txt')

effect <- effect_loci |> inner_join(select(effect_genes, GenomicLocus, ensg, symbol), by='GenomicLocus')

het <- het_loci |> inner_join(select(het_genes, GenomicLocus, ensg, symbol), by='GenomicLocus')

loci <-
bind_rows(effect, het) |>
left_join(sumstats, by=c('rsID'='SNP')) |>
select(SNP=rsID, CHR, BP, A1, A2, MAF, Beta=est, SE=se_c, P=Pval_Estimate, Q, Q_df, Q_pval, gene=symbol, ensg, locus.left=start, locus.right=end)
