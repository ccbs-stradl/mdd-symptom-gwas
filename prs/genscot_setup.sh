# link Generation Scotland data
# ln -s /path/to/genscot genscot

# prepare phenotype file
R
library(readr)
library(dplyr)

# open SCID
scid <- read_csv('genscot/phenotypes/SCID_QC_201113_GWASids.csv')

# open fam file
fam <- read_table('genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA_1.fam', col_names=c('FID', 'IID', 'father', 'mother', 'sex', 'pheno'))

# get affect/unaffected for depression
scid_pheno <- 
scid %>% select(gwas, SCID_Diagnosis) %>%
mutate(MDD=case_when(SCID_Diagnosis %in% 0 ~ 0L,
                     SCID_Diagnosis %in% 1:2 ~ 1L,
                     TRUE ~ NA_integer_
)) %>%
filter(!is.na(MDD)) %>%
inner_join(fam, by=c('gwas'='IID')) %>%
select(FID, IID=gwas, MDD)

write_tsv(scid_pheno, 'scid.pheno.txt')

