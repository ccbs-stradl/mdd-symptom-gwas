---
title: GenomicSEM of MDD symptom factors with other phenotypes setup
author: Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel Nivard 
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

Test genetic relationship between [symptom factors](mdd-symptom-gsem-model.md) and a selection of other phenotypes that are genetically correlated with MDD, as well as to MDD itself. Phenotypes to examine:

- Major depressive disorder: Clinical cohorts from [Wray et al](https://www.nature.com/articles/s41588-018-0090-3%5C) and all cohorts from [Howard et al](https://www.nature.com/articles/s41593-018-%200326-7). Download from [PGC](https://www.med.unc.edu/pgc/download-results/) and obtain via [data access](https://www.med.unc.edu/pgc/shared-methods/how-to/).
- bipolar disorder: [Mullins et al](https://pubmed.ncbi.nlm.nih.gov/34002096/). Download from [PGC](https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594)
- anxiety: [meta-analysis](https://drive.google.com/drive/folders/1fguHvz7l2G45sbMI9h_veQun4aXNTy1v) of [UKBB, iPSYCH](https://www.nature.com/articles/s41380-019-0559-1), and [ANGST](https://pubmed.ncbi.nlm.nih.gov/26754954/), from [Grotzinger et al medRxiv](https://www.medrxiv.org/content/10.1101/2020.09.22.20196089v1.full)
- PTSD: [Nievergelt et al](https://pubmed.ncbi.nlm.nih.gov/31594949/). Download from [PGC](https://figshare.com/articles/dataset/ptsd2019/14672133)
- tobacco use. Cigarettes per day [Liu et al](https://www.nature.com/articles/s41588-018-0307-5). Download from [UofM](https://conservancy.umn.edu/handle/11299/201564)
- alcohol dependence. [Walters et al](https://www.nature.com/articles/s41593-018-0275-1). Download from [PGC](https://doi.org/10.6084/m9.figshare.14672187)
- educational attainment. [Okbay et. al](https://www.nature.com/articles/s41588-022-01016-z). Download from [SSGAC Data Portal](https://thessgac.com).
- BMI, sex combined [Pulit et al](https://academic.oup.com/hmg/article/28/1/166/5098227). Download from [GIANT/Broad](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#2018_GIANT_and_UK_BioBank_Meta-analysis).
- neuroticism: [Nagel et al](https://www.nature.com/articles/s41588-018-0151-7). Download from [CNCR](https://ctg.cncr.nl/software/summary_statistics)
- pain: multisite chronic pain [Johnston et al](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008164). Download from [UoG](https://researchdata.gla.ac.uk/822/)
- chronotype: long sleep duration. [Dashti et al](https://www.ncbi.nlm.nih.gov/pubmed/30846698). Download from [SDKP](https://sleep.hugeamp.org/downloads.html).

# Setup

## R packages

R version


```r
R.version
```

```
##                _                           
## platform       x86_64-apple-darwin17.0     
## arch           x86_64                      
## os             darwin17.0                  
## system         x86_64, darwin17.0          
## status                                     
## major          4                           
## minor          2.1                         
## year           2022                        
## month          06                          
## day            23                          
## svn rev        82513                       
## language       R                           
## version.string R version 4.2.1 (2022-06-23)
## nickname       Funny-Looking Kid
```

Package installation


```r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

library(devtools)

if(!require(GenomicSEM)) install_github("MichelNivard/GenomicSEM")

if(!require(tidySEM)) install_github("cjvanlissa/tidySEM")
```

GenomicSEM version

```r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(corrplot)
require(GenomicSEM)

packageVersion("GenomicSEM")
```

```
## [1] '0.0.5'
```

# Process external sumstats

## Reformat

Format the external sumstats for reading by GenomicSEM with columns `SNP`, `A1` (effect allele), `A2` (non-effect allele), `BETA`/`OR`, `P`, `INFO`, and `N`.


```r
# Major depressive disorder
mdd <- read_table('sumstats/PGC_UKB_23andMe_depression_genome-wide_info_N.txt.gz')
mdd_sumstats <- mdd %>%
    mutate(Nca=UKB_Ncases+PGC_Ncases+X23andMe_Ncases,
           Nco=UKB_Ncontrols+PGC_Ncontrols+X23andMe_Ncontrols) %>%
    transmute(SNP=MarkerName, A1=toupper(Allele1), A2=toupper(Allele2),
              BETA=Effect, SE=StdErr, P=P.value, 
              FREQ=Freq1, N=4*Nca*Nco/(Nca+Nco))
write_tsv(mdd_sumstats, 'sumstats/MD.txt')

mdd_clin <- read_tsv('sumstats/daner_MDD29.0515a_mds6.0316.gz')
mdd_clin_sumstats <- mdd_clin %>%
    transmute(SNP, A1, A2, BETA=log(OR), SE, FREQ=FRQ_U_25632, INFO, P,
              N=4*Nca*Nco/(Nca+Nco))
write_tsv(mdd_clin_sumstats, 'sumstats/MDD.txt')

# Bipolar disorder
bip <- read_tsv('sumstats/pgc-bip2021-all.vcf.tsv.gz', comment='##')
bip_sumstats <- bip %>%
    filter(IMPINFO >= 0.6) %>%
    transmute(SNP=ID, A1, A2, BETA, SE, P=PVAL, INFO=IMPINFO, N=2*NEFFDIV2)
write_tsv(bip_sumstats, 'sumstats/BIP.txt')

# alcohol dependence
alcdep <- read_table('sumstats/pgc_alcdep.eur_discovery.aug2018_release.txt.gz')
alcdep_sumstats <- alcdep %>%
mutate(SNP=str_split_fixed(SNP, pattern=":", n=2)[,1]) %>%
select(SNP, A1, A2, BETA=Z, P, N=Weight)
write_tsv(alcdep_sumstats, 'sumstats/AlcDep.txt')

# body-mass index
bmi <- read_table('sumstats/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz')
bmi_sumstats <- bmi %>%
    separate(SNP, into=c('SNP', 'A1', 'A2'), sep=':') %>%
    filter(INFO >= 0.6) %>%
    select(SNP, A1=Tested_Allele, A2=Other_Allele, BETA, SE, P, INFO, N)
write_tsv(bmi_sumstats, 'sumstats/BMI.txt')

# educational attainment
ea <- read_tsv('sumstats/EA4_additive_excl_23andMe.txt.gz')
ea_sumstats <-  ea %>%
    transmute(SNP=rsID, A1=Effect_allele, A2=Other_allele, FREQ=EAF_HRC, BETA=Beta, SE_unadj, P=P_unadj, N=765283)
write_tsv(ea_sumstats, 'sumstats/EA.txt')

# multisite chronic pain
pain <- read_tsv('sumstats/chronic_pain-bgen.stats.gz')
pain_sumstats <- pain %>%
    filter(INFO >= 0.6) %>%
    transmute(SNP, A1=ALLELE1, A2=ALLELE0, BETA, SE, P=P_BOLT_LMM_INF, INFO, N=387649) 
write_tsv(pain_sumstats, 'sumstats/Pain.txt')

# tabacco use
smoking <- read_table('sumstats/CigarettesPerDay.txt.gz')
smoking_sumstats <- smoking %>%
    select(SNP=RSID, A1=ALT, A2=REF, BETA, SE, P=PVALUE, N)
write_tsv(smoking_sumstats, 'sumstats/Smoking.txt')

# chronotype (long sleep duration)
sleep <- read_table('sumstats/longsumstats.txt.gz')
sleep_sumstats <- sleep %>%
    filter(INFO >= 0.6) %>%
    transmute(SNP, A1=ALLELE1, A2=ALLELE0, BETA=BETA_LONGSLEEP, SE=SE_LONGSLEEP, P=P_LONGSLEEP, INFO, N=4*34184*(305742-34184)/305742)
write_tsv(sleep_sumstats, 'sumstats/Sleep.txt')

# anxiety disorder
anxiety <- read_table('sumstats/META_UKBB_iPSYCH_ANGST_wNcol.sumstats.gz')
anxiety_sumstats <- anxiety %>%
    mutate(Nca=25453+12655+7016, Nco=58113+19225+14745) %>%
    transmute(SNP, A1=Allele1, A2=Allele2, BETA=Effect, SE=StdErr, P, N=4*Nca*Nco/(Nca+Nco))
write_tsv(anxiety_sumstats, 'sumstats/Anxiety.txt')    

# post-traumatic stress disorder
ptsd <- read_tsv('sumstats/pts_eur_freeze2_overall.results.gz')
ptsd_sumstats <- ptsd %>%
    transmute(SNP,  A1, A2, OR, SE, P, INFO, N=4*Nca*Nco/(Nca+Nco))
write_tsv(ptsd_sumstats, 'sumstats/PTSD.txt')
    
# neuroticism
neu <- read_tsv('sumstats/sumstats_neuroticism_ctg_format.txt.gz')
neu_sumstats <- neu %>%
    filter(!is.na(INFO_UKB)) %>%
    select(SNP=RSID, A1, A2, BETA=Z, P, INFO=INFO_UKB, N)
write_tsv(neu_sumstats, 'sumstats/Neu.txt')
```

## Munge


```r
ext_traits <- c('AlcDep'=0.159, 'Anxiety'=0.16, 'BIP'=0.01,
                'BMI'=NA, 'EA'=NA, 'MD'=0.3, 'MDD'=0.15, 'Neu'=NA,
                'PTSD'=0.3, 'Pain'=NA, 'Sleep'=0.11, 'Smoking'=NA)
ext_trait_names <- names(ext_traits)
```


```r
munge(file.path('sumstats', paste(ext_trait_names, 'txt', sep='.')),
trait.names=ext_trait_names,
hm3=here::here("sumstats/reference/w_hm3.snplist"),
info.filter = 0.9, maf.filter = 0.01)
```

# Symptom labels

MDD DSM symptoms are numbered 1-9:


```r
# plot labels

dsm_mdd_symptoms_labels <-
read_delim("
MDD1;Mood;Mood;Dep
MDD2;Interest;Interest;Anh
MDD3;Weight⇅;Weight⇆;App
MDD3a;Weight⇊;Weight⇇;AppDec
MDD3b;Weight⇈;Weight⇉;AppInc
MDD4;Sleep⇅;Sleep⇆;Sle
MDD4a;Sleep⇊;Sleep⇇;SleDec
MDD4b;Sleep⇈;Sleep⇉;SleInc
MDD5;Motor⇅;Motor⇆;Moto
MDD5a;Motor⇈;Motor⇉;MotoInc
MDD5b;Motor⇊;Motor⇇;MotoDec
MDD6;Fatigue;Fatigue;Fatig
MDD7;Guilt;Guilt;Guilt
MDD8;Concentrate;Concentrate;Conc
MDD9;Suicidality;Suicidality;Sui
", col_names=c('ref', 'h', 'v', 'abbv'), delim=';')
```

```
## Rows: 15 Columns: 4
## ── Column specification ────────────────────────────────────────────────────────────────────
## Delimiter: ";"
## chr (4): ref, h, v, abbv
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
dsm_mdd_symptoms_reference <-
read_delim("
MDD1;Depressed mood most of the day, nearly every day
MDD2;Markedly diminished interest or pleasure in all, or almost all, activities most of the day, nearly every day
MDD3;Significant change in weight or appetite
MDD3a;Significant weight loss or decrease in appetite
MDD3b;Significant weight gain or increase in appetite
MDD4;Sleeping too much or not sleeping enough
MDD4a;Insomnia nearly every day
MDD4b;Hypersomnia nearly every day
MDD5;Changes in speed/amount of moving or speaking
MDD5a;Psychomotor agitation nearly every day
MDD5b;Psychomotor slowing nearly every day
MDD6;Fatigue or loss of energy nearly every day
MDD7;Feelings of worthlessness or excessive or inappropriate guilt
MDD8;Diminished ability to think or concentrate, or indecisiveness
MDD9;Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for attempting suicide
", col_names=c('Reference', 'Description'), delim=';')
```

```
## Rows: 15 Columns: 2
## ── Column specification ────────────────────────────────────────────────────────────────────
## Delimiter: ";"
## chr (2): Reference, Description
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


# Symptom prevalences

Load previously calculated symptom prevalences:


```r
all_covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
all_sumstats_prevs <- read_tsv(here::here('ldsc', paste(all_covstruct_prefix, 'prevs', 'txt', sep='.'))) 
```

```
## Rows: 38 Columns: 9
## ── Column specification ────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (5): cohorts, symptom, sumstats, filename, trait_name
## dbl (4): Nca, Nco, samp_prev, pop_prev
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


# Multivariable LDSC estimation

Calculate LDSC covariance structure for syptoms used in the combined structural model and the external phenotypes.


```r
covstruct_prefix <- 'symptoms.external.covstruct'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

if(!file.exists(covstruct_r)) {
    
  symptoms_sumstats_prevs <- all_sumstats_prevs %>%
  left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
  mutate(samp_prev=0.5,
         cohort=case_when(str_detect(filename, 'AGDS_PGC') ~ 'Clin',
                          str_detect(filename, 'ALSPAC_UKB') ~ 'Comm',
                          str_detect(filename, 'UKBt') ~ 'Ukb')) %>%
  transmute(filename=here::here(filename), samp_prev, pop_prev, trait_name=paste0(cohort, abbv))

  # external files, prevalences and trait names
  external_sumstats_prevs <-
    tibble(filename=paste(ext_trait_names, 'sumstats', 'gz', sep='.'),
           trait_name=ext_trait_names,
           pop_prev=ext_traits) %>%
    mutate(samp_prev=if_else(!is.na(pop_prev), true=0.5, false=NA_real_))
    
  sumstats_prevs <- bind_rows(symptoms_sumstats_prevs, external_sumstats_prevs)
    

  symptoms_covstruct <- ldsc(traits=sumstats_prevs$filename,
                             sample.prev=sumstats_prevs$samp_prev,
                             population.prev=sumstats_prevs$pop_prev,
                             ld=here::here('sumstats/reference/eur_w_ld_chr/'),
                             wld=here::here('sumstats/reference/eur_w_ld_chr/'),
                             trait.names=sumstats_prevs$trait_name)

  dput(symptoms_covstruct, covstruct_r, control=c('exact'))
  saveRDS(symptoms_covstruct, covstruct_rds)
  
  # check for exact match of deparsed object
  identical(dget(covstruct_r), symptoms_covstruct)

} else {

  symptoms_covstruct <- dget(covstruct_r)

}
```

# Models

## Symptom factors

Base model of symptom factors, with Melancholic and Atypical. The mood symptoms factor was highly correlated with the Melancholic factor, so those symptoms have been combined. 


```r
model <- "
DEP =~ NA*CommAnh + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommGuilt + CommDep + UkbDep + ClinSui + CommConc + CommSui 
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

DEP ~~ 1*DEP
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*DEP + 0*ATY
"
fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##   1.441 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355035 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.63713097734563 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

## External phenotypes

Compare symptom factors against each external phenotype. Single regression of each external phenotype on each symptom or symptom cluster.


```r
ext.glue <- "
DEP =~ NA*CommAnh + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommGuilt + CommDep + UkbDep + ClinSui + CommConc + CommSui 
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

DEP ~~ 1*DEP
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*DEP + 0*ATY

AlcDep ~ {symptom}
Anxiety ~  {symptom}
BIP ~ {symptom}
BMI ~ {symptom}
EA ~ {symptom}
MD ~ {symptom}
MDD ~ {symptom}
Neu ~ {symptom}
PTSD ~ {symptom}
Pain ~ {symptom}
Sleep ~ {symptom}
Smoking ~ {symptom}
"

ext.model_list <- lapply(c('DEP', 'ATY', 'GATE'), function(symptom) str_glue_data(list(symptom=symptom), ext.glue))

ext.fit_list <- lapply(ext.model_list, function(model) usermodel(symptoms_covstruct, estimation='DWLS', model=model))
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  23.706 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058033 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  18.418 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058033 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  20.863 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058033 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model): A
## difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```


```r
ext_full <-
bind_rows(lapply(ext.fit_list, function(fit) mutate(fit$results, p_value=as.character(p_value)))) %>%
select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
filter(lhs %in% ext_trait_names, rhs %in% c('DEP', 'ATY', 'GATE')) %>%
mutate(Beta='Total', Factor=rhs, Phenotype=lhs, p_value=as.numeric(p_value))
```

Multiple regression of each phenotype on symptom factors, to estimate relationship after condition on each of the other factors. 



```r
ext_mult.model <- "
DEP =~ NA*CommAnh + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommGuilt + CommDep + UkbDep + ClinSui + CommConc + CommSui 
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

DEP ~~ 1*DEP
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*DEP + 0*ATY

AlcDep ~ DEP + ATY + GATE
Anxiety ~  DEP + ATY + GATE
BIP ~ DEP + ATY + GATE
BMI ~ DEP + ATY + GATE
EA ~ DEP + ATY + GATE
MD ~ DEP + ATY + GATE
MDD ~ DEP + ATY + GATE
Neu ~ DEP + ATY + GATE
PTSD ~ DEP + ATY + GATE
Pain ~ DEP + ATY + GATE
Sleep ~ DEP + ATY + GATE
Smoking ~ DEP + ATY + GATE
"
ext_mult.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=ext_mult.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  23.814 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058033 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_mult.model): A difference greater than .025 was observed pre- and
## post-smoothing in the genetic covariance matrix. This reflects a large
## difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_mult.model): A difference greater than .025 was observed pre- and
## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
## a large difference and results should be interpreted with caution!! This can
## often result from including low powered traits, and you might consider removing
## those traits from the model. If you are going to run a multivariate GWAS we
## strongly recommend setting the smooth_check argument to true to check smoothing
## for each SNP.
```


```r
ext_partial <-
bind_rows(lapply(list(ext_mult.fit),
                 function(fit) fit$results)) %>%
  select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
  filter(lhs %in% ext_trait_names, rhs %in% c('DEP', 'ATY', 'GATE')) %>%
  mutate(Beta='Specific', Factor=rhs, Phenotype=lhs)
```


```r
bind_rows(ext_full, ext_partial) %>%
select(Factor, Phenotype, Model=Beta, STD_Genotype, STD_Genotype_SE, p_value) %>%
arrange(Factor, Phenotype, Model)
```

<div class="kable-table">

|Factor |Phenotype |Model    | STD_Genotype|STD_Genotype_SE    |   p_value|
|:------|:---------|:--------|------------:|:------------------|---------:|
|ATY    |AlcDep    |Specific |    0.3640767|0.254228487814086  | 0.1521144|
|ATY    |AlcDep    |Total    |    0.5005632|0.0801365998233833 | 0.0000000|
|ATY    |Anxiety   |Specific |   -0.1002595|0.15789788053738   | 0.5254618|
|ATY    |Anxiety   |Total    |    0.7461638|0.0610498972634396 | 0.0000000|
|ATY    |BIP       |Specific |   -0.1925331|0.151262513376248  | 0.2030751|
|ATY    |BIP       |Total    |    0.5253496|0.0433068909754881 | 0.0000000|
|ATY    |BMI       |Specific |    1.1769822|0.36900876973096   | 0.0014241|
|ATY    |BMI       |Total    |    0.2236476|0.0293760349283482 | 0.0000000|
|ATY    |EA        |Specific |   -0.4022932|0.13930888661245   | 0.0038787|
|ATY    |EA        |Total    |   -0.1059923|0.0303955945887964 | 0.0004883|
|ATY    |MD        |Specific |    0.1362367|0.118657902329921  | 0.2508930|
|ATY    |MD        |Total    |    0.9337764|0.0553271353460975 | 0.0000000|
|ATY    |MDD       |Specific |    0.1932694|0.224332429760117  | 0.3889392|
|ATY    |MDD       |Total    |    0.8343727|0.0844201722693633 | 0.0000000|
|ATY    |Neu       |Specific |   -0.1438345|0.157369633342377  | 0.3607263|
|ATY    |Neu       |Total    |    0.6951237|0.0475181114949475 | 0.0000000|
|ATY    |PTSD      |Specific |    0.1306751|0.204823292306298  | 0.5234703|
|ATY    |PTSD      |Total    |    0.7948335|0.0761178908892611 | 0.0000000|
|ATY    |Pain      |Specific |    0.6255118|0.193529566742197  | 0.0012282|
|ATY    |Pain      |Total    |    0.6028382|0.0446756382218185 | 0.0000000|
|ATY    |Sleep     |Specific |    0.6465696|0.243383591610749  | 0.0078918|
|ATY    |Sleep     |Total    |    0.3062883|0.0495371839928435 | 0.0000000|
|ATY    |Smoking   |Specific |    0.6431627|0.23347309145044   | 0.0058719|
|ATY    |Smoking   |Total    |    0.2991301|0.039228791444688  | 0.0000000|
|DEP    |AlcDep    |Specific |    0.1885050|0.231210188105455  | 0.4148641|
|DEP    |AlcDep    |Total    |    0.4737221|0.0729119638494336 | 0.0000000|
|DEP    |Anxiety   |Specific |    0.6166142|0.144183584948087  | 0.0000190|
|DEP    |Anxiety   |Total    |    0.7085221|0.0504659271394381 | 0.0000000|
|DEP    |BIP       |Specific |    0.5534770|0.140238803314301  | 0.0000792|
|DEP    |BIP       |Total    |    0.4993681|0.0364951812168613 | 0.0000000|
|DEP    |BMI       |Specific |   -0.5753820|0.359698684262476  | 0.1096769|
|DEP    |BMI       |Total    |    0.2096220|0.027165906311088  | 0.0000000|
|DEP    |EA        |Specific |    0.0031493|0.133246052921943  | 0.9812037|
|DEP    |EA        |Total    |   -0.0994529|0.0286656998615448 | 0.0005216|
|DEP    |MD        |Specific |    0.5660656|0.109409228135781  | 0.0000002|
|DEP    |MD        |Total    |    0.8861225|0.0390501161886385 | 0.0000000|
|DEP    |MDD       |Specific |    0.4055509|0.209299460427074  | 0.0526553|
|DEP    |MDD       |Total    |    0.7912590|0.072414808169817  | 0.0000000|
|DEP    |Neu       |Specific |    0.7589834|0.142804349095257  | 0.0000001|
|DEP    |Neu       |Total    |    0.6604487|0.0370095039478862 | 0.0000000|
|DEP    |PTSD      |Specific |    0.5261640|0.179137179821488  | 0.0033107|
|DEP    |PTSD      |Total    |    0.7537116|0.064209039992394  | 0.0000000|
|DEP    |Pain      |Specific |    0.0876166|0.193299649274597  | 0.6502778|
|DEP    |Pain      |Total    |    0.5702388|0.0366128094171443 | 0.0000000|
|DEP    |Sleep     |Specific |   -0.2822782|0.236718610214829  | 0.2330832|
|DEP    |Sleep     |Total    |    0.2880514|0.0472773039715776 | 0.0000000|
|DEP    |Smoking   |Specific |   -0.2684021|0.224016527027553  | 0.2308656|
|DEP    |Smoking   |Total    |    0.2820754|0.0363413353422488 | 0.0000000|
|GATE   |AlcDep    |Specific |    0.0719889|0.134107865360387  | 0.5914075|
|GATE   |AlcDep    |Total    |    0.5775986|0.137025033250661  | 0.0000249|
|GATE   |Anxiety   |Specific |    0.3895213|0.0865468839695091 | 0.0000068|
|GATE   |Anxiety   |Total    |    1.0186390|0.141087258397134  | 0.0000000|
|GATE   |BIP       |Specific |    0.2203518|0.0760515302325504 | 0.0037629|
|GATE   |BIP       |Total    |    0.6907365|0.104545146671682  | 0.0000000|
|GATE   |BMI       |Specific |   -0.1092501|0.0579927773313639 | 0.0595845|
|GATE   |BMI       |Total    |    0.1629996|0.0465424467425572 | 0.0004615|
|GATE   |EA        |Specific |    0.3189161|0.074797933480592  | 0.0000201|
|GATE   |EA        |Total    |    0.0076982|0.0423127806802913 | 0.8556408|
|GATE   |MD        |Specific |    0.5070460|0.0743490245529634 | 0.0000000|
|GATE   |MD        |Total    |    1.2682565|0.167656586562129  | 0.0000000|
|GATE   |MDD       |Specific |    0.5092440|0.134643957082049  | 0.0001555|
|GATE   |MDD       |Total    |    1.1340963|0.177589744720942  | 0.0000000|
|GATE   |Neu       |Specific |    0.0924421|0.0899956413225501 | 0.3043413|
|GATE   |Neu       |Total    |    0.8423175|0.120011235300734  | 0.0000000|
|GATE   |PTSD      |Specific |    0.3205619|0.117818043669866  | 0.0065124|
|GATE   |PTSD      |Total    |    1.0304451|0.160874779070875  | 0.0000000|
|GATE   |Pain      |Specific |    0.0900431|0.0758069477957605 | 0.2349187|
|GATE   |Pain      |Total    |    0.6810387|0.0993562767980237 | 0.0000000|
|GATE   |Sleep     |Specific |    0.1469529|0.088841538295064  | 0.0981105|
|GATE   |Sleep     |Total    |    0.3400320|0.0826835343196338 | 0.0000391|
|GATE   |Smoking   |Specific |    0.1616602|0.0635913333922893 | 0.0110170|
|GATE   |Smoking   |Total    |    0.3702149|0.0694564471208399 | 0.0000001|

</div>


```r
ggplot(bind_rows(ext_full, ext_partial),
       aes(x=factor(Factor, levels=c('DEP', 'ATY', 'GATE')),
           y=STD_Genotype,
           color=factor(Beta, levels=c('Specific', 'Total')),
           shape=factor(Beta, levels=c('Specific', 'Total')),
          ymin=qnorm(0.025, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)),
          ymax=qnorm(0.975, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)))) +
geom_hline(yintercept=0, col='gray') +
geom_pointrange(position=position_dodge(width=0.5)) +
facet_wrap(~Phenotype) +
scale_x_discrete('Factor') +
scale_y_continuous(expression(r[g]), breaks=c(-1, 0, 1)) +
scale_color_discrete('Correlation: ') +
coord_flip(ylim=c(-1, 1)) +
theme_bw() +
theme(axis.text.y=element_text(size=16),
      strip.text=element_text(size=16),
      legend.title=element_text(size=12),
      legend.text=element_text(size=14),
      legend.position='top') +
labs(color  = "Correlation: ", shape = "Correlation: ")
```

![](mdd-symptom-gsem-ext_files/figure-html/ex_plot-1.png)<!-- -->

Test for attenuation of correlations


```r
ext_factor_wide <-
bind_rows(ext_full, ext_partial) %>%
as_tibble() %>%
select(Factor, Phenotype, STD_Genotype, STD_Genotype_SE, p_value, Beta) %>%
pivot_wider(id_cols=c(Phenotype, Factor),
            names_from=Beta, 
            values_from=c(STD_Genotype, STD_Genotype_SE, p_value))

ext_attenuation <-     
ext_factor_wide %>%
mutate(rg_Total=STD_Genotype_Total, rg_Specific=STD_Genotype_Specific,
       se_Total=as.numeric(STD_Genotype_SE_Total),
       se_Specific=as.numeric(STD_Genotype_SE_Specific)) %>%
mutate(attenuation_z=(rg_Total-rg_Specific)/sqrt(se_Total^2+se_Specific^2)) %>%
mutate(attenuation_p=2*pnorm(abs(attenuation_z), lower.tail=F)) %>%
select(Phenotype, Factor, p_value_Total, p_value_Specific, attenuation_p)

knitr::kable(ext_attenuation)
```



|Phenotype |Factor | p_value_Total| p_value_Specific| attenuation_p|
|:---------|:------|-------------:|----------------:|-------------:|
|AlcDep    |DEP    |     0.0000000|        0.4148641|     0.2394062|
|Anxiety   |DEP    |     0.0000000|        0.0000190|     0.5474086|
|BIP       |DEP    |     0.0000000|        0.0000792|     0.7088530|
|BMI       |DEP    |     0.0000000|        0.1096769|     0.0295406|
|EA        |DEP    |     0.0005216|        0.9812037|     0.4515719|
|MD        |DEP    |     0.0000000|        0.0000002|     0.0058676|
|MDD       |DEP    |     0.0000000|        0.0526553|     0.0815854|
|Neu       |DEP    |     0.0000000|        0.0000001|     0.5041776|
|PTSD      |DEP    |     0.0000000|        0.0033107|     0.2317940|
|Pain      |DEP    |     0.0000000|        0.6502778|     0.0141615|
|Sleep     |DEP    |     0.0000000|        0.2330832|     0.0181446|
|Smoking   |DEP    |     0.0000000|        0.2308656|     0.0152832|
|AlcDep    |ATY    |     0.0000000|        0.1521144|     0.6086301|
|Anxiety   |ATY    |     0.0000000|        0.5254618|     0.0000006|
|BIP       |ATY    |     0.0000000|        0.2030751|     0.0000051|
|BMI       |ATY    |     0.0000000|        0.0014241|     0.0100138|
|EA        |ATY    |     0.0004883|        0.0038787|     0.0377052|
|MD        |ATY    |     0.0000000|        0.2508930|     0.0000000|
|MDD       |ATY    |     0.0000000|        0.3889392|     0.0074795|
|Neu       |ATY    |     0.0000000|        0.3607263|     0.0000003|
|PTSD      |ATY    |     0.0000000|        0.5234703|     0.0023698|
|Pain      |ATY    |     0.0000000|        0.0012282|     0.9091142|
|Sleep     |ATY    |     0.0000000|        0.0078918|     0.1706751|
|Smoking   |ATY    |     0.0000000|        0.0058719|     0.1461757|
|AlcDep    |GATE   |     0.0000249|        0.5914075|     0.0083624|
|Anxiety   |GATE   |     0.0000000|        0.0000068|     0.0001442|
|BIP       |GATE   |     0.0000000|        0.0037629|     0.0002743|
|BMI       |GATE   |     0.0004615|        0.0595845|     0.0002510|
|EA        |GATE   |     0.8556408|        0.0000201|     0.0002929|
|MD        |GATE   |     0.0000000|        0.0000000|     0.0000332|
|MDD       |GATE   |     0.0000000|        0.0001555|     0.0050508|
|Neu       |GATE   |     0.0000000|        0.3043413|     0.0000006|
|PTSD      |GATE   |     0.0000000|        0.0065124|     0.0003708|
|Pain      |GATE   |     0.0000000|        0.2349187|     0.0000023|
|Sleep     |GATE   |     0.0000391|        0.0981105|     0.1116319|
|Smoking   |GATE   |     0.0000001|        0.0110170|     0.0267838|

