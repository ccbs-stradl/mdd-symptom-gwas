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
## platform       aarch64-apple-darwin20      
## arch           aarch64                     
## os             darwin20                    
## system         aarch64, darwin20           
## status                                     
## major          4                           
## minor          2.2                         
## year           2022                        
## month          10                          
## day            31                          
## svn rev        83211                       
## language       R                           
## version.string R version 4.2.2 (2022-10-31)
## nickname       Innocent and Trusting
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
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
MDD9;Recurrent thoughts of death or suicide or a suicide attempt or a Multiple plan for attempting suicide
", col_names=c('Reference', 'Description'), delim=';')
```

```
## Rows: 15 Columns: 2
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
## ── Column specification ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=model, CFIcalc=TRUE)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##   1.201 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355036 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.63713097734562 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model, :
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model, :
## A difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

## Sample factors

Base model of Clinical and Community sample factors. 


```r
model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
GATE ~~ 0*COMM + 0*CLIN
"
fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=model, CFIcalc=TRUE)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##   1.425 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355036 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.63713097734562 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model, :
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model, :
## A difference greater than .025 was observed pre- and post-smoothing for
## Z-statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

## External phenotypes

### Symptom factors

Compare symptom factors against each external phenotype. Single regression of each external phenotype on each symptom or symptom cluster.


```r
ext_symp.glue <- "
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

ext_symp.model_list <- lapply(c('DEP', 'ATY', 'GATE'), function(symptom) str_glue_data(list(symptom=symptom), ext_symp.glue))

ext_symp.fit_list <- lapply(ext_symp.model_list, function(model) usermodel(symptoms_covstruct, estimation='DWLS', model=model))
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  22.319 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
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
##  17.236 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
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
##  19.836 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
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

### Sample factors


```r
ext_samp.glue <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
GATE ~~ 0*COMM + 0*CLIN

cc < 0.99
CLIN ~~ cc*COMM

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

ext_samp.model_list <- lapply(c('CLIN', 'COMM'), function(symptom) str_glue_data(list(symptom=symptom), ext_samp.glue))

ext_samp.fit_list <- lapply(ext_samp.model_list, function(model) usermodel(symptoms_covstruct, estimation='DWLS', model=model))
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  91.023 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
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
##  76.325 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
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
ext_single <-
bind_rows(lapply(c(ext_symp.fit_list, ext_samp.fit_list), function(fit) mutate(fit$results, p_value=as.character(p_value)))) %>%
select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
filter(lhs %in% ext_trait_names, rhs %in% c('DEP', 'ATY', 'GATE', 'CLIN', 'COMM')) %>%
mutate(Beta='Single', Factor=rhs, Phenotype=lhs, p_value=as.numeric(p_value))
```

### Symptom factors

Multiple regression of each phenotype on symptom factors, to estimate relationship after condition on each of the other factors. 



```r
ext_mult_symp.model <- "
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
ext_mult_symp.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=ext_mult_symp.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  22.932 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_mult_symp.model): A difference greater than .025 was observed pre- and
## post-smoothing in the genetic covariance matrix. This reflects a large
## difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_mult_symp.model): A difference greater than .025 was observed pre- and
## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
## a large difference and results should be interpreted with caution!! This can
## often result from including low powered traits, and you might consider removing
## those traits from the model. If you are going to run a multivariate GWAS we
## strongly recommend setting the smooth_check argument to true to check smoothing
## for each SNP.
```

### Sample factors



```r
ext_samp_mult.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
GATE ~~ 0*COMM + 0*CLIN

AlcDep ~ CLIN + COMM 
Anxiety ~  CLIN + COMM
BIP ~ CLIN + COMM 
BMI ~ CLIN + COMM 
EA ~ CLIN + COMM 
MD ~ CLIN + COMM
MDD ~ CLIN + COMM
Neu ~ CLIN + COMM
PTSD ~ CLIN + COMM
Pain ~ CLIN + COMM
Sleep ~ CLIN + COMM
Smoking ~ CLIN + COMM
"
ext_samp_mult.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=ext_samp_mult.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  41.639 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0499086254058032 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.97466132311611 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_samp_mult.model): A difference greater than .025 was observed pre- and
## post-smoothing in the genetic covariance matrix. This reflects a large
## difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## ext_samp_mult.model): A difference greater than .025 was observed pre- and
## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
## a large difference and results should be interpreted with caution!! This can
## often result from including low powered traits, and you might consider removing
## those traits from the model. If you are going to run a multivariate GWAS we
## strongly recommend setting the smooth_check argument to true to check smoothing
## for each SNP.
```


```r
ext_multiple <-
bind_rows(lapply(list(ext_mult_symp.fit, ext_samp_mult.fit),
                 function(fit) fit$results)) %>%
  select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
  filter(lhs %in% ext_trait_names, rhs %in% c('DEP', 'ATY', 'GATE', 'CLIN', 'COMM')) %>%
  mutate(Beta='Multiple', Factor=rhs, Phenotype=lhs)
```


```r
bind_rows(ext_single, ext_multiple) %>%
select(Factor, Phenotype, Model=Beta, STD_Genotype, STD_Genotype_SE, p_value) %>%
arrange(Factor, Phenotype, Model)
```

<div class="kable-table">

|Factor |Phenotype |Model    | STD_Genotype|STD_Genotype_SE    |   p_value|
|:------|:---------|:--------|------------:|:------------------|---------:|
|ATY    |AlcDep    |Multiple |    0.3640712|0.254223785791238  | 0.1521133|
|ATY    |AlcDep    |Single   |    0.5005631|0.0801366140442318 | 0.0000000|
|ATY    |Anxiety   |Multiple |   -0.1002574|0.157895205544124  | 0.5254606|
|ATY    |Anxiety   |Single   |    0.7461647|0.061049930775689  | 0.0000000|
|ATY    |BIP       |Multiple |   -0.1925292|0.151259371408457  | 0.2030753|
|ATY    |BIP       |Single   |    0.5253500|0.0433069091428086 | 0.0000000|
|ATY    |BMI       |Multiple |    1.1769646|0.368995219662611  | 0.0014241|
|ATY    |BMI       |Single   |    0.2236477|0.0293760425533776 | 0.0000000|
|ATY    |EA        |Multiple |   -0.4022878|0.139304956948227  | 0.0038787|
|ATY    |EA        |Single   |   -0.1059925|0.0303956020755315 | 0.0004883|
|ATY    |MD        |Multiple |    0.1362361|0.118656132658375  | 0.2508968|
|ATY    |MD        |Single   |    0.9337767|0.0553271604614671 | 0.0000000|
|ATY    |MDD       |Multiple |    0.1932670|0.224328961589488  | 0.3889393|
|ATY    |MDD       |Single   |    0.8343739|0.084420219045379  | 0.0000000|
|ATY    |Neu       |Multiple |   -0.1438306|0.157366347646865  | 0.3607231|
|ATY    |Neu       |Single   |    0.6951241|0.04751813414257   | 0.0000000|
|ATY    |PTSD      |Multiple |    0.1306745|0.204820296419978  | 0.5234706|
|ATY    |PTSD      |Single   |    0.7948341|0.0761179204369635 | 0.0000000|
|ATY    |Pain      |Multiple |    0.6255035|0.193523306355692  | 0.0012282|
|ATY    |Pain      |Single   |    0.6028379|0.044675637984252  | 0.0000000|
|ATY    |Sleep     |Multiple |    0.6465586|0.243375939792957  | 0.0078918|
|ATY    |Sleep     |Single   |    0.3062879|0.0495371868573943 | 0.0000000|
|ATY    |Smoking   |Multiple |    0.6431527|0.233465617055737  | 0.0058718|
|ATY    |Smoking   |Single   |    0.2991304|0.039228804765846  | 0.0000000|
|CLIN   |AlcDep    |Multiple |    0.5539985|0.451475966420007  | 0.2197182|
|CLIN   |AlcDep    |Single   |    0.4770973|0.0822308805109108 | 0.0000000|
|CLIN   |Anxiety   |Multiple |   -0.4473043|0.395788775975659  | 0.2583376|
|CLIN   |Anxiety   |Single   |    0.7129587|0.0791675441466363 | 0.0000000|
|CLIN   |BIP       |Multiple |   -0.4511859|0.377105431259951  | 0.2314438|
|CLIN   |BIP       |Single   |    0.5024973|0.0558807401700145 | 0.0000000|
|CLIN   |BMI       |Multiple |    1.1754744|0.7937046512339    | 0.1385062|
|CLIN   |BMI       |Single   |    0.2113343|0.0327028634246196 | 0.0000000|
|CLIN   |EA        |Multiple |   -0.8634334|0.585694284005175  | 0.1403267|
|CLIN   |EA        |Single   |   -0.1003859|0.0291786839273907 | 0.0005809|
|CLIN   |MD        |Multiple |   -0.3147681|0.320838492341552  | 0.3265085|
|CLIN   |MD        |Single   |    0.8916358|0.0842683837826528 | 0.0000000|
|CLIN   |MDD       |Multiple |   -0.5264491|0.511747314532823  | 0.3035501|
|CLIN   |MDD       |Single   |    0.7963982|0.0961948348056698 | 0.0000000|
|CLIN   |Neu       |Multiple |   -0.0881142|0.195992511987339  | 0.6530585|
|CLIN   |Neu       |Single   |    0.6645534|0.0671021281930817 | 0.0000000|
|CLIN   |PTSD      |Multiple |   -0.2234120|0.34801383869207   | 0.5209005|
|CLIN   |PTSD      |Single   |    0.7587889|0.0971739423320513 | 0.0000000|
|CLIN   |Pain      |Multiple |    0.3942759|0.294699439656623  | 0.1808442|
|CLIN   |Pain      |Single   |    0.5742093|0.0608871788983916 | 0.0000000|
|CLIN   |Sleep     |Multiple |    0.1801937|0.233665669666763  | 0.4405708|
|CLIN   |Sleep     |Single   |    0.2903798|0.0547082088136701 | 0.0000001|
|CLIN   |Smoking   |Multiple |    0.6438998|0.447983662807281  | 0.1505291|
|CLIN   |Smoking   |Single   |    0.2841731|0.0440562581722677 | 0.0000000|
|COMM   |AlcDep    |Multiple |    0.1208426|0.370476961320492  | 0.7440459|
|COMM   |AlcDep    |Single   |    0.4727444|0.072715189272104  | 0.0000000|
|COMM   |Anxiety   |Multiple |    0.9872020|0.334212763993674  | 0.0031289|
|COMM   |Anxiety   |Single   |    0.7067274|0.0502600320551331 | 0.0000000|
|COMM   |BIP       |Multiple |    0.7762264|0.324816025541961  | 0.0168257|
|COMM   |BIP       |Single   |    0.4980530|0.0363035201204751 | 0.0000000|
|COMM   |BMI       |Multiple |   -0.5128830|0.722403664125655  | 0.4776757|
|COMM   |BMI       |Single   |    0.2094419|0.027101582747042  | 0.0000000|
|COMM   |EA        |Multiple |    0.4417557|0.534646529525153  | 0.4085868|
|COMM   |EA        |Single   |   -0.0993831|0.028591720855334  | 0.0005091|
|COMM   |MD        |Multiple |    1.0876584|0.263989857042839  | 0.0000376|
|COMM   |MD        |Single   |    0.8838678|0.0387631464950772 | 0.0000000|
|COMM   |MDD       |Multiple |    1.1230484|0.425872159719764  | 0.0083437|
|COMM   |MDD       |Single   |    0.7894503|0.0721919040779266 | 0.0000000|
|COMM   |Neu       |Multiple |    0.7164269|0.141036769670343  | 0.0000004|
|COMM   |Neu       |Single   |    0.6587390|0.036818313492786  | 0.0000000|
|COMM   |PTSD      |Multiple |    0.8913235|0.252573458582583  | 0.0004157|
|COMM   |PTSD      |Single   |    0.7520557|0.0640070272834876 | 0.0000000|
|COMM   |Pain      |Multiple |    0.3269622|0.250509815204117  | 0.1915965|
|COMM   |Pain      |Single   |    0.5691436|0.0364114685612383 | 0.0000000|
|COMM   |Sleep     |Multiple |    0.1733625|0.171980850568423  | 0.3132514|
|COMM   |Sleep     |Single   |    0.2877816|0.0470993019805111 | 0.0000000|
|COMM   |Smoking   |Multiple |   -0.1162823|0.401478185650852  | 0.7721715|
|COMM   |Smoking   |Single   |    0.2816045|0.0362003139181727 | 0.0000000|
|DEP    |AlcDep    |Multiple |    0.1885105|0.231205774371916  | 0.4148665|
|DEP    |AlcDep    |Single   |    0.4737222|0.0729119786423913 | 0.0000000|
|DEP    |Anxiety   |Multiple |    0.6166120|0.144181255958035  | 0.0000190|
|DEP    |Anxiety   |Single   |    0.7085225|0.0504659378980927 | 0.0000000|
|DEP    |BIP       |Multiple |    0.5534740|0.140235803335444  | 0.0000792|
|DEP    |BIP       |Single   |    0.4993684|0.0364951938484106 | 0.0000000|
|DEP    |BMI       |Multiple |   -0.5753649|0.359685270470339  | 0.1096768|
|DEP    |BMI       |Single   |    0.2096219|0.027165910520389  | 0.0000000|
|DEP    |EA        |Multiple |    0.0031440|0.13324209758816   | 0.9812037|
|DEP    |EA        |Single   |   -0.0994526|0.0286657057514955 | 0.0005216|
|DEP    |MD        |Multiple |    0.5660664|0.109407664287715  | 0.0000002|
|DEP    |MD        |Single   |    0.8861226|0.0390501287043192 | 0.0000000|
|DEP    |MDD       |Multiple |    0.4055531|0.209296351426157  | 0.0526557|
|DEP    |MDD       |Single   |    0.7912586|0.0724148184158745 | 0.0000000|
|DEP    |Neu       |Multiple |    0.7589806|0.1428013375112    | 0.0000001|
|DEP    |Neu       |Single   |    0.6604488|0.0370095138762591 | 0.0000000|
|DEP    |PTSD      |Multiple |    0.5261650|0.179134534455853  | 0.0033107|
|DEP    |PTSD      |Single   |    0.7537116|0.0642090506830633 | 0.0000000|
|DEP    |Pain      |Multiple |    0.0876251|0.193293420557446  | 0.6502747|
|DEP    |Pain      |Single   |    0.5702388|0.0366128165925743 | 0.0000000|
|DEP    |Sleep     |Multiple |   -0.2822677|0.236711060933285  | 0.2330837|
|DEP    |Sleep     |Single   |    0.2880514|0.047277313201574  | 0.0000000|
|DEP    |Smoking   |Multiple |   -0.2683926|0.224009138006028  | 0.2308655|
|DEP    |Smoking   |Single   |    0.2820754|0.0363413425996212 | 0.0000000|
|GATE   |AlcDep    |Multiple |    0.0719891|0.134107967929391  | 0.5914044|
|GATE   |AlcDep    |Single   |    0.5776047|0.13702762426172   | 0.0000249|
|GATE   |Anxiety   |Multiple |    0.3895213|0.086546921738715  | 0.0000068|
|GATE   |Anxiety   |Single   |    1.0186501|0.141091949002172  | 0.0000000|
|GATE   |BIP       |Multiple |    0.2203511|0.0760516143009265 | 0.0037628|
|GATE   |BIP       |Single   |    0.6907434|0.104548297042388  | 0.0000000|
|GATE   |BMI       |Multiple |   -0.1092503|0.0579928842170064 | 0.0595848|
|GATE   |BMI       |Single   |    0.1630011|0.046543171653155  | 0.0004615|
|GATE   |EA        |Multiple |    0.3189160|0.0747980148618434 | 0.0000201|
|GATE   |EA        |Single   |    0.0076979|0.0423132101001359 | 0.8556397|
|GATE   |MD        |Multiple |    0.5070457|0.0743490668763843 | 0.0000000|
|GATE   |MD        |Single   |    1.2682699|0.167662730177194  | 0.0000000|
|GATE   |MDD       |Multiple |    0.5092444|0.134644047312628  | 0.0001555|
|GATE   |MDD       |Single   |    1.1341098|0.177595017100915  | 0.0000000|
|GATE   |Neu       |Multiple |    0.0924412|0.0899958325231919 | 0.3043390|
|GATE   |Neu       |Single   |    0.8423266|0.120015254657171  | 0.0000000|
|GATE   |PTSD      |Multiple |    0.3205616|0.117818130546624  | 0.0065124|
|GATE   |PTSD      |Single   |    1.0304568|0.160879568321433  | 0.0000000|
|GATE   |Pain      |Multiple |    0.0900426|0.0758070755318798 | 0.2349191|
|GATE   |Pain      |Single   |    0.6810453|0.0993593432370933 | 0.0000000|
|GATE   |Sleep     |Multiple |    0.1469520|0.0888416030424469 | 0.0981094|
|GATE   |Sleep     |Single   |    0.3400354|0.0826849166197801 | 0.0000392|
|GATE   |Smoking   |Multiple |    0.1616598|0.0635913767023067 | 0.0110170|
|GATE   |Smoking   |Single   |    0.3702171|0.0694578418999869 | 0.0000001|

</div>


```r
ggplot(bind_rows(ext_multiple, ext_single),
       aes(x=factor(Factor,
	   				levels=c('DEP', 'ATY', 'GATE', 'COMM', 'CLIN'),
					   labels=c("Depression", "Atypical", "Gating", "Community", "Clinical")),
           y=STD_Genotype,
           color=factor(Beta, levels=c('Single', 'Multiple')),
           shape=factor(Beta, levels=c('Single', 'Multiple')),
          ymin=qnorm(0.025, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)),
          ymax=qnorm(0.975, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)))) +
geom_hline(yintercept=0, col='gray') +
geom_pointrange(position=position_dodge(width=-0.5)) +
facet_wrap(~Phenotype) +
scale_x_discrete('Factor') +
scale_y_continuous('Beta') +
scale_color_discrete('Regression: ') +
coord_flip() +
theme_bw() +
theme(axis.text.y=element_text(size=16),
      strip.text=element_text(size=16),
      legend.title=element_text(size=12),
      legend.text=element_text(size=14),
      legend.position='top') +
labs(color  = "Regression: ", shape = "Regression: ")
```

```
## Warning: `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
## `position_dodge()` requires non-overlapping x intervals
```

![](mdd-symptom-gsem-ext_files/figure-html/ex_plot-1.png)<!-- -->

Test for attenuation of correlations


```r
ext_factor_wide <-
bind_rows(ext_single, ext_multiple) %>%
as_tibble() %>%
select(Factor, Phenotype, STD_Genotype, STD_Genotype_SE, p_value, Beta) %>%
pivot_wider(id_cols=c(Phenotype, Factor),
            names_from=Beta, 
            values_from=c(STD_Genotype, STD_Genotype_SE, p_value))

ext_attenuation <-     
ext_factor_wide %>%
mutate(rg_Single=STD_Genotype_Single, rg_Multiple=STD_Genotype_Multiple,
       se_Single=as.numeric(STD_Genotype_SE_Single),
       se_Multiple=as.numeric(STD_Genotype_SE_Multiple)) %>%
mutate(attenuation_z=(rg_Single-rg_Multiple)/sqrt(se_Single^2+se_Multiple^2)) %>%
mutate(attenuation_p=2*pnorm(abs(attenuation_z), lower.tail=F)) %>%
select(Phenotype, Factor, p_value_Single, p_value_Multiple, attenuation_p)

knitr::kable(ext_attenuation)
```



|Phenotype |Factor | p_value_Single| p_value_Multiple| attenuation_p|
|:---------|:------|--------------:|----------------:|-------------:|
|AlcDep    |DEP    |      0.0000000|        0.4148665|     0.2394069|
|Anxiety   |DEP    |      0.0000000|        0.0000190|     0.5473917|
|BIP       |DEP    |      0.0000000|        0.0000792|     0.7088645|
|BMI       |DEP    |      0.0000000|        0.1096768|     0.0295381|
|EA        |DEP    |      0.0005216|        0.9812037|     0.4515836|
|MD        |DEP    |      0.0000000|        0.0000002|     0.0058671|
|MDD       |DEP    |      0.0000000|        0.0526557|     0.0815835|
|Neu       |DEP    |      0.0000000|        0.0000001|     0.5041817|
|PTSD      |DEP    |      0.0000000|        0.0033107|     0.2317900|
|Pain      |DEP    |      0.0000000|        0.6502747|     0.0141602|
|Sleep     |DEP    |      0.0000000|        0.2330837|     0.0181432|
|Smoking   |DEP    |      0.0000000|        0.2308655|     0.0152817|
|AlcDep    |ATY    |      0.0000000|        0.1521133|     0.6086098|
|Anxiety   |ATY    |      0.0000000|        0.5254606|     0.0000006|
|BIP       |ATY    |      0.0000000|        0.2030753|     0.0000051|
|BMI       |ATY    |      0.0000000|        0.0014241|     0.0100124|
|EA        |ATY    |      0.0004883|        0.0038787|     0.0377036|
|MD        |ATY    |      0.0000000|        0.2508968|     0.0000000|
|MDD       |ATY    |      0.0000000|        0.3889393|     0.0074783|
|Neu       |ATY    |      0.0000000|        0.3607231|     0.0000003|
|PTSD      |ATY    |      0.0000000|        0.5234706|     0.0023694|
|Pain      |ATY    |      0.0000000|        0.0012282|     0.9091432|
|Sleep     |ATY    |      0.0000000|        0.0078918|     0.1706756|
|Smoking   |ATY    |      0.0000000|        0.0058718|     0.1461753|
|AlcDep    |GATE   |      0.0000249|        0.5914044|     0.0083622|
|Anxiety   |GATE   |      0.0000000|        0.0000068|     0.0001442|
|BIP       |GATE   |      0.0000000|        0.0037628|     0.0002743|
|BMI       |GATE   |      0.0004615|        0.0595848|     0.0002510|
|EA        |GATE   |      0.8556397|        0.0000201|     0.0002929|
|MD        |GATE   |      0.0000000|        0.0000000|     0.0000332|
|MDD       |GATE   |      0.0000000|        0.0001555|     0.0050508|
|Neu       |GATE   |      0.0000000|        0.3043390|     0.0000006|
|PTSD      |GATE   |      0.0000000|        0.0065124|     0.0003708|
|Pain      |GATE   |      0.0000000|        0.2349191|     0.0000023|
|Sleep     |GATE   |      0.0000392|        0.0981094|     0.1116271|
|Smoking   |GATE   |      0.0000001|        0.0110170|     0.0267837|
|AlcDep    |CLIN   |      0.0000000|        0.2197182|     0.8669169|
|Anxiety   |CLIN   |      0.0000000|        0.2583376|     0.0040457|
|BIP       |CLIN   |      0.0000000|        0.2314438|     0.0123620|
|BMI       |CLIN   |      0.0000000|        0.1385062|     0.2248607|
|EA        |CLIN   |      0.0005809|        0.1403267|     0.1931918|
|MD        |CLIN   |      0.0000000|        0.3265085|     0.0002760|
|MDD       |CLIN   |      0.0000000|        0.3035501|     0.0110704|
|Neu       |CLIN   |      0.0000000|        0.6530585|     0.0002799|
|PTSD      |CLIN   |      0.0000000|        0.5209005|     0.0065614|
|Pain      |CLIN   |      0.0000000|        0.1808442|     0.5498820|
|Sleep     |CLIN   |      0.0000001|        0.4405708|     0.6461351|
|Smoking   |CLIN   |      0.0000000|        0.1505291|     0.4242118|
|AlcDep    |COMM   |      0.0000000|        0.7440459|     0.3512964|
|Anxiety   |COMM   |      0.0000000|        0.0031289|     0.4066075|
|BIP       |COMM   |      0.0000000|        0.0168257|     0.3947117|
|BMI       |COMM   |      0.0000000|        0.4776757|     0.3177036|
|EA        |COMM   |      0.0005091|        0.4085868|     0.3121606|
|MD        |COMM   |      0.0000000|        0.0000376|     0.4450019|
|MDD       |COMM   |      0.0000000|        0.0083437|     0.4399300|
|Neu       |COMM   |      0.0000000|        0.0000004|     0.6922790|
|PTSD      |COMM   |      0.0000000|        0.0004157|     0.5929963|
|Pain      |COMM   |      0.0000000|        0.1915965|     0.3387181|
|Sleep     |COMM   |      0.0000000|        0.3132514|     0.5210853|
|Smoking   |COMM   |      0.0000000|        0.7721715|     0.3236180|

