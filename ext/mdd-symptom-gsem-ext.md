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
## minor          1.0                         
## year           2021                        
## month          05                          
## day            18                          
## svn rev        80317                       
## language       R                           
## version.string R version 4.1.0 (2021-05-18)
## nickname       Camp Pontanezen
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
## [1] '0.0.3'
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
## Rows: 15 Columns: 4── Column specification ────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ";"
## chr (4): ref, h, v, abbv
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
## Rows: 15 Columns: 2── Column specification ────────────────────────────────────────────────────────────────────────────────────
## Delimiter: ";"
## chr (2): Reference, Description
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```


# Symptom prevalences

Load previously calculated symptom prevalences:


```r
symptoms_sample_prev <- read_tsv(here::here('meta/symptoms_prev.txt'))
```

```
## Rows: 24 Columns: 6── Column specification ────────────────────────────────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr (3): cohorts, symptom, sumstats
## dbl (3): Nca, Nco, samp_prev
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

```r
pop_prevs_w <-
symptoms_sample_prev %>%
mutate(w=case_when(cohorts == 'AGDS_PGC' ~ 0.15,
                   symptom %in% c('MDD1', 'MDD2') ~ 1.0,
                   TRUE ~ 0.57)) %>%
mutate(pop_prev=samp_prev*w) %>%
select(symptom, cohorts, pop_prev) %>%
pivot_wider(names_from=cohorts, values_from=pop_prev) %>%
group_by(symptom) %>%
mutate(pop_prev=mean(c(AGDS_PGC, ALSPAC_UKB))) %>%
select(symptom, pop_prev)
```


# Multivariable LDSC estimation

Calculate LDSC covariance structure for syptoms used in the combined structural model and the external phenotypes.


```r
covstruct_prefix <- 'symptoms.external.covstruct'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

if(!file.exists(covstruct_r)) {

  # list sumstats distribution directories
  sumstats_files <- list.files(here::here('meta', 'munged'), '.gz', full.names=TRUE)

  # pull out which cohorts and symptom 'x' this is from the filename (COHORTS_MDDx_*)
  cohorts_symptoms <- str_match(basename(sumstats_files), '([A-Z_]+).(MDD[:digit:](a|b)?)')[,1]

  sumstats_paths <- data.frame(filename=sumstats_files, sumstats=str_remove(basename(sumstats_files), '.sumstats.gz'))

  # symptom files, prevalences, and trait names
  symptoms_sumstats_prevs <- 
    symptoms_sample_prev %>%
    left_join(sumstats_paths, by='sumstats') %>%
    left_join(pop_prevs_w, by='symptom') %>%
    left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
    mutate(Sample=if_else(cohorts=="AGDS_PGC", true='Clin', false='Pop')) %>%
    mutate(trait_name=paste0(Sample, abbv)) %>%
    filter(trait_name %in% c('ClinAppDec', 'ClinAppInc',
                             'ClinSleDec', 'ClinSleInc',
                             'ClinMotoInc', 'ClinSui',
                             'PopDep', 'PopGuilt', 'PopSui',
                             'PopAnh', 'PopAppInc', 'PopAppDec',
                             'PopSleInc', 'PopSleDec', 'PopFatig',
                             'PopConc'))
    
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

  dput(symptoms_covstruct, covstruct_r, control=c('all', 'digits17'))
  saveRDS(symptoms_covstruct, covstruct_rds)
  
  # check for exact match of deparsed object
  identical(dget(covstruct_r), symptoms_covstruct)

} else {

  symptoms_covstruct <- dget(covstruct_r)

}
```

# Models

## Clinical and population factors

Base model of symptom factors


```r
clin_pop.model <- "
ClinSoma =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
ClinAPPDEC =~ ClinAppDec
ClinSUI =~ ClinSui
PopAffect =~ NA*PopDep + PopGuilt + PopSui
PopNeuroveg =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
ClinSoma ~~ 1*ClinSoma
ClinAPPDEC ~~ 1*ClinAPPDEC
ClinSUI ~~ 1*ClinSUI
PopAffect ~~ 1*PopAffect
PopNeuroveg ~~ 1*PopNeuroveg
PopDep ~~ PopAnh
ClinAppDec ~~ 0*ClinAppDec
ClinSui ~~ 0*ClinSui
"
clin_pop.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## clin_pop.model): CFI estimates below 0 should not be trusted, and indicate that
## the other model fit estimates should be interpreted with caution. A negative CFI
## estimates typically appears due to negative residual variances.
```

```
## elapsed 
##   0.561 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493224 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.40833704981454 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## clin_pop.model): A difference greater than .025 was observed pre- and post-
## smoothing in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## clin_pop.model): A difference greater than .025 was observed pre- and post-
## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
## large difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```

```r
clin_pop.fit$modelfit
```

<div class="kable-table">

|   |    chisq| df| p_chisq|      AIC|       CFI|     SRMR|
|:--|--------:|--:|-------:|--------:|---------:|--------:|
|df | 3448.373| 97|       0| 3526.373| -0.115681| 1.102515|

</div>

```r
clin_pop.fit$results[c(1,2,3,6,7,9)]
```

<div class="kable-table">

|   |lhs         |op |rhs         | STD_Genotype|STD_Genotype_SE    |   p_value|
|:--|:-----------|:--|:-----------|------------:|:------------------|---------:|
|11 |ClinSoma    |=~ |ClinAppInc  |    0.6172062|0.234390128455356  | 0.0084563|
|13 |ClinSoma    |=~ |ClinSleDec  |    0.3285757|0.242923623806449  | 0.1761654|
|14 |ClinSoma    |=~ |ClinSleInc  |    0.5015887|0.309176037184633  | 0.1047204|
|12 |ClinSoma    |=~ |ClinMotoInc |    0.5274582|0.226754280285284  | 0.0200104|
|25 |PopAffect   |=~ |PopDep      |    0.7283865|0.0732018717045886 | 0.0000000|
|26 |PopAffect   |=~ |PopGuilt    |    0.7787916|0.0971568272755929 | 0.0000000|
|27 |PopAffect   |=~ |PopSui      |    0.6749296|0.100534114227696  | 0.0000000|
|38 |PopNeuroveg |=~ |PopAnh      |    0.8028936|0.0670154414621765 | 0.0000000|
|40 |PopNeuroveg |=~ |PopAppInc   |    0.4533701|0.0752779490778464 | 0.0000000|
|39 |PopNeuroveg |=~ |PopAppDec   |    0.1924650|0.0838677628122933 | 0.0217413|
|44 |PopNeuroveg |=~ |PopSleInc   |    0.4974910|0.0965615899155911 | 0.0000003|
|43 |PopNeuroveg |=~ |PopSleDec   |    0.6510688|0.105169851039874  | 0.0000000|
|42 |PopNeuroveg |=~ |PopFatig    |    0.7454795|0.0913072363565786 | 0.0000000|
|41 |PopNeuroveg |=~ |PopConc     |    0.7882212|0.0948884391549751 | 0.0000000|
|34 |PopDep      |~~ |PopAnh      |    0.4209173|0.0996560173775461 | 0.0000240|
|7  |ClinAppInc  |~~ |ClinAppInc  |    0.6190542|0.401558790372386  | 0.1231666|
|9  |ClinSleDec  |~~ |ClinSleDec  |    0.8920397|0.541607943124004  | 0.0995651|
|10 |ClinSleInc  |~~ |ClinSleInc  |    0.7484088|0.834442725160734  | 0.3697854|
|8  |ClinMotoInc |~~ |ClinMotoInc |    0.7217901|0.516345053208815  | 0.1621531|
|35 |PopDep      |~~ |PopDep      |    0.4694540|0.132257516933774  | 0.0003860|
|37 |PopGuilt    |~~ |PopGuilt    |    0.3934857|0.191936555717315  | 0.0403605|
|48 |PopSui      |~~ |PopSui      |    0.5444704|0.260221006990255  | 0.0364101|
|30 |PopAnh      |~~ |PopAnh      |    0.3553602|0.117065889942085  | 0.0024009|
|32 |PopAppInc   |~~ |PopAppInc   |    0.7944563|0.155436064272031  | 0.0000003|
|31 |PopAppDec   |~~ |PopAppDec   |    0.9629531|0.231638254044186  | 0.0000322|
|47 |PopSleInc   |~~ |PopSleInc   |    0.7525028|0.262559880735371  | 0.0041561|
|46 |PopSleDec   |~~ |PopSleDec   |    0.5761092|0.311828552354797  | 0.0646733|
|36 |PopFatig    |~~ |PopFatig    |    0.4442614|0.293619849236642  | 0.1302702|
|33 |PopConc     |~~ |PopConc     |    0.3787081|0.25274231558892   | 0.1340361|
|15 |ClinSoma    |~~ |ClinAPPDEC  |    0.0477607|0.231902170703992  | 0.8368009|
|17 |ClinSoma    |~~ |ClinSUI     |    0.1357271|0.272372423144733  | 0.6182715|
|18 |ClinSoma    |~~ |PopAffect   |    0.0951326|0.154827174950651  | 0.5389527|
|19 |ClinSoma    |~~ |PopNeuroveg |    0.5096870|0.201608268600431  | 0.0114669|
|4  |ClinAPPDEC  |~~ |ClinSUI     |   -0.0556810|0.189730900415623  | 0.7691499|
|5  |ClinAPPDEC  |~~ |PopAffect   |   -0.1106276|0.100359337548617  | 0.2703286|
|6  |ClinAPPDEC  |~~ |PopNeuroveg |   -0.0939028|0.10294879486182   | 0.3617013|
|23 |ClinSUI     |~~ |PopAffect   |    0.7341822|0.132496233531843  | 0.0000000|
|24 |ClinSUI     |~~ |PopNeuroveg |    0.6615762|0.12073058532012   | 0.0000000|
|29 |PopAffect   |~~ |PopNeuroveg |    0.8228287|0.0610313347858923 | 0.0000000|
|1  |ClinAPPDEC  |=~ |ClinAppDec  |    1.0000000|                   |        NA|
|20 |ClinSUI     |=~ |ClinSui     |    1.0000000|                   |        NA|
|16 |ClinSoma    |~~ |ClinSoma    |    1.0000000|                   |        NA|
|3  |ClinAPPDEC  |~~ |ClinAPPDEC  |    1.0000000|                   |        NA|
|22 |ClinSUI     |~~ |ClinSUI     |    1.0000000|                   |        NA|
|28 |PopAffect   |~~ |PopAffect   |    1.0000000|                   |        NA|
|45 |PopNeuroveg |~~ |PopNeuroveg |    1.0000000|                   |        NA|
|2  |ClinAppDec  |~~ |ClinAppDec  |    0.0000000|                   |        NA|
|21 |ClinSui     |~~ |ClinSui     |    0.0000000|                   |        NA|

</div>

## External phenotypes

Compare symptom factors against each external phenotype. Single regression of each external phenotype on each symptom or symptom cluster.


```r
pop_ext.glue <- "
ClinSoma =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
PopAffect =~ NA*PopDep + PopGuilt + PopSui
PopNeuroveg =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
ClinSoma ~~ 1*ClinSoma
PopAffect ~~ 1*PopAffect
PopNeuroveg ~~ 1*PopNeuroveg
PopDep ~~ PopAnh
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

ext.model_list <- lapply(c('ClinSoma', 'ClinAppDec', 'ClinSui', 'PopAffect', 'PopNeuroveg'), function(symptom) str_glue_data(list(symptom=symptom), pop_ext.glue))

ext.fit_list <- lapply(ext.model_list, function(model) usermodel(symptoms_covstruct, estimation='DWLS', model=model))
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##  13.196 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0344079562870582 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.72652759701003 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing for Z-
## statistics in the genetic covariance matrix. This reflects a large difference
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
##  14.318 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.038381882103906 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.7596833326812 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing for Z-
## statistics in the genetic covariance matrix. This reflects a large difference
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
##  14.056 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.034730283398439 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.74014666206664 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing for Z-
## statistics in the genetic covariance matrix. This reflects a large difference
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
##  12.458 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0344079562870582 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.72652759701003 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing for Z-
## statistics in the genetic covariance matrix. This reflects a large difference
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
##   12.21 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0344079562870582 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.72652759701003 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing in the
## genetic covariance matrix. This reflects a large difference and results should
## be interpreted with caution!! This can often result from including low powered
## traits, and you might consider removing those traits from the model. If you are
## going to run a multivariate GWAS we strongly recommend setting the smooth_check
## argument to true to check smoothing for each SNP.

## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = model):
## A difference greater than .025 was observed pre- and post-smoothing for Z-
## statistics in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```


```r
clin_pop_ext_full <-
bind_rows(lapply(ext.fit_list, function(fit) fit$results)) %>%
select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
filter(lhs %in% ext_trait_names, rhs %in% c('ClinSoma', 'ClinAppDec', 'ClinSui', 'PopAffect', 'PopNeuroveg')) %>%
mutate(Beta='Full', Factor=rhs, Phenotype=lhs)
```

Multiple regression of each phenotype on the clinical/population symptom factors, to estimate relationship after condition on each of the other factors. 


```r
clin_ext_mult.model <- "
ClinSoma =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
ClinSoma ~~ 1*ClinSoma
AlcDep ~ ClinSoma + ClinAppDec + ClinSui
Anxiety ~ ClinSoma + ClinAppDec + ClinSui
BIP ~ ClinSoma + ClinAppDec + ClinSui
BMI ~ ClinSoma + ClinAppDec + ClinSui
EA ~ ClinSoma + ClinAppDec + ClinSui
MD ~ ClinSoma + ClinAppDec + ClinSui
MDD ~ ClinSoma + ClinAppDec + ClinSui
Neu ~ ClinSoma + ClinAppDec + ClinSui
PTSD ~ ClinSoma + ClinAppDec + ClinSui
Pain ~ ClinSoma + ClinAppDec + ClinSui
Sleep ~ ClinSoma + ClinAppDec + ClinSui
Smoking ~ ClinSoma + ClinAppDec + ClinSui
"

clin_ext_mult.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_ext_mult.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##   4.967 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0393892077130448 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.76209824965711 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## clin_ext_mult.model): A difference greater than .025 was observed pre- and post-
## smoothing in the genetic covariance matrix. This reflects a large difference
## and results should be interpreted with caution!! This can often result from
## including low powered traits, and you might consider removing those traits from
## the model. If you are going to run a multivariate GWAS we strongly recommend
## setting the smooth_check argument to true to check smoothing for each SNP.
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## clin_ext_mult.model): A difference greater than .025 was observed pre- and post-
## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
## large difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```



```r
pop_ext_mult.model <- "
PopAffect =~ NA*PopDep + PopGuilt + PopSui
PopNeuroveg =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
PopAffect ~~ 1*PopAffect
PopNeuroveg ~~ 1*PopNeuroveg
PopDep ~~ PopAnh
AlcDep ~ PopAffect + PopNeuroveg
Anxiety ~ PopAffect + PopNeuroveg
BIP ~ PopAffect + PopNeuroveg
BMI ~ PopAffect + PopNeuroveg
EA ~ PopAffect + PopNeuroveg
MD ~ PopAffect + PopNeuroveg
MDD ~ PopAffect + PopNeuroveg
Neu ~ PopAffect + PopNeuroveg
PTSD ~ PopAffect + PopNeuroveg
Pain ~ PopAffect + PopNeuroveg
Sleep ~ PopAffect + PopNeuroveg
Smoking ~ PopAffect + PopNeuroveg
"
pop_ext_mult.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_ext_mult.model)
```

```
## [1] "Running primary model"
## [1] "Calculating CFI"
## [1] "Calculating Standardized Results"
## [1] "Calculating SRMR"
## elapsed 
##   6.599 
## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00517308597252561 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.516833421474208 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."
```

```
## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
## pop_ext_mult.model): A difference greater than .025 was observed pre- and post-
## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
## large difference and results should be interpreted with caution!! This can often
## result from including low powered traits, and you might consider removing those
## traits from the model. If you are going to run a multivariate GWAS we strongly
## recommend setting the smooth_check argument to true to check smoothing for each
## SNP.
```


```r
clin_pop_ext_partial <-
bind_rows(lapply(list(clin_ext_mult.fit, pop_ext_mult.fit),
                 function(fit) fit$results)) %>%
  select(lhs, op, rhs, STD_Genotype, STD_Genotype_SE, p_value) %>%
  filter(lhs %in% ext_trait_names, rhs %in% c('ClinSoma', 'ClinAppDec', 'ClinSui', 'PopAffect', 'PopNeuroveg')) %>%
  mutate(Beta='Partial', Factor=rhs, Phenotype=lhs)
```


```r
ggplot(bind_rows(clin_pop_ext_full, clin_pop_ext_partial),
       aes(x=factor(Factor, levels=c('ClinSui', 'ClinAppDec', 'ClinSoma', 'PopNeuroveg', 'PopAffect')),
           y=STD_Genotype,
           color=factor(Beta, levels=c('Partial', 'Full')),
           shape=factor(Beta, levels=c('Partial', 'Full')),
          ymin=qnorm(0.025, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)),
          ymax=qnorm(0.975, mean=STD_Genotype, sd=as.numeric(STD_Genotype_SE)))) +
geom_hline(yintercept=0, col='gray') +
geom_pointrange(position=position_dodge(width=0.5)) +
facet_wrap(~Phenotype) +
scale_x_discrete('Symptom/Factor') +
scale_y_continuous(expression(r[g]), breaks=c(-1, 0, 1)) +
scale_color_discrete('Model: ') +
coord_flip(ylim=c(-1, 1)) +
theme_bw() +
theme(axis.text.y=element_text(size=16),
      strip.text=element_text(size=16),
      legend.title=element_text(size=12),
      legend.text=element_text(size=14),
      legend.position='top') +
labs(color  = "Model: ", shape = "Model: ")
```

![](mdd-symptom-gsem-ext_files/figure-html/clin_pop_ex_plot-1.png)<!-- -->

Test for attenuation of correlations


```r
clin_pop_ext_factor_wide <-
bind_rows(clin_pop_ext_full, clin_pop_ext_partial) %>%
as_tibble() %>%
select(Factor, Phenotype, STD_Genotype, STD_Genotype_SE, p_value, Beta) %>%
pivot_wider(id_cols=c(Phenotype, Factor),
            names_from=Beta, 
            values_from=c(STD_Genotype, STD_Genotype_SE, p_value))

clin_pop_ext_attenuation <-     
clin_pop_ext_factor_wide %>%
mutate(rg_full=STD_Genotype_Full, rg_partial=STD_Genotype_Partial,
       se_full=as.numeric(STD_Genotype_SE_Full),
       se_partial=as.numeric(STD_Genotype_SE_Partial)) %>%
mutate(attenuation_z=(rg_full-rg_partial)/sqrt(se_full^2+se_partial^2)) %>%
mutate(attenuation_p=2*pnorm(abs(attenuation_z), lower.tail=F)) %>%
select(Phenotype, Factor, p_value_Full, p_value_Partial, attenuation_p)
```

