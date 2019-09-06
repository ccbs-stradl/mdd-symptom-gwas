---
title: GenomicSEM of MDD symptoms
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

# Sumstats munging

## Setup

The LDSC support files first need to be downloaded and unpacked


```bash
# LD Score reference files
mkdir -p sumstats/reference
curl https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2 > sumstats/reference/eur_w_ld_chr.tar.bz2
curl https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2 > sumstats/reference/w_hm3.snplist.bz2

tar -xjf sumstats/reference/eur_w_ld_chr.tar.bz2 -C sumstats/reference
rm sumstats/reference/eur_w_ld_chr.tar.bz2
bunzip2 sumstats/reference/w_hm3.snplist.bz2
```

## PGC

PGC sumstats are output in the Ricopoli [daner](https://docs.google.com/document/d/1TWIhr8-qpCXB13WCXcU1_HDio8lC_MeWoAg2jlggrtU/edit) ("**D**osage **An**alyz**er**) format. These can be munged with `munge_sumstats.py` from the [ldsc](https://github.com/bulik/ldsc) program. We find all daner files in the PGC directory and loop them through the munge step.


```bash
# Munge sumstats for all cohorts symptom GWASs

for sumstats in $(ls sumstats/PGC/CasesAllCohorts/daner_MDD*.meta.gz); do

        prefix=$(basename $sumstats .gz)

        munge_sumstats.py --daner-n \
        --sumstats $sumstats \
        --merge-alleles sumstats/reference/w_hm3.snplist \
        --out sumstats/PGC/CasesAllCohorts/${prefix}.ldsc

done

```

## UKB

UKB CIDI sumstats are in a format easily digestible by `munge_sumstats.py`.


```r
# Munge sumstats for UKB CIDI

for sumstats in $(ls sumstats/UKB/CIDI/UKB_CIDI_MDD*.gz); do

        prefix=$(basename $sumstats .gz)

        munge_sumstats.py \
        --sumstats $sumstats \
        --N-cas-col Nca \
        --N-con-col Nco \
        --signed-sumstats OR,1 \
        --p P \
        --merge-alleles sumstats/reference/w_hm3.snplist \
        --out sumstats/UKB/CIDI/${prefix}.ldsc

done
```

# LDSC estimation


### R packages

R version


```r
R.version
```

```
##                _                           
## platform       x86_64-redhat-linux-gnu     
## arch           x86_64                      
## os             linux-gnu                   
## system         x86_64, linux-gnu           
## status                                     
## major          3                           
## minor          6.0                         
## year           2019                        
## month          04                          
## day            26                          
## svn rev        76424                       
## language       R                           
## version.string R version 3.6.0 (2019-04-26)
## nickname       Planting of a Tree
```

Package installation


```r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

library(devtools)

if(!require(GenomicSEM)) install_github("MichelNivard/GenomicSEM")
```

GenomicSEM version

```r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(GenomicSEM)

packageVersion("GenomicSEM")
```

```
## [1] '0.0.2'
```

## Symptom labels

MDD DSM symptoms are numbered 1-9:


```r
dsm_mdd_symptoms_reference <-
read_delim("
MDD1;	Depressed mood most of the day, nearly every day
MDD2;	Markedly diminished interest or pleasure in all, or almost all, activities most of the day, nearly every day
MDD3a;	Significant weight loss or decrease in appetite
MDD3b;	Significant weight gain or increase in appetite
MDD4a;	Insomnia nearly every day
MDD4b;	Hypersomnia nearly every day
MDD5a;	Psychomotor agitation nearly every day
MDD5b;	Psychomotor retardation nearly every day
MDD6;	Fatigue or loss of energy nearly every day
MDD7;	Feelings of worthlessness or excessive or inappropriate guilt
MDD8;	Diminished ability to think or concentrate, or indecisiveness
MDD9;	Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for committing suicide
", col_names=c('Reference', 'Description'), delim=';')

dsm_mdd_symptoms_reference
```

<div class="kable-table">

Reference   Description                                                                                                  
----------  -------------------------------------------------------------------------------------------------------------
MDD1        Depressed mood most of the day, nearly every day                                                             
MDD2        Markedly diminished interest or pleasure in all, or almost all, activities most of the day, nearly every day 
MDD3a       Significant weight loss or decrease in appetite                                                              
MDD3b       Significant weight gain or increase in appetite                                                              
MDD4a       Insomnia nearly every day                                                                                    
MDD4b       Hypersomnia nearly every day                                                                                 
MDD5a       Psychomotor agitation nearly every day                                                                       
MDD5b       Psychomotor retardation nearly every day                                                                     
MDD6        Fatigue or loss of energy nearly every day                                                                   
MDD7        Feelings of worthlessness or excessive or inappropriate guilt                                                
MDD8        Diminished ability to think or concentrate, or indecisiveness                                                
MDD9        Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for committing suicide        

</div>


## Symptom prevalences

Running [multivariable LDSC](https://github.com/MichelNivard/GenomicSEM/wiki/3.-Models-without-Individual-SNP-effects) requires knowing the sample prevalences and population prevalences of each symptom. Sample prevalences can be calculated from the GWAS summary statistics output put population prevalences have to be estimated.

### Population prevalences

We include a table of counts of symptom presence and absence for PGC cohorts


```r
pgc_symptom_counts <- read_table2('sumstats/PGC/CasesAllCohorts/pgc_dsm_symptom_counts.txt')

pgc_symptom_counts %>%
  spread(Status, N) %>%
  unite(AbsentPresent, Absent, Present, sep=':') %>%
  mutate(MDD_counts=paste(MDD, '(absent:present)')) %>%
  select(-MDD) %>%
  spread(MDD_counts, AbsentPresent)
```

<div class="kable-table">

Symptom   Case (absent:present)   Control (absent:present) 
--------  ----------------------  -------------------------
MDD1      914:12689               2809:364                 
MDD2      1451:11671              2591:234                 
MDD3a     6370:7060               2742:83                  
MDD3b     9421:2930               2804:21                  
MDD4a     3340:10209              2639:186                 
MDD4b     7859:3440               2793:32                  
MDD5a     5533:5510               2738:87                  
MDD5b     6405:5815               2756:69                  
MDD6      1736:11833              2602:223                 
MDD7      3072:10113              2703:122                 
MDD8      1501:11209              2719:106                 
MDD9      6194:7221               2779:46                  

</div>

Calculate symptom prevalences separately for cases and controls:


```r
pgc_symptom_prevalences <- 
pgc_symptom_counts %>%
  spread(Status, N) %>%
  mutate(prev=Present / (Present + Absent)) %>%
  select(-Absent, -Present) %>%
  spread(MDD, prev)

pgc_symptom_prevalences
```

<div class="kable-table">

Symptom         Case     Control
--------  ----------  ----------
MDD1       0.9328089   0.1147179
MDD2       0.8894223   0.0828319
MDD3a      0.5256888   0.0293805
MDD3b      0.2372278   0.0074336
MDD4a      0.7534873   0.0658407
MDD4b      0.3044517   0.0113274
MDD5a      0.4989586   0.0307965
MDD5b      0.4758592   0.0244248
MDD6       0.8720613   0.0789381
MDD7       0.7670080   0.0431858
MDD8       0.8819040   0.0375221
MDD9       0.5382780   0.0162832

</div>

```r
pgc_symptom_sample_sizes <- 
pgc_symptom_counts %>%
group_by(Symptom) %>%
summarize(Ntotal=sum(N))
```

Estimation of population prevalence based on average case/control estimates depends on the prevalence of MDD (e.g., [15% in high income countries](https://www.annualreviews.org/doi/10.1146/annurev-publhealth-031912-114409))


```r
pgc_symptom_prev_size <- pgc_symptom_prevalences %>% 
  left_join(pgc_symptom_sample_sizes, by='Symptom')

case_control_prev_lm <- 
lm(Case ~ Control, data=pgc_symptom_prev_size, weights=Ntotal)

summary(case_control_prev_lm)
```

```
## 
## Call:
## lm(formula = Case ~ Control, data = pgc_symptom_prev_size, weights = Ntotal)
## 
## Weighted Residuals:
##     Min      1Q  Median      3Q     Max 
## -21.833  -8.900  -2.081   4.932  35.727 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.36979    0.06711   5.510 0.000258 ***
## Control      6.00913    1.18889   5.054 0.000496 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 16.6 on 10 degrees of freedom
## Multiple R-squared:  0.7187,	Adjusted R-squared:  0.6906 
## F-statistic: 25.55 on 1 and 10 DF,  p-value: 0.000496
```

```r
ggplot(pgc_symptom_prev_size, aes(x=Control, y=Case, weight=Ntotal)) +
  geom_point() +
  geom_text(hjust=-0.1, aes(label=Symptom)) +
  stat_smooth(method='lm') +
  geom_abline() +
  scale_x_continuous('Sympotom prevalence among contols') +
  scale_y_continuous('Symptom prevalence among cases')
```

![](mdd-symptom-gsem_files/figure-html/pgc_symptom_prev-1.png)<!-- -->

Symptoms are are more likely to be present in MDD cases also have higher prevalence in MDD controls. Calculate symptom population prevalences weighted by MDD prevalence 


```r
k <- 0.15
pgc_symptoms_pop_prev <- 
pgc_symptom_prevalences %>%
  transmute(Symptom, pop_prev=k*Case + (1-k)*Control)
 
pgc_symptoms_pop_prev
```

<div class="kable-table">

Symptom     pop_prev
--------  ----------
MDD1       0.2374316
MDD2       0.2038204
MDD3a      0.1038268
MDD3b      0.0419027
MDD4a      0.1689877
MDD4b      0.0552961
MDD5a      0.1010208
MDD5b      0.0921399
MDD6       0.1979065
MDD7       0.1517592
MDD8       0.1641794
MDD9       0.0945824

</div>

### Sample prevalences

Read in headers from PGC daner files. The daner format contains headers for the frequency of the referenec allele in cases (A=Affected) and controls (U=Unaffected) where the column name includes the sample size (`FRQ_A_NNNN`, `FRQ_U_MMMM`)


```r
pgc_daner_meta_gz <- list.files('sumstats/PGC/CasesAllCohorts', pattern='meta.gz', full.names=TRUE)

names(pgc_daner_meta_gz) <- sapply(str_split(basename(pgc_daner_meta_gz), '_'), function(x) x[2])

pgc_daner_meta_frq_cols <- 
bind_rows(
lapply(pgc_daner_meta_gz, function(daner) {
        daner_header <-read.table(daner, nrows=1, stringsAsFactors=F)
        return(data.frame(frq_a_col=daner_header$V6, frq_u_col=daner_header$V7))
}), .id='Symptom')

pgc_daner_meta_frq_cols %>%
gather(key='col', value='frq', frq_a_col:frq_u_col) %>%
select(-col) %>%
separate(frq, into=c('frq', 'status', 'Count')) %>%
mutate(presence=recode(status, 'A'='Present', 'U'='Absent'),
       N=as.integer(Count)) %>%
select(-frq, -status, -Count) %>%
spread(presence, N) %>%
mutate(samp_prev=Present / (Present + Absent))
```

<div class="kable-table">

Symptom    Absent   Present   samp_prev
--------  -------  --------  ----------
MDD1         5221      2538   0.3271040
MDD2         5410      2349   0.3027452
MDD3a        5873      6520   0.5261034
MDD3b        8688      2684   0.2360183
MDD4a        3141      9332   0.7481761
MDD4b        7163      3210   0.3094572
MDD5a        5032      5072   0.5019794
MDD5b        5911      5295   0.4725147
MDD6         1579     10913   0.8735991
MDD7         2862      9363   0.7658896
MDD8         1421     10281   0.8785678
MDD9         5631      6721   0.5441224

</div>
