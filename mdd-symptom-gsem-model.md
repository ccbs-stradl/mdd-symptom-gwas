GenomicSEM of MDD symptoms
================
Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel
Nivard

# Setup

## R packages

R version

``` r
R.version
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

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

if(!require(GenomicSEM)) remotes::install_github("MichelNivard/GenomicSEM")
```

GenomicSEM version

``` r
require(readr)
```

    ## Loading required package: readr

``` r
require(tidyr)
```

    ## Loading required package: tidyr

``` r
require(stringr)
```

    ## Loading required package: stringr

``` r
require(dplyr)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
require(ggplot2)
```

    ## Loading required package: ggplot2

``` r
require(corrplot)
```

    ## Loading required package: corrplot

    ## corrplot 0.92 loaded

``` r
require(GenomicSEM)
```

    ## Loading required package: GenomicSEM

``` r
packageVersion("GenomicSEM")
```

    ## [1] '0.0.5'

# Symptom labels

MDD DSM symptoms are numbered 1-9:

``` r
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

    ## Rows: 15 Columns: 4
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ";"
    ## chr (4): ref, h, v, abbv
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
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

    ## Rows: 15 Columns: 2
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: ";"
    ## chr (2): Reference, Description
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
dsm_mdd_symptoms_reference %>%
left_join(dsm_mdd_symptoms_labels, by=c('Reference'='ref')) %>%
select(Reference, Abbreviation=abbv, Label=h, Description)
```

    ## # A tibble: 15 × 4
    ##    Reference Abbreviation Label       Description                               
    ##    <chr>     <chr>        <chr>       <chr>                                     
    ##  1 MDD1      Dep          Mood        Depressed mood most of the day, nearly ev…
    ##  2 MDD2      Anh          Interest    Markedly diminished interest or pleasure …
    ##  3 MDD3      App          Weight⇅     Significant change in weight or appetite  
    ##  4 MDD3a     AppDec       Weight⇊     Significant weight loss or decrease in ap…
    ##  5 MDD3b     AppInc       Weight⇈     Significant weight gain or increase in ap…
    ##  6 MDD4      Sle          Sleep⇅      Sleeping too much or not sleeping enough  
    ##  7 MDD4a     SleDec       Sleep⇊      Insomnia nearly every day                 
    ##  8 MDD4b     SleInc       Sleep⇈      Hypersomnia nearly every day              
    ##  9 MDD5      Moto         Motor⇅      Changes in speed/amount of moving or spea…
    ## 10 MDD5a     MotoInc      Motor⇈      Psychomotor agitation nearly every day    
    ## 11 MDD5b     MotoDec      Motor⇊      Psychomotor slowing nearly every day      
    ## 12 MDD6      Fatig        Fatigue     Fatigue or loss of energy nearly every day
    ## 13 MDD7      Guilt        Guilt       Feelings of worthlessness or excessive or…
    ## 14 MDD8      Conc         Concentrate Diminished ability to think or concentrat…
    ## 15 MDD9      Sui          Suicidality Recurrent thoughts of death or suicide or…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 38 Columns: 9
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): cohorts, symptom, sumstats, filename, trait_name
    ## dbl (4): Nca, Nco, samp_prev, pop_prev
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Rename samples: AGDS/PGC is the **Clin**ical meta-analysis sample
(`Clin`) and ALSPAC/UKB is the **Comm**unity meta-analysis sample
(`Comm`); and rename symptoms numbers (`MDD1`, `MDD2`) to abbreviations
(`Dep`, `Anh`). There are also extra measures of `MDD1` and `MDD2` from
**UKB** Baseline data.

``` r
cohorts_sample_symptoms <-
sumstats_prevs %>%
left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
select(cohorts, symptom, trait_name, abbv) %>%
mutate(Sample=case_when(cohorts %in% 'AGDS_PGC' ~ 'Clin',
                        cohorts %in% 'ALSPAC_UKB' ~ 'Comm',
                        cohorts %in% 'UKBt' ~ 'Ukb',
                        TRUE ~ NA_character_)) %>%
mutate(sample_symptom=paste0(Sample, abbv))

sample_symptoms <- cohorts_sample_symptoms$sample_symptom
names(sample_symptoms) <- cohorts_sample_symptoms$trait_name

# rename traits in covstruct
dimnames(symptoms_covstruct$S)[[2]] <-
as.vector(sample_symptoms[dimnames(symptoms_covstruct$S)[[2]]])
```

# Structural models

Symptoms with positive variances

``` r
symptoms_S_var <- diag(symptoms_covstruct$S)
names(symptoms_S_var) <- dimnames(symptoms_covstruct$S)[[2]]

symptoms_S_var[which(symptoms_S_var > 0)]
```

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ## 0.128871019 0.051543416 0.013663704 0.013924961 0.039896861 0.002375007 
    ##     ClinSui     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec 
    ## 0.103892914 0.082198484 0.087343514 0.035741853 0.079200630 0.043111949 
    ##  CommSleInc   CommFatig   CommGuilt    CommConc     CommSui      UkbDep 
    ## 0.047553540 0.057539570 0.060012455 0.056887041 0.033166683 0.042463388 
    ##      UkbAnh 
    ## 0.047031487

## Common factor

Common factor across symptoms from both cohorts, as a general MDD factor

``` r
commonfactor.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

MDD ~~ 1*MDD
"

commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.722 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
commonfactor.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC      CFI      SRMR
    ## df 1336.652 152 1.909226e-188 1412.652 0.950653 0.1751983

``` r
commonfactor.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs == 'MDD') %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1  MDD =~  ClinAppDec        0.031           0.075
    ## 2  MDD =~  ClinAppInc       -0.056           0.104
    ## 3  MDD =~  ClinSleDec       -0.095           0.127
    ## 4  MDD =~  ClinSleInc       -0.135           0.143
    ## 5  MDD =~ ClinMotoInc       -0.178           0.107
    ## 6  MDD =~ ClinMotoDec       -0.169           0.165
    ## 7  MDD =~     ClinSui       -0.483           0.094
    ## 8  MDD =~     CommDep       -0.903           0.047
    ## 9  MDD =~     CommAnh       -0.981           0.039
    ## 10 MDD =~  CommAppDec       -0.161           0.085
    ## 11 MDD =~  CommAppInc       -0.321           0.076
    ## 12 MDD =~  CommSleDec       -0.530           0.092
    ## 13 MDD =~  CommSleInc       -0.463           0.102
    ## 14 MDD =~   CommFatig       -0.575           0.091
    ## 15 MDD =~   CommGuilt       -0.584           0.080
    ## 16 MDD =~    CommConc       -0.682           0.085
    ## 17 MDD =~     CommSui       -0.550           0.089
    ## 18 MDD =~      UkbDep       -0.861           0.063
    ## 19 MDD =~      UkbAnh       -0.889           0.056
    ## 20 MDD ~~         MDD        1.000              NA

## Directional symptoms

Common factor with residual correlations among paired directional
symptoms.

``` r
common_dir.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
ClinMotoInc ~~ ClinMotoDec
CommAppDec ~~ CommAppInc
CommSleDec ~~ CommSleInc

MDD ~~ 1*MDD
"

common_dir.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=common_dir.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.653 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## common_dir.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## common_dir.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
common_dir.fit$modelfit
```

    ##     chisq  df       p_chisq    AIC       CFI      SRMR
    ## df 1222.6 147 7.177833e-169 1308.6 0.9551956 0.1664624

``` r
common_dir.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs == 'MDD' | (lhs != 'MDD' & lhs != rhs)) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##            lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1          MDD =~  ClinAppDec        0.029           0.075
    ## 2          MDD =~  ClinAppInc       -0.056           0.104
    ## 3          MDD =~  ClinSleDec       -0.095           0.126
    ## 4          MDD =~  ClinSleInc       -0.134           0.143
    ## 5          MDD =~ ClinMotoInc       -0.180           0.107
    ## 6          MDD =~ ClinMotoDec       -0.173           0.165
    ## 7          MDD =~     ClinSui       -0.483           0.094
    ## 8          MDD =~     CommDep       -0.902           0.047
    ## 9          MDD =~     CommAnh       -0.980           0.038
    ## 10         MDD =~  CommAppDec       -0.171           0.085
    ## 11         MDD =~  CommAppInc       -0.325           0.076
    ## 12         MDD =~  CommSleDec       -0.542           0.092
    ## 13         MDD =~  CommSleInc       -0.477           0.102
    ## 14         MDD =~   CommFatig       -0.575           0.091
    ## 15         MDD =~   CommGuilt       -0.583           0.080
    ## 16         MDD =~    CommConc       -0.682           0.085
    ## 17         MDD =~     CommSui       -0.550           0.089
    ## 18         MDD =~      UkbDep       -0.860           0.063
    ## 19         MDD =~      UkbAnh       -0.888           0.056
    ## 20  ClinAppDec ~~  ClinAppInc       -0.423           0.185
    ## 21  ClinSleDec ~~  ClinSleInc        0.383           0.459
    ## 22 ClinMotoInc ~~ ClinMotoDec       -0.340           0.424
    ## 23  CommAppDec ~~  CommAppInc       -0.191           0.134
    ## 24  CommSleDec ~~  CommSleInc       -0.275           0.185
    ## 25         MDD ~~         MDD        1.000              NA

## Ascertainment-specific factors

### Clinical-Community

Symptoms were assessed in Clinical and Community samples. Make factors
representing sample type.

``` r
clin_comm.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
COMM ~~ 1*COMM

"

clin_comm.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_comm.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.821 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
clin_comm.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1327.451 151 3.811313e-187 1405.451 0.9509946 0.1749241

``` r
clin_comm.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('CLIN', 'COMM'), rhs %in% c('CLIN', 'COMM'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 CLIN ~~ COMM         0.72            0.37   0.052

### Gating symptoms

Then consider the community cardinal symptoms (Depression and
Anhedonia), which are gating items to the rest of the community-sample
symptoms.

``` r
gate.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

MDD ~~ 1*MDD
GATE ~~ 1*GATE
MDD ~~ 0*GATE
"

gate.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=gate.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.634 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## gate.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## gate.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
gate.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI     SRMR
    ## df 1150.775 148 9.954511e-155 1234.775 0.9582291 0.162555

``` r
gate.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('GATE')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op     rhs STD_Genotype STD_Genotype_SE
    ## 1 GATE =~ CommDep         0.65           0.090
    ## 2 GATE =~ CommAnh         0.54           0.089
    ## 3 GATE =~  UkbDep         0.86           0.112
    ## 4 GATE =~  UkbAnh         0.59           0.110
    ## 5 GATE ~~    GATE         1.00              NA

### Gate-Community-Clinical (Spectrum)

Comparison of ascertainment and measurement. Gating symptoms from
community sample (distinguish controls from subthreshold), symptoms from
community sample symptoms (subthreshold from cases), and clinical cohort
symptoms (distinguish cases from each other)

``` r
measure.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
GATE ~~ 0*COMM + 0*CLIN
"

measure.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=measure.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.891 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
measure.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1142.614 147 1.254613e-153 1228.614 0.9585274 0.1625625

``` r
measure.fit$results[c(1,2,3,6,7, 9)] %>%
     filter(lhs %in% c('CLIN', 'COMM', 'GATE'), rhs %in% c('CLIN', 'COMM', 'GATE'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 CLIN ~~ COMM         0.88            0.41   0.031
    ## 2 COMM ~~ GATE         0.00              NA      NA
    ## 3 CLIN ~~ GATE         0.00              NA      NA

### MDD Subtypes

Group symptoms from the clinical cohorts (that had positive genetic
variance) together with the same symptoms from the community cohorts (to
create a dimension of MDD subtypes), then a separate factor for all
other community cohort symptoms (that separates cases from controls)

``` r
subtype.model <- "
SUBTYPE =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommSui
MDD =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommFatig + CommGuilt + CommConc 

SUBTYPE ~~ 1*SUBTYPE
MDD ~~ 1*MDD
"

subtype.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=subtype.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.759 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## subtype.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## subtype.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
subtype.fit$modelfit
```

    ##      chisq  df      p_chisq     AIC       CFI      SRMR
    ## df 1324.77 151 1.25306e-186 1402.77 0.9511063 0.1739938

``` r
subtype.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('SUBTYPE', 'MDD'), rhs %in% c('SUBTYPE', 'MDD'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##       lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 SUBTYPE ~~ MDD        -0.86           0.099

## Two-factor models

[Elhai Psychiat Res
2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Krause Rehab Psychol
2008](https://psycnet.apa.org/record/2008-17022-011), [Krause Arch Psys
Med Rehab
2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> slowing as a fourth somatic item

``` r
psych_soma.model <- "
PSYCH =~ NA*ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SOMA =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
SOMA ~~ 1*SOMA
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*SOMA
"

psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_soma.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.906 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
psych_soma.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1136.492 147 1.815648e-152 1222.492 0.9587824 0.1625163

``` r
psych_soma.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'SOMA'), rhs %in% c('PSYCH', 'SOMA'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ SOMA        -0.98            0.13

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
psych_veg.model <- "
PSYCH =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*VEG
"
psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_veg.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.792 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
psych_veg.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 1081.985 147 3.55372e-142 1167.985 0.9610529 0.1614363

``` r
psych_veg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'VEG'), rhs %in% c('PSYCH', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ VEG         -0.9            0.11

### Affective-Neurovegetative (Elhai Model 2c)

``` r
affect_neuroveg.model <- "
AFFECT =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommSui + UkbDep + UkbAnh
NEUROVEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommConc

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

AFFECT ~~ 1*AFFECT
NEUROVEG ~~ 1*NEUROVEG
GATE ~~ 1*GATE
GATE ~~ 0*AFFECT + 0*NEUROVEG
"
affect_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=affect_neuroveg.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.795 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
affect_neuroveg.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1148.944 147 7.896611e-155 1234.944 0.9582637 0.1608248

``` r
affect_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('AFFECT', 'NEUROVEG'), rhs %in% c('AFFECT', 'NEUROVEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##      lhs op      rhs STD_Genotype STD_Genotype_SE
    ## 1 AFFECT ~~ NEUROVEG        -0.89             0.1

## Three factor models

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
cog_mood_neuroveg.model <- "
COG =~ NA*ClinSui + CommGuilt + CommConc + CommSui
MOOD =~ NA*CommDep + CommAnh + CommGuilt + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

COG ~~ 1*COG
MOOD ~~ 1*MOOD
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*COG + 0*MOOD + 0*VEG
"
cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_mood_neuroveg.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.973 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
cog_mood_neuroveg.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1088.313 144 1.104593e-144 1180.313 0.9606644 0.1615334

``` r
cog_mood_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('COG', 'MOOD', 'VEG'), rhs %in% c('COG', 'MOOD', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1  COG ~~ MOOD         0.86            0.12
    ## 2  COG ~~  VEG        -0.87            0.12
    ## 3 MOOD ~~  VEG        -0.79            0.13

### Cognitive-Appetite-Vegetative (van Loo)

``` r
cog_app_veg.model <- "
COG =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + CommGuilt + ClinSui + CommSui
APP =~ NA*ClinAppInc + ClinAppDec + CommAppDec + app_co3b*CommAppInc
VEG =~ NA*ClinSleInc + CommSleInc + ClinMotoDec + CommFatig + CommConc
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

app_co3b < 1.0

COG ~~ 1*COG
APP ~~ 1*APP
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*COG + 0*APP + 0*VEG

u1 > 0.001
UkbDep ~~ u1*UkbDep
co3b > 0.001
CommAppInc ~~ co3b*CommAppInc
"
cog_app_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_app_veg.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   8.003 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_app_veg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_app_veg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
cog_app_veg.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 862.0416 145 2.498191e-102 952.0416 0.9701314 0.1525697

``` r
cog_app_veg.fit$results[c(1,2,3,6,7, 9)] %>%
     filter(lhs %in% c('COG', 'APP', 'VEG'), rhs %in% c('COG', 'APP', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##   lhs op rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 COG ~~ APP         0.29           0.097 2.9e-03
    ## 2 COG ~~ VEG         0.88           0.110 1.5e-15
    ## 3 APP ~~ VEG         0.35           0.143 1.6e-02

## Melancholic and atypical

Account for directionality of symptoms using melancholic and atypical
classifications, plus remaining affective/cognitive symptoms

``` r
mel_aty_afc.model <- "
MEL =~ NA*CommAnh + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommGuilt
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
AFC =~ NA*CommDep + UkbDep + ClinSui + CommConc + CommSui 
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

MEL ~~ 1*MEL
ATY ~~ 1*ATY
AFC ~~ 1*AFC
GATE ~~ 1*GATE
GATE ~~ 0*AFC + 0*MEL + 0*ATY
"

mel_aty_afc.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=mel_aty_afc.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.856 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty_afc.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty_afc.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
mel_aty_afc.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1044.033 145 6.482037e-136 1134.033 0.9625505 0.1588098

``` r
mel_aty_afc.fit$results[c(1,2,3,6,7,9)] %>%
     filter(rhs %in% c('MEL', 'ATY', 'AFC'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##   lhs op rhs STD_Genotype STD_Genotype_SE  p_value
    ## 1 MEL ~~ ATY         0.84           0.125  2.2e-11
    ## 2 MEL ~~ AFC         0.99           0.039 2.5e-141
    ## 3 ATY ~~ AFC         0.69           0.129  9.4e-08

### Model comparisons

``` r
model_list=
list("A"=list(name="Common", model=commonfactor.fit),
     "B"=list(name="Clinical-Community", model=clin_comm.fit),
     "C"=list(name="Common-Gating", model=gate.fit),
     "D"=list(name="Clinical-Community-Gating", model=measure.fit),
     "E"=list(name="Psych-Somatic", model=psych_soma.fit),
     "F"=list(name="Psych-Neuroveg", model=psych_veg.fit),
     "G"=list(name="Affect-Neuroveg", model=affect_neuroveg.fit),
     "H"=list(name="Cog-Mood-Neuroveg", model=cog_mood_neuroveg.fit),
     "I"=list(name="Mood-Appetite-Vegetative", model=cog_app_veg.fit),
     "J"=list(name="Depression-Melancholic-Atypical", model=mel_aty_afc.fit)
     )

model_fits <- 
data.frame(Model=names(model_list),
           Name=sapply(model_list, function(m) m$name)) %>%
bind_cols(
bind_rows(
lapply(model_list, function(m) m$model$modelfit)
))
rownames(model_fits) <- NULL

knitr::kable(
model_fits %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~signif(., 3))
)
```

| Model | Name                            | chisq |  df | p_chisq |  AIC |   CFI |  SRMR | dAIC |
|:------|:--------------------------------|------:|----:|--------:|-----:|------:|------:|-----:|
| A     | Common                          |  1340 | 152 |       0 | 1410 | 0.951 | 0.175 |  461 |
| B     | Clinical-Community              |  1330 | 151 |       0 | 1410 | 0.951 | 0.175 |  453 |
| C     | Common-Gating                   |  1150 | 148 |       0 | 1230 | 0.958 | 0.163 |  283 |
| D     | Clinical-Community-Gating       |  1140 | 147 |       0 | 1230 | 0.959 | 0.163 |  277 |
| E     | Psych-Somatic                   |  1140 | 147 |       0 | 1220 | 0.959 | 0.163 |  270 |
| F     | Psych-Neuroveg                  |  1080 | 147 |       0 | 1170 | 0.961 | 0.161 |  216 |
| G     | Affect-Neuroveg                 |  1150 | 147 |       0 | 1230 | 0.958 | 0.161 |  283 |
| H     | Cog-Mood-Neuroveg               |  1090 | 144 |       0 | 1180 | 0.961 | 0.162 |  228 |
| I     | Mood-Appetite-Vegetative        |   862 | 145 |       0 |  952 | 0.970 | 0.153 |    0 |
| J     | Depression-Melancholic-Atypical |  1040 | 145 |       0 | 1130 | 0.963 | 0.159 |  182 |

The Mood-Appetite-Vegetative model is the best, SRMR is high across all
the models, indicating that there are high residual correlations.

## Modifications

Add residual correlations betwee same-item symptoms across cohorts.

``` r
cog_app_veg_mod.model <- "
COG =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + CommGuilt + ClinSui + CommSui
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + app_co3b*CommAppInc
VEG =~ NA*ClinSleInc + CommSleInc + ClinMotoDec + CommFatig + CommConc
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

app_co3b > -1.0

COG ~~ 1*COG
APP ~~ 1*APP
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*COG + 0*APP + 0*VEG

u1 > 0.001
UkbDep ~~ u1*UkbDep
co3b > 0.001
CommAppInc ~~ co3b*CommAppInc

ClinAppDec ~~ CommAppDec
ClinAppInc ~~ CommAppInc
ClinSleDec ~~ CommSleDec
ClinSleInc ~~ CommSleInc
ClinSui ~~ CommSui
"

cog_app_veg_mod.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_app_veg_mod.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   8.687 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355012 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782737 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_app_veg_mod.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_app_veg_mod.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
cog_app_veg_mod.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC      CFI      SRMR
    ## df 1143.425 140 5.971991e-157 1243.425 0.958202 0.1429939

``` r
cog_app_veg_mod.fit$results[c(1,2,3,6,7,9)] %>%
     filter(rhs %in% c('COG', 'APP', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##   lhs op rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 COG ~~ APP        -0.34            0.12 6.9e-03
    ## 2 COG ~~ VEG         0.89            0.11 1.8e-15
    ## 3 APP ~~ VEG        -0.39            0.17 2.7e-02

#### Model coefficients

Add modified model to list

``` r
mod_list <- model_list
mod_list[["M"]] <- list(name="Modified", model=cog_app_veg_mod.fit)
```

Check that all symptoms have been included in each model

``` r
sapply(mod_list, function(m) all(names(symptoms_S_var) %in% m$model$results$rhs))
```

    ##    A    B    C    D    E    F    G    H    I    J    M 
    ## TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
model_coefs <- bind_rows(lapply(mod_list, function(m) m$model$results), .id='model')
```

#### Model Implied and residual genetic correlations

Plot the implied genetic and residual correlations, but scaled by the
total genetic correlation.

``` r
fit_cor <- function(fit) {
    # get the implied and residual covariances
    cov_imp <-fit$resid_cov$`Model Implied Covariance Matrix`
    cov_res <- fit$resid_cov$`Residual Covariance Matrix: Calculated as Observed Cov - Model Implied Cov`
    
    # get the residual variances from the model
    symp <- fit$results %>% filter(lhs == rhs) %>% pull(lhs)
    res <- fit$results %>% filter(lhs == rhs) %>% pull(Unstand_Est)
    names(res) <- symp
    
    # convert the implied covariance matrix to a correlation matrix
    cor_imp_var <- cov2cor(cov_imp)
    
    # subtract proportion of residual variance from the
    # correlation matrix
    highpass1 <- function(x) if_else(x > 1, true=1, false=x)
    cov_var <- diag(cov_imp)
    diag(cor_imp_var) <- highpass1((cov_var - res[names(cov_var)]) / cov_var)
    
    # do same for residual matrix
    cov_res_var <- cov_res
    diag(cov_res_var) <- diag(cov_imp)
    
    # convert to correlation matrix
    cor_res_var <- cov2cor(cov_res_var)
    
    # replace diagonals
    diag(cor_res_var) <- res[names(cov_var)] / cov_var
    
    return(list(imp_cor=cor_imp_var, res_cor=cor_res_var))
}
```

``` r
d_fit_cor <- fit_cor(measure.fit)

corrplot(d_fit_cor$imp_cor, is.corr=FALSE, col.lim=c(-1, 1))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptoms_d_imp-1.png)<!-- -->

``` r
corrplot(d_fit_cor$res_cor, is.corr=FALSE, col.lim=c(-1, 1))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptoms_d_resid-1.png)<!-- -->

``` r
i_fit_cor <- fit_cor(cog_app_veg.fit)

corrplot(i_fit_cor$imp_cor, is.corr=FALSE, col.lim=c(-1, 1))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptoms_i_imp-1.png)<!-- -->

``` r
corrplot(i_fit_cor$res_cor, is.corr=FALSE, col.lim=c(-1, 1))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptoms_i_resid-1.png)<!-- -->

# Exploratory factor analysis

Get the genetic covariance matrix for symptoms with a positive
heritability

``` r
symptoms_cov <- symptoms_covstruct$S
k <- nrow(symptoms_cov)
symptoms_se <- matrix(0, k, k)
symptoms_se[lower.tri(symptoms_se, diag=TRUE)] <- sqrt(diag(symptoms_covstruct$V))

symptoms_se[upper.tri(symptoms_se)] <- t(symptoms_se)[upper.tri(symptoms_se)]

symptoms_cov_keep <- which(diag(symptoms_cov > 0))

symptoms_cov_pos <- symptoms_cov[symptoms_cov_keep,symptoms_cov_keep]
```

Smooth the genetic covariance matrix so that it is positive definite

``` r
# smooth the covariance matrix
symptoms_cov_pd <- as.matrix(Matrix::nearPD(symptoms_covstruct$S, corr=FALSE)$mat)

corrplot(cov2cor(symptoms_cov_pd))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_pd-1.png)<!-- -->

Check eigen values of the correlation matrix

``` r
symptoms_eigen <- eigen(cov2cor(symptoms_cov_pd)) 

plot(symptoms_eigen$values, ylab='Eigenvalue')
lines(symptoms_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_eigen-1.png)<!-- -->

``` r
symptoms_efa <- factanal(covmat=symptoms_cov_pd, factors=3, rotation='varimax')
print(symptoms_efa, cut=0.4)
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ##       0.082       0.069       0.005       0.831       0.895       0.968 
    ##     ClinSui     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec 
    ##       0.798       0.112       0.050       0.869       0.290       0.440 
    ##  CommSleInc   CommFatig   CommGuilt    CommConc     CommSui      UkbDep 
    ##       0.757       0.646       0.585       0.583       0.712       0.122 
    ##      UkbAnh 
    ##       0.168 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec           0.567  -0.753 
    ## ClinAppInc           0.408   0.874 
    ## ClinSleDec           0.990         
    ## ClinSleInc           0.402         
    ## ClinMotoInc                        
    ## ClinMotoDec                        
    ## ClinSui      0.439                 
    ## CommDep      0.910                 
    ## CommAnh      0.946                 
    ## CommAppDec                         
    ## CommAppInc                   0.738 
    ## CommSleDec           0.641         
    ## CommSleInc   0.475                 
    ## CommFatig    0.572                 
    ## CommGuilt    0.577                 
    ## CommConc     0.632                 
    ## CommSui      0.535                 
    ## UkbDep       0.906                 
    ## UkbAnh       0.886                 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      5.423   2.371   2.227
    ## Proportion Var   0.285   0.125   0.117
    ## Cumulative Var   0.285   0.410   0.527
    ## 
    ## The degrees of freedom for the model is 117 and the fit was 76.6301
