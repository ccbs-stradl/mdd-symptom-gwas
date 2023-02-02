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
    ## platform       x86_64-pc-linux-gnu         
    ## arch           x86_64                      
    ## os             linux-gnu                   
    ## system         x86_64, linux-gnu           
    ## status                                     
    ## major          4                           
    ## minor          1.2                         
    ## year           2021                        
    ## month          11                          
    ## day            01                          
    ## svn rev        81115                       
    ## language       R                           
    ## version.string R version 4.1.2 (2021-11-01)
    ## nickname       Bird Hippie

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

if(!require(GenomicSEM)) remotes::install_github("MichelNivard/GenomicSEM")
```

GenomicSEM version

``` r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(corrplot)
require(GenomicSEM)

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
    ## ── Column specification ──────────────────────────────────────────────────────────
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
    ## ── Column specification ──────────────────────────────────────────────────────────
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
    ##  1 MDD1      Dep          Mood        Depressed mood most of the day, nearly ever…
    ##  2 MDD2      Anh          Interest    Markedly diminished interest or pleasure in…
    ##  3 MDD3      App          Weight⇅     Significant change in weight or appetite    
    ##  4 MDD3a     AppDec       Weight⇊     Significant weight loss or decrease in appe…
    ##  5 MDD3b     AppInc       Weight⇈     Significant weight gain or increase in appe…
    ##  6 MDD4      Sle          Sleep⇅      Sleeping too much or not sleeping enough    
    ##  7 MDD4a     SleDec       Sleep⇊      Insomnia nearly every day                   
    ##  8 MDD4b     SleInc       Sleep⇈      Hypersomnia nearly every day                
    ##  9 MDD5      Moto         Motor⇅      Changes in speed/amount of moving or speaki…
    ## 10 MDD5a     MotoInc      Motor⇈      Psychomotor agitation nearly every day      
    ## 11 MDD5b     MotoDec      Motor⇊      Psychomotor slowing nearly every day        
    ## 12 MDD6      Fatig        Fatigue     Fatigue or loss of energy nearly every day  
    ## 13 MDD7      Guilt        Guilt       Feelings of worthlessness or excessive or i…
    ## 14 MDD8      Conc         Concentrate Diminished ability to think or concentrate,…
    ## 15 MDD9      Sui          Suicidality Recurrent thoughts of death or suicide or a…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 19 Columns: 9
    ## ── Column specification ──────────────────────────────────────────────────────────
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

commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.331 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

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

## Correlation among directional symptoms

Directional symptoms that can increase or decrease (appetite/weight
changes, sleep, psychomotor) will likely be negatively correlated, so
add an indepdendent factor to capture this.

``` r
commonfactor_app.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc

MDD ~~ 1*MDD
APP ~~ 1*APP
MDD ~~ 0*APP
"

commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.486 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor_app.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1191.863 148 1.534121e-162 1275.863 0.9565176 0.1578159

``` r
commonfactor_app.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('MDD', 'APP')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE  p_value
    ## 1  MDD =~  ClinAppDec       0.0079           0.075  9.2e-01
    ## 2  MDD =~  ClinAppInc      -0.0279           0.104  7.9e-01
    ## 3  MDD =~  ClinSleDec      -0.0957           0.127  4.5e-01
    ## 4  MDD =~  ClinSleInc      -0.1347           0.143  3.5e-01
    ## 5  MDD =~ ClinMotoInc      -0.1785           0.107  9.5e-02
    ## 6  MDD =~ ClinMotoDec      -0.1691           0.165  3.1e-01
    ## 7  MDD =~     ClinSui      -0.4833           0.094  3.1e-07
    ## 8  MDD =~     CommDep      -0.9044           0.047  4.6e-82
    ## 9  MDD =~     CommAnh      -0.9805           0.039 2.0e-142
    ## 10 MDD =~  CommAppDec      -0.1773           0.085  3.8e-02
    ## 11 MDD =~  CommAppInc      -0.3211           0.076  2.7e-05
    ## 12 MDD =~  CommSleDec      -0.5303           0.092  7.2e-09
    ## 13 MDD =~  CommSleInc      -0.4626           0.102  5.8e-06
    ## 14 MDD =~   CommFatig      -0.5744           0.091  3.0e-10
    ## 15 MDD =~   CommGuilt      -0.5820           0.080  2.9e-13
    ## 16 MDD =~    CommConc      -0.6817           0.085  8.3e-16
    ## 17 MDD =~     CommSui      -0.5502           0.090  7.9e-10
    ## 18 MDD =~      UkbDep      -0.8622           0.063  1.3e-42
    ## 19 MDD =~      UkbAnh      -0.8890           0.056  4.4e-57
    ## 20 APP =~  ClinAppDec       0.6207           0.126  8.9e-07
    ## 21 APP =~  ClinAppInc      -0.8232           0.174  2.3e-06
    ## 22 APP =~  CommAppDec       0.3437           0.128  7.2e-03
    ## 23 APP =~  CommAppInc      -0.7992           0.180  9.4e-06
    ## 24 MDD ~~         MDD       1.0000              NA       NA
    ## 25 APP ~~         APP       1.0000              NA       NA
    ## 26 MDD ~~         APP       0.0000              NA       NA

``` r
commonfactor_sle.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

MDD ~~ 1*MDD
SLE ~~ 1*SLE
MDD ~~ 0*SLE

c4a > 0
CommSleDec ~~ c4a*CommSleDec
"

commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   9.163 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_sle.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_sle.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor_sle.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 1301.647 148 1.36484e-183 1385.647 0.9519445 0.1697398

``` r
commonfactor_sle.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('MDD', 'SLE')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1  MDD =~  ClinAppDec        0.031           0.075
    ## 2  MDD =~  ClinAppInc       -0.056           0.104
    ## 3  MDD =~  ClinSleDec       -0.076           0.128
    ## 4  MDD =~  ClinSleInc       -0.134           0.143
    ## 5  MDD =~ ClinMotoInc       -0.178           0.107
    ## 6  MDD =~ ClinMotoDec       -0.169           0.165
    ## 7  MDD =~     ClinSui       -0.483           0.094
    ## 8  MDD =~     CommDep       -0.903           0.047
    ## 9  MDD =~     CommAnh       -0.981           0.038
    ## 10 MDD =~  CommAppDec       -0.161           0.085
    ## 11 MDD =~  CommAppInc       -0.321           0.076
    ## 12 MDD =~  CommSleDec       -0.535           0.092
    ## 13 MDD =~  CommSleInc       -0.474           0.102
    ## 14 MDD =~   CommFatig       -0.575           0.091
    ## 15 MDD =~   CommGuilt       -0.584           0.080
    ## 16 MDD =~    CommConc       -0.682           0.085
    ## 17 MDD =~     CommSui       -0.550           0.089
    ## 18 MDD =~      UkbDep       -0.861           0.063
    ## 19 MDD =~      UkbAnh       -0.888           0.056
    ## 20 SLE =~  ClinSleDec        0.596           0.637
    ## 21 SLE =~  ClinSleInc        0.044           0.320
    ## 22 SLE =~  CommSleDec        0.869           0.936
    ## 23 SLE =~  CommSleInc       -0.240           0.260
    ## 24 MDD ~~         MDD        1.000              NA
    ## 25 SLE ~~         SLE        1.000              NA
    ## 26 MDD ~~         SLE        0.000              NA

## Ascertainment-specific factors

Symptoms were assessed in Clinical and Community samples. Make factors
representing sample type.

``` r
clin_comm.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
COMM ~~ 1*COMM

"

clin_comm.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_comm.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.531 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_comm.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1327.448 151 3.816421e-187 1405.448 0.9509947 0.1749241

``` r
clin_comm.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('CLIN', 'COMM'), rhs %in% c('CLIN', 'COMM'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 CLIN ~~ COMM         0.72            0.37   0.052

Add the sample factors as a bifactor model to the general MDD factor

``` r
clin_comm_bif.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

MDD ~~ 1*MDD
CLIN ~~ 1*CLIN
COMM ~~ 1*COMM

MDD ~~ 0*CLIN + 0*COMM
"

clin_comm_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_comm_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.116 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_comm_bif.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 784.4446 132 2.521846e-93 900.4446 0.9728222 0.1475522

``` r
clin_comm_bif.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CLIN', 'COMM'), rhs %in% c('CLIN', 'COMM'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 CLIN ~~ COMM         0.64            0.33

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

gate.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=gate.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.384 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = gate.model):
    ## A difference greater than .025 was observed pre- and post-smoothing in the
    ## genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered
    ## traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check
    ## argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = gate.model):
    ## A difference greater than .025 was observed pre- and post-smoothing for
    ## Z-statistics in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

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

Comparison of ascertainment and measurement. Gating symptoms from
community sample (distinguish controls from subthreshold), other
symptoms from community sample symptoms (subthreshold from cases), and
clinical cohort symptoms (distinguish cases from each other)

``` r
measure.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
COMM =~ NA*CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui 

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
cl_co < 1
CLIN ~~ cl_co*COMM
"

measure.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=measure.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   9.436 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
measure.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1480.445 149 7.669439e-218 1562.445 0.9445383 0.1595514

``` r
measure.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CLIN', 'COMM', 'GATE'), rhs %in% c('CLIN', 'COMM', 'GATE'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 CLIN ~~ COMM         1.00           0.348
    ## 2 CLIN ~~ GATE         0.44           0.177
    ## 3 GATE ~~ COMM         0.79           0.073

Group symptoms from the clinical cohorts (that had positive genetic
variance) together with the same symptoms from the community cohorts,
then a separate factor for all other community cohort symptoms

``` r
divergence.model <- "
CASE =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommSui
CASECON =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommFatig + CommGuilt + CommConc 

CASE ~~ 1*CASE
CASECON ~~ 1*CASECON
"

divergence.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=divergence.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.561 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## divergence.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## divergence.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
divergence.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1324.772 151 1.252082e-186 1402.772 0.9511062 0.1739938

``` r
divergence.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CASE', 'CASECON'), rhs %in% c('CASE', 'CASECON'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op     rhs STD_Genotype STD_Genotype_SE
    ## 1 CASE ~~ CASECON        -0.86           0.099

## Baseline gating and directional symptoms factors

``` r
base.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

MDD ~~ 1*MDD
GATE ~~ 1*GATE
APP ~~ 1*APP
SLE ~~ 1*SLE
MDD ~~ 0*GATE + 0*APP + 0*SLE
GATE ~~ 0*APP + 0*SLE
APP ~~ 0*SLE

c4a > 0
CommSleDec ~~ c4a*CommSleDec
u1 > 0
UkbDep ~~ u1*UkbDep
"

base.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=base.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  12.207 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = base.model):
    ## A difference greater than .025 was observed pre- and post-smoothing in the
    ## genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered
    ## traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check
    ## argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = base.model):
    ## A difference greater than .025 was observed pre- and post-smoothing for
    ## Z-statistics in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

``` r
base.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 549.4165 140 7.404182e-50 649.4165 0.9829456 0.1379429

## Two-factor models

[Elhai Psychiat
Res 2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Krause Rehab
Psychol 2008](https://psycnet.apa.org/record/2008-17022-011), [Krause
Arch Psys Med
Rehab 2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> slowing as a fourth somatic item

``` r
psych_soma.model <- "
PSYCH =~ NA*ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SOMA =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

PSYCH ~~ 1*PSYCH
SOMA ~~ 1*SOMA
GATE ~~ 1*GATE
APP ~~ 1*APP
SLE ~~ 1*SLE
PSYCH ~~ 0*GATE + 0*APP + 0*SLE
SOMA ~~ 0*GATE + 0*APP + 0*SLE
GATE ~~ 0*APP + 0*SLE
APP ~~ 0*SLE
"

psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  12.679 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
psych_soma.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 532.6555 139 1.957202e-47 634.6555 0.9836022 0.1361536

``` r
psych_soma.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'SOMA'), rhs %in% c('PSYCH', 'SOMA'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ SOMA        -0.92            0.12

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
psych_veg.model <- "
PSYCH =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec +  + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

PSYCH ~~ 1*PSYCH
VEG ~~ 1*VEG
GATE ~~ 1*GATE
APP ~~ 1*APP
SLE ~~ 1*SLE
PSYCH ~~ 0*GATE + 0*APP + 0*SLE
VEG ~~ 0*GATE + 0*APP + 0*SLE
GATE ~~ 0*APP + 0*SLE
APP ~~ 0*SLE
"
psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.167 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
psych_veg.fit$modelfit
```

    ##       chisq  df     p_chisq      AIC       CFI      SRMR
    ## df 533.5023 139 1.42816e-47 635.5023 0.9835669 0.1342801

``` r
psych_veg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'VEG'), rhs %in% c('PSYCH', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ VEG        -0.86           0.098

### Affective-Neurovegetative (Elhai Model 2c)

``` r
affect_neuroveg.model <- "
AFFECT =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommSui + UkbDep + UkbAnh
NEUROVEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec +  + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommConc

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

AFFECT ~~ 1*AFFECT
NEUROVEG ~~ 1*NEUROVEG
GATE ~~ 1*GATE
APP ~~ 1*APP
SLE ~~ 1*SLE
AFFECT ~~ 0*GATE + 0*APP + 0*SLE
NEUROVEG ~~ 0*GATE + 0*APP + 0*SLE
GATE ~~ 0*APP + 0*SLE
APP ~~ 0*SLE
"
affect_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=affect_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  10.619 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
affect_neuroveg.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 557.8332 139 1.554835e-51 659.8332 0.9825534 0.1339122

``` r
affect_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('AFFECT', 'NEUROVEG'), rhs %in% c('AFFECT', 'NEUROVEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##      lhs op      rhs STD_Genotype STD_Genotype_SE
    ## 1 AFFECT ~~ NEUROVEG        -0.87           0.096

## Three factor models

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
cog_mood_neuroveg.model <- "
COG =~ NA*ClinSui + CommGuilt + CommConc + CommSui
MOOD =~ NA*CommDep + CommAnh + CommGuilt + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec +  + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

COG ~~ 1*COG
MOOD ~~ 1*MOOD
VEG ~~ 1*VEG
APP ~~ 1*APP
SLE ~~ 1*SLE
COG ~~ 0*APP + 0*SLE
MOOD ~~ 0*APP + 0*SLE
VEG ~~ 0*APP + 0*SLE
APP ~~ 0*SLE
"
cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.936 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
cog_mood_neuroveg.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 786.5594 140 1.217301e-90 886.5594 0.9730674 0.1363686

``` r
cog_mood_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('COG', 'MOOD', 'VEG'), rhs %in% c('COG', 'MOOD', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1  COG ~~ MOOD         0.74           0.091
    ## 2  COG ~~  VEG        -0.82           0.112
    ## 3 MOOD ~~  VEG        -0.65           0.086

### Model comparisons

``` r
model_list=
list("1a"=list(name="Common", model=commonfactor.fit),
     "1b"=list(name="Common (App)", model=commonfactor_app.fit),
     "1c"=list(name="Common (Sle)", model=commonfactor_sle.fit),
     "1d"=list(name="Common (gating)", model=gate.fit),
     "1e"=list(name="Ascertainment", model=clin_comm.fit),
     "1f"=list(name="Measurement", model=measure.fit),
     "1f"=list(name="Divergence", model=divergence.fit),
     "1g"=list(name="Baseline", model=base.fit),
     "2a"=list(name="Psych-Somatic", model=psych_soma.fit),
     "2b"=list(name="Psych-Neuroveg", model=psych_veg.fit),
     "2c"=list(name="Affect-Neuroveg", model=affect_neuroveg.fit),
     "3"=list(name="Cog-Mood-Neuroveg", model=cog_mood_neuroveg.fit)
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

| Model | Name              | chisq |  df | p\_chisq |  AIC |   CFI |  SRMR |    dAIC |
| :---- | :---------------- | ----: | --: | -------: | ---: | ----: | ----: | ------: |
| 1a    | Common            |  1340 | 152 |        0 | 1410 | 0.951 | 0.175 | 778.000 |
| 1b    | Common (App)      |  1190 | 148 |        0 | 1280 | 0.957 | 0.158 | 641.000 |
| 1c    | Common (Sle)      |  1300 | 148 |        0 | 1390 | 0.952 | 0.170 | 751.000 |
| 1d    | Common (gating)   |  1150 | 148 |        0 | 1230 | 0.958 | 0.163 | 600.000 |
| 1e    | Ascertainment     |  1330 | 151 |        0 | 1410 | 0.951 | 0.175 | 771.000 |
| 1f    | Measurement       |  1480 | 149 |        0 | 1560 | 0.945 | 0.160 | 928.000 |
| 1f    | Divergence        |  1320 | 151 |        0 | 1400 | 0.951 | 0.174 | 768.000 |
| 1g    | Baseline          |   549 | 140 |        0 |  649 | 0.983 | 0.138 |  14.800 |
| 2a    | Psych-Somatic     |   533 | 139 |        0 |  635 | 0.984 | 0.136 |   0.000 |
| 2b    | Psych-Neuroveg    |   534 | 139 |        0 |  636 | 0.984 | 0.134 |   0.847 |
| 2c    | Affect-Neuroveg   |   558 | 139 |        0 |  660 | 0.983 | 0.134 |  25.200 |
| 3     | Cog-Mood-Neuroveg |   787 | 140 |        0 |  887 | 0.973 | 0.136 | 252.000 |

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
