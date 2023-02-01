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

Common factor across symptoms from both cohorts

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
    ##   1.726 
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
commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 18         MDD =~  ClinAppDec   0.03087295  0.075424956510523
    ## 19         MDD =~  ClinAppInc  -0.05648025  0.104221102483004
    ## 22         MDD =~  ClinSleDec  -0.09548149  0.126520788178042
    ## 23         MDD =~  ClinSleInc  -0.13452986   0.14288723079386
    ## 21         MDD =~ ClinMotoInc  -0.17848264  0.106888684166469
    ## 20         MDD =~ ClinMotoDec  -0.16890134  0.165079547090001
    ## 24         MDD =~     ClinSui  -0.48346865 0.0944034155754929
    ## 29         MDD =~     CommDep  -0.90349685 0.0471059012941759
    ## 25         MDD =~     CommAnh  -0.98119884 0.0385543877169788
    ## 26         MDD =~  CommAppDec  -0.16094698 0.0846229910108325
    ## 27         MDD =~  CommAppInc  -0.32149011 0.0764573829795638
    ## 32         MDD =~  CommSleDec  -0.52989405 0.0917467992980665
    ## 33         MDD =~  CommSleInc  -0.46312540   0.10197792538428
    ## 30         MDD =~   CommFatig  -0.57496480 0.0911192605362049
    ## 31         MDD =~   CommGuilt  -0.58375220 0.0797733005190659
    ## 28         MDD =~    CommConc  -0.68184269 0.0847235197175214
    ## 34         MDD =~     CommSui  -0.55048448 0.0894730634356845
    ## 36         MDD =~      UkbDep  -0.86145790 0.0629920917984675
    ## 35         MDD =~      UkbAnh  -0.88850556 0.0557914448769199
    ## 1   ClinAppDec ~~  ClinAppDec   0.99904927    0.2186295553047
    ## 2   ClinAppInc ~~  ClinAppInc   0.99680855  0.314823426474471
    ## 5   ClinSleDec ~~  ClinSleDec   0.99089135  0.584327069042853
    ## 6   ClinSleInc ~~  ClinSleInc   0.98190706  0.667353923131451
    ## 4  ClinMotoInc ~~ ClinMotoInc   0.96814428  0.393246555771617
    ## 3  ClinMotoDec ~~ ClinMotoDec   0.97143513   1.06541442972338
    ## 7      ClinSui ~~     ClinSui   0.76625689  0.314179499681655
    ## 12     CommDep ~~     CommDep   0.18369295 0.0590847119029301
    ## 8      CommAnh ~~     CommAnh   0.03724902 0.0484852438802214
    ## 9   CommAppDec ~~  CommAppDec   0.97409508  0.239268461363474
    ## 10  CommAppInc ~~  CommAppInc   0.89664759  0.141851052316786
    ## 15  CommSleDec ~~  CommSleDec   0.71921169  0.298824919936334
    ## 16  CommSleInc ~~  CommSleInc   0.78551602   0.23416861051656
    ## 13   CommFatig ~~   CommFatig   0.66941495  0.294638082011445
    ## 14   CommGuilt ~~   CommGuilt   0.65923735  0.157052606498463
    ## 11    CommConc ~~    CommConc   0.53509011  0.239315034547649
    ## 17     CommSui ~~     CommSui   0.69696655  0.228809840679776
    ## 39      UkbDep ~~      UkbDep   0.25788873  0.109417270543876
    ## 38      UkbAnh ~~      UkbAnh   0.21055631 0.0979218365898245
    ## 37         MDD ~~         MDD   1.00000000

## Correlation among directional symptoms

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
    ##   1.636 
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
commonfactor_app.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE       p_value
    ## 23         MDD =~  ClinAppDec  0.007903225 0.0753956169509163  9.165139e-01
    ## 24         MDD =~  ClinAppInc -0.027856253  0.104330416495867  7.894691e-01
    ## 27         MDD =~  ClinSleDec -0.095653697  0.126565019307597  4.497891e-01
    ## 28         MDD =~  ClinSleInc -0.134724300  0.142923605443394  3.458694e-01
    ## 26         MDD =~ ClinMotoInc -0.178486125  0.106779847575673  9.461651e-02
    ## 25         MDD =~ ClinMotoDec -0.169126495  0.165104662229313  3.056651e-01
    ## 29         MDD =~     ClinSui -0.483289481 0.0944170375395253  3.077112e-07
    ## 34         MDD =~     CommDep -0.904444418 0.0471334889445929  4.574495e-82
    ## 30         MDD =~     CommAnh -0.980515609  0.038590261707413 2.038532e-142
    ## 31         MDD =~  CommAppDec -0.177302374 0.0852914224918241  3.763731e-02
    ## 32         MDD =~  CommAppInc -0.321128985 0.0764452986784956  2.660131e-05
    ## 37         MDD =~  CommSleDec -0.530283476 0.0916442258059075  7.194284e-09
    ## 38         MDD =~  CommSleInc -0.462590118  0.101996245579616  5.750402e-06
    ## 35         MDD =~   CommFatig -0.574352912 0.0912073026870912  3.030155e-10
    ## 36         MDD =~   CommGuilt -0.581981356 0.0797253514832343  2.881324e-13
    ## 33         MDD =~    CommConc -0.681747444 0.0846942118967457  8.312176e-16
    ## 39         MDD =~     CommSui -0.550228234  0.089524371620989  7.939938e-10
    ## 41         MDD =~      UkbDep -0.862187251 0.0630225497238655  1.325647e-42
    ## 40         MDD =~      UkbAnh -0.889015188 0.0558310092853427  4.366030e-57
    ## 1          APP =~  ClinAppDec  0.620693147  0.126317535515851  8.935816e-07
    ## 2          APP =~  ClinAppInc -0.823223427    0.1742135344484  2.297099e-06
    ## 3          APP =~  CommAppDec  0.343677693  0.127792502124557  7.159301e-03
    ## 4          APP =~  CommAppInc -0.799213220   0.18035776040401  9.368509e-06
    ## 6   ClinAppDec ~~  ClinAppDec  0.614676780  0.250717535111211  1.421816e-02
    ## 7   ClinAppInc ~~  ClinAppInc  0.321529081  0.355379851257604  3.656034e-01
    ## 10  ClinSleDec ~~  ClinSleDec  0.990841576  0.584283574879615  8.991746e-02
    ## 11  ClinSleInc ~~  ClinSleInc  0.981830113  0.667393166034087  1.412449e-01
    ## 9  ClinMotoInc ~~ ClinMotoInc  0.968140677  0.393155705808944  1.379749e-02
    ## 8  ClinMotoDec ~~ ClinMotoDec  0.971430135   1.06543161734121  3.619067e-01
    ## 12     ClinSui ~~     ClinSui  0.766431405  0.314267075745392  1.473631e-02
    ## 17     CommDep ~~     CommDep  0.181980409 0.0591455647562796  2.092170e-03
    ## 13     CommAnh ~~     CommAnh  0.038589531 0.0486159314627347  4.273385e-01
    ## 14  CommAppDec ~~  CommAppDec  0.850448755  0.236841264387029  3.296670e-04
    ## 15  CommAppInc ~~  CommAppInc  0.258134614  0.294849134806033  3.813189e-01
    ## 20  CommSleDec ~~  CommSleDec  0.718799183  0.298803623514558  1.614647e-02
    ## 21  CommSleInc ~~  CommSleInc  0.786009511  0.234149872050473  7.883069e-04
    ## 18   CommFatig ~~   CommFatig  0.670118309  0.294708746425304  2.297602e-02
    ## 19   CommGuilt ~~   CommGuilt  0.661298337  0.157070434289938  2.551468e-05
    ## 16    CommConc ~~    CommConc  0.535220682   0.23915975934024  2.522623e-02
    ## 22     CommSui ~~     CommSui  0.697249173   0.22882205360819  2.310418e-03
    ## 45      UkbDep ~~      UkbDep  0.256634068  0.109449792316928  1.903959e-02
    ## 44      UkbAnh ~~      UkbAnh  0.209651957 0.0979630654927365  3.234539e-02
    ## 43         MDD ~~         MDD  1.000000000                               NA
    ## 5          APP ~~         APP  1.000000000                               NA
    ## 42         MDD ~~         APP  0.000000000                               NA

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
    ##   9.315 
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
commonfactor_sle.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE       p_value
    ## 18         MDD =~  ClinAppDec  3.105632e-02 0.0754163702775368  6.804884e-01
    ## 19         MDD =~  ClinAppInc -5.634641e-02  0.104210442833294  5.887150e-01
    ## 22         MDD =~  ClinSleDec -7.580523e-02  0.128254437349517  5.544789e-01
    ## 23         MDD =~  ClinSleInc -1.335384e-01   0.14280816813168  3.497416e-01
    ## 21         MDD =~ ClinMotoInc -1.783848e-01  0.106877595022174  9.510649e-02
    ## 20         MDD =~ ClinMotoDec -1.692327e-01  0.165021010912492  3.051167e-01
    ## 24         MDD =~     ClinSui -4.833931e-01 0.0943623497303377  3.011434e-07
    ## 29         MDD =~     CommDep -9.027163e-01 0.0470495292047369  4.798316e-82
    ## 25         MDD =~     CommAnh -9.807173e-01 0.0384682114938047 2.284468e-143
    ## 26         MDD =~  CommAppDec -1.608740e-01 0.0846419510359896  5.734921e-02
    ## 27         MDD =~  CommAppInc -3.213593e-01 0.0764596331251346  2.634145e-05
    ## 32         MDD =~  CommSleDec -5.346713e-01 0.0922620971123642  6.826634e-09
    ## 33         MDD =~  CommSleInc -4.740095e-01   0.10240369668969  3.677247e-06
    ## 30         MDD =~   CommFatig -5.753984e-01 0.0911018911962135  2.684483e-10
    ## 31         MDD =~   CommGuilt -5.837368e-01  0.079725651870396  2.446336e-13
    ## 28         MDD =~    CommConc -6.818719e-01 0.0846821103057028  8.135958e-16
    ## 34         MDD =~     CommSui -5.504059e-01 0.0894186843737062  7.489785e-10
    ## 36         MDD =~      UkbDep -8.610323e-01 0.0629678630189304  1.448845e-42
    ## 35         MDD =~      UkbAnh -8.883374e-01 0.0557753080812947  4.111514e-57
    ## 39         SLE =~  ClinSleDec  5.955909e-01  0.637041755732376  3.498387e-01
    ## 40         SLE =~  ClinSleInc  4.431292e-02  0.319575363471062  8.897342e-01
    ## 41         SLE =~  CommSleDec  8.691263e-01  0.935508611238659  3.528830e-01
    ## 42         SLE =~  CommSleInc -2.399897e-01  0.259700017789543  3.554465e-01
    ## 15  CommSleDec ~~  CommSleDec  2.582653e-08   1.67981891989989  9.999998e-01
    ## 1   ClinAppDec ~~  ClinAppDec  9.990355e-01  0.218628388700028  4.887496e-06
    ## 2   ClinAppInc ~~  ClinAppInc  9.968251e-01  0.314822083694479  1.543849e-03
    ## 5   ClinSleDec ~~  ClinSleDec  6.394295e-01  0.951906220126985  5.016385e-01
    ## 6   ClinSleInc ~~  ClinSleInc  9.801955e-01  0.670374818628035  1.436933e-01
    ## 4  ClinMotoInc ~~ ClinMotoInc  9.681789e-01  0.393254241822881  1.381774e-02
    ## 3  ClinMotoDec ~~ ClinMotoDec  9.713763e-01   1.06546169078945  3.619371e-01
    ## 7      ClinSui ~~     ClinSui  7.663314e-01  0.314128983848222  1.470582e-02
    ## 12     CommDep ~~     CommDep  1.851034e-01 0.0589312233421475  1.683684e-03
    ## 8      CommAnh ~~     CommAnh  3.819358e-02 0.0483431604167105  4.294994e-01
    ## 9   CommAppDec ~~  CommAppDec  9.741196e-01  0.239276792639746  4.679158e-05
    ## 10  CommAppInc ~~  CommAppInc  8.967285e-01  0.141848532656417  2.586466e-10
    ## 16  CommSleInc ~~  CommSleInc  7.177199e-01  0.245522794967052  3.464715e-03
    ## 13   CommFatig ~~   CommFatig  6.689170e-01  0.294644555482452  2.319238e-02
    ## 14   CommGuilt ~~   CommGuilt  6.592515e-01  0.156995886373114  2.679005e-05
    ## 11    CommConc ~~    CommConc  5.350504e-01  0.239305584219658  2.536182e-02
    ## 17     CommSui ~~     CommSui  6.970533e-01  0.228800071184429  2.314739e-03
    ## 45      UkbDep ~~      UkbDep  2.586234e-01  0.109384556091017  1.806187e-02
    ## 44      UkbAnh ~~      UkbAnh  2.108565e-01 0.0979396487158097  3.132488e-02
    ## 37         MDD ~~         MDD  1.000000e+00                               NA
    ## 43         SLE ~~         SLE  1.000000e+00                               NA
    ## 38         MDD ~~         SLE  0.000000e+00                               NA

## Ascertainment-specific factors

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
    ##   1.569 
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
clin_comm.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE       p_value
    ## 1         CLIN =~  ClinAppDec  -0.04258709  0.103295607943293  6.801256e-01
    ## 2         CLIN =~  ClinAppInc   0.08587407  0.143895919537491  5.506612e-01
    ## 5         CLIN =~  ClinSleDec   0.12785444  0.171315438777727  4.554820e-01
    ## 6         CLIN =~  ClinSleInc   0.19301351  0.204056125108648  3.442150e-01
    ## 4         CLIN =~ ClinMotoInc   0.24117067  0.174489866052847  1.669300e-01
    ## 3         CLIN =~ ClinMotoDec   0.22686025  0.230972871037498  3.260067e-01
    ## 7         CLIN =~     ClinSui   0.67182617   0.33218786047748  4.313615e-02
    ## 21        COMM =~     CommDep   0.90359949 0.0470894905974054  4.575227e-82
    ## 17        COMM =~     CommAnh   0.98144684 0.0385785630900186 9.059976e-143
    ## 18        COMM =~  CommAppDec   0.16067209  0.084653893137147  5.769727e-02
    ## 19        COMM =~  CommAppInc   0.32159989 0.0764572808746309  2.596246e-05
    ## 24        COMM =~  CommSleDec   0.52978607 0.0917451058984792  7.716507e-09
    ## 25        COMM =~  CommSleInc   0.46321246  0.101982724850881  5.570600e-06
    ## 22        COMM =~   CommFatig   0.57491092 0.0911159149815166  2.796330e-10
    ## 23        COMM =~   CommGuilt   0.58385890 0.0797793854972886  2.509118e-13
    ## 20        COMM =~    CommConc   0.68177911 0.0847280767453358  8.507301e-16
    ## 26        COMM =~     CommSui   0.55054899 0.0894817682653463  7.620273e-10
    ## 28        COMM =~      UkbDep   0.86155265 0.0629787644345039  1.336140e-42
    ## 27        COMM =~      UkbAnh   0.88850413 0.0557955020364048  4.297986e-57
    ## 10  ClinAppDec ~~  ClinAppDec   0.99818778  0.218943419854232  5.137117e-06
    ## 11  ClinAppInc ~~  ClinAppInc   0.99262123  0.315281430896365  1.641787e-03
    ## 14  ClinSleDec ~~  ClinSleDec   0.98365265  0.581359307329133  9.065010e-02
    ## 15  ClinSleInc ~~  ClinSleInc   0.96274576  0.669215942187155  1.502566e-01
    ## 13 ClinMotoInc ~~ ClinMotoInc   0.94183879  0.389435933693039  1.558620e-02
    ## 12 ClinMotoDec ~~ ClinMotoDec   0.94854238   1.07523346923892  3.776830e-01
    ## 16     ClinSui ~~     ClinSui   0.54865569  0.493520378544786  2.662466e-01
    ## 34     CommDep ~~     CommDep   0.18350749 0.0590361240655969  1.881017e-03
    ## 30     CommAnh ~~     CommAnh   0.03676170 0.0485243168783239  4.486969e-01
    ## 31  CommAppDec ~~  CommAppDec   0.97418506     0.239270379975  4.671460e-05
    ## 32  CommAppInc ~~  CommAppInc   0.89657356  0.141856207038843  2.610470e-10
    ## 37  CommSleDec ~~  CommSleDec   0.71932489  0.298836286772321  1.607940e-02
    ## 38  CommSleInc ~~  CommSleInc   0.78543412  0.234162179257341  7.958363e-04
    ## 35   CommFatig ~~   CommFatig   0.66947573  0.294663787479242  2.308663e-02
    ## 36   CommGuilt ~~   CommGuilt   0.65911062  0.157030796983279  2.700880e-05
    ## 33    CommConc ~~    CommConc   0.53517746  0.239303838864662  2.532633e-02
    ## 39     CommSui ~~     CommSui   0.69689618  0.228822612568779  2.322418e-03
    ## 41      UkbDep ~~      UkbDep   0.25772463  0.109398158498847  1.848010e-02
    ## 40      UkbAnh ~~      UkbAnh   0.21055869 0.0979042439958563  3.150257e-02
    ## 9         CLIN ~~        COMM   0.71654461  0.368272154113315  5.169597e-02
    ## 8         CLIN ~~        CLIN   1.00000000                               NA
    ## 29        COMM ~~        COMM   1.00000000                               NA

Community gating symptom factors (Depression and Anhedonia)

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
    ##   1.417 
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
gate.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 23         MDD =~  ClinAppDec   0.06301077 0.0862796172566476
    ## 24         MDD =~  ClinAppInc  -0.11435466  0.118161345805648
    ## 27         MDD =~  ClinSleDec  -0.11456742  0.137745732394995
    ## 28         MDD =~  ClinSleInc  -0.18316849  0.161609612865441
    ## 26         MDD =~ ClinMotoInc  -0.20980074  0.120098508018152
    ## 25         MDD =~ ClinMotoDec  -0.23599392  0.189406340005686
    ## 29         MDD =~     ClinSui  -0.58707979  0.106694333523352
    ## 34         MDD =~     CommDep  -0.67784590 0.0670091099792996
    ## 30         MDD =~     CommAnh  -0.80621775 0.0644741497837596
    ## 31         MDD =~  CommAppDec  -0.18730404 0.0929877225854634
    ## 32         MDD =~  CommAppInc  -0.39244995 0.0856228642127589
    ## 37         MDD =~  CommSleDec  -0.62750309  0.107198744968471
    ## 38         MDD =~  CommSleInc  -0.51333087  0.114366879382452
    ## 35         MDD =~   CommFatig  -0.69268770   0.10101861025261
    ## 36         MDD =~   CommGuilt  -0.68977702 0.0895215139766066
    ## 33         MDD =~    CommConc  -0.81116771  0.098555478518087
    ## 39         MDD =~     CommSui  -0.62489641 0.0988832960246222
    ## 41         MDD =~      UkbDep  -0.53490989 0.0928888484823965
    ## 40         MDD =~      UkbAnh  -0.67853104 0.0873324002032185
    ## 19        GATE =~     CommDep   0.64807811 0.0904538901668868
    ## 18        GATE =~     CommAnh   0.53572570 0.0888120837093477
    ## 21        GATE =~      UkbDep   0.86195131  0.111651266938297
    ## 20        GATE =~      UkbAnh   0.59319773  0.109710352474906
    ## 1   ClinAppDec ~~  ClinAppDec   0.99602912  0.218792956013309
    ## 2   ClinAppInc ~~  ClinAppInc   0.98692176  0.315710579965843
    ## 5   ClinSleDec ~~  ClinSleDec   0.98687544  0.583078299814053
    ## 6   ClinSleInc ~~  ClinSleInc   0.96645210  0.667152647041503
    ## 4  ClinMotoInc ~~ ClinMotoInc   0.95598499  0.392430380061258
    ## 3  ClinMotoDec ~~ ClinMotoDec   0.94429605   1.07277497017868
    ## 7      ClinSui ~~     ClinSui   0.65533761  0.322936430300814
    ## 12     CommDep ~~     CommDep   0.12051942 0.0527649422689811
    ## 8      CommAnh ~~     CommAnh   0.06301044 0.0434619297811766
    ## 9   CommAppDec ~~  CommAppDec   0.96491402  0.238790544776165
    ## 10  CommAppInc ~~  CommAppInc   0.84598363  0.143816585322348
    ## 15  CommSleDec ~~  CommSleDec   0.60624000  0.297071554896014
    ## 16  CommSleInc ~~  CommSleInc   0.73649257  0.236569320944274
    ## 13   CommFatig ~~   CommFatig   0.52018426   0.29542775949499
    ## 14   CommGuilt ~~   CommGuilt   0.52420743  0.166160239347965
    ## 11    CommConc ~~    CommConc   0.34200508  0.242744930666403
    ## 17     CommSui ~~     CommSui   0.60950508  0.232531462593675
    ## 45      UkbDep ~~      UkbDep  -0.02908761  0.123768775908716
    ## 44      UkbAnh ~~      UkbAnh   0.18771147 0.0843467524937709
    ## 43         MDD ~~         MDD   1.00000000                   
    ## 22        GATE ~~        GATE   1.00000000                   
    ## 42         MDD ~~        GATE   0.00000000

Gating symptoms, case and control symptoms, and case-only symptoms

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
    ##   9.656 
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
measure.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE       p_value
    ## 1         CLIN =~  ClinAppDec  -0.11470422  0.109301466805721  2.939816e-01
    ## 2         CLIN =~  ClinAppInc   0.19700902  0.144314643334302  1.722117e-01
    ## 5         CLIN =~  ClinSleDec   0.13610893  0.164439013763595  4.078310e-01
    ## 6         CLIN =~  ClinSleInc   0.25309949  0.203626911695486  2.138864e-01
    ## 4         CLIN =~ ClinMotoInc   0.25163612  0.158420111064599  1.121946e-01
    ## 3         CLIN =~ ClinMotoDec   0.32702418  0.235507271814263  1.649608e-01
    ## 7         CLIN =~     ClinSui   0.77865868  0.255022076574024  2.263596e-03
    ## 38        GATE =~     CommDep   0.92718282  0.048268439337077  3.118824e-82
    ## 37        GATE =~     CommAnh   1.00910846 0.0393277234192472 3.358178e-145
    ## 40        GATE =~      UkbDep   0.88421665 0.0636077042935733  6.238258e-44
    ## 39        GATE =~      UkbAnh   0.90007449 0.0558607291573742  2.072924e-58
    ## 18        COMM =~  CommAppDec   0.17601501  0.088601003232122  4.696666e-02
    ## 19        COMM =~  CommAppInc   0.37606979 0.0825883701482292  5.274902e-06
    ## 23        COMM =~  CommSleDec   0.59663925  0.102789581515603  6.456861e-09
    ## 24        COMM =~  CommSleInc   0.49092940  0.111336662928006  1.036500e-05
    ## 21        COMM =~   CommFatig   0.66194861 0.0994643961075617  2.830438e-11
    ## 22        COMM =~   CommGuilt   0.66384908 0.0872647530127947  2.798753e-14
    ## 20        COMM =~    CommConc   0.76988540 0.0965676568981218  1.554938e-15
    ## 25        COMM =~     CommSui   0.60244747 0.0973595075704742  6.098390e-10
    ## 9         CLIN ~~        COMM   1.00000036  0.347667481922324  4.023805e-03
    ## 11  ClinAppDec ~~  ClinAppDec   0.98684273  0.220453877526029  7.590703e-06
    ## 12  ClinAppInc ~~  ClinAppInc   0.96118536  0.318803577259212  2.569934e-03
    ## 15  ClinSleDec ~~  ClinSleDec   0.98147010  0.581960398190021  9.170015e-02
    ## 16  ClinSleInc ~~  ClinSleInc   0.93590686  0.666569836082144  1.602844e-01
    ## 14 ClinMotoInc ~~ ClinMotoInc   0.93667409  0.392670050283812  1.706008e-02
    ## 13 ClinMotoDec ~~ ClinMotoDec   0.89283303   1.09146331363914  4.132311e-01
    ## 17     ClinSui ~~     ClinSui   0.39365221   0.45365539670254  3.854810e-01
    ## 31     CommDep ~~     CommDep   0.14033172 0.0588849191382344  1.716513e-02
    ## 27     CommAnh ~~     CommAnh  -0.01830002 0.0505235295664466  7.171930e-01
    ## 44      UkbDep ~~      UkbDep   0.21816123  0.106211803180794  3.997393e-02
    ## 43      UkbAnh ~~      UkbAnh   0.18986541  0.096207518745071  4.843870e-02
    ## 28  CommAppDec ~~  CommAppDec   0.96901872   0.23934430766459  5.151634e-05
    ## 29  CommAppInc ~~  CommAppInc   0.85857157  0.142569686645977  1.721517e-09
    ## 34  CommSleDec ~~  CommSleDec   0.64402135  0.297439838934825  3.037128e-02
    ## 35  CommSleInc ~~  CommSleInc   0.75898835  0.233764802779492  1.167103e-03
    ## 32   CommFatig ~~   CommFatig   0.56182373  0.291481957980129  5.392074e-02
    ## 33   CommGuilt ~~   CommGuilt   0.55930421  0.163292031576354  6.144091e-04
    ## 30    CommConc ~~    CommConc   0.40727630  0.241546362216286  9.177330e-02
    ## 36     CommSui ~~     CommSui   0.63705691   0.23223449726845  6.085082e-03
    ## 10        CLIN ~~        GATE   0.43622752   0.17713039689377  1.378827e-02
    ## 41        GATE ~~        COMM   0.79497026 0.0734396126631191  2.625852e-27
    ## 8         CLIN ~~        CLIN   1.00000000                               NA
    ## 42        GATE ~~        GATE   1.00000000                               NA
    ## 26        COMM ~~        COMM   1.00000000                               NA

Symptoms with positive genetic variance between cases vs others

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
    ##   1.594 
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
divergence.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE       p_value
    ## 1         CASE =~  ClinAppDec   0.03992252 0.0846394984030654  6.371481e-01
    ## 2         CASE =~  ClinAppInc  -0.08044936  0.116751598286193  4.907821e-01
    ## 5         CASE =~  ClinSleDec  -0.11795962  0.138566241846208  3.946059e-01
    ## 6         CASE =~  ClinSleInc  -0.15840993  0.159123926302155  3.194841e-01
    ## 4         CASE =~ ClinMotoInc  -0.20367475  0.121780641277375  9.442861e-02
    ## 3         CASE =~ ClinMotoDec  -0.19938254   0.18226738166505  2.739967e-01
    ## 7         CASE =~     ClinSui  -0.54368606   0.11313615815691  1.542790e-06
    ## 8         CASE =~  CommAppDec  -0.18144657 0.0916020927360044  4.761267e-02
    ## 9         CASE =~  CommAppInc  -0.36652952  0.089293307985168  4.046542e-05
    ## 10        CASE =~  CommSleDec  -0.59851179  0.111960209971917  9.003104e-08
    ## 11        CASE =~  CommSleInc  -0.51038871  0.121984653440016  2.863501e-05
    ## 12        CASE =~     CommSui  -0.61810262  0.106381884424632  6.238486e-09
    ## 17     CASECON =~     CommDep   0.90634336 0.0473207367541238  9.102521e-82
    ## 15     CASECON =~     CommAnh   0.98491145  0.038386436939519 3.459219e-145
    ## 21     CASECON =~      UkbDep   0.86417730  0.062870768001688  5.434765e-43
    ## 20     CASECON =~      UkbAnh   0.89001702 0.0557073717461034  1.858791e-57
    ## 18     CASECON =~   CommFatig   0.57416401 0.0912125523178972  3.077917e-10
    ## 19     CASECON =~   CommGuilt   0.58386714 0.0799392771987675  2.795808e-13
    ## 16     CASECON =~    CommConc   0.68151956 0.0848733118104004  9.759953e-16
    ## 23  ClinAppDec ~~  ClinAppDec   0.99840964  0.218718780821048  4.999945e-06
    ## 24  ClinAppInc ~~  ClinAppInc   0.99352682    0.3150388343945  1.612364e-03
    ## 27  ClinSleDec ~~  ClinSleDec   0.98608692  0.582723896000264  9.060876e-02
    ## 28  ClinSleInc ~~  ClinSleInc   0.97490923  0.667880501266091  1.443722e-01
    ## 26 ClinMotoInc ~~ ClinMotoInc   0.95851654  0.392071926393186  1.449559e-02
    ## 25 ClinMotoDec ~~ ClinMotoDec   0.96024234   1.06845570100574  3.687987e-01
    ## 29     ClinSui ~~     ClinSui   0.70440132  0.321570154889791  2.848745e-02
    ## 31  CommAppDec ~~  CommAppDec   0.96707426  0.239630175318302  5.443733e-05
    ## 32  CommAppInc ~~  CommAppInc   0.86565546  0.145093821098716  2.428920e-09
    ## 37  CommSleDec ~~  CommSleDec   0.64177901  0.299164455492572  3.193322e-02
    ## 38  CommSleInc ~~  CommSleInc   0.73950259   0.24309289441317  2.349788e-03
    ## 39     CommSui ~~     CommSui   0.61794944  0.248925832982639  1.304786e-02
    ## 34     CommDep ~~     CommDep   0.17854270 0.0592005979110395  2.562222e-03
    ## 30     CommAnh ~~     CommAnh   0.02994833 0.0477741887071577  5.307396e-01
    ## 41      UkbDep ~~      UkbDep   0.25319827  0.108807003187553  1.996375e-02
    ## 40      UkbAnh ~~      UkbAnh   0.20786993 0.0979091693255113  3.374747e-02
    ## 35   CommFatig ~~   CommFatig   0.67033333  0.294706389373076  2.293113e-02
    ## 36   CommGuilt ~~   CommGuilt   0.65909898   0.15709878602602  2.723278e-05
    ## 33    CommConc ~~    CommConc   0.53553217  0.239505832963301  2.535299e-02
    ## 14        CASE ~~     CASECON  -0.86110496 0.0985262113833212  2.333173e-18
    ## 13        CASE ~~        CASE   1.00000000                               NA
    ## 22     CASECON ~~     CASECON   1.00000000                               NA

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
    ##  12.426 
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

``` r
base.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 31         MDD =~  ClinAppDec  2.400068e-02 0.0864732737761684
    ## 32         MDD =~  ClinAppInc -6.777572e-02  0.118782418523133
    ## 35         MDD =~  ClinSleDec -9.463892e-02  0.141568722485598
    ## 36         MDD =~  ClinSleInc -1.859500e-01  0.161586077392536
    ## 34         MDD =~ ClinMotoInc -2.096724e-01  0.119853836735535
    ## 33         MDD =~ ClinMotoDec -2.369060e-01  0.189238380320869
    ## 37         MDD =~     ClinSui -5.857554e-01   0.10658247736388
    ## 42         MDD =~     CommDep -6.785403e-01 0.0666081806460785
    ## 38         MDD =~     CommAnh -8.041504e-01 0.0639453263673254
    ## 39         MDD =~  CommAppDec -2.170140e-01 0.0945431782626874
    ## 40         MDD =~  CommAppInc -3.861764e-01 0.0853398317871467
    ## 45         MDD =~  CommSleDec -6.404276e-01  0.108133494892181
    ## 46         MDD =~  CommSleInc -5.370399e-01  0.115844244887878
    ## 43         MDD =~   CommFatig -6.917903e-01  0.101284837985503
    ## 44         MDD =~   CommGuilt -6.846651e-01 0.0890415489322731
    ## 41         MDD =~    CommConc -8.102995e-01 0.0983059233796283
    ## 47         MDD =~     CommSui -6.235091e-01 0.0988628815622949
    ## 49         MDD =~      UkbDep -5.398409e-01 0.0926459193608735
    ## 48         MDD =~      UkbAnh -6.799796e-01  0.087306463689032
    ## 25        GATE =~     CommDep  6.499977e-01 0.0899567720384444
    ## 24        GATE =~     CommAnh  5.384938e-01 0.0884080473416342
    ## 27        GATE =~      UkbDep  8.471919e-01   0.11220816918733
    ## 26        GATE =~      UkbAnh  5.923275e-01  0.110035159769216
    ## 1          APP =~  ClinAppDec  6.220021e-01  0.126553168740481
    ## 2          APP =~  ClinAppInc -8.107827e-01  0.172221587510311
    ## 3          APP =~  CommAppDec  3.688119e-01  0.129980752876614
    ## 4          APP =~  CommAppInc -7.896364e-01  0.177078726112111
    ## 54         SLE =~  ClinSleDec  5.322908e-01  0.431774728136065
    ## 55         SLE =~  ClinSleInc -6.175547e-02  0.323305715951815
    ## 56         SLE =~  CommSleDec  8.164809e-01  0.684911209097567
    ## 57         SLE =~  CommSleInc -3.770270e-01  0.312124535717793
    ## 21  CommSleDec ~~  CommSleDec  1.210468e-07   1.18206419461726
    ## 60      UkbDep ~~      UkbDep  1.152302e-07  0.120629704030633
    ## 7   ClinAppDec ~~  ClinAppDec  6.125373e-01  0.251471738194413
    ## 8   ClinAppInc ~~  ClinAppInc  3.380380e-01  0.351181499317523
    ## 11  ClinSleDec ~~  ClinSleDec  7.076890e-01  0.726164936019098
    ## 12  ClinSleInc ~~  ClinSleInc  9.616073e-01  0.663941167196385
    ## 10 ClinMotoInc ~~ ClinMotoInc  9.560375e-01  0.392183983187364
    ## 9  ClinMotoDec ~~ ClinMotoDec  9.438651e-01   1.07290450728844
    ## 13     ClinSui ~~     ClinSui  6.568901e-01  0.322775100185114
    ## 18     CommDep ~~     CommDep  1.170853e-01 0.0528169608654347
    ## 14     CommAnh ~~     CommAnh  6.336664e-02 0.0434192988207158
    ## 15  CommAppDec ~~  CommAppDec  8.168826e-01  0.238068746298768
    ## 16  CommAppInc ~~  CommAppInc  2.273421e-01  0.287469010168785
    ## 22  CommSleInc ~~  CommSleInc  5.694376e-01  0.308078297024581
    ## 19   CommFatig ~~   CommFatig  5.214255e-01  0.295704881700229
    ## 20   CommGuilt ~~   CommGuilt  5.312334e-01  0.165677359061429
    ## 17    CommConc ~~    CommConc  3.434136e-01   0.24174815764722
    ## 23     CommSui ~~     CommSui  6.112360e-01  0.232475299925383
    ## 59      UkbAnh ~~      UkbAnh  1.867756e-01 0.0844687013511417
    ## 52         MDD ~~         MDD  1.000000e+00                   
    ## 29        GATE ~~        GATE  1.000000e+00                   
    ## 5          APP ~~         APP  1.000000e+00                   
    ## 58         SLE ~~         SLE  1.000000e+00                   
    ## 51         MDD ~~        GATE  0.000000e+00                   
    ## 50         MDD ~~         APP  0.000000e+00                   
    ## 53         MDD ~~         SLE  0.000000e+00                   
    ## 28        GATE ~~         APP  0.000000e+00                   
    ## 30        GATE ~~         SLE  0.000000e+00                   
    ## 6          APP ~~         SLE  0.000000e+00

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
    ##  12.709 
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
psych_soma.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 32       PSYCH =~ ClinMotoInc  2.084517e-01  0.119650388070368
    ## 31       PSYCH =~ ClinMotoDec  2.347327e-01   0.18880250847073
    ## 33       PSYCH =~     ClinSui  5.841014e-01  0.106566487818872
    ## 36       PSYCH =~     CommDep  6.854636e-01 0.0681319686382665
    ## 34       PSYCH =~     CommAnh  8.102457e-01 0.0655284251543423
    ## 37       PSYCH =~   CommGuilt  6.840831e-01  0.089205699420619
    ## 35       PSYCH =~    CommConc  8.079536e-01 0.0987083276209315
    ## 38       PSYCH =~     CommSui  6.236881e-01 0.0989292978553161
    ## 40       PSYCH =~      UkbDep  5.422966e-01 0.0935838245493764
    ## 39       PSYCH =~      UkbAnh  6.858665e-01 0.0883663972959003
    ## 51        SOMA =~  ClinAppDec  2.268940e-02  0.091491099334393
    ## 52        SOMA =~  ClinAppInc -7.591392e-02  0.125018897292158
    ## 53        SOMA =~  ClinSleDec -8.043465e-02  0.149670895725705
    ## 54        SOMA =~  ClinSleInc -1.998658e-01  0.171611115183355
    ## 55        SOMA =~  CommAppDec -2.294578e-01 0.0984669316261766
    ## 56        SOMA =~  CommAppInc -4.077254e-01 0.0934650517441628
    ## 58        SOMA =~  CommSleDec -6.881446e-01  0.121951257680068
    ## 59        SOMA =~  CommSleInc -5.763906e-01  0.133070218519894
    ## 57        SOMA =~   CommFatig -7.316242e-01  0.117748567607308
    ## 25        GATE =~     CommDep  6.400216e-01  0.092987280964847
    ## 24        GATE =~     CommAnh  5.283738e-01 0.0919270513241979
    ## 27        GATE =~      UkbDep  8.585669e-01  0.113856616135817
    ## 26        GATE =~      UkbAnh  5.849863e-01  0.113420435261398
    ## 1          APP =~  ClinAppDec  6.227219e-01  0.126351557524694
    ## 2          APP =~  ClinAppInc -8.057128e-01  0.171497108741214
    ## 3          APP =~  CommAppDec  3.763046e-01  0.130384994959979
    ## 4          APP =~  CommAppInc -7.894496e-01  0.176717616638823
    ## 46         SLE =~  ClinSleDec  2.260639e-02   8.18706046793016
    ## 47         SLE =~  ClinSleInc -5.062027e-04  0.184204814504322
    ## 48         SLE =~  CommSleDec  2.385098e+01   8638.15261380878
    ## 49         SLE =~  CommSleInc -1.732637e-02   6.27541475704098
    ## 10 ClinMotoInc ~~ ClinMotoInc  9.565480e-01  0.392190358608507
    ## 9  ClinMotoDec ~~ ClinMotoDec  9.449001e-01   1.07263183162876
    ## 13     ClinSui ~~     ClinSui  6.588255e-01  0.322692610857926
    ## 18     CommDep ~~     CommDep  1.205120e-01 0.0527580221436207
    ## 14     CommAnh ~~     CommAnh  6.432303e-02 0.0436930596830622
    ## 20   CommGuilt ~~   CommGuilt  5.320304e-01  0.165743073024942
    ## 17    CommConc ~~    CommConc  3.472112e-01  0.241044106348979
    ## 23     CommSui ~~     CommSui  6.110132e-01  0.232372791860449
    ## 65      UkbDep ~~      UkbDep -3.122257e-02  0.126346444888124
    ## 64      UkbAnh ~~      UkbAnh  1.873782e-01 0.0843719057183106
    ## 7   ClinAppDec ~~  ClinAppDec  6.117025e-01  0.251668474405079
    ## 8   ClinAppInc ~~  ClinAppInc  3.450639e-01  0.349094191797315
    ## 11  ClinSleDec ~~  ClinSleDec  9.930192e-01  0.678006870822972
    ## 12  ClinSleInc ~~  ClinSleInc  9.600528e-01  0.666508490724797
    ## 15  CommAppDec ~~  CommAppDec  8.057440e-01  0.237897751882888
    ## 16  CommAppInc ~~  CommAppInc  2.105294e-01  0.286398365825205
    ## 21  CommSleDec ~~  CommSleDec -5.683429e+02   412056.882270266
    ## 22  CommSleInc ~~  CommSleInc  6.674735e-01   0.30944044703178
    ## 19   CommFatig ~~   CommFatig  4.647260e-01  0.295714882086838
    ## 45       PSYCH ~~        SOMA -9.240060e-01  0.115064473100508
    ## 43       PSYCH ~~       PSYCH  1.000000e+00                   
    ## 63        SOMA ~~        SOMA  1.000000e+00                   
    ## 29        GATE ~~        GATE  1.000000e+00                   
    ## 5          APP ~~         APP  1.000000e+00                   
    ## 50         SLE ~~         SLE  1.000000e+00                   
    ## 42       PSYCH ~~        GATE  0.000000e+00                   
    ## 41       PSYCH ~~         APP  0.000000e+00                   
    ## 44       PSYCH ~~         SLE  0.000000e+00                   
    ## 61        SOMA ~~        GATE  0.000000e+00                   
    ## 60        SOMA ~~         APP  0.000000e+00                   
    ## 62        SOMA ~~         SLE  0.000000e+00                   
    ## 28        GATE ~~         APP  0.000000e+00                   
    ## 30        GATE ~~         SLE  0.000000e+00                   
    ## 6          APP ~~         SLE  0.000000e+00

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
    ##   6.297 
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
psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 31       PSYCH =~     ClinSui    0.58175613  0.106408512176395
    ## 34       PSYCH =~     CommDep    0.69446205 0.0696577576090295
    ## 32       PSYCH =~     CommAnh    0.81871690 0.0664908484357224
    ## 35       PSYCH =~   CommGuilt    0.68293698 0.0894787770564488
    ## 33       PSYCH =~    CommConc    0.80476244  0.099082817707634
    ## 36       PSYCH =~     CommSui    0.62378908 0.0990981965298188
    ## 38       PSYCH =~      UkbDep    0.54896770 0.0947575352040295
    ## 37       PSYCH =~      UkbAnh    0.69315991 0.0894562164219068
    ## 51         VEG =~  ClinAppDec    0.01984847 0.0956309119937124
    ## 52         VEG =~  ClinAppInc   -0.08569187  0.130334200391164
    ## 55         VEG =~  ClinSleDec   -0.08813394     0.156079318894
    ## 56         VEG =~  ClinSleInc   -0.21535040  0.179561618701526
    ## 54         VEG =~ ClinMotoInc   -0.24222051  0.132340022571323
    ## 53         VEG =~ ClinMotoDec   -0.27149511  0.206356936092365
    ## 57         VEG =~  CommAppDec   -0.24161384  0.101801918761671
    ## 58         VEG =~  CommAppInc   -0.42783562 0.0958184607928428
    ## 60         VEG =~  CommSleDec   -0.72570463  0.122644555712771
    ## 61         VEG =~  CommSleInc   -0.60942775  0.135777402050544
    ## 59         VEG =~   CommFatig   -0.76686929  0.117199150581111
    ## 25        GATE =~     CommDep    0.62960530 0.0965247715319635
    ## 24        GATE =~     CommAnh    0.51637388 0.0955363705962653
    ## 27        GATE =~      UkbDep    0.85807032  0.116600266537722
    ## 26        GATE =~      UkbAnh    0.57612056  0.116579703278469
    ## 1          APP =~  ClinAppDec    0.62408425  0.126494614565211
    ## 2          APP =~  ClinAppInc   -0.79982355  0.170544697794542
    ## 3          APP =~  CommAppDec    0.38392662  0.131761200277741
    ## 4          APP =~  CommAppInc   -0.78929421  0.176672718964991
    ## 44         SLE =~  ClinSleDec    0.04539645   3.65849717967441
    ## 45         SLE =~  ClinSleInc   -0.00267862   0.21837545272079
    ## 46         SLE =~  CommSleDec   11.68669371   942.096158902801
    ## 47         SLE =~  CommSleInc   -0.03926610   3.16593973530122
    ## 13     ClinSui ~~     ClinSui    0.66155961  0.322601192497847
    ## 18     CommDep ~~     CommDep    0.12132005 0.0529366574864233
    ## 14     CommAnh ~~     CommAnh    0.06306136 0.0438424386462029
    ## 20   CommGuilt ~~   CommGuilt    0.53359768  0.165880419755608
    ## 17    CommConc ~~    CommConc    0.35235802  0.241526079662097
    ## 23     CommSui ~~     CommSui    0.61088769  0.232365469947806
    ## 50      UkbDep ~~      UkbDep   -0.03765081  0.130039486670347
    ## 49      UkbAnh ~~      UkbAnh    0.18761432 0.0843292750856142
    ## 7   ClinAppDec ~~  ClinAppDec    0.61012609  0.251673905698184
    ## 8   ClinAppInc ~~  ClinAppInc    0.35293750  0.347233279396766
    ## 11  ClinSleDec ~~  ClinSleDec    0.99017320  0.657011115299653
    ## 12  ClinSleInc ~~  ClinSleInc    0.95361531  0.667162859595402
    ## 10 ClinMotoInc ~~ ClinMotoInc    0.94133080  0.392445179699618
    ## 9  ClinMotoDec ~~ ClinMotoDec    0.92629027   1.07943651847992
    ## 15  CommAppDec ~~  CommAppDec    0.79422257  0.239868525673042
    ## 16  CommAppInc ~~  CommAppInc    0.19397167  0.287589075742453
    ## 21  CommSleDec ~~  CommSleDec -136.10544150   22020.0052188246
    ## 22  CommSleInc ~~  CommSleInc    0.62705573  0.336230561398465
    ## 19   CommFatig ~~   CommFatig    0.41191043  0.299372730505467
    ## 43       PSYCH ~~         VEG   -0.85610564 0.0978104956103518
    ## 41       PSYCH ~~       PSYCH    1.00000000                   
    ## 65         VEG ~~         VEG    1.00000000                   
    ## 29        GATE ~~        GATE    1.00000000                   
    ## 5          APP ~~         APP    1.00000000                   
    ## 48         SLE ~~         SLE    1.00000000                   
    ## 40       PSYCH ~~        GATE    0.00000000                   
    ## 39       PSYCH ~~         APP    0.00000000                   
    ## 42       PSYCH ~~         SLE    0.00000000                   
    ## 63         VEG ~~        GATE    0.00000000                   
    ## 62         VEG ~~         APP    0.00000000                   
    ## 64         VEG ~~         SLE    0.00000000                   
    ## 28        GATE ~~         APP    0.00000000                   
    ## 30        GATE ~~         SLE    0.00000000                   
    ## 6          APP ~~         SLE    0.00000000

### Affective-Neurovegetative (Elhai Model 2c)

``` r
affect_neuroveg.model <- "
PSYCH =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommSui + UkbDep + UkbAnh
NEUROVEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec +  + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommConc

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

PSYCH ~~ 1*PSYCH
NEUROVEG ~~ 1*NEUROVEG
GATE ~~ 1*GATE
APP ~~ 1*APP
SLE ~~ 1*SLE
PSYCH ~~ 0*GATE + 0*APP + 0*SLE
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
    ##  10.646 
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
affect_neuroveg.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 47       PSYCH =~     ClinSui  5.768017e-01  0.106619629953544
    ## 49       PSYCH =~     CommDep  7.110583e-01 0.0775435933957163
    ## 48       PSYCH =~     CommAnh  8.365208e-01 0.0764986896616674
    ## 50       PSYCH =~   CommGuilt  6.797795e-01  0.089548861361782
    ## 51       PSYCH =~     CommSui  6.231153e-01 0.0988002500286814
    ## 53       PSYCH =~      UkbDep  5.604969e-01 0.0993477157840231
    ## 52       PSYCH =~      UkbAnh  7.078245e-01  0.096198618738819
    ## 31    NEUROVEG =~  ClinAppDec  1.943200e-02  0.091947672229563
    ## 32    NEUROVEG =~  ClinAppInc -7.735987e-02   0.12609941620127
    ## 35    NEUROVEG =~  ClinSleDec -8.428226e-02  0.150479757713121
    ## 36    NEUROVEG =~  ClinSleInc -2.054559e-01  0.172719772229238
    ## 34    NEUROVEG =~ ClinMotoInc -2.299386e-01  0.127784237414298
    ## 33    NEUROVEG =~ ClinMotoDec -2.644900e-01  0.200262740738586
    ## 37    NEUROVEG =~  CommAppDec -2.318555e-01 0.0991148257800385
    ## 38    NEUROVEG =~  CommAppInc -4.158165e-01 0.0938334101803996
    ## 41    NEUROVEG =~  CommSleDec -6.984393e-01   0.11295078145849
    ## 42    NEUROVEG =~  CommSleInc -5.803194e-01  0.124126620257122
    ## 40    NEUROVEG =~   CommFatig -7.378543e-01  0.110752453564182
    ## 39    NEUROVEG =~    CommConc -8.700379e-01  0.102109265570664
    ## 25        GATE =~     CommDep  6.090544e-01  0.110889523314977
    ## 24        GATE =~     CommAnh  4.910007e-01  0.114352722778641
    ## 27        GATE =~      UkbDep  8.601007e-01  0.122087309818489
    ## 26        GATE =~      UkbAnh  5.574465e-01   0.12819509295066
    ## 1          APP =~  ClinAppDec  6.233709e-01  0.126391079150481
    ## 2          APP =~  ClinAppInc -8.037523e-01  0.171442706927283
    ## 3          APP =~  CommAppDec  3.778443e-01  0.130395156411876
    ## 4          APP =~  CommAppInc -7.904912e-01  0.177040018529995
    ## 59         SLE =~  ClinSleDec  2.679652e-02   6.72789289277692
    ## 60         SLE =~  ClinSleInc -9.050982e-04  0.228146826783487
    ## 61         SLE =~  CommSleDec  1.999004e+01   5019.24823349776
    ## 62         SLE =~  CommSleInc -2.110702e-02   5.29983331753858
    ## 13     ClinSui ~~     ClinSui  6.672999e-01  0.322871189533751
    ## 18     CommDep ~~     CommDep  1.234489e-01 0.0538416003481985
    ## 14     CommAnh ~~     CommAnh  5.915129e-02 0.0453612515886917
    ## 20   CommGuilt ~~   CommGuilt  5.378998e-01  0.166332673404001
    ## 23     CommSui ~~     CommSui  6.117273e-01  0.232155094587155
    ## 65      UkbDep ~~      UkbDep -5.392999e-02  0.142074233709883
    ## 64      UkbAnh ~~      UkbAnh  1.882378e-01 0.0843390077570837
    ## 7   ClinAppDec ~~  ClinAppDec  6.110309e-01  0.251454542778138
    ## 8   ClinAppInc ~~  ClinAppInc  3.479979e-01  0.348375110740686
    ## 11  ClinSleDec ~~  ClinSleDec  9.921774e-01  0.671875022816404
    ## 12  ClinSleInc ~~  ClinSleInc  9.577874e-01   0.66699737172517
    ## 10 ClinMotoInc ~~ ClinMotoInc  9.471274e-01   0.39207598594071
    ## 9  ClinMotoDec ~~ ClinMotoDec  9.300449e-01   1.07774673285778
    ## 15  CommAppDec ~~  CommAppDec  8.034771e-01  0.238766638345283
    ## 16  CommAppInc ~~  CommAppInc  2.022203e-01  0.289602338779471
    ## 21  CommSleDec ~~  CommSleDec -3.990893e+02   200669.931433786
    ## 22  CommSleInc ~~  CommSleInc  6.627839e-01  0.308930617685936
    ## 19   CommFatig ~~   CommFatig  4.555711e-01  0.294946065937522
    ## 17    CommConc ~~    CommConc  2.430338e-01  0.250119172340718
    ## 56       PSYCH ~~    NEUROVEG -8.716159e-01 0.0963162072447003
    ## 57       PSYCH ~~       PSYCH  1.000000e+00                   
    ## 45    NEUROVEG ~~    NEUROVEG  1.000000e+00                   
    ## 29        GATE ~~        GATE  1.000000e+00                   
    ## 5          APP ~~         APP  1.000000e+00                   
    ## 63         SLE ~~         SLE  1.000000e+00                   
    ## 55       PSYCH ~~        GATE  0.000000e+00                   
    ## 54       PSYCH ~~         APP  0.000000e+00                   
    ## 58       PSYCH ~~         SLE  0.000000e+00                   
    ## 44    NEUROVEG ~~        GATE  0.000000e+00                   
    ## 43    NEUROVEG ~~         APP  0.000000e+00                   
    ## 46    NEUROVEG ~~         SLE  0.000000e+00                   
    ## 28        GATE ~~         APP  0.000000e+00                   
    ## 30        GATE ~~         SLE  0.000000e+00                   
    ## 6          APP ~~         SLE  0.000000e+00

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
    ##   7.078 
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
cog_mood_neuroveg.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 14         COG =~     ClinSui  5.935415e-01  0.112726409271895
    ## 16         COG =~   CommGuilt  7.202635e-01  0.274455215718593
    ## 15         COG =~    CommConc  8.189770e-01  0.110637354339811
    ## 17         COG =~     CommSui  6.370582e-01   0.10377120723321
    ## 34        MOOD =~     CommDep  9.304694e-01 0.0483742827172213
    ## 33        MOOD =~     CommAnh  1.005890e+00 0.0393782875413622
    ## 35        MOOD =~   CommGuilt -2.046430e-02  0.262088529175322
    ## 37        MOOD =~      UkbDep  8.845797e-01 0.0636038860729348
    ## 36        MOOD =~      UkbAnh  9.003483e-01  0.055777682556211
    ## 49         VEG =~  ClinAppDec  1.682393e-02 0.0967530334234969
    ## 50         VEG =~  ClinAppInc -7.869966e-02  0.131388032460375
    ## 53         VEG =~  ClinSleDec -8.680981e-02  0.158813147594721
    ## 54         VEG =~  ClinSleInc -2.106561e-01  0.181349797051195
    ## 52         VEG =~ ClinMotoInc -2.423954e-01   0.13391459336253
    ## 51         VEG =~ ClinMotoDec -2.651344e-01  0.207566158688483
    ## 55         VEG =~  CommAppDec -2.414672e-01  0.102941892932068
    ## 56         VEG =~  CommAppInc -4.248922e-01 0.0965631568040698
    ## 58         VEG =~  CommSleDec -7.297881e-01  0.123610481885544
    ## 59         VEG =~  CommSleInc -6.125996e-01   0.13722133945583
    ## 57         VEG =~   CommFatig -7.649797e-01  0.118131502647133
    ## 1          APP =~  ClinAppDec  6.234312e-01  0.126208893706376
    ## 2          APP =~  ClinAppInc -8.011289e-01  0.170665923169305
    ## 3          APP =~  CommAppDec  3.816154e-01  0.131649861511371
    ## 4          APP =~  CommAppInc -7.923287e-01  0.177388464939124
    ## 42         SLE =~  ClinSleDec  3.945978e-02   4.17218851143763
    ## 43         SLE =~  ClinSleInc -2.126389e-03  0.226736268743274
    ## 44         SLE =~  CommSleDec  1.346070e+01   1423.56739075866
    ## 45         SLE =~  CommSleInc -3.444768e-02   3.64360593829778
    ## 13     ClinSui ~~     ClinSui  6.477085e-01  0.329572428937458
    ## 29   CommGuilt ~~   CommGuilt  5.026794e-01  0.207370911155568
    ## 26    CommConc ~~    CommConc  3.292767e-01  0.256522508922798
    ## 32     CommSui ~~     CommSui  5.941569e-01    0.2329927337406
    ## 27     CommDep ~~     CommDep  1.342266e-01   0.05950816673173
    ## 23     CommAnh ~~     CommAnh -1.181476e-02 0.0506752586916915
    ## 48      UkbDep ~~      UkbDep  2.175189e-01  0.106314799882181
    ## 47      UkbAnh ~~      UkbAnh  1.893729e-01 0.0964506016517212
    ## 7   ClinAppDec ~~  ClinAppDec  6.110506e-01  0.251175301785318
    ## 8   ClinAppInc ~~  ClinAppInc  3.519992e-01  0.347647652137485
    ## 11  ClinSleDec ~~  ClinSleDec  9.909069e-01  0.655441836175985
    ## 12  ClinSleInc ~~  ClinSleInc  9.556195e-01  0.667307368654419
    ## 10 ClinMotoInc ~~ ClinMotoInc  9.412446e-01    0.3922726078707
    ## 9  ClinMotoDec ~~ ClinMotoDec  9.297039e-01   1.07896396495165
    ## 24  CommAppDec ~~  CommAppDec  7.960632e-01  0.240110246397063
    ## 25  CommAppInc ~~  CommAppInc  1.916820e-01  0.290124852045905
    ## 30  CommSleDec ~~  CommSleDec -1.807230e+02   38324.4489504604
    ## 31  CommSleInc ~~  CommSleInc  6.235350e-01  0.339428416475126
    ## 28   CommFatig ~~   CommFatig  4.148061e-01  0.300633992732814
    ## 20         COG ~~        MOOD  7.421412e-01 0.0910344207381464
    ## 22         COG ~~         VEG -8.168233e-01  0.111566493647679
    ## 41        MOOD ~~         VEG -6.492587e-01 0.0863499256625542
    ## 19         COG ~~         COG  1.000000e+00                   
    ## 39        MOOD ~~        MOOD  1.000000e+00                   
    ## 62         VEG ~~         VEG  1.000000e+00                   
    ## 5          APP ~~         APP  1.000000e+00                   
    ## 46         SLE ~~         SLE  1.000000e+00                   
    ## 18         COG ~~         APP  0.000000e+00                   
    ## 21         COG ~~         SLE  0.000000e+00                   
    ## 38        MOOD ~~         APP  0.000000e+00                   
    ## 40        MOOD ~~         SLE  0.000000e+00                   
    ## 60         VEG ~~         APP  0.000000e+00                   
    ## 61         VEG ~~         SLE  0.000000e+00                   
    ## 6          APP ~~         SLE  0.000000e+00

### Model comparisons

``` r
model_list=
list("1a"=list(name="Common", model=commonfactor.fit),
     "1b"=list(name="Common (App)", model=commonfactor_app.fit),
     "1c"=list(name="Common (Sle)", model=commonfactor_sle.fit),
     "1d"=list(name="Common (gating)", model=gate.fit),
     "1e"=list(name="Ascertainment", model=clin_comm.fit),
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
