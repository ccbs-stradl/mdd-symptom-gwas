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
    ## minor          1.0                         
    ## year           2021                        
    ## month          05                          
    ## day            18                          
    ## svn rev        80317                       
    ## language       R                           
    ## version.string R version 4.1.0 (2021-05-18)
    ## nickname       Camp Pontanezen

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot', 'mvtnorm')
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

    ## [1] '0.0.3'

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
MDD5;Motor⇅;Motor⇆;Psyc
MDD5a;Motor⇈;Motor⇉;PsycInc
MDD5b;Motor⇊;Motor⇇;PsycDec
MDD6;Fatigue;Fatigue;Fatig
MDD7;Guilt;Guilt;Guilt
MDD8;Concentrate;Concentrate;Conc
MDD9;Suicidality;Suicidality;Sui
", col_names=c('ref', 'h', 'v', 'abbv'), delim=';')

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
MDD5b;Psychomotor retardation nearly every day
MDD6;Fatigue or loss of energy nearly every day
MDD7;Feelings of worthlessness or excessive or inappropriate guilt
MDD8;Diminished ability to think or concentrate, or indecisiveness
MDD9;Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for attempting suicide
", col_names=c('Reference', 'Description'), delim=';')

dsm_mdd_symptoms_reference %>%
left_join(dsm_mdd_symptoms_labels, by=c('Reference'='ref')) %>%
select(Reference, Abbreviation=abbv, Label=h, Description)
```

    ## # A tibble: 15 x 4
    ##    Reference Abbreviation Label     Description                                 
    ##    <chr>     <chr>        <chr>     <chr>                                       
    ##  1 MDD1      Dep          Mood      Depressed mood most of the day, nearly ever…
    ##  2 MDD2      Anh          Interest  Markedly diminished interest or pleasure in…
    ##  3 MDD3      App          Weight⇅   Significant change in weight or appetite    
    ##  4 MDD3a     AppDec       Weight⇊   Significant weight loss or decrease in appe…
    ##  5 MDD3b     AppInc       Weight⇈   Significant weight gain or increase in appe…
    ##  6 MDD4      Sle          Sleep⇅    Sleeping too much or not sleeping enough    
    ##  7 MDD4a     SleDec       Sleep⇊    Insomnia nearly every day                   
    ##  8 MDD4b     SleInc       Sleep⇈    Hypersomnia nearly every day                
    ##  9 MDD5      Psyc         Motor⇅    Changes in speed/amount of moving or speaki…
    ## 10 MDD5a     PsycInc      Motor⇈    Psychomotor agitation nearly every day      
    ## 11 MDD5b     PsycDec      Motor⇊    Psychomotor retardation nearly every day    
    ## 12 MDD6      Fatig        Fatigue   Fatigue or loss of energy nearly every day  
    ## 13 MDD7      Guilt        Guilt     Feelings of worthlessness or excessive or i…
    ## 14 MDD8      Conc         Concentr… Diminished ability to think or concentrate,…
    ## 15 MDD9      Sui          Suicidal… Recurrent thoughts of death or suicide or a…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'hdl.covstruct', 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'hdl.covstruct', 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'covstruct', 'prevs', 'txt', sep='.')))
```

    ## 
    ## ── Column specification ───────────────────────────────────────────────────────────────
    ## cols(
    ##   cohorts = col_character(),
    ##   symptom = col_character(),
    ##   sumstats = col_character(),
    ##   Nca = col_double(),
    ##   Nco = col_double(),
    ##   samp_prev = col_double(),
    ##   filename = col_character(),
    ##   pop_prev = col_double(),
    ##   trait_name = col_character()
    ## )

Rename samples: AGDS/PGC is the **Clin**ical sample (`Clin`) and
ALSPAC/UKB is the **Pop**ulation sample (`Pop`); and rename symptoms
numbers (`MDD1`, `MDD2`) to abbreviations (`Dep`, `Anh`)

``` r
cohorts_sample_symptoms <-
sumstats_prevs %>%
left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
select(cohorts, symptom, trait_name, abbv) %>%
mutate(Sample=case_when(cohorts %in% 'AGDS_PGC' ~ 'Clin',
                        cohorts %in% 'ALSPAC_UKB' ~ 'Pop',
                        TRUE ~ NA_character_)) %>%
mutate(sample_symptom=paste0(Sample, abbv))

sample_symptoms <- cohorts_sample_symptoms$sample_symptom
names(sample_symptoms) <- cohorts_sample_symptoms$trait_name

# rename traits in covstruct
dimnames(symptoms_covstruct$S)[[2]] <-
as.vector(sample_symptoms[dimnames(symptoms_covstruct$S)[[2]]])
```

# Structural models

## ADGS-PGC

### Common factor

Common factor model. Allow residual negative correlation between
directional symptoms

``` r
pgc_commonfactor.model <- "
A1 =~ NA*ClinDep + ClinAnh + ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinPsycDec + ClinFatig + ClinGuilt + ClinConc + ClinSui
A1 ~~ 1*A1
"
pgc_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pgc_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.913 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.42189770725478 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

    ## [1] "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  1.68614501877486e-10 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

``` r
pgc_commonfactor.fit$modelfit
```

    ##       chisq df p_chisq      AIC       CFI      SRMR
    ## df 38626871 54       0 38626919 0.4890952 0.2197145

``` r
pgc_commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 5           A1 =~     ClinDep   0.10772008   0.13145991597515
    ## 1           A1 =~     ClinAnh   0.28703562  0.113291047642627
    ## 2           A1 =~  ClinAppDec  -0.18687351  0.151003216293719
    ## 3           A1 =~  ClinAppInc   0.71470363  0.111063838503204
    ## 10          A1 =~  ClinSleDec  -0.22700067  0.160199861461656
    ## 11          A1 =~  ClinSleInc   0.69343278  0.110461913940812
    ## 9           A1 =~ ClinPsycInc   0.03478905  0.200030350852541
    ## 8           A1 =~ ClinPsycDec   0.75910840  0.133109541501726
    ## 6           A1 =~   ClinFatig   0.77438545  0.182155523753146
    ## 7           A1 =~   ClinGuilt   0.47142969  0.140828427843276
    ## 4           A1 =~    ClinConc   0.89811135  0.199716226258574
    ## 12          A1 =~     ClinSui   0.62538339  0.116823341379837
    ## 18     ClinDep ~~     ClinDep   0.98839606 0.0620626508567864
    ## 14     ClinAnh ~~     ClinAnh   0.91761048 0.0710719097240964
    ## 15  ClinAppDec ~~  ClinAppDec   0.96507875  0.144087191913337
    ## 16  ClinAppInc ~~  ClinAppInc   0.48919843  0.161907127482506
    ## 23  ClinSleDec ~~  ClinSleDec   0.94847071  0.154475177721915
    ## 24  ClinSleInc ~~  ClinSleInc   0.51915046   0.18202907446396
    ## 22 ClinPsycInc ~~ ClinPsycInc   0.99879213  0.283816983116301
    ## 21 ClinPsycDec ~~ ClinPsycDec   0.42375454  0.205051998168403
    ## 19   ClinFatig ~~   ClinFatig   0.40032774  0.301414558866454
    ## 20   ClinGuilt ~~   ClinGuilt   0.77775395  0.223071972695511
    ## 17    ClinConc ~~    ClinConc   0.19339599  0.361144441224099
    ## 25     ClinSui ~~     ClinSui   0.60889566  0.184808243871065
    ## 13          A1 ~~          A1   1.00000000

Add negative correlations for directional symptoms

``` r
pgc_commonfactor_dir.model <- "
A1 =~ NA*ClinDep + ClinAnh + ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinPsycDec + ClinFatig + ClinGuilt + ClinConc + ClinSui
A1 ~~ 1*A1
ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
"
pgc_commonfactor_dir.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pgc_commonfactor_dir.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.444 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.42189770725478 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor_dir.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor_dir.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

    ## [1] "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  1.68614501877486e-10 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

``` r
pgc_commonfactor_dir.fit$modelfit
```

    ##       chisq df p_chisq      AIC      CFI      SRMR
    ## df 33158398 52       0 33158450 0.561425 0.2119889

``` r
pgc_commonfactor_dir.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 5           A1 =~     ClinDep  0.124528403   0.13460966204486 3.549128e-01
    ## 1           A1 =~     ClinAnh  0.293514640  0.111682949794449 8.586038e-03
    ## 2           A1 =~  ClinAppDec -0.008169196  0.162051031280854 9.597974e-01
    ## 3           A1 =~  ClinAppInc  0.694075213  0.108709760210574 1.717729e-10
    ## 10          A1 =~  ClinSleDec -0.121327635  0.181859670780267 5.046725e-01
    ## 11          A1 =~  ClinSleInc  0.689188399  0.111721551186912 6.880520e-10
    ## 9           A1 =~ ClinPsycInc  0.071393093  0.204386232652756 7.268515e-01
    ## 8           A1 =~ ClinPsycDec  0.791210451  0.132732380564229 2.507818e-09
    ## 6           A1 =~   ClinFatig  0.761729764  0.182914383108398 3.121777e-05
    ## 7           A1 =~   ClinGuilt  0.465558294   0.14124871674545 9.806205e-04
    ## 4           A1 =~    ClinConc  0.901957593  0.197438629519415 4.916792e-06
    ## 12          A1 =~     ClinSui  0.636074833  0.120028324533571 1.161946e-07
    ## 16  ClinAppDec ~~  ClinAppInc -0.575536618   0.19181717486919 2.695968e-03
    ## 25  ClinSleDec ~~  ClinSleInc -0.314044816   0.24252058241091 1.953480e-01
    ## 19     ClinDep ~~     ClinDep  0.984492465 0.0634553827283735 2.756618e-54
    ## 14     ClinAnh ~~     ClinAnh  0.913849302 0.0715258652811533 2.219170e-37
    ## 15  ClinAppDec ~~  ClinAppDec  0.999933663  0.122064777429076 2.572910e-16
    ## 17  ClinAppInc ~~  ClinAppInc  0.518259940  0.155698461140796 8.727777e-04
    ## 24  ClinSleDec ~~  ClinSleDec  0.985279541  0.147077621322504 2.097862e-11
    ## 26  ClinSleInc ~~  ClinSleInc  0.525018564  0.185583541685131 4.669206e-03
    ## 23 ClinPsycInc ~~ ClinPsycInc  0.994902602   0.28794342437423 5.498787e-04
    ## 22 ClinPsycDec ~~ ClinPsycDec  0.373984596  0.211165905922801 7.655113e-02
    ## 20   ClinFatig ~~   ClinFatig  0.419767636  0.298714753183607 1.599429e-01
    ## 21   ClinGuilt ~~   ClinGuilt  0.783256635   0.22274930042823 4.376109e-04
    ## 18    ClinConc ~~    ClinConc  0.186472159  0.358605714970313 6.030777e-01
    ## 27     ClinSui ~~     ClinSui  0.595408979  0.188636764347223 1.597502e-03
    ## 13          A1 ~~          A1  1.000000000                              NA

### Two-factor models

[Elhai Psychiat Res
2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Kruse Rehab Psychol
2008](https://psycnet.apa.org/record/2008-17022-011), [Kruse Arch Psys
Med Rehab
2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> retardation as a fourth somatic item

``` r
clin_psych_soma.model <- "
A1 =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui + ClinPsycDec + ClinPsycInc
A2 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
a12 > -1
a12 < 1
A1 ~~ a12*A2
"
clin_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_psych_soma.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##            lhs op         rhs Unstandardized_Estimate          SE
    ## 1           A1 =~     ClinDep             0.085707463 0.092620532
    ## 2           A1 =~     ClinAnh             0.195201558 0.073450946
    ## 3           A1 =~   ClinGuilt             0.144754560 0.044523206
    ## 4           A1 =~    ClinConc             0.327052163 0.071992936
    ## 5           A1 =~     ClinSui             0.169717414 0.034193966
    ## 6           A1 =~ ClinPsycDec             0.170590869 0.029092217
    ## 7           A1 =~ ClinPsycInc             0.008911857 0.025516895
    ## 8           A2 =~  ClinAppDec             0.002080627 0.041277209
    ## 9           A2 =~  ClinAppInc            -0.244341981 0.043644379
    ## 10          A2 =~  ClinSleDec             0.032124696 0.048232139
    ## 11          A2 =~  ClinSleInc            -0.172024564 0.028611869
    ## 12          A2 =~   ClinFatig            -0.294971935 0.070838118
    ## 15  ClinAppDec ~~  ClinAppInc            -0.051602643 0.017196571
    ## 16  ClinSleDec ~~  ClinSleInc            -0.020755227 0.015992263
    ## 17          A1 ~~          A2            -1.000000057 0.180430003
    ## 18     ClinDep ~~     ClinDep             0.466368435 0.030068636
    ## 19     ClinAnh ~~     ClinAnh             0.404188983 0.031444197
    ## 20   ClinGuilt ~~   ClinGuilt             0.075721246 0.021812872
    ## 21    ClinConc ~~    ClinConc             0.024516007 0.047185781
    ## 22     ClinSui ~~     ClinSui             0.042388754 0.013549763
    ## 23 ClinPsycDec ~~ ClinPsycDec             0.017385346 0.010254641
    ## 24 ClinPsycInc ~~ ClinPsycInc             0.015502792 0.004488919
    ## 25  ClinAppDec ~~  ClinAppDec             0.064861528 0.007917880
    ## 26  ClinAppInc ~~  ClinAppInc             0.064228972 0.021526432
    ## 27  ClinSleDec ~~  ClinSleDec             0.069075035 0.010319988
    ## 28  ClinSleInc ~~  ClinSleInc             0.032710091 0.011616897
    ## 29   ClinFatig ~~   ClinFatig             0.062945493 0.044866499

``` r
clin_psych_soma.fit$modelfit
```

    ## NULL

``` r
clin_psych_soma.fit$results[c(1, 2, 3, 6, 7, 9)]
```

    ## NULL

Bifactor model

``` r
clin_psych_soma_bif.model <- "
A =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui + ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycDec + ClinPsycInc + ClinFatig
A1 =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui + ClinPsycDec + ClinPsycInc
A2 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A ~~ 1*A
A ~~ 0*A1 + 0*A2
A1 ~~ 0*A2
c2 > 0.001
ClinAnh ~~ c2*ClinAnh
c3a > 0.001
ClinAppDec ~~ c3a*ClinAppDec
"
clin_psych_soma_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_psych_soma_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.815 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.42189770725478 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_soma_bif.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_soma_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

    ## [1] "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  1.68614501877486e-10 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

``` r
clin_psych_soma_bif.fit$modelfit
```

    ##       chisq df p_chisq      AIC       CFI      SRMR
    ## df 29822148 42       0 29822220 0.6055524 0.1755414

``` r
clin_psych_soma_bif.fit$results[c(1, 2, 3, 6, 7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 5            A =~     ClinDep  0.261243167 0.199597954322161 1.911013e-01
    ## 1            A =~     ClinAnh  0.479802986 0.138733714495713 5.497690e-04
    ## 7            A =~   ClinGuilt  0.459214099 0.139289924025651 9.761203e-04
    ## 4            A =~    ClinConc  0.895997789 0.209224440396602 1.854278e-05
    ## 12           A =~     ClinSui  0.645479099 0.130232011650598 7.185017e-07
    ## 2            A =~  ClinAppDec  0.100187950 0.162345103860735 5.410357e-01
    ## 3            A =~  ClinAppInc  0.619794426 0.113741104968002 5.189943e-08
    ## 10           A =~  ClinSleDec -0.191293272 0.176899004289822 2.807452e-01
    ## 11           A =~  ClinSleInc  0.670295703 0.122657310149091 4.903257e-08
    ## 8            A =~ ClinPsycDec  0.813276811 0.134190441516677 1.368807e-09
    ## 9            A =~ ClinPsycInc  0.095136255 0.208549353486906 6.486268e-01
    ## 6            A =~   ClinFatig  0.636771714 0.199826188074397 1.449618e-03
    ## 18          A1 =~     ClinDep  0.946339079 0.355492890268405 7.795314e-03
    ## 16          A1 =~     ClinAnh  0.877129150 0.358695337979806 1.448911e-02
    ## 19          A1 =~   ClinGuilt -0.228244441  0.19310097769283 2.375865e-01
    ## 17          A1 =~    ClinConc -0.421156674  0.27884950596739 1.312344e-01
    ## 22          A1 =~     ClinSui -0.207352657 0.157545051798138 1.886539e-01
    ## 20          A1 =~ ClinPsycDec -0.245405801 0.180794851682091 1.750739e-01
    ## 21          A1 =~ ClinPsycInc  0.157386749 0.209095097033484 4.514479e-01
    ## 25          A2 =~  ClinAppDec  1.002757975 0.454705503827001 2.724441e-02
    ## 26          A2 =~  ClinAppInc -0.532828312 0.279374549819039 5.619780e-02
    ## 28          A2 =~  ClinSleDec  0.257876902 0.222263152380692 2.468481e-01
    ## 29          A2 =~  ClinSleInc -0.127704311 0.199898920152302 5.243900e-01
    ## 27          A2 =~   ClinFatig -0.410398479 0.252421634979461 1.041867e-01
    ## 31     ClinAnh ~~     ClinAnh  0.001000079 0.640690402266755 9.971816e-01
    ## 32  ClinAppDec ~~  ClinAppDec  0.001000108 0.909747678548447 9.862638e-01
    ## 35     ClinDep ~~     ClinDep  0.036193895 0.713295199500809 9.608617e-01
    ## 37   ClinGuilt ~~   ClinGuilt  0.737026388  0.24001089773347 2.140459e-03
    ## 34    ClinConc ~~    ClinConc  0.019815032  0.41880475702145 9.626582e-01
    ## 42     ClinSui ~~     ClinSui  0.540361460 0.217563800900229 1.301150e-02
    ## 33  ClinAppInc ~~  ClinAppInc  0.331948318 0.298534550133766 2.713191e-01
    ## 40  ClinSleDec ~~  ClinSleDec  0.896905627 0.170046959485865 1.357319e-07
    ## 41  ClinSleInc ~~  ClinSleInc  0.534394916 0.179426514992822 2.906891e-03
    ## 38 ClinPsycDec ~~ ClinPsycDec  0.278356722 0.251727271040841 2.691783e-01
    ## 39 ClinPsycInc ~~ ClinPsycInc  0.966177650 0.299793007361215 1.269094e-03
    ## 36   ClinFatig ~~   ClinFatig  0.426094752 0.293818379651391 1.470549e-01
    ## 23          A1 ~~          A1  1.000000000                             NA
    ## 30          A2 ~~          A2  1.000000000                             NA
    ## 13           A ~~           A  1.000000000                             NA
    ## 14           A ~~          A1  0.000000000                             NA
    ## 15           A ~~          A2  0.000000000                             NA
    ## 24          A1 ~~          A2  0.000000000                             NA

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
clin_psych_veg.model <- "
A1 =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui 
A2 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycDec + ClinPsycInc + ClinFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
a12 > -1
a12 < 1
A1 ~~ a12*A2
"
clin_psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_psych_veg.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.391 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.42189770725478 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_veg.model, : A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_veg.model, : A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

    ## [1] "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  1.68614501877486e-10 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

``` r
clin_psych_veg.fit$modelfit
```

    ##       chisq df p_chisq      AIC       CFI     SRMR
    ## df 33158461 51       0 33158515 0.5614241 0.211989

``` r
clin_psych_veg.fit$results[c(1, 2, 3, 6, 7, 9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 3           A1 =~     ClinDep  0.124526879  0.136454334384617 3.614530e-01
    ## 1           A1 =~     ClinAnh  0.293513563  0.114162506438981 1.013987e-02
    ## 4           A1 =~   ClinGuilt  0.465560472   0.14254120101459 1.090275e-03
    ## 2           A1 =~    ClinConc  0.901971860  0.242927069169378 2.048998e-04
    ## 5           A1 =~     ClinSui  0.636076853  0.151437406612846 2.666346e-05
    ## 8           A2 =~  ClinAppDec  0.008169072  0.162036410795969 9.597948e-01
    ## 9           A2 =~  ClinAppInc -0.694073746  0.110087846422189 2.887187e-10
    ## 13          A2 =~  ClinSleDec  0.121326475  0.181810097784469 5.045591e-01
    ## 14          A2 =~  ClinSleInc -0.689186963  0.110459722951935 4.396636e-10
    ## 11          A2 =~ ClinPsycDec -0.791210743  0.130655040157318 1.397987e-09
    ## 12          A2 =~ ClinPsycInc -0.071392821  0.204456218890843 7.269505e-01
    ## 10          A2 =~   ClinFatig -0.761729606  0.182842715549095 3.099452e-05
    ## 18  ClinAppDec ~~  ClinAppInc -0.575535493  0.191762354489261 2.688323e-03
    ## 27  ClinSleDec ~~  ClinSleInc -0.314047508  0.242651496315845 1.955874e-01
    ## 7           A1 ~~          A2 -1.000000533  0.199438382793882 5.329650e-07
    ## 21     ClinDep ~~     ClinDep  0.984493103 0.0632113081840991 1.083658e-54
    ## 16     ClinAnh ~~     ClinAnh  0.913849908 0.0725014038544077 1.994455e-36
    ## 23   ClinGuilt ~~   ClinGuilt  0.783251518  0.219756441000349 3.649714e-04
    ## 20    ClinConc ~~    ClinConc  0.186446647  0.435724747884466 6.686949e-01
    ## 29     ClinSui ~~     ClinSui  0.595403702  0.219445907711995 6.662993e-03
    ## 17  ClinAppDec ~~  ClinAppDec  0.999933240  0.122063901827155 2.571636e-16
    ## 19  ClinAppInc ~~  ClinAppInc  0.518261500  0.154603311556676 8.017835e-04
    ## 26  ClinSleDec ~~  ClinSleDec  0.985279984  0.147112583867526 2.120873e-11
    ## 28  ClinSleInc ~~  ClinSleInc  0.525020973  0.185282628558927 4.602519e-03
    ## 24 ClinPsycDec ~~ ClinPsycDec  0.373984746  0.211519543577415 7.704455e-02
    ## 25 ClinPsycInc ~~ ClinPsycInc  0.994903171  0.287968597503438 5.504977e-04
    ## 22   ClinFatig ~~   ClinFatig  0.419767997  0.298413463132903 1.595262e-01
    ## 6           A1 ~~          A1  1.000000000                              NA
    ## 15          A2 ~~          A2  1.000000000                              NA

Bifactor model

``` r
clin_psych_veg_bif.model <- "
A =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui + ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycDec + ClinPsycInc + ClinFatig
A1 =~ NA*ClinDep + ClinAnh + ClinGuilt + ClinConc + ClinSui 
A2 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycDec + ClinPsycInc + ClinFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A ~~ 1*A
A ~~ 0*A1 + 0*A2
A1 ~~ 0*A2
c2 > 0.001
ClinAnh ~~ c2*ClinAnh
c3a > 0.001
ClinAppDec ~~ c3a*ClinAppDec
"
clin_psych_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_psych_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.893 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.42189770725478 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_veg_bif.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_psych_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

    ## [1] "The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  1.68614501877486e-10 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  4.0842542727239 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

``` r
clin_psych_veg_bif.fit$modelfit
```

    ##       chisq df p_chisq      AIC       CFI      SRMR
    ## df 31149197 42       0 31149269 0.5879999 0.1555306

``` r
clin_psych_veg_bif.fit$results[c(1, 2, 3, 6, 7, 9)]
```

    ##            lhs op         rhs  STD_Genotype   STD_Genotype_SE      p_value
    ## 5            A =~     ClinDep  0.1148439700 0.154119840582366 4.575677e-01
    ## 1            A =~     ClinAnh  0.3406413769 0.110571446008195 2.082780e-03
    ## 7            A =~   ClinGuilt  0.4856483784 0.142670224319411 6.636223e-04
    ## 4            A =~    ClinConc  0.9338161843 0.208680176675807 7.668842e-06
    ## 12           A =~     ClinSui  0.6618669191 0.128549760341861 2.639343e-07
    ## 2            A =~  ClinAppDec -0.0013181410 0.180704595272179 9.891214e-01
    ## 3            A =~  ClinAppInc  0.6573794014 0.118046201001687 2.611930e-08
    ## 10           A =~  ClinSleDec -0.1390705516 0.189252312277256 4.651586e-01
    ## 11           A =~  ClinSleInc  0.6694349668 0.119719462853956 2.358026e-08
    ## 8            A =~ ClinPsycDec  0.8090050710 0.135412734819458 2.328716e-09
    ## 9            A =~ ClinPsycInc  0.1984795613  0.25049266153501 4.266322e-01
    ## 6            A =~   ClinFatig  0.6888497157 0.194679518534647 4.055047e-04
    ## 18          A1 =~     ClinDep  0.9705635058 0.541499794291101 7.342517e-02
    ## 16          A1 =~     ClinAnh  0.9398070987 0.521038644369331 7.160779e-02
    ## 19          A1 =~   ClinGuilt -0.1576589930 0.192570188634872 4.136454e-01
    ## 17          A1 =~    ClinConc -0.2966395332  0.26051833577976 2.552490e-01
    ## 20          A1 =~     ClinSui -0.1091936910 0.132406429094967 4.106600e-01
    ## 23          A2 =~  ClinAppDec  1.0050632740 0.351571189013089 4.101048e-03
    ## 24          A2 =~  ClinAppInc -0.4345720014 0.219835315913071 4.805683e-02
    ## 28          A2 =~  ClinSleDec  0.3323753076 0.227476287781106 1.425181e-01
    ## 29          A2 =~  ClinSleInc -0.1162837502 0.184121795363081 5.248825e-01
    ## 26          A2 =~ ClinPsycDec  0.0920573247 0.195176357574437 6.347477e-01
    ## 27          A2 =~ ClinPsycInc  0.5154025017 0.350554105795028 1.402606e-01
    ## 25          A2 =~   ClinFatig -0.3730371737 0.271830855121083 1.700814e-01
    ## 31     ClinAnh ~~     ClinAnh  0.0009999963 0.985005788930005 9.981687e-01
    ## 32  ClinAppDec ~~  ClinAppDec  0.0009999932 0.697308934490187 9.820365e-01
    ## 35     ClinDep ~~     ClinDep  0.0448165208  1.05568395498905 9.671590e-01
    ## 37   ClinGuilt ~~   ClinGuilt  0.7392895937 0.240432482738937 2.106785e-03
    ## 34    ClinConc ~~    ClinConc  0.0399924189 0.427852908524838 9.257014e-01
    ## 42     ClinSui ~~     ClinSui  0.5500089867 0.214389527586824 1.027730e-02
    ## 33  ClinAppInc ~~  ClinAppInc  0.3789996280 0.232470562905145 1.024854e-01
    ## 40  ClinSleDec ~~  ClinSleDec  0.8701865027 0.188526843241306 4.521602e-06
    ## 41  ClinSleInc ~~  ClinSleInc  0.5383345941 0.181991829250475 3.071601e-03
    ## 38 ClinPsycDec ~~ ClinPsycDec  0.3370361870 0.225205445290888 1.370045e-01
    ## 39 ClinPsycInc ~~ ClinPsycInc  0.6949569516 0.553911974127128 2.143097e-01
    ## 36   ClinFatig ~~   ClinFatig  0.3863289318 0.310350615169514 2.150418e-01
    ## 21          A1 ~~          A1  1.0000000000                             NA
    ## 30          A2 ~~          A2  1.0000000000                             NA
    ## 13           A ~~           A  1.0000000000                             NA
    ## 14           A ~~          A1  0.0000000000                             NA
    ## 15           A ~~          A2  0.0000000000                             NA
    ## 22          A1 ~~          A2  0.0000000000                             NA

## ALSPAC-UKB (Population)

### Common factor

Common factor model

``` r
pop_commonfactor.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
"
pop_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.184 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_commonfactor.fit$modelfit
```

    ##     chisq df      p_chisq    AIC       CFI      SRMR
    ## df 245.29 35 2.128852e-33 285.29 0.8227773 0.1111411

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep   0.73459939 0.0740541788887663
    ## 1         A1 =~    PopAnh   0.95016705 0.0663930755586976
    ## 2         A1 =~ PopAppDec   0.03904363 0.0974078842154991
    ## 3         A1 =~ PopAppInc   0.41273431 0.0952632993430382
    ## 8         A1 =~ PopSleDec   0.51548333  0.117153315316911
    ## 9         A1 =~ PopSleInc   0.49496516 0.0863292719868153
    ## 6         A1 =~  PopFatig   0.61717106  0.120138231313964
    ## 7         A1 =~  PopGuilt   0.57586175 0.0850086537778713
    ## 4         A1 =~   PopConc   0.79005037  0.103433356446807
    ## 10        A1 =~    PopSui   0.41087866  0.128534540529035
    ## 16    PopDep ~~    PopDep   0.46036450 0.0972048661276899
    ## 12    PopAnh ~~    PopAnh   0.09718340   0.11218421702619
    ## 13 PopAppDec ~~ PopAppDec   0.99847587   0.17091830590454
    ## 14 PopAppInc ~~ PopAppInc   0.82965050  0.148711128277044
    ## 19 PopSleDec ~~ PopSleDec   0.73427721  0.217346704206379
    ## 20 PopSleInc ~~ PopSleInc   0.75500932  0.177379711076552
    ## 17  PopFatig ~~  PopFatig   0.61910017  0.177071556988041
    ## 18  PopGuilt ~~  PopGuilt   0.66838399  0.156451870798405
    ## 15   PopConc ~~   PopConc   0.37582099  0.201340106328164
    ## 21    PopSui ~~    PopSui   0.83117828  0.174781904705785
    ## 11        A1 ~~        A1   1.00000000

Remove common variance shared between the gating items (Mood:
`UKB_CIDI1`, Interest: `UKB_CIDI2`) that is uncorrelated with the common
factor variance, to recover the genetic structure among gated items

``` r
pop_commonfactor_gating.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_commonfactor_gating.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_gating.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.712 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_gating.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_gating.fit$modelfit
```

    ##       chisq df      p_chisq      AIC      CFI      SRMR
    ## df 188.8045 34 2.292352e-23 230.8045 0.869538 0.1106203

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs  STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep  0.8110685526  0.123875507468551
    ## 1         A1 =~    PopAnh  1.0046181056  0.111460589243591
    ## 2         A1 =~ PopAppDec  0.0388258774  0.093896558358791
    ## 3         A1 =~ PopAppInc  0.4001994398 0.0969012898255562
    ## 8         A1 =~ PopSleDec  0.4971465151  0.111028583451577
    ## 9         A1 =~ PopSleInc  0.4801487333 0.0857451879407244
    ## 6         A1 =~  PopFatig  0.5908054776  0.116755495964051
    ## 7         A1 =~  PopGuilt  0.5586478780 0.0830007595164537
    ## 4         A1 =~   PopConc  0.7551502387  0.102716104214085
    ## 10        A1 =~    PopSui  0.4128035490  0.125482223454412
    ## 16    PopDep ~~    PopAnh -0.1887574005  0.195672103359941
    ## 12    PopAnh ~~    PopAnh  0.0009994275  0.210456644296615
    ## 17    PopDep ~~    PopDep  0.3421691082  0.193306115452014
    ## 13 PopAppDec ~~ PopAppDec  0.9984925476  0.170935976179391
    ## 14 PopAppInc ~~ PopAppInc  0.8398402971  0.148311059874262
    ## 20 PopSleDec ~~ PopSleDec  0.7528442480  0.214380188621573
    ## 21 PopSleInc ~~ PopSleInc  0.7694568499  0.176723027345913
    ## 18  PopFatig ~~  PopFatig  0.6509480567  0.173088415463523
    ## 19  PopGuilt ~~  PopGuilt  0.6879119109  0.152601931210766
    ## 15   PopConc ~~   PopConc  0.4297469270  0.194173772179395
    ## 22    PopSui ~~    PopSui  0.8295932363  0.172633316672397
    ## 11        A1 ~~        A1  1.0000000000

Check if model is improved by allowing residual correlations between the
directional symptoms.

``` r
pop_commonfactor_app.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    0.74 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_app.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_app.fit$modelfit
```

    ##       chisq df     p_chisq      AIC       CFI      SRMR
    ## df 181.7319 33 1.80947e-22 225.7319 0.8746557 0.1096971

``` r
pop_commonfactor_sle.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.939 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_sle.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_sle.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 174.3416 33 3.859704e-21 218.3416 0.8808839 0.1050545

``` r
pop_commonfactor_app_sle.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_commonfactor_app_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_app_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.029 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_app_sle.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_app_sle.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 168.0815 32 2.170078e-20 214.0815 0.8853168 0.1040781

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
pop_cog_mood_neuroveg.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
"
pop_cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##          lhs op       rhs Unstandardized_Estimate          SE
    ## 1         A1 =~  PopGuilt            -0.002554036 0.044489190
    ## 2         A1 =~   PopConc             0.121490232 0.060372132
    ## 3         A1 =~    PopSui             0.051530276 0.029553591
    ## 4         A2 =~    PopDep             0.220252924 0.035731148
    ## 5         A2 =~    PopAnh             0.289266668 0.036589809
    ## 6         A2 =~  PopGuilt             0.126891260 0.071048781
    ## 7         A3 =~ PopSleDec             0.115049082 0.025860794
    ## 8         A3 =~ PopSleInc             0.134802195 0.026724122
    ## 9         A3 =~  PopFatig             0.148061577 0.032913093
    ## 10        A3 =~ PopAppDec             0.009787183 0.019741767
    ## 11        A3 =~ PopAppInc             0.114204417 0.028240290
    ## 12    PopDep ~~    PopAnh            -0.023290086 0.017672160
    ## 13 PopAppDec ~~ PopAppInc            -0.004994815 0.007958542
    ## 14 PopSleDec ~~ PopSleInc            -0.014133540 0.006796112
    ## 18  PopGuilt ~~  PopGuilt             0.037495933 0.008292746
    ## 19   PopConc ~~   PopConc             0.039028434 0.016673053
    ## 20    PopSui ~~    PopSui             0.028319681 0.005905529
    ## 21    PopDep ~~    PopDep             0.012878378 0.015432946
    ## 22    PopAnh ~~    PopAnh            -0.015769195 0.020887170
    ## 23 PopSleDec ~~ PopSleDec             0.026120438 0.009239081
    ## 24 PopSleInc ~~ PopSleInc             0.038940899 0.011066008
    ## 25  PopFatig ~~  PopFatig             0.031946340 0.010240294
    ## 26 PopAppDec ~~ PopAppDec             0.035450290 0.006046625
    ## 27 PopAppInc ~~ PopAppInc             0.052754002 0.010873314
    ## 28        A1 ~~        A2             1.317222167 0.678452164
    ## 29        A1 ~~        A3             1.536295349 0.731953631
    ## 30        A2 ~~        A3             0.788587435 0.157746926

Add constraints to prevent variances from being negative and
correlations from going out of bounds.

``` r
pop_cog_mood_neuroveg_constr.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
a13 < 1
a13 > -1.0
A1 ~~ a13*A3
a23 < 1
A2 ~~ a23*A3
c7 > 0.001
PopGuilt ~~ c7*PopGuilt
"

pop_cog_mood_neuroveg_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.993 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_cog_mood_neuroveg_constr.model, : A difference greater than .025 was
    ## observed pre- and post-smoothing for Z-statistics in the genetic covariance
    ## matrix. This reflects a large difference and results should be interpreted
    ## with caution!! This can often result from including low powered traits, and you
    ## might consider removing those traits from the model. If you are going to run
    ## a multivariate GWAS we strongly recommend setting the smooth_check argument to
    ## true to check smoothing for each SNP.

``` r
pop_cog_mood_neuroveg_constr.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI       SRMR
    ## df 238.3433 28 3.090656e-35 292.3433 0.8227324 0.09914118

``` r
pop_cog_mood_neuroveg_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype   STD_Genotype_SE
    ## 2         A1 =~  PopGuilt  -0.04432851 0.757612088496248
    ## 1         A1 =~   PopConc   0.70180094  0.20492476859453
    ## 3         A1 =~    PopSui   0.38503363 0.141365566678758
    ## 8         A2 =~    PopDep   0.88268544   0.1826409888757
    ## 7         A2 =~    PopAnh   1.10328705  0.20459820275158
    ## 9         A2 =~  PopGuilt   0.58529764 0.856445552856894
    ## 15        A3 =~ PopSleDec   0.60286071 0.132638993449343
    ## 16        A3 =~ PopSleInc   0.59284724 0.113312587245224
    ## 14        A3 =~  PopFatig   0.66652800 0.142700747011123
    ## 12        A3 =~ PopAppDec   0.05398970 0.108025151903329
    ## 13        A3 =~ PopAppInc   0.46014834 0.110358526692984
    ## 23    PopDep ~~    PopAnh  -0.34779962 0.379650872693125
    ## 20 PopAppDec ~~ PopAppInc  -0.10501217 0.165144071807478
    ## 28 PopSleDec ~~ PopSleInc  -0.32839556 0.145657161410977
    ## 6         A1 ~~        A3   0.99999975  0.29707003269183
    ## 11        A2 ~~        A3   0.76111330 0.192657095678813
    ## 26  PopGuilt ~~  PopGuilt   0.70693008 0.178356333921437
    ## 22   PopConc ~~   PopConc   0.50747320 0.316535431267782
    ## 30    PopSui ~~    PopSui   0.85174816 0.190932821096428
    ## 24    PopDep ~~    PopDep   0.22086635 0.324100591122958
    ## 18    PopAnh ~~    PopAnh  -0.21724239 0.451629860270724
    ## 27 PopSleDec ~~ PopSleDec   0.63655798 0.238272345162784
    ## 29 PopSleInc ~~ PopSleInc   0.64853208 0.198715321863806
    ## 25  PopFatig ~~  PopFatig   0.55574030 0.196016163317985
    ## 19 PopAppDec ~~ PopAppDec   0.99708506 0.170050548458619
    ## 21 PopAppInc ~~ PopAppInc   0.78826366 0.166375102324416
    ## 5         A1 ~~        A2   0.99186370 0.315204210723144
    ## 4         A1 ~~        A1   1.00000000                  
    ## 10        A2 ~~        A2   1.00000000                  
    ## 17        A3 ~~        A3   1.00000000

### Two-factor models

[Elhai Psychiat Res
2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Kruse Rehab Psychol
2008](https://psycnet.apa.org/record/2008-17022-011), [Kruse Arch Psys
Med Rehab
2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> retardation as a fourth somatic item

``` r
pop_psych_soma.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    1.14 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_soma.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 195.5946 31 8.528719e-26 243.5946 0.8612873 0.1029518

``` r
pop_psych_soma.fit$results[c(1, 2, 3, 6, 7, 9)]
```

    ##          lhs op       rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 3         A1 =~    PopDep  0.814883163 0.120527093045237 1.339589e-11
    ## 1         A1 =~    PopAnh  1.004490920 0.108032629053411 1.632983e-20
    ## 4         A1 =~  PopGuilt  0.558576279 0.082633336831527 1.328651e-11
    ## 2         A1 =~   PopConc  0.754609421 0.102013611266956 1.242977e-13
    ## 5         A1 =~    PopSui  0.413989404 0.125411784372342 9.704275e-04
    ## 8         A2 =~ PopAppDec  0.051922138 0.104225143871416 6.182737e-01
    ## 9         A2 =~ PopAppInc  0.438823852 0.108874864507072 5.583432e-05
    ## 11        A2 =~ PopSleDec  0.581273877 0.131158316256823 9.228572e-06
    ## 12        A2 =~ PopSleInc  0.560054202 0.111892876780786 5.596998e-07
    ## 10        A2 =~  PopFatig  0.647380947 0.142899624166892 5.920558e-06
    ## 19    PopDep ~~    PopAnh -0.192486909 0.192191607115801 3.305793e-01
    ## 16 PopAppDec ~~ PopAppInc -0.102953739 0.164069576377373 5.303309e-01
    ## 24 PopSleDec ~~ PopSleInc -0.296536174 0.143462710982098 3.853850e-02
    ## 14    PopAnh ~~    PopAnh  0.001000421 0.204361265469114 9.420025e-01
    ## 20    PopDep ~~    PopDep  0.335964955 0.189953341475806 7.235377e-02
    ## 22  PopGuilt ~~  PopGuilt  0.687992225 0.152584178720675 6.941489e-06
    ## 18   PopConc ~~   PopConc  0.430562810 0.195195341484521 2.882530e-02
    ## 26    PopSui ~~    PopSui  0.828612685 0.172516389508268 1.584569e-06
    ## 15 PopAppDec ~~ PopAppDec  0.997304089 0.170177209489666 4.615086e-09
    ## 17 PopAppInc ~~ PopAppInc  0.807433677 0.163163322404025 7.462385e-07
    ## 23 PopSleDec ~~ PopSleDec  0.662120450 0.234274943477595 4.771944e-03
    ## 25 PopSleInc ~~ PopSleInc  0.686339239 0.193081725324797 3.769509e-04
    ## 21  PopFatig ~~  PopFatig  0.580897891 0.192527772952175 2.528030e-03
    ## 7         A1 ~~        A2  0.880203209 0.128064705299898 5.942898e-12
    ## 6         A1 ~~        A1  1.000000000                             NA
    ## 13        A2 ~~        A2  1.000000000                             NA

Bifactor model

``` r
pop_psych_soma_bif.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A  =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A  ~~ 0*A1
A  ~~ 0*A2
A1 ~~ 0*A2
"
pop_psych_soma_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.236 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_soma_bif.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI       SRMR
    ## df 157.3717 25 3.641802e-21 217.3717 0.8884433 0.07667289

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 14        A1 =~    PopDep   0.49211077  0.208349782485547
    ## 12        A1 =~    PopAnh   0.26150814  0.221700214446572
    ## 15        A1 =~  PopGuilt   0.44381967  0.156247262658778
    ## 13        A1 =~   PopConc  -0.17956274   0.25007432626088
    ## 16        A1 =~    PopSui   0.57880456  0.270044706099227
    ## 20        A2 =~ PopAppDec   0.39533203  0.225719894500856
    ## 21        A2 =~ PopAppInc   0.02374671  0.179078586035038
    ## 23        A2 =~ PopSleDec   0.73584946  0.416654330907976
    ## 24        A2 =~ PopSleInc  -0.44412030  0.270372100113865
    ## 22        A2 =~  PopFatig  -0.11999426  0.211396610394473
    ## 5          A =~    PopDep   0.63242509  0.124273508169526
    ## 1          A =~    PopAnh   0.87316185  0.111959918435586
    ## 7          A =~  PopGuilt   0.46029723  0.102009782843046
    ## 4          A =~   PopConc   0.93831706  0.125458760750997
    ## 10         A =~    PopSui   0.29160742  0.147810565898414
    ## 2          A =~ PopAppDec   0.04666501  0.109346109671257
    ## 3          A =~ PopAppInc   0.43906532  0.104845165656758
    ## 8          A =~ PopSleDec   0.58468212  0.111998773374315
    ## 9          A =~ PopSleInc   0.57527702 0.0956861211859184
    ## 6          A =~  PopFatig   0.66411825  0.132655671636109
    ## 31    PopDep ~~    PopDep   0.35786596  0.131254009837816
    ## 27    PopAnh ~~    PopAnh   0.16920218  0.113478423645472
    ## 33  PopGuilt ~~  PopGuilt   0.59115073  0.192721192151209
    ## 30   PopConc ~~   PopConc   0.08731827  0.291068307566445
    ## 36    PopSui ~~    PopSui   0.57995095  0.301353447099617
    ## 28 PopAppDec ~~ PopAppDec   0.84153456  0.218316594503329
    ## 29 PopAppInc ~~ PopAppInc   0.80665785  0.157019432952041
    ## 34 PopSleDec ~~ PopSleDec   0.11667276  0.646363683014892
    ## 35 PopSleInc ~~ PopSleInc   0.47181355  0.282623645313942
    ## 32  PopFatig ~~  PopFatig   0.54454831  0.190049247387685
    ## 18        A1 ~~        A1   1.00000000                   
    ## 26        A2 ~~        A2   1.00000000                   
    ## 11         A ~~         A   1.00000000                   
    ## 17        A1 ~~         A   0.00000000                   
    ## 25        A2 ~~         A   0.00000000                   
    ## 19        A1 ~~        A2   0.00000000

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
pop_psych_veg.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
A1 ~~ 1*A1
A2 ~~ 1*A2
"
pop_psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.249 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
pop_psych_veg.fit$modelfit
```

    ##      chisq df     p_chisq     AIC       CFI       SRMR
    ## df 244.311 31 5.47383e-35 292.311 0.8202314 0.09396672

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~    PopDep   0.94197789  0.141560029517489
    ## 1         A1 =~    PopAnh   1.17303097  0.136166379063755
    ## 3         A1 =~  PopGuilt   0.51981158 0.0814625075000047
    ## 4         A1 =~    PopSui   0.41599079  0.123572961949974
    ## 7         A2 =~ PopAppDec   0.05094247  0.107591267538842
    ## 8         A2 =~ PopAppInc   0.45523098  0.104887626484827
    ## 11        A2 =~ PopSleDec   0.58652299    0.1179633126525
    ## 12        A2 =~ PopSleInc   0.59268093  0.100957936849494
    ## 10        A2 =~  PopFatig   0.66198896  0.131395746288747
    ## 9         A2 =~   PopConc   0.86699386  0.116176096615581
    ## 19    PopDep ~~    PopAnh  -0.47891310  0.272111212571393
    ## 16 PopAppDec ~~ PopAppInc  -0.10335908  0.165177479535517
    ## 24 PopSleDec ~~ PopSleInc  -0.31861151  0.145569151459495
    ## 20    PopDep ~~    PopDep   0.11267803   0.25369763156479
    ## 14    PopAnh ~~    PopAnh  -0.37600127  0.310290185538584
    ## 22  PopGuilt ~~  PopGuilt   0.72979641  0.151055440415343
    ## 26    PopSui ~~    PopSui   0.82695183  0.170180359217284
    ## 15 PopAppDec ~~ PopAppDec   0.99740501  0.170202854664273
    ## 17 PopAppInc ~~ PopAppInc   0.79276512  0.159054321455887
    ## 23 PopSleDec ~~ PopSleDec   0.65598992  0.224334153098137
    ## 25 PopSleInc ~~ PopSleInc   0.64872894   0.18723261919431
    ## 21  PopFatig ~~  PopFatig   0.56176937  0.189389962368915
    ## 18   PopConc ~~   PopConc   0.24832234  0.227358469748024
    ## 6         A1 ~~        A2   0.72030207  0.110177951909754
    ## 5         A1 ~~        A1   1.00000000                   
    ## 13        A2 ~~        A2   1.00000000

``` r
pop_psych_veg_bif.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
"
pop_psych_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.238 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_veg_bif.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI       SRMR
    ## df 185.2048 25 2.089639e-26 245.2048 0.8649868 0.07579896

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype   STD_Genotype_SE
    ## 13        A1 =~    PopDep   0.51805448 0.181111321410942
    ## 12        A1 =~    PopAnh   0.31713591 0.195212216225347
    ## 14        A1 =~  PopGuilt   0.47537814 0.156612822064148
    ## 15        A1 =~    PopSui   0.60182061 0.282689357930054
    ## 19        A2 =~ PopAppDec   0.38244594 0.214728153293995
    ## 20        A2 =~ PopAppInc   0.01497302 0.182137765440824
    ## 23        A2 =~ PopSleDec   0.69327495 0.394453671945075
    ## 24        A2 =~ PopSleInc  -0.50693446 0.286420383879953
    ## 22        A2 =~  PopFatig  -0.13428788 0.213585530486709
    ## 21        A2 =~   PopConc  -0.15950868 0.233348493142459
    ## 5          A =~    PopDep   0.60403391 0.113576775223284
    ## 1          A =~    PopAnh   0.85542393 0.119451220642437
    ## 7          A =~  PopGuilt   0.43436967  0.10208642100046
    ## 10         A =~    PopSui   0.25117375 0.130258257958491
    ## 2          A =~ PopAppDec   0.06232718 0.113727135106465
    ## 3          A =~ PopAppInc   0.45026244 0.103600926190511
    ## 8          A =~ PopSleDec   0.61782263 0.118948913014807
    ## 9          A =~ PopSleInc   0.57967017 0.110033491691525
    ## 6          A =~  PopFatig   0.68000416 0.136097840942603
    ## 4          A =~   PopConc   0.87875817 0.121149015731455
    ## 31    PopDep ~~    PopDep   0.36676255 0.128681382418366
    ## 27    PopAnh ~~    PopAnh   0.16767504 0.116381722365414
    ## 33  PopGuilt ~~  PopGuilt   0.58533935 0.195382221817964
    ## 36    PopSui ~~    PopSui   0.57472321  0.32252994655306
    ## 28 PopAppDec ~~ PopAppDec   0.84984974 0.209343811478884
    ## 29 PopAppInc ~~ PopAppInc   0.79703973 0.159166476072502
    ## 34 PopSleDec ~~ PopSleDec   0.13766461 0.603162576391705
    ## 35 PopSleInc ~~ PopSleInc   0.40699999 0.316791321597809
    ## 32  PopFatig ~~  PopFatig   0.51956161 0.191751258995465
    ## 30   PopConc ~~   PopConc   0.20234092 0.236243631240835
    ## 17        A1 ~~        A1   1.00000000                  
    ## 26        A2 ~~        A2   1.00000000                  
    ## 11         A ~~         A   1.00000000                  
    ## 16        A1 ~~         A   0.00000000                  
    ## 25        A2 ~~         A   0.00000000                  
    ## 18        A1 ~~        A2   0.00000000

### Affective-Neurovegetative (Elhai Model 2c)

``` r
pop_affect_veg.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
A1 ~~ 1*A1
A2 ~~ 1*A2
"
pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.249 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_affect_veg.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI       SRMR
    ## df 161.4343 31 1.431445e-19 209.4343 0.8900761 0.09673964

``` r
pop_affect_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep   0.87867914   0.12021069284951
    ## 2         A1 =~  PopGuilt   0.61955658 0.0979883446902635
    ## 3         A1 =~    PopSui   0.46968805  0.140632456292199
    ## 6         A2 =~    PopAnh   0.99249623  0.109222646509391
    ## 8         A2 =~ PopAppInc   0.41904127  0.100509086143725
    ## 7         A2 =~ PopAppDec   0.04589373 0.0992847351055375
    ## 12        A2 =~ PopSleInc   0.53373809 0.0882204432828656
    ## 11        A2 =~ PopSleDec   0.53965511  0.112255648081532
    ## 10        A2 =~  PopFatig   0.62720105  0.120810360415341
    ## 9         A2 =~   PopConc   0.79857949  0.106167690448579
    ## 19    PopDep ~~    PopAnh  -0.06230039  0.185085634537978
    ## 16 PopAppInc ~~ PopAppDec  -0.09940013  0.162907914212608
    ## 24 PopSleInc ~~ PopSleDec  -0.25902616   0.14333326932045
    ## 20    PopDep ~~    PopDep   0.22792278   0.20549664417154
    ## 22  PopGuilt ~~  PopGuilt   0.61614941  0.168658322105171
    ## 26    PopSui ~~    PopSui   0.77939390  0.177347496506463
    ## 14    PopAnh ~~    PopAnh   0.01495152  0.204263476717302
    ## 17 PopAppInc ~~ PopAppInc   0.82440539  0.151365766743462
    ## 15 PopAppDec ~~ PopAppDec   0.99789441  0.170453117043963
    ## 25 PopSleInc ~~ PopSleInc   0.71512415  0.178462548536937
    ## 23 PopSleDec ~~ PopSleDec   0.70877143  0.217672336188014
    ## 21  PopFatig ~~  PopFatig   0.60661888  0.177747362540625
    ## 18   PopConc ~~   PopConc   0.36227097  0.204901825353906
    ## 5         A1 ~~        A2   0.78932205 0.0697927688160416
    ## 4         A1 ~~        A1   1.00000000                   
    ## 13        A2 ~~        A2   1.00000000

Bifactor model

``` r
pop_affect_veg_bif.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A =~ NA*PopDep + PopGuilt + PopSui + PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
"
pop_affect_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.256 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00320242979228456 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.37139331670346 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_affect_veg_bif.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI       SRMR
    ## df 61.07617 25 7.418933e-05 121.0762 0.9695967 0.06869808

``` r
pop_affect_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 12        A1 =~    PopDep   0.69574366  0.805590687898018
    ## 13        A1 =~  PopGuilt   0.20053677  0.253671901109028
    ## 14        A1 =~    PopSui   0.35692072  0.427796844378391
    ## 18        A2 =~    PopAnh   0.35475414   0.14017191087658
    ## 20        A2 =~ PopAppInc   0.08182344  0.179169920769544
    ## 19        A2 =~ PopAppDec  -0.30865414  0.207639955675287
    ## 24        A2 =~ PopSleInc   0.56418920  0.223162671388468
    ## 23        A2 =~ PopSleDec  -0.54248712   0.32023084120434
    ## 22        A2 =~  PopFatig   0.36640284  0.212965458198938
    ## 21        A2 =~   PopConc   0.38864891  0.230603303261995
    ## 5          A =~    PopDep   0.70377123 0.0882088285715635
    ## 7          A =~  PopGuilt   0.54434268  0.102323226546384
    ## 10         A =~    PopSui   0.37422215  0.125449411699762
    ## 1          A =~    PopAnh   0.90275518 0.0892916751258092
    ## 3          A =~ PopAppInc   0.41028571 0.0978852662171092
    ## 2          A =~ PopAppDec   0.13386271  0.116782439773122
    ## 9          A =~ PopSleInc   0.41203749  0.123734070312046
    ## 8          A =~ PopSleDec   0.73285768  0.154803692195913
    ## 6          A =~  PopFatig   0.56827506  0.145636684280155
    ## 4          A =~   PopConc   0.73472929  0.126063761490523
    ## 31    PopDep ~~    PopDep   0.02064676   1.12733095703791
    ## 33  PopGuilt ~~  PopGuilt   0.66347615  0.180645430518929
    ## 36    PopSui ~~    PopSui   0.73256550  0.314357657277161
    ## 27    PopAnh ~~    PopAnh   0.05918273  0.132655265003761
    ## 29 PopAppInc ~~ PopAppInc   0.82497044  0.150209298869333
    ## 28 PopAppDec ~~ PopAppDec   0.88681309  0.197888898081329
    ## 35 PopSleInc ~~ PopSleInc   0.51191552  0.251185975932993
    ## 34 PopSleDec ~~ PopSleDec   0.16862684  0.513988382933592
    ## 32  PopFatig ~~  PopFatig   0.54281247  0.190297637612518
    ## 30   PopConc ~~   PopConc   0.30912518  0.212222833590392
    ## 16        A1 ~~        A1   1.00000000                   
    ## 26        A2 ~~        A2   1.00000000                   
    ## 11         A ~~         A   1.00000000                   
    ## 15        A1 ~~         A   0.00000000                   
    ## 25        A2 ~~         A   0.00000000                   
    ## 17        A1 ~~        A2   0.00000000

### Model comparisons

``` r
model_fits <- 
data.frame(Model=c('1a', '1b', '1c', '1d', '1e',
                   '2a(i)', '2a(ii)', '2b(i)', '2b(ii)',
                   '2c(i)', '2c(ii)',
                   '3(i)'),
       Name=c('Common',
              'Common (gating)',
              'Common (App)',
              'Common (Sle)',
              'Common(App,Sle)',
              'Psych-Somatic',
              'Psych-Somatic (BiF)',
              'Psych-Neuroveg',
              'Psych-Neuroveg (BiF)',
              'Affect-Neuroveg',
              'Affect-Neuroveg (BiF)',
              'Cog-Mood-Neuroveg'
              )) %>%
bind_cols(
bind_rows(
lapply(list(pop_commonfactor.fit,
            pop_commonfactor_gating.fit,
            pop_commonfactor_app.fit,
            pop_commonfactor_sle.fit,
            pop_commonfactor_app_sle.fit,
            pop_psych_soma.fit,
            pop_psych_soma_bif.fit,
            pop_psych_veg.fit,
            pop_psych_veg_bif.fit,
            pop_affect_veg.fit,
            pop_affect_veg_bif.fit,
            pop_cog_mood_neuroveg_constr.fit),
       function(fit) fit$modelfit)
))
rownames(model_fits) <- NULL

model_fits %>%
select(-chisq, -df) %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~signif(., 4))
```

    ##     Model                  Name   p_chisq   AIC    CFI    SRMR   dAIC
    ## 1      1a                Common 2.129e-33 285.3 0.8228 0.11110 164.20
    ## 2      1b       Common (gating) 2.292e-23 230.8 0.8695 0.11060 109.70
    ## 3      1c          Common (App) 1.809e-22 225.7 0.8747 0.10970 104.70
    ## 4      1d          Common (Sle) 3.860e-21 218.3 0.8809 0.10510  97.27
    ## 5      1e       Common(App,Sle) 2.170e-20 214.1 0.8853 0.10410  93.01
    ## 6   2a(i)         Psych-Somatic 8.529e-26 243.6 0.8613 0.10300 122.50
    ## 7  2a(ii)   Psych-Somatic (BiF) 3.642e-21 217.4 0.8884 0.07667  96.30
    ## 8   2b(i)        Psych-Neuroveg 5.474e-35 292.3 0.8202 0.09397 171.20
    ## 9  2b(ii)  Psych-Neuroveg (BiF) 2.090e-26 245.2 0.8650 0.07580 124.10
    ## 10  2c(i)       Affect-Neuroveg 1.431e-19 209.4 0.8901 0.09674  88.36
    ## 11 2c(ii) Affect-Neuroveg (BiF) 7.419e-05 121.1 0.9696 0.06870   0.00
    ## 12   3(i)     Cog-Mood-Neuroveg 3.091e-35 292.3 0.8227 0.09914 171.30

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
symptoms_cov_pd <- as.matrix(Matrix::nearPD(symptoms_cov_pos, corr=FALSE)$mat)

corrplot(cov2cor(symptoms_cov_pd))
```

![](mdd-symptom-gsem-model-hdl_files/figure-gfm/mdd_symptom_gsem_efa_pd_hdl-1.png)<!-- -->

## PGC/AGDS

Check eigen values of the correlation matrix

``` r
symptoms_clin_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Clin'))

symptoms_clin_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx])) 

plot(symptoms_clin_eigen$values, ylab='Eigenvalue')
lines(symptoms_clin_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model-hdl_files/figure-gfm/mdd_symptom_gsem_clin_efa_eigen_hdl-1.png)<!-- -->

``` r
symptoms_clin_efa <- factanal(covmat=symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx], factors=3, rotation='varimax')
symptoms_clin_efa
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd[symptoms_clin_idx,     symptoms_clin_idx], rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     ClinDep     ClinAnh  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc 
    ##       0.005       0.005       0.005       0.189       0.642       0.570 
    ## ClinPsycInc ClinPsycDec   ClinFatig   ClinGuilt    ClinConc     ClinSui 
    ##       0.806       0.303       0.645       0.871       0.023       0.640 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinDep              0.982   0.150 
    ## ClinAnh      0.168   0.976   0.119 
    ## ClinAppDec           0.183   0.976 
    ## ClinAppInc   0.641   0.276  -0.569 
    ## ClinSleDec          -0.427   0.418 
    ## ClinSleInc   0.640   0.141         
    ## ClinPsycInc  0.157           0.400 
    ## ClinPsycDec  0.815           0.183 
    ## ClinFatig    0.534          -0.262 
    ## ClinGuilt    0.206          -0.293 
    ## ClinConc     0.974  -0.140         
    ## ClinSui      0.598                 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      3.189   2.261   1.849
    ## Proportion Var   0.266   0.188   0.154
    ## Cumulative Var   0.266   0.454   0.608
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 6.5329

## ALSPAC/UKB

Check eigen values of the correlation matrix

``` r
symptoms_pop_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Pop'))

symptoms_pop_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_pop_idx,symptoms_pop_idx])) 

plot(symptoms_pop_eigen$values, ylab='Eigenvalue')
lines(symptoms_pop_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model-hdl_files/figure-gfm/mdd_symptom_gsem_pop_efa_eigen_hdl-1.png)<!-- -->

``` r
symptoms_pop_efa <- factanal(covmat=symptoms_cov_pd[symptoms_pop_idx,symptoms_pop_idx], factors=3, rotation='varimax')
symptoms_pop_efa
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd[symptoms_pop_idx,     symptoms_pop_idx], rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     PopDep     PopAnh  PopAppDec  PopAppInc  PopSleDec  PopSleInc PopPsycInc 
    ##      0.469      0.092      0.963      0.786      0.696      0.667      0.005 
    ## PopPsycDec   PopFatig   PopGuilt    PopConc     PopSui 
    ##      0.005      0.503      0.667      0.347      0.733 
    ## 
    ## Loadings:
    ##            Factor1 Factor2 Factor3
    ## PopDep      0.511   0.498   0.147 
    ## PopAnh      0.806   0.437   0.258 
    ## PopAppDec                   0.187 
    ## PopAppInc   0.385   0.124  -0.226 
    ## PopSleDec   0.513  -0.139  -0.147 
    ## PopSleInc   0.569                 
    ## PopPsycInc         -0.153   0.983 
    ## PopPsycDec -0.240   0.952  -0.176 
    ## PopFatig    0.688          -0.152 
    ## PopGuilt    0.430   0.321   0.211 
    ## PopConc     0.781           0.188 
    ## PopSui      0.258   0.348   0.281 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      3.044   1.638   1.384
    ## Proportion Var   0.254   0.136   0.115
    ## Cumulative Var   0.254   0.390   0.506
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 1.3684

## All symptoms

Check eigen values of the correlation matrix

``` r
symptoms_cov_pd.eigen <- eigen(cov2cor(symptoms_cov_pd))

signif(symptoms_cov_pd.eigen$values, 3)
```

    ##  [1] 5.53e+00 3.03e+00 2.61e+00 2.23e+00 2.07e+00 1.66e+00 1.50e+00 1.13e+00
    ##  [9] 1.03e+00 7.09e-01 6.19e-01 5.75e-01 4.34e-01 3.56e-01 2.99e-01 2.07e-01
    ## [17] 2.69e-07 2.33e-07 1.94e-07 1.72e-07 1.67e-07 1.12e-07 4.69e-08 3.06e-08

``` r
plot(eigen(cov2cor(symptoms_cov_pd))$values, ylab='Eigenvalue')
lines(eigen(cov2cor(symptoms_cov_pd))$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model-hdl_files/figure-gfm/mdd_symptom_gsem_efa_eigen_hdl-1.png)<!-- -->

Simulate uncertainty in **S** using the **V** matrix

``` r
# replicates
m <- 100

# simulate with S as the mean and V as the variance
S_lowertri_sim = mvtnorm::rmvnorm(m,
    mean=symptoms_covstruct$S[lower.tri(symptoms_covstruct$S, diag=TRUE)],
    sigma=symptoms_covstruct$V)
    
# reshape into a 24x24xm array
S_sim <- plyr::aaply(S_lowertri_sim, 1, function(x, k=nrow(symptoms_covstruct$S)){
    S <- matrix(NA, ncol=k, nrow=k)
    S[lower.tri(S, diag=T)] <- x
    S[upper.tri(S, diag=T)] <- t(S)[upper.tri(S, diag=T)]
    return(S)
})
dimnames(S_sim) <- list(1:m, colnames(symptoms_cov), colnames(symptoms_cov))

# find which symptoms have positive variances across all replicates
sim_cov_keep <- which(colSums(plyr::aaply(S_sim, 1, diag) > 0) == m)

eigen(cov2cor(symptoms_covstruct$S[sim_cov_keep,sim_cov_keep]))$values
```

    ##  [1]  7.5052818  5.1712519  3.7425629  2.6386342  2.3730357  1.6913698
    ##  [7]  1.4492533  1.1899598  0.9419901  0.7532313  0.6077189  0.5502851
    ## [13]  0.4101393  0.3049773  0.2529266 -0.0350118 -0.1127862 -0.6430436
    ## [19] -1.2420131 -1.5318487 -5.0179145

``` r
S_sim_pos <- S_sim[,sim_cov_keep,sim_cov_keep]

S_sim_pos_ev <- plyr::aaply(S_sim_pos, 1, function(x) eigen(cov2cor(x))$values)

summary(S_sim_pos_ev)
```

    ##        1                2               3               4        
    ##  Min.   : 6.508   Min.   :4.750   Min.   :3.407   Min.   :2.695  
    ##  1st Qu.: 7.898   1st Qu.:5.618   1st Qu.:4.116   1st Qu.:3.166  
    ##  Median : 8.675   Median :6.134   Median :4.643   Median :3.425  
    ##  Mean   : 8.912   Mean   :6.366   Mean   :4.666   Mean   :3.459  
    ##  3rd Qu.: 9.478   3rd Qu.:6.972   3rd Qu.:5.150   3rd Qu.:3.738  
    ##  Max.   :15.738   Max.   :9.755   Max.   :6.630   Max.   :4.521  
    ##        5               6               7               8        
    ##  Min.   :2.166   Min.   :1.716   Min.   :1.199   Min.   :1.029  
    ##  1st Qu.:2.606   1st Qu.:1.997   1st Qu.:1.581   1st Qu.:1.282  
    ##  Median :2.749   Median :2.211   Median :1.740   Median :1.409  
    ##  Mean   :2.836   Mean   :2.222   Mean   :1.745   Mean   :1.404  
    ##  3rd Qu.:3.079   3rd Qu.:2.404   3rd Qu.:1.900   3rd Qu.:1.509  
    ##  Max.   :3.619   Max.   :3.101   Max.   :2.375   Max.   :1.954  
    ##        9                10               11               12        
    ##  Min.   :0.7591   Min.   :0.5356   Min.   :0.3569   Min.   :0.2367  
    ##  1st Qu.:0.9944   1st Qu.:0.7945   1st Qu.:0.5885   1st Qu.:0.3891  
    ##  Median :1.1075   Median :0.8711   Median :0.6644   Median :0.4542  
    ##  Mean   :1.1117   Mean   :0.8693   Mean   :0.6648   Mean   :0.4626  
    ##  3rd Qu.:1.2174   3rd Qu.:0.9628   3rd Qu.:0.7345   3rd Qu.:0.5300  
    ##  Max.   :1.5115   Max.   :1.1621   Max.   :0.9492   Max.   :0.7806  
    ##        13                14                  15                 16          
    ##  Min.   :0.01697   Min.   :-0.253745   Min.   :-0.66408   Min.   :-1.00461  
    ##  1st Qu.:0.20104   1st Qu.: 0.005646   1st Qu.:-0.19554   1st Qu.:-0.55929  
    ##  Median :0.26536   Median : 0.088012   Median :-0.09583   Median :-0.39819  
    ##  Mean   :0.26486   Mean   : 0.087327   Mean   :-0.12747   Mean   :-0.43317  
    ##  3rd Qu.:0.33224   3rd Qu.: 0.159616   3rd Qu.:-0.04322   3rd Qu.:-0.28762  
    ##  Max.   :0.48736   Max.   : 0.322849   Max.   : 0.11118   Max.   :-0.03694  
    ##        17                18                19                20        
    ##  Min.   :-1.6929   Min.   :-2.3755   Min.   :-5.4061   Min.   :-7.413  
    ##  1st Qu.:-0.9855   1st Qu.:-1.4765   1st Qu.:-2.2810   1st Qu.:-3.772  
    ##  Median :-0.7302   Median :-1.2537   Median :-1.9170   Median :-2.898  
    ##  Mean   :-0.7776   Mean   :-1.2892   Mean   :-2.0340   Mean   :-3.176  
    ##  3rd Qu.:-0.5914   3rd Qu.:-1.0355   3rd Qu.:-1.6065   3rd Qu.:-2.338  
    ##  Max.   :-0.3061   Max.   :-0.6951   Max.   :-0.9767   Max.   :-1.509  
    ##        21         
    ##  Min.   :-11.898  
    ##  1st Qu.: -7.097  
    ##  Median : -6.062  
    ##  Mean   : -6.235  
    ##  3rd Qu.: -5.147  
    ##  Max.   : -4.208

``` r
symptoms_efa3 <- factanal(covmat=symptoms_cov_pd, factors=3, rotation='varimax')

symptoms_efa3
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     ClinDep     ClinAnh  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc 
    ##       0.005       0.005       0.868       0.494       0.857       0.566 
    ## ClinPsycInc ClinPsycDec   ClinFatig   ClinGuilt    ClinConc     ClinSui 
    ##       0.896       0.376       0.669       0.716       0.017       0.316 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.453       0.231       0.951       0.901       0.855       0.753 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.889       0.957       0.695       0.564       0.491       0.712 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinDep      0.213           0.973 
    ## ClinAnh      0.227   0.209   0.949 
    ## ClinAppDec          -0.161   0.321 
    ## ClinAppInc   0.188   0.673   0.129 
    ## ClinSleDec                  -0.373 
    ## ClinSleInc           0.647         
    ## ClinPsycInc  0.294                 
    ## ClinPsycDec  0.197   0.764         
    ## ClinFatig            0.564         
    ## ClinGuilt    0.497   0.142  -0.135 
    ## ClinConc     0.101   0.962  -0.218 
    ## ClinSui      0.644   0.500  -0.139 
    ## PopDep       0.731   0.116         
    ## PopAnh       0.869   0.117         
    ## PopAppDec            0.209         
    ## PopAppInc    0.312                 
    ## PopSleDec    0.377                 
    ## PopSleInc    0.467   0.120   0.119 
    ## PopPsycInc   0.140   0.101  -0.284 
    ## PopPsycDec   0.123          -0.153 
    ## PopFatig     0.539  -0.107         
    ## PopGuilt     0.616   0.173  -0.165 
    ## PopConc      0.709                 
    ## PopSui       0.524           0.118 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      4.173   3.195   2.389
    ## Proportion Var   0.174   0.133   0.100
    ## Cumulative Var   0.174   0.307   0.407
    ## 
    ## The degrees of freedom for the model is 207 and the fit was 114.9038

``` r
symptoms_efa4 <- factanal(covmat=symptoms_cov_pd, factors=4, rotation='varimax')

symptoms_efa4
```

    ## 
    ## Call:
    ## factanal(factors = 4, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     ClinDep     ClinAnh  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc 
    ##       0.005       0.005       0.005       0.185       0.631       0.585 
    ## ClinPsycInc ClinPsycDec   ClinFatig   ClinGuilt    ClinConc     ClinSui 
    ##       0.735       0.319       0.638       0.604       0.005       0.295 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.459       0.260       0.843       0.645       0.846       0.731 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.881       0.933       0.682       0.560       0.498       0.710 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4
    ## ClinDep      0.191           0.971         
    ## ClinAnh      0.206   0.186   0.951   0.115 
    ## ClinAppDec          -0.254   0.238   0.934 
    ## ClinAppInc   0.248   0.720   0.184  -0.449 
    ## ClinSleDec                  -0.417   0.434 
    ## ClinSleInc           0.628   0.104         
    ## ClinPsycInc  0.272                   0.431 
    ## ClinPsycDec  0.201   0.728           0.326 
    ## ClinFatig            0.569          -0.173 
    ## ClinGuilt    0.554   0.156  -0.108  -0.231 
    ## ClinConc     0.108   0.966  -0.206         
    ## ClinSui      0.649   0.482  -0.134   0.183 
    ## PopDep       0.720   0.115                 
    ## PopAnh       0.849   0.111                 
    ## PopAppDec            0.172           0.349 
    ## PopAppInc    0.368                  -0.465 
    ## PopSleDec    0.381                         
    ## PopSleInc    0.472   0.111   0.144  -0.114 
    ## PopPsycInc   0.146   0.104  -0.287         
    ## PopPsycDec   0.101          -0.163   0.161 
    ## PopFatig     0.543  -0.123                 
    ## PopGuilt     0.623   0.171  -0.151         
    ## PopConc      0.697                   0.116 
    ## PopSui       0.523           0.126         
    ## 
    ##                Factor1 Factor2 Factor3 Factor4
    ## SS loadings      4.225   3.198   2.387   2.130
    ## Proportion Var   0.176   0.133   0.099   0.089
    ## Cumulative Var   0.176   0.309   0.409   0.498
    ## 
    ## The degrees of freedom for the model is 186 and the fit was 112.5211

``` r
symptoms_efa6 <- factanal(covmat=symptoms_cov_pd, factors=6, rotation='varimax')

symptoms_efa6
```

    ## 
    ## Call:
    ## factanal(factors = 6, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     ClinDep     ClinAnh  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc 
    ##       0.005       0.005       0.264       0.005       0.671       0.335 
    ## ClinPsycInc ClinPsycDec   ClinFatig   ClinGuilt    ClinConc     ClinSui 
    ##       0.738       0.369       0.005       0.442       0.005       0.321 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.393       0.175       0.758       0.337       0.785       0.483 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.317       0.285       0.615       0.379       0.628 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## ClinDep     -0.137           0.165   0.973                 
    ## ClinAnh      0.134           0.148   0.971                 
    ## ClinAppDec  -0.194  -0.101   0.328   0.262  -0.714         
    ## ClinAppInc   0.649                   0.236   0.704   0.148 
    ## ClinSleDec          -0.243   0.313  -0.405                 
    ## ClinSleInc   0.694   0.188  -0.155   0.189  -0.167   0.247 
    ## ClinPsycInc                  0.490          -0.103         
    ## ClinPsycDec  0.767   0.123   0.136                         
    ## ClinFatig    0.582   0.283  -0.332           0.109  -0.671 
    ## ClinGuilt    0.158   0.266   0.234           0.518   0.370 
    ## ClinConc     0.970           0.164  -0.140                 
    ## ClinSui      0.534   0.385   0.464                   0.168 
    ## PopDep       0.103   0.357   0.670           0.126         
    ## PopAnh       0.139   0.626   0.626                   0.115 
    ## PopAppDec    0.247                          -0.401         
    ## PopAppInc            0.253   0.133           0.735  -0.196 
    ## PopSleDec            0.406                   0.128  -0.134 
    ## PopSleInc    0.180   0.658           0.200           0.104 
    ## PopPsycInc   0.177                  -0.183           0.956 
    ## PopPsycDec          -0.230   0.624  -0.233          -0.421 
    ## PopFatig             0.829                          -0.149 
    ## PopGuilt     0.204   0.385   0.393           0.159         
    ## PopConc      0.105   0.709   0.299                   0.126 
    ## PopSui               0.161   0.524   0.155           0.217 
    ## 
    ##                Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## SS loadings      3.377   3.041   2.725   2.424   2.124   1.992
    ## Proportion Var   0.141   0.127   0.114   0.101   0.088   0.083
    ## Cumulative Var   0.141   0.267   0.381   0.482   0.570   0.653
    ## 
    ## The degrees of freedom for the model is 147 and the fit was 106.2657
