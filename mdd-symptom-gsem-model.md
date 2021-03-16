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
    ## platform       x86_64-conda-linux-gnu      
    ## arch           x86_64                      
    ## os             linux-gnu                   
    ## system         x86_64, linux-gnu           
    ## status                                     
    ## major          4                           
    ## minor          0.3                         
    ## year           2020                        
    ## month          10                          
    ## day            10                          
    ## svn rev        79318                       
    ## language       R                           
    ## version.string R version 4.0.3 (2020-10-10)
    ## nickname       Bunny-Wunnies Freak Out

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

if(!require(GenomicSEM)) remotes::install_github("MichelNivard/GenomicSEM")

if(!require(tidySEM)) remotes::install_github("cjvanlissa/tidySEM")
```

GenomicSEM version

``` r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(corrplot)
require(tidySEM)
require(GenomicSEM)

packageVersion("GenomicSEM")
```

    ## [1] '0.0.2'

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
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────
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
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
A1 ~~ 1*A1
c3a3b > -1
ClinAppDec ~~ c3a3b*ClinAppInc
"
pgc_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pgc_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.111 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0367078980491331 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  Inf . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

``` r
pgc_commonfactor.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 2.604373  8 0.9566861 28.60437   1 0.1435089

``` r
pgc_commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE
    ## 1           A1 =~  ClinAppDec  0.671742966 0.433122026645965
    ## 2           A1 =~  ClinAppInc  0.729689563 0.521524441752345
    ## 3           A1 =~  ClinSleDec  0.746888149 0.474570264504223
    ## 4           A1 =~  ClinSleInc  0.612271642 0.461874360406521
    ## 5           A1 =~ ClinPsycInc  0.279997292 0.237173888205013
    ## 6           A1 =~     ClinSui  0.007296054 0.178669809982275
    ## 7           A1 ~~          A1  1.000000000                  
    ## 8   ClinAppDec ~~  ClinAppDec  0.548761368 0.574094487130573
    ## 9   ClinAppDec ~~  ClinAppInc -0.861169631 0.618301627013953
    ## 10  ClinAppInc ~~  ClinAppInc  0.467553163  0.79129051639139
    ## 11  ClinSleDec ~~  ClinSleDec  0.442157966 0.900065138452166
    ## 12  ClinSleInc ~~  ClinSleInc  0.625124988  1.61291792720481
    ## 13 ClinPsycInc ~~ ClinPsycInc  0.921601537 0.470933828537936
    ## 14     ClinSui ~~     ClinSui  0.999946802 0.363500519041763

``` r
fit_graph <- function(results, ...) {

  results_sort <- results %>% arrange(lhs, rhs)

  node_names <- unique(c(results_sort$lhs, results_sort$rhs))
  
  node_idx <- seq_along(node_names)
  names(node_idx) <- node_names
  
  graph <- create_graph(
    nodes_df=create_node_df(n=length(node_names),
                            label=node_labels[node_names], 
                            shape='oval', width=1,
                            fillcolor=node_colors[node_names],
                            fontcolor='black'),
    edges_df=create_edge_df(from=node_idx[results_sort$lhs],
                            to=node_idx[results_sort$rhs],
                            label=round(results_sort$STD_Genotype, 2),
                            penwidth=0.3+abs(2*results_sort$STD_Genotype),
                            dir=edge_dir[results_sort$op]),
    attr_theme="tb")
  
  return(graph)

}

# render_fit <- function(results) render_graph(fit_graph(results))

# render_fit(pgc_commonfactor.fit$results)

# pgc_commonfactor.graph <- fit_graph(pgc_commonfactor.fit$results)

add_rank_same <- function(gv, top, bottom) {

  # add in block to specify node ranks
  gv.list <- str_split(gv, '\n\n')[[1]]
  
  # move the nodes/edges element to the end
  gv.list[6] <- gv.list[5]
  
  # add manually made ranks
  gv.list[5] <- paste("{rank=same",
  paste(paste0("'", top, "'"), collapse=' '),
  "}\n{rank=same",
  paste(paste0("'", bottom, "'"), collapse=' '),
  "}")
  
  rank.gv <- paste(gv.list, collapse='\n\n')

}


# pgc_commonfactor.gv <- add_rank_same(generate_dot(pgc_commonfactor.graph), 1, 2:5)
# grViz(pgc_commonfactor.gv)

# # output as a GraphViz dot file. Replace single quotes with double quotes 
# # as that's what the command line utility expects
# cat(str_replace_all(pgc_commonfactor_rank.gv, "'", '"'), file='mdd-symptom-gsem_files/pgc_commonfactor.gv')
```

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
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.801 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##      chisq df     p_chisq     AIC       CFI      SRMR
    ## df 59.8087 35 0.005595033 99.8087 0.9928165 0.1257256

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep  0.863027381 0.0466869327427978
    ## 2         A1 =~    PopSui  0.587162338 0.0983552009075276
    ## 3         A1 =~    PopAnh  0.997079793 0.0429873353388014
    ## 4         A1 =~ PopAppDec  0.170614118  0.088175575191152
    ## 5         A1 =~ PopAppInc  0.406193508 0.0809098960160941
    ## 6         A1 =~ PopSleDec  0.576324497  0.094833161999191
    ## 7         A1 =~ PopSleInc  0.520220088  0.101746585361935
    ## 8         A1 =~  PopFatig  0.671995433 0.0937553413527298
    ## 9         A1 =~  PopGuilt  0.627685840  0.080178149750251
    ## 10        A1 =~   PopConc  0.701639753 0.0903611152361196
    ## 11        A1 ~~        A1  1.000000000                   
    ## 12    PopDep ~~    PopDep  0.255183639 0.0738078886665451
    ## 13    PopSui ~~    PopSui  0.655234771  0.253105388023696
    ## 14    PopAnh ~~    PopAnh  0.005831204 0.0676895104525014
    ## 15 PopAppDec ~~ PopAppDec  0.970892075  0.251433695546643
    ## 16 PopAppInc ~~ PopAppInc  0.835007396  0.166782777541358
    ## 17 PopSleDec ~~ PopSleDec  0.667846685  0.331077724393484
    ## 18 PopSleInc ~~ PopSleInc  0.729373604  0.284423004694919
    ## 19  PopFatig ~~  PopFatig  0.548420644  0.322517150516019
    ## 20  PopGuilt ~~  PopGuilt  0.606011274  0.166002077115738
    ## 21   PopConc ~~   PopConc  0.507701935  0.268448841475589

Remove common variance shared between the gating items (Mood:
`UKB_CIDI1`, Interest: `UKB_CIDI2`) that is uncorrelated with the common
factor variance, to recover the genetic structure among gated items

``` r
pop_commonfactor_gating.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
"
pop_commonfactor_gating.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_gating.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.862 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI     SRMR
    ## df 40.57809 34 0.2029523 82.57809 0.9980953 0.117905

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7357950 0.0783402516990925
    ## 2         A1 =~    PopSui    0.6213223  0.105725480677441
    ## 3         A1 =~    PopAnh    0.8736869 0.0734377965485583
    ## 4         A1 =~ PopAppDec    0.1846453 0.0936391380069975
    ## 5         A1 =~ PopAppInc    0.4363183 0.0877965382345287
    ## 6         A1 =~ PopSleDec    0.6278670  0.103964416706934
    ## 7         A1 =~ PopSleInc    0.5482410  0.109301879710098
    ## 8         A1 =~  PopFatig    0.7276567  0.101483938810399
    ## 9         A1 =~  PopGuilt    0.6808364 0.0887935292254339
    ## 10        A1 =~   PopConc    0.7679949 0.0985288654752275
    ## 11        A1 ~~        A1    1.0000000                   
    ## 12    PopDep ~~    PopDep    0.4586053   0.13038309661583
    ## 13    PopDep ~~    PopAnh    0.2809443  0.119166070898074
    ## 14    PopSui ~~    PopSui    0.6139594  0.258378889036349
    ## 15    PopAnh ~~    PopAnh    0.2366716  0.130111342800632
    ## 16 PopAppDec ~~ PopAppDec    0.9659057  0.251348401913751
    ## 17 PopAppInc ~~ PopAppInc    0.8096253  0.171640370800988
    ## 18 PopSleDec ~~ PopSleDec    0.6057813  0.333924066420959
    ## 19 PopSleInc ~~ PopSleInc    0.6994331  0.285330009876885
    ## 20  PopFatig ~~  PopFatig    0.4705154  0.324391546022596
    ## 21  PopGuilt ~~  PopGuilt    0.5364621  0.171168778682847
    ## 22   PopConc ~~   PopConc    0.4101842  0.270813269624111

Check if model is improved by allowing residual correlations between the
directional symptoms.

``` r
pop_commonfactor_app.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
"
pop_commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.045 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##      chisq df   p_chisq     AIC       CFI      SRMR
    ## df 39.8628 33 0.1912615 83.8628 0.9980128 0.1155375

``` r
pop_commonfactor_sle.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
"
pop_commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.021 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df p_chisq      AIC       CFI      SRMR
    ## df 38.33266 33 0.24035 82.33266 0.9984559 0.1063218

### Kendler Neale model

``` r
pop_kendler_neale.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
"
pop_kendler_neale.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_kendler_neale.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_kendler_neale.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##           lhs op       rhs Unstandardized_Estimate
    ## 1          A1 =~  PopGuilt           -1.587317e-01
    ## 2          A1 =~   PopConc           -1.687143e-01
    ## 3          A1 =~    PopSui           -1.074497e-01
    ## 4          A2 =~    PopDep            3.257335e+00
    ## 5          A2 =~    PopAnh            4.014837e+00
    ## 6          A2 =~  PopGuilt           -2.624842e-04
    ## 7          A3 =~ PopSleDec            1.320644e-01
    ## 8          A3 =~ PopSleInc            1.194367e-01
    ## 9          A3 =~  PopFatig            1.680700e-01
    ## 10         A3 =~ PopAppDec            3.698864e-02
    ## 11         A3 =~ PopAppInc            1.212150e-01
    ## 15     PopDep ~~    PopAnh           -1.300040e+01
    ## 16  PopSleDec ~~ PopSleInc           -1.635657e-02
    ## 17  PopAppDec ~~ PopAppInc           -8.481291e-03
    ## 133    PopDep ~~    PopDep           -1.052974e+01
    ## 134    PopAnh ~~    PopAnh           -1.603200e+01
    ## 135 PopAppDec ~~ PopAppDec            3.381163e-02
    ## 136 PopAppInc ~~ PopAppInc            6.345856e-02
    ## 137 PopSleDec ~~ PopSleDec            2.456743e-02
    ## 138 PopSleInc ~~ PopSleInc            3.060018e-02
    ## 139  PopFatig ~~  PopFatig            2.728230e-02
    ## 140  PopGuilt ~~  PopGuilt            3.214604e-02
    ## 141   PopConc ~~   PopConc            2.304104e-02
    ## 142    PopSui ~~    PopSui            2.077681e-02
    ## 195        A1 ~~        A2           -6.630526e-02
    ## 196        A1 ~~        A3           -1.096331e+00
    ## 197        A2 ~~        A3            6.439752e-02

Add constraints to prevent variances from being negative and
correlations from going out of bounds.

``` r
pop_kendler_neale_constr.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
a13 > -1.0
A1 ~~ a13*A3
"
pop_kendler_neale_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_kendler_neale_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   7.678 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_kendler_neale_constr.model, : A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_kendler_neale_constr.fit$modelfit
```

    ##     chisq df    p_chisq    AIC      CFI      SRMR
    ## df 38.485 28 0.08956883 92.485 0.996964 0.1034545

``` r
pop_kendler_neale_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopSui   0.40812631  0.208103708831721
    ## 2         A1 =~  PopGuilt  -0.07345475  0.211342376306899
    ## 3         A1 =~   PopConc   0.50277825  0.251333320462741
    ## 4         A1 ~~        A1   1.00000000                   
    ## 5         A1 ~~        A2   1.53797968  0.735737883976811
    ## 6         A1 ~~        A3   1.71222908  0.784418260902913
    ## 7         A2 =~    PopDep   0.72962243 0.0844352429616578
    ## 8         A2 =~    PopAnh   0.86516791 0.0858289439411759
    ## 9         A2 =~  PopGuilt   0.78746635  0.296277503826445
    ## 10        A2 ~~        A2   1.00000000                   
    ## 11        A2 ~~        A3   1.01406796  0.160411495645851
    ## 12        A3 =~ PopAppDec   0.19657622 0.0946055813127982
    ## 13        A3 =~ PopAppInc   0.43542653 0.0956821673874205
    ## 14        A3 =~ PopSleDec   0.64316133  0.121428998511701
    ## 15        A3 =~ PopSleInc   0.56507097  0.117485895533465
    ## 16        A3 =~  PopFatig   0.71093275  0.123634929781074
    ## 17        A3 ~~        A3   1.00000000                   
    ## 18    PopDep ~~    PopDep   0.46765113  0.141645139230262
    ## 19    PopDep ~~    PopAnh   0.29255256  0.136434315852559
    ## 20    PopSui ~~    PopSui   0.83343046  0.294706962559891
    ## 21    PopAnh ~~    PopAnh   0.25148441  0.153844883182964
    ## 22 PopAppDec ~~ PopAppDec   0.96135682  0.249622147199745
    ## 23 PopAppDec ~~ PopAppInc  -0.16183684  0.140945795040004
    ## 24 PopAppInc ~~ PopAppInc   0.81040419  0.170094672488897
    ## 25 PopSleDec ~~ PopSleDec   0.58634096  0.329259854076718
    ## 26 PopSleDec ~~ PopSleInc  -0.37686382  0.217024607044167
    ## 27 PopSleInc ~~ PopSleInc   0.68069864  0.287671473332995
    ## 28  PopFatig ~~  PopFatig   0.49457474   0.32138896180926
    ## 29  PopGuilt ~~  PopGuilt   0.55242426  0.182793159275404
    ## 30   PopConc ~~   PopConc   0.74721633  0.328998385787683

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
"
pop_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##           lhs op       rhs Unstandardized_Estimate          SE
    ## 1          A1 =~    PopDep              0.20414183 0.021805609
    ## 2          A1 =~    PopAnh              0.25222890 0.021822434
    ## 3          A1 =~  PopGuilt              0.16343160 0.021214806
    ## 4          A1 =~   PopConc              0.17576042 0.022291012
    ## 5          A1 =~    PopSui              0.11160311 0.018944332
    ## 6          A2 =~ PopAppDec              0.03170017 0.016348249
    ## 7          A2 =~ PopAppInc              0.11236656 0.027040822
    ## 8          A2 =~ PopSleDec              0.11837431 0.025122227
    ## 9          A2 =~ PopSleInc              0.10665399 0.024434674
    ## 10         A2 =~  PopFatig              0.15586198 0.030523387
    ## 13     PopDep ~~    PopAnh              0.02578165 0.009731849
    ## 109    PopDep ~~    PopDep              0.03882323 0.010197936
    ## 110    PopAnh ~~    PopAnh              0.02329877 0.011290945
    ## 111 PopAppDec ~~ PopAppDec              0.03417495 0.008844264
    ## 112 PopAppInc ~~ PopAppInc              0.06552526 0.013258847
    ## 113 PopSleDec ~~ PopSleDec              0.02799581 0.013850274
    ## 114 PopSleInc ~~ PopSleInc              0.03349036 0.013092085
    ## 115  PopFatig ~~  PopFatig              0.03123707 0.017848820
    ## 116  PopGuilt ~~  PopGuilt              0.03062658 0.009824874
    ## 117   PopConc ~~   PopConc              0.02061383 0.013897763
    ## 118    PopSui ~~    PopSui              0.01986703 0.008349318
    ## 173        A1 ~~        A2              1.13705974 0.181377065

``` r
pop_psych_soma_constr.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
a12 < 1
A1 ~~ a12*A2
"
pop_psych_soma_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_constr.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   7.197 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma_constr.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_psych_soma_constr.fit$modelfit
```

    ##      chisq df   p_chisq     AIC       CFI     SRMR
    ## df 40.5781 33 0.1709602 84.5781 0.9978057 0.117905

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7357950 0.0792508192209112
    ## 2         A1 =~    PopSui    0.6213223    0.1057176430222
    ## 3         A1 =~    PopAnh    0.8736869 0.0765156687368674
    ## 4         A1 =~  PopGuilt    0.6808365  0.089011984983211
    ## 5         A1 =~   PopConc    0.7679945 0.0986678216361639
    ## 6         A1 ~~        A1    1.0000000                   
    ## 7         A1 ~~        A2    1.0000003  0.141132499767876
    ## 8         A2 =~ PopAppDec    0.1846453 0.0939232566718942
    ## 9         A2 =~ PopAppInc    0.4363189 0.0977850757401601
    ## 10        A2 =~ PopSleDec    0.6278683  0.122381880737089
    ## 11        A2 =~ PopSleInc    0.5482415  0.117355780853881
    ## 12        A2 =~  PopFatig    0.7276579  0.126881601662239
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.4586059  0.131677951348167
    ## 15    PopDep ~~    PopAnh    0.2809442  0.121956938326632
    ## 16    PopSui ~~    PopSui    0.6139587  0.258567172795995
    ## 17    PopAnh ~~    PopAnh    0.2366713  0.136399205876294
    ## 18 PopAppDec ~~ PopAppDec    0.9659057  0.251477286523358
    ## 19 PopAppInc ~~ PopAppInc    0.8096255   0.17290631705096
    ## 20 PopSleDec ~~ PopSleDec    0.6057727  0.330123572288157
    ## 21 PopSleInc ~~ PopSleInc    0.6994260  0.292864378544478
    ## 22  PopFatig ~~  PopFatig    0.4705010  0.323305867238873
    ## 23  PopGuilt ~~  PopGuilt    0.5364622  0.171536381474381
    ## 24   PopConc ~~   PopConc    0.4101863    0.2702224547922

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
PopDep ~~ PopAnh
"
pop_psych_soma_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.064 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_soma_bif.fit$model
```

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 19.30104 24 0.7357501 81.30104   1 0.06566416

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep   0.65273490 0.0879947493287824
    ## 2          A =~    PopSui   0.48424415  0.145463090665366
    ## 3          A =~    PopAnh   0.82379047 0.0707052850702519
    ## 4          A =~ PopAppDec   0.19156961 0.0967300993241066
    ## 5          A =~ PopAppInc   0.45823800 0.0910152721602024
    ## 6          A =~ PopSleDec   0.70126734  0.110643490808361
    ## 7          A =~ PopSleInc   0.59392057  0.113075771149899
    ## 8          A =~  PopFatig   0.77320679  0.108032711745068
    ## 9          A =~  PopGuilt   0.63722296  0.101702366888652
    ## 10         A =~   PopConc   0.90658882  0.119618485817347
    ## 11         A ~~         A   1.00000000                   
    ## 12        A1 =~    PopDep   0.40288551  0.206117064013125
    ## 13        A1 =~    PopSui   0.75735768  0.370827526630514
    ## 14        A1 =~    PopAnh   0.14774125  0.171426700924149
    ## 15        A1 =~  PopGuilt   0.26332226  0.180894678184873
    ## 16        A1 =~   PopConc  -0.29735944  0.230311426119716
    ## 17        A1 ~~        A1   1.00000000                   
    ## 18        A2 =~ PopAppDec   0.11914916  0.200027466516649
    ## 19        A2 =~ PopAppInc   0.17686410  0.175194239286978
    ## 20        A2 =~ PopSleDec   0.69647515   0.35831829159435
    ## 21        A2 =~ PopSleInc  -0.63559173  0.305039379114667
    ## 22        A2 =~  PopFatig  -0.36873940  0.232125680420696
    ## 23        A2 ~~        A2   1.00000000                   
    ## 24    PopDep ~~    PopDep   0.41162024   0.18120245528551
    ## 25    PopDep ~~    PopAnh   0.32655898  0.122413606059505
    ## 26    PopSui ~~    PopSui   0.19191696  0.571442627841785
    ## 27    PopAnh ~~    PopAnh   0.29954175  0.117151176702797
    ## 28 PopAppDec ~~ PopAppDec   0.94910425  0.251104080725651
    ## 29 PopAppInc ~~ PopAppInc   0.75873679  0.195977352175425
    ## 30 PopSleDec ~~ PopSleDec   0.02314588  0.602650982750194
    ## 31 PopSleInc ~~ PopSleInc   0.24328158  0.467301264901229
    ## 32  PopFatig ~~  PopFatig   0.26618140  0.359908764484935
    ## 33  PopGuilt ~~  PopGuilt   0.52460849  0.176553657162092
    ## 34   PopConc ~~   PopConc   0.08967470   0.33243748099494

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
pop_psych_veg.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
"
pop_psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    1.15 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 41.47244 33 0.1478822 85.47244 0.9975468 0.1149995

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7649529 0.0877460088812335
    ## 2         A1 =~    PopSui    0.6202713  0.106120502177647
    ## 3         A1 =~    PopAnh    0.9061497 0.0904406188290761
    ## 4         A1 =~  PopGuilt    0.6746057 0.0894314818186037
    ## 5         A1 ~~        A1    1.0000000                   
    ## 6         A1 ~~        A2    0.9053703   0.10797320766407
    ## 7         A2 =~ PopAppDec    0.1909812 0.0970905126852466
    ## 8         A2 =~ PopAppInc    0.4556835 0.0946454494590175
    ## 9         A2 =~ PopSleDec    0.6552291  0.113786029273202
    ## 10        A2 =~ PopSleInc    0.5689440  0.114775020116461
    ## 11        A2 =~  PopFatig    0.7614467   0.11332843985583
    ## 12        A2 =~   PopConc    0.8091802  0.108009763400932
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.4148471  0.150292567851095
    ## 15    PopDep ~~    PopAnh    0.2306365  0.147433800177904
    ## 16    PopSui ~~    PopSui    0.6152621  0.258901771525553
    ## 17    PopAnh ~~    PopAnh    0.1788925  0.168660973672598
    ## 18 PopAppDec ~~ PopAppDec    0.9635267  0.251249718756764
    ## 19 PopAppInc ~~ PopAppInc    0.7923524  0.174299623045794
    ## 20 PopSleDec ~~ PopSleDec    0.5706749  0.332694013267791
    ## 21 PopSleInc ~~ PopSleInc    0.6763023  0.287503353472929
    ## 22  PopFatig ~~  PopFatig    0.4201995  0.320052552601852
    ## 23  PopGuilt ~~  PopGuilt    0.5449067  0.172137817822894
    ## 24   PopConc ~~   PopConc    0.3452274  0.277230141692826

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
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
"
pop_psych_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##     1.4 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 16.23926 23 0.8448367 80.23926   1 0.07370005

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep    0.6044215 0.0927994004537886
    ## 2          A =~    PopSui    0.4222853  0.143066476019575
    ## 3          A =~    PopAnh    0.8108810 0.0833271537563961
    ## 4          A =~ PopAppDec    0.1984547 0.0999354666394873
    ## 5          A =~ PopAppInc    0.4775928  0.101548708709006
    ## 6          A =~ PopSleDec    0.6960469  0.115740465627768
    ## 7          A =~ PopSleInc    0.6418009  0.121180930554654
    ## 8          A =~  PopFatig    0.8235210   0.12227599625143
    ## 9          A =~  PopGuilt    0.6028376  0.105367035761329
    ## 10         A =~   PopConc    0.8399459  0.108928700908508
    ## 11         A ~~         A    1.0000000                   
    ## 12        A1 =~    PopDep    0.4780921   0.22271931764991
    ## 13        A1 =~    PopSui    0.7683683  0.367447436910188
    ## 14        A1 =~    PopAnh    0.2018479  0.182338090101011
    ## 15        A1 =~  PopGuilt    0.3399819  0.175772028898417
    ## 16        A1 ~~        A1    1.0000000                   
    ## 17        A2 =~ PopAppDec   -0.1397059  0.204427023225692
    ## 18        A2 =~ PopAppInc    0.5471581  0.349724078788114
    ## 19        A2 =~ PopSleDec    0.2431100  0.252900071929552
    ## 20        A2 =~ PopSleInc   -0.3203945  0.231981597806744
    ## 21        A2 =~  PopFatig   -0.4029902   0.31844720874543
    ## 22        A2 =~   PopConc    0.2683290  0.237398141130081
    ## 23        A2 ~~        A2    1.0000000                   
    ## 24    PopDep ~~    PopDep    0.4061026  0.201642634216033
    ## 25    PopDep ~~    PopAnh    0.3371827  0.123843070877741
    ## 26    PopSui ~~    PopSui    0.2312853  0.565285197828883
    ## 27    PopAnh ~~    PopAnh    0.3017294  0.121030631777169
    ## 28 PopAppDec ~~ PopAppDec    0.9410979  0.252177695734157
    ## 29 PopAppInc ~~ PopAppInc    0.4725234  0.429664389786041
    ## 30 PopSleDec ~~ PopSleDec    0.4564163  0.343348715211575
    ## 31 PopSleDec ~~ PopSleInc   -0.3822669  0.245234658454021
    ## 32 PopSleInc ~~ PopSleInc    0.4854390  0.329096298891754
    ## 33  PopFatig ~~  PopFatig    0.1594123  0.438515215181782
    ## 34  PopGuilt ~~  PopGuilt    0.5209993  0.180172863184273
    ## 35   PopConc ~~   PopConc    0.2224904    0.3197174543961

### Affective-Neurovegetative (Elhai Model 2c)

``` r
pop_affect_veg.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
"
pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.337 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 30.35265 33 0.5995863 74.35265   1 0.1123385

``` r
pop_affect_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7816756 0.0817350006188186
    ## 2         A1 =~    PopSui    0.6837606  0.116024418012113
    ## 3         A1 =~  PopGuilt    0.7371762    0.0983352342247
    ## 4         A1 ~~        A1    1.0000000                   
    ## 5         A1 ~~        A2    0.8289466 0.0648763855973035
    ## 6         A2 =~    PopAnh    0.8801195 0.0732212439140747
    ## 7         A2 =~ PopAppDec    0.1885545 0.0966207293963238
    ## 8         A2 =~ PopAppInc    0.4522554  0.090247603271216
    ## 9         A2 =~ PopSleDec    0.6438603  0.107671194832654
    ## 10        A2 =~ PopSleInc    0.5659896  0.111444754853955
    ## 11        A2 =~  PopFatig    0.7562459  0.104205049519063
    ## 12        A2 =~   PopConc    0.7982981  0.100661116967664
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.3889830  0.143675580153658
    ## 15    PopDep ~~    PopAnh    0.3535095  0.108409839891768
    ## 16    PopSui ~~    PopSui    0.5324713  0.265339067835176
    ## 17    PopAnh ~~    PopAnh    0.2253890  0.130825510573078
    ## 18 PopAppDec ~~ PopAppDec    0.9644476  0.251317571310415
    ## 19 PopAppInc ~~ PopAppInc    0.7954645  0.173231705659927
    ## 20 PopSleDec ~~ PopSleDec    0.5854420  0.332410748034128
    ## 21 PopSleInc ~~ PopSleInc    0.6796555  0.285924240977905
    ## 22  PopFatig ~~  PopFatig    0.4280916   0.32106464107293
    ## 23  PopGuilt ~~  PopGuilt    0.4565710  0.182906234773029
    ## 24   PopConc ~~   PopConc    0.3627193  0.269143524901785

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
PopDep ~~ PopAnh
"
pop_affect_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): CFI estimates below 0 should not be trusted, and
    ## indicate that the other model fit estimates should be interpreted with caution.
    ## A negative CFI estimates typically appears due to negative residual variances.

    ## elapsed 
    ##   1.102 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_affect_veg_bif_constr.model <- "
a2 < 1
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A =~ NA*PopDep + PopGuilt + PopSui + a2*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
PopDep ~~ PopAnh
c4a > 0.001
PopSleDec ~~ c4a*PopSleDec
c4b > 0.001
PopSleInc ~~ c4b*PopSleInc
"
pop_affect_veg_bif_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_bif_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.127 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif_constr.model, : A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_affect_veg_bif_constr.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 17.14824 24 0.8422709 79.14824   1 0.06815131

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep  0.627332287 0.0863700007367264
    ## 2          A =~    PopSui  0.510751966  0.108481053294169
    ## 3          A =~    PopAnh  0.856559079 0.0813452432490731
    ## 4          A =~ PopAppDec  0.203986003 0.0994409801322619
    ## 5          A =~ PopAppInc  0.467005284 0.0914896337698255
    ## 6          A =~ PopSleDec  0.764568865  0.125881594340571
    ## 7          A =~ PopSleInc  0.534946357  0.127271456512818
    ## 8          A =~  PopFatig  0.738520614  0.112053975395942
    ## 9          A =~  PopGuilt  0.641209445 0.0944833041121248
    ## 10         A =~   PopConc  0.823003795  0.105970218528199
    ## 11         A ~~         A  1.000000000                   
    ## 12        A1 =~    PopDep  0.451178030  0.265063627659287
    ## 13        A1 =~    PopSui  0.663759684   0.36788525267666
    ## 14        A1 =~  PopGuilt  0.283690739  0.165843630971553
    ## 15        A1 ~~        A1  1.000000000                   
    ## 16        A2 =~    PopAnh  0.210043179  0.111698031732063
    ## 17        A2 =~ PopAppDec -0.122773299  0.191515956721645
    ## 18        A2 =~ PopAppInc -0.112424095  0.175249331694881
    ## 19        A2 =~ PopSleDec -0.692900191  0.348172500388403
    ## 20        A2 =~ PopSleInc  0.642253672  0.256241685091402
    ## 21        A2 =~  PopFatig  0.378446165  0.224545052708205
    ## 22        A2 =~   PopConc -0.083851048  0.190901268363022
    ## 23        A2 ~~        A2  1.000000000                   
    ## 24    PopDep ~~    PopDep  0.402892601  0.227471265316279
    ## 25    PopDep ~~    PopAnh  0.386451141  0.111718449484741
    ## 26    PopSui ~~    PopSui  0.298555201  0.514634431948263
    ## 27    PopAnh ~~    PopAnh  0.222188339  0.140480896247316
    ## 28 PopAppDec ~~ PopAppDec  0.943316728  0.250727900584799
    ## 29 PopAppInc ~~ PopAppInc  0.769266355   0.19003535685593
    ## 30 PopSleDec ~~ PopSleDec  0.001000239  0.597977430761541
    ## 31 PopSleInc ~~ PopSleInc  0.301339046  0.387421642262966
    ## 32  PopFatig ~~  PopFatig  0.311365930  0.340687522675193
    ## 33  PopGuilt ~~  PopGuilt  0.508370094    0.1788505330553
    ## 34   PopConc ~~   PopConc  0.315633448  0.276353067832643

### Model comparisons

``` r
model_fits <- 
data.frame(Model=c('1a', '1b', '2a', '2a(ii)', '2b(i)', '2b(ii)', '2c(i)', '2c(ii)', '3'),
       Name=c('Common',
              'Common (gating)',
              'Psych-Somatic',
              'Psych-Somatic (BiF)',
              'Psych-Neuroveg',
              'Psych-Neuroveg (BiF)',
              'Affect-Neuroveg',
              'Affect-Neuroveg (BiF)',
              'Cog-Mood-Neuroveg')) %>%
bind_cols(
bind_rows(
lapply(list(pop_commonfactor.fit,
            pop_commonfactor_gating.fit,
            pop_psych_soma_constr.fit,
            pop_psych_soma_bif.fit,
            pop_psych_veg.fit,
            pop_psych_veg_bif.fit,
            pop_affect_veg.fit,
            pop_affect_veg_bif_constr.fit,
            pop_kendler_neale_constr.fit),
       function(fit) signif(fit$modelfit))
))
rownames(model_fits) <- NULL
model_fits
```

    ##    Model                  Name   chisq df    p_chisq     AIC      CFI      SRMR
    ## 1     1a                Common 59.8087 35 0.00559503 99.8087 0.992817 0.1257260
    ## 2     1b       Common (gating) 40.5781 34 0.20295200 82.5781 0.998095 0.1179050
    ## 3     2a         Psych-Somatic 40.5781 33 0.17096000 84.5781 0.997806 0.1179050
    ## 4 2a(ii)   Psych-Somatic (BiF) 19.3010 24 0.73575000 81.3010 1.000000 0.0656642
    ## 5  2b(i)        Psych-Neuroveg 41.4724 33 0.14788200 85.4724 0.997547 0.1149990
    ## 6 2b(ii)  Psych-Neuroveg (BiF) 16.2393 23 0.84483700 80.2393 1.000000 0.0737000
    ## 7  2c(i)       Affect-Neuroveg 30.3526 33 0.59958600 74.3526 1.000000 0.1123390
    ## 8 2c(ii) Affect-Neuroveg (BiF) 17.1482 24 0.84227100 79.1482 1.000000 0.0681513
    ## 9      3     Cog-Mood-Neuroveg 38.4850 28 0.08956880 92.4850 0.996964 0.1034540

# Regression models

Regress clinical symptoms on the population symptom factors

``` r
pop_psych_veg_bif_clin_reg_a.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  31.765 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493139 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  Inf . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif_clin_reg_a.model): A difference greater than .025 was observed
    ## pre- and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif_clin_reg_a.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_psych_veg_bif_clin_reg_a.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 325.4651 92 8.014794e-28 413.4651 0.9739853 0.1697943

``` r
pop_psych_veg_bif_clin_reg_a.fit$results[c(1:3,6,7,9)] %>%
filter(op == '~')
```

    ##   lhs op         rhs STD_Genotype   STD_Genotype_SE   p_value
    ## 1   A  ~  ClinAppDec   -0.1572080 0.220790433662912 0.4765732
    ## 2   A  ~  ClinAppInc    0.4392811  0.40460264386533 0.2776314
    ## 3   A  ~  ClinSleDec    0.1257140 0.260378340917268 0.6292028
    ## 4   A  ~  ClinSleInc    0.3921574 0.576628130467399 0.4964694
    ## 5   A  ~ ClinPsycInc    0.4384564   0.4219321254146 0.2987260
    ## 6   A  ~     ClinSui    1.3947971  1.22376987172374 0.2543979

``` r
pop_psych_veg_bif_clin_reg_a1.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A1 ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a1.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a1.model)
```

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_psych_veg_bif_clin_reg_a1.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

``` r
pop_psych_veg_bif_clin_reg_a2.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A2 ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a2.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a2.model)
```

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_psych_veg_bif_clin_reg_a2.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

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

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_pd-1.png)<!-- -->

## PGC/AGDS

Check eigen values of the correlation matrix

``` r
symptoms_clin_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Clin'))

symptoms_clin_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx])) 

plot(symptoms_clin_eigen$values, ylab='Eigenvalue')
lines(symptoms_clin_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_clin_efa_eigen-1.png)<!-- -->

``` r
symptoms_clin_efa <- factanal(covmat=symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx], factors=3, rotation='varimax')
symptoms_clin_efa
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd[symptoms_clin_idx,     symptoms_clin_idx], rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.005       0.173       0.005       0.470       0.918       0.005 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec           0.989  -0.130 
    ## ClinAppInc   0.767  -0.450   0.191 
    ## ClinSleDec   0.940   0.285  -0.174 
    ## ClinSleInc   0.556   0.104   0.459 
    ## ClinPsycInc          0.253   0.130 
    ## ClinSui                      0.997 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      1.782   1.337   1.305
    ## Proportion Var   0.297   0.223   0.217
    ## Cumulative Var   0.297   0.520   0.737
    ## 
    ## The degrees of freedom for the model is 0 and the fit was 0.582

## ALSPAC/UKB

Check eigen values of the correlation matrix

``` r
symptoms_pop_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Pop'))

symptoms_pop_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_pop_idx,symptoms_pop_idx])) 

plot(symptoms_pop_eigen$values, ylab='Eigenvalue')
lines(symptoms_pop_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_pop_efa_eigen-1.png)<!-- -->

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
    ##      0.005      0.005      0.811      0.490      0.773      0.437      0.005 
    ## PopPsycDec   PopFatig   PopGuilt    PopConc     PopSui 
    ##      0.416      0.507      0.614      0.497      0.571 
    ## 
    ## Loadings:
    ##            Factor1 Factor2 Factor3
    ## PopDep      0.944   0.320         
    ## PopAnh      0.746   0.639  -0.173 
    ## PopAppDec   0.170   0.104   0.386 
    ## PopAppInc   0.194   0.216  -0.652 
    ## PopSleDec   0.297   0.361         
    ## PopSleInc           0.745         
    ## PopPsycInc          0.208   0.975 
    ## PopPsycDec  0.354           0.675 
    ## PopFatig    0.175   0.658  -0.171 
    ## PopGuilt    0.421   0.408   0.204 
    ## PopConc     0.208   0.672         
    ## PopSui      0.629           0.161 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      2.381   2.359   2.128
    ## Proportion Var   0.198   0.197   0.177
    ## Cumulative Var   0.198   0.395   0.572
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 8.1788

## All symptoms

Check eigen values of the correlation matrix

``` r
symptoms_cov_pd.eigen <- eigen(cov2cor(symptoms_cov_pd))

signif(symptoms_cov_pd.eigen$values, 3)
```

    ##  [1] 4.98e+00 3.14e+00 2.60e+00 1.90e+00 1.33e+00 1.21e+00 9.80e-01 4.83e-01
    ##  [9] 4.58e-01 3.89e-01 2.73e-01 2.58e-01 1.27e-07 1.24e-07 8.63e-08 7.59e-08
    ## [17] 6.21e-08 3.37e-08

``` r
plot(eigen(cov2cor(symptoms_cov_pd))$values, ylab='Eigenvalue')
lines(eigen(cov2cor(symptoms_cov_pd))$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_eigen-1.png)<!-- -->

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

    ##  [1]  5.141297642  1.829910381  1.213248010  1.018928407  0.859131851
    ##  [6]  0.702233999  0.388599832  0.293480958  0.007868141 -0.028475750
    ## [11] -0.426223470

``` r
S_sim_pos <- S_sim[,sim_cov_keep,sim_cov_keep]

S_sim_pos_ev <- plyr::aaply(S_sim_pos, 1, function(x) eigen(cov2cor(x))$values)

summary(S_sim_pos_ev)
```

    ##        1               2               3               4         
    ##  Min.   :4.298   Min.   :1.564   Min.   :1.127   Min.   :0.8009  
    ##  1st Qu.:4.937   1st Qu.:1.920   1st Qu.:1.349   1st Qu.:1.0508  
    ##  Median :5.289   Median :2.080   Median :1.480   Median :1.1451  
    ##  Mean   :5.490   Mean   :2.104   Mean   :1.477   Mean   :1.1533  
    ##  3rd Qu.:5.843   3rd Qu.:2.253   3rd Qu.:1.596   3rd Qu.:1.2494  
    ##  Max.   :9.541   Max.   :3.431   Max.   :1.867   Max.   :1.5283  
    ##        5                6                7                8           
    ##  Min.   :0.6220   Min.   :0.3968   Min.   :0.1447   Min.   :-0.03045  
    ##  1st Qu.:0.8046   1st Qu.:0.5608   1st Qu.:0.3170   1st Qu.: 0.07870  
    ##  Median :0.8677   Median :0.6339   Median :0.4178   Median : 0.17195  
    ##  Mean   :0.8762   Mean   :0.6320   Mean   :0.4063   Mean   : 0.16820  
    ##  3rd Qu.:0.9343   3rd Qu.:0.7219   3rd Qu.:0.4771   3rd Qu.: 0.24652  
    ##  Max.   :1.1353   Max.   :1.0520   Max.   :0.6923   Max.   : 0.46514  
    ##        9                  10                  11         
    ##  Min.   :-0.27138   Min.   :-0.914782   Min.   :-5.2784  
    ##  1st Qu.:-0.06645   1st Qu.:-0.416622   1st Qu.:-1.2814  
    ##  Median :-0.01120   Median :-0.193797   Median :-0.8558  
    ##  Mean   :-0.02187   Mean   :-0.280802   Mean   :-1.0052  
    ##  3rd Qu.: 0.03430   3rd Qu.:-0.116528   3rd Qu.:-0.5088  
    ##  Max.   : 0.17940   Max.   :-0.003535   Max.   :-0.2225

``` r
symptoms_efa3 <- factanal(covmat=symptoms_cov_pd, factors=3, rotation='varimax')

symptoms_efa3
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.895       0.438       0.545       0.271       0.950       0.226 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.114       0.016       0.567       0.494       0.581       0.673 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.005       0.609       0.600       0.527       0.661 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec                  -0.323 
    ## ClinAppInc   0.140   0.152   0.720 
    ## ClinSleDec                   0.671 
    ## ClinSleInc   0.128   0.690   0.486 
    ## ClinPsycInc  0.163          -0.154 
    ## ClinSui      0.491   0.726         
    ## PopDep       0.903          -0.253 
    ## PopAnh       0.977  -0.161         
    ## PopAppDec    0.169   0.390  -0.503 
    ## PopAppInc    0.316  -0.630         
    ## PopSleDec    0.495           0.405 
    ## PopSleInc    0.537   0.165   0.106 
    ## PopPsycInc           0.993         
    ## PopPsycDec   0.178   0.637  -0.747 
    ## PopFatig     0.586           0.196 
    ## PopGuilt     0.596   0.211         
    ## PopConc      0.622   0.157   0.247 
    ## PopSui       0.526   0.107  -0.228 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      4.142   3.128   2.553
    ## Proportion Var   0.230   0.174   0.142
    ## Cumulative Var   0.230   0.404   0.546
    ## 
    ## The degrees of freedom for the model is 102 and the fit was 87.385

``` r
symptoms_efa4 <- factanal(covmat=symptoms_cov_pd, factors=4, rotation='varimax')

symptoms_efa4
```

    ## 
    ## Call:
    ## factanal(factors = 4, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.327       0.413       0.005       0.253       0.901       0.157 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.089       0.047       0.451       0.467       0.408       0.712 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.005       0.494       0.506       0.497       0.628 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4
    ## ClinAppDec  -0.233           0.121   0.777 
    ## ClinAppInc   0.144           0.681  -0.310 
    ## ClinSleDec  -0.148           0.958   0.230 
    ## ClinSleInc   0.106   0.634   0.567  -0.108 
    ## ClinPsycInc  0.103                   0.297 
    ## ClinSui      0.586   0.707                 
    ## PopDep       0.829  -0.126           0.456 
    ## PopAnh       0.915  -0.215           0.254 
    ## PopAppDec            0.396  -0.173   0.596 
    ## PopAppInc    0.326  -0.649                 
    ## PopSleDec    0.409           0.621   0.198 
    ## PopSleInc    0.507   0.138                 
    ## PopPsycInc   0.137   0.978   0.138         
    ## PopPsycDec   0.194   0.678  -0.531   0.464 
    ## PopFatig     0.654  -0.126          -0.243 
    ## PopGuilt     0.677   0.182                 
    ## PopConc      0.655   0.106   0.225  -0.110 
    ## PopSui       0.548          -0.107   0.229 
    ## 
    ##                Factor1 Factor2 Factor3 Factor4
    ## SS loadings      4.196   3.055   2.516   1.867
    ## Proportion Var   0.233   0.170   0.140   0.104
    ## Cumulative Var   0.233   0.403   0.543   0.646
    ## 
    ## The degrees of freedom for the model is 87 and the fit was 84.9414

``` r
symptoms_efa6 <- factanal(covmat=symptoms_cov_pd, factors=6, rotation='varimax')

symptoms_efa6
```

    ## 
    ## Call:
    ## factanal(factors = 6, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.005       0.176       0.005       0.005       0.864       0.005 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.090       0.005       0.388       0.473       0.246       0.177 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.052       0.005       0.163       0.408       0.357       0.574 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## ClinAppDec          -0.104                   0.986         
    ## ClinAppInc                   0.752          -0.493         
    ## ClinSleDec          -0.108   0.935           0.262   0.154 
    ## ClinSleInc   0.556  -0.183   0.470   0.213           0.620 
    ## ClinPsycInc          0.149           0.105   0.262  -0.165 
    ## ClinSui      0.779   0.304           0.532                 
    ## PopDep               0.888           0.262   0.226         
    ## PopAnh      -0.141   0.865           0.422   0.111   0.181 
    ## PopAppDec    0.496   0.305          -0.130   0.418  -0.277 
    ## PopAppInc   -0.516   0.340           0.205  -0.251  -0.176 
    ## PopSleDec            0.328   0.683   0.317   0.167  -0.224 
    ## PopSleInc            0.375           0.171           0.799 
    ## PopPsycInc   0.943           0.168                   0.155 
    ## PopPsycDec   0.776   0.313  -0.458           0.273         
    ## PopFatig    -0.101   0.201           0.872           0.159 
    ## PopGuilt     0.250   0.518           0.425  -0.176  -0.210 
    ## PopConc      0.112   0.290   0.185   0.710                 
    ## PopSui       0.225   0.585                  -0.129   0.109 
    ## 
    ##                Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## SS loadings      3.104   3.001   2.435   2.253   1.824   1.385
    ## Proportion Var   0.172   0.167   0.135   0.125   0.101   0.077
    ## Cumulative Var   0.172   0.339   0.474   0.600   0.701   0.778
    ## 
    ## The degrees of freedom for the model is 60 and the fit was 79.9357
