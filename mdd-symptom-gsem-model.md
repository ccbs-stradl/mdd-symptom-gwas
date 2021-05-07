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
    ## platform       x86_64-apple-darwin17.0     
    ## arch           x86_64                      
    ## os             darwin17.0                  
    ## system         x86_64, darwin17.0          
    ## status                                     
    ## major          4                           
    ## minor          0.5                         
    ## year           2021                        
    ## month          03                          
    ## day            31                          
    ## svn rev        80133                       
    ## language       R                           
    ## version.string R version 4.0.5 (2021-03-31)
    ## nickname       Shake and Throw

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
    ## ── Column specification ────────────────────────────────────────────────────────
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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## Warning in lav_object_post_check(object): lavaan WARNING: the covariance matrix of the residuals of the observed
    ##                 variables (theta) is not positive definite;
    ##                 use lavInspect(fit, "theta") to investigate.

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   4.736 
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
    ## 1           A1 =~  ClinAppDec  0.671742966 0.433122026645966
    ## 2           A1 =~  ClinAppInc  0.729689563 0.521524441752345
    ## 3           A1 =~  ClinSleDec  0.746888149 0.474570264504222
    ## 4           A1 =~  ClinSleInc  0.612271642 0.461874360406522
    ## 5           A1 =~ ClinPsycInc  0.279997292 0.237173888205013
    ## 6           A1 =~     ClinSui  0.007296054 0.178669809982275
    ## 7           A1 ~~          A1  1.000000000                  
    ## 8   ClinAppDec ~~  ClinAppDec  0.548761368 0.574094487130575
    ## 9   ClinAppDec ~~  ClinAppInc -0.861169631 0.618301627013954
    ## 10  ClinAppInc ~~  ClinAppInc  0.467553163 0.791290516391392
    ## 11  ClinSleDec ~~  ClinSleDec  0.442157966 0.900065138452165
    ## 12  ClinSleInc ~~  ClinSleInc  0.625124988  1.61291792720481
    ## 13 ClinPsycInc ~~ ClinPsycInc  0.921601537 0.470933828537937
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
    ##   1.834 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 1         A1 =~    PopDep  0.863027381 0.0466869327427861
    ## 2         A1 =~    PopSui  0.587162338 0.0983552009075191
    ## 3         A1 =~    PopAnh  0.997079793 0.0429873353388082
    ## 4         A1 =~ PopAppDec  0.170614118 0.0881755751911459
    ## 5         A1 =~ PopAppInc  0.406193508 0.0809098960160876
    ## 6         A1 =~ PopSleDec  0.576324497 0.0948331619991829
    ## 7         A1 =~ PopSleInc  0.520220088  0.101746585361924
    ## 8         A1 =~  PopFatig  0.671995433 0.0937553413527222
    ## 9         A1 =~  PopGuilt  0.627685840  0.080178149750249
    ## 10        A1 =~   PopConc  0.701639753 0.0903611152361146
    ## 11        A1 ~~        A1  1.000000000                   
    ## 12    PopSui ~~    PopSui  0.655234771  0.253105388023692
    ## 13    PopDep ~~    PopDep  0.255183639 0.0738078886665163
    ## 14    PopAnh ~~    PopAnh  0.005831204  0.067689510452544
    ## 15 PopAppDec ~~ PopAppDec  0.970892075  0.251433695546643
    ## 16 PopAppInc ~~ PopAppInc  0.835007396   0.16678277754136
    ## 17 PopSleDec ~~ PopSleDec  0.667846685  0.331077724393483
    ## 18 PopSleInc ~~ PopSleInc  0.729373604  0.284423004694911
    ## 19  PopFatig ~~  PopFatig  0.548420644  0.322517150516019
    ## 20  PopGuilt ~~  PopGuilt  0.606011274   0.16600207711574
    ## 21   PopConc ~~   PopConc  0.507701935  0.268448841475581

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.912 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 1         A1 =~    PopDep    0.7357950 0.0783402516990839
    ## 2         A1 =~    PopSui    0.6213223  0.105725480677445
    ## 3         A1 =~    PopAnh    0.8736869 0.0734377965485503
    ## 4         A1 =~ PopAppDec    0.1846453 0.0936391380070021
    ## 5         A1 =~ PopAppInc    0.4363183 0.0877965382345317
    ## 6         A1 =~ PopSleDec    0.6278670  0.103964416706937
    ## 7         A1 =~ PopSleInc    0.5482410  0.109301879710101
    ## 8         A1 =~  PopFatig    0.7276567    0.1014839388104
    ## 9         A1 =~  PopGuilt    0.6808364 0.0887935292254337
    ## 10        A1 =~   PopConc    0.7679949 0.0985288654752324
    ## 11        A1 ~~        A1    1.0000000                   
    ## 12    PopSui ~~    PopSui    0.6139594  0.258378889036357
    ## 13    PopDep ~~    PopDep    0.4586053  0.130383096615807
    ## 14    PopDep ~~    PopAnh    0.2809443   0.11916607089805
    ## 15    PopAnh ~~    PopAnh    0.2366716  0.130111342800608
    ## 16 PopAppDec ~~ PopAppDec    0.9659057  0.251348401913751
    ## 17 PopAppInc ~~ PopAppInc    0.8096253  0.171640370800988
    ## 18 PopSleDec ~~ PopSleDec    0.6057813  0.333924066420954
    ## 19 PopSleInc ~~ PopSleInc    0.6994331  0.285330009876883
    ## 20  PopFatig ~~  PopFatig    0.4705154  0.324391546022596
    ## 21  PopGuilt ~~  PopGuilt    0.5364621  0.171168778682845
    ## 22   PopConc ~~   PopConc    0.4101842  0.270813269624119

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2
    ##      V3 ~~ V4

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.992 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2
    ##      V5 ~~ V6

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.707 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2
    ##      V3 ~~ V4
    ##      V5 ~~ V6

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   8.496 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 38.15579 32 0.2097686 84.15579 0.9982176 0.1036837

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
pop_cog_mood_neuroveg.model <- "
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
pop_cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. \n            The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note\n            that these results should not be interpreted."
    ##           lhs op       rhs Unstandardized_Estimate
    ## 1          A1 =~  PopGuilt           -1.587319e-01
    ## 2          A1 =~   PopConc           -1.687144e-01
    ## 3          A1 =~    PopSui           -1.074497e-01
    ## 4          A2 =~    PopDep            3.228731e+00
    ## 5          A2 =~    PopAnh            3.979588e+00
    ## 6          A2 =~  PopGuilt           -2.648084e-04
    ## 7          A3 =~ PopSleDec            1.320644e-01
    ## 8          A3 =~ PopSleInc            1.194367e-01
    ## 9          A3 =~  PopFatig            1.680701e-01
    ## 10         A3 =~ PopAppDec            3.698864e-02
    ## 11         A3 =~ PopAppInc            1.212151e-01
    ## 15     PopDep ~~    PopAnh           -1.277175e+01
    ## 16  PopSleDec ~~ PopSleInc           -1.635657e-02
    ## 17  PopAppDec ~~ PopAppInc           -8.481290e-03
    ## 133    PopDep ~~    PopDep           -1.034421e+01
    ## 134    PopAnh ~~    PopAnh           -1.575020e+01
    ## 135 PopAppDec ~~ PopAppDec            3.381163e-02
    ## 136 PopAppInc ~~ PopAppInc            6.345856e-02
    ## 137 PopSleDec ~~ PopSleDec            2.456743e-02
    ## 138 PopSleInc ~~ PopSleInc            3.060018e-02
    ## 139  PopFatig ~~  PopFatig            2.728229e-02
    ## 140  PopGuilt ~~  PopGuilt            3.214609e-02
    ## 141   PopConc ~~   PopConc            2.304101e-02
    ## 142    PopSui ~~    PopSui            2.077680e-02
    ## 195        A1 ~~        A2           -6.689256e-02
    ## 196        A1 ~~        A3           -1.096331e+00
    ## 197        A2 ~~        A3            6.496797e-02

Add constraints to prevent variances from being negative and
correlations from going out of bounds.

``` r
pop_cog_mood_neuroveg_constr.model <- "
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
pop_cog_mood_neuroveg_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##           lhs op       rhs Unstandardized_Estimate           SE
    ## 1          A1 =~  PopGuilt           -0.1789654529  2.945093455
    ## 2          A1 =~   PopConc           -0.1728821919  0.039931740
    ## 3          A1 =~    PopSui           -0.1105456607  0.026597695
    ## 4          A2 =~    PopDep            0.2378590271  5.706871077
    ## 5          A2 =~    PopAnh            0.2931343534  7.033628625
    ## 6          A2 =~  PopGuilt           -0.0170870858  2.982141222
    ## 7          A3 =~ PopSleDec            0.1359365347  0.024928134
    ## 8          A3 =~ PopSleInc            0.1235493294  0.025118011
    ## 9          A3 =~  PopFatig            0.1732064215  0.029057216
    ## 10         A3 =~ PopAppDec            0.0382316494  0.018124152
    ## 11         A3 =~ PopAppInc            0.1242794193  0.026775043
    ## 15     PopDep ~~    PopAnh            0.0075474437  3.346286169
    ## 16  PopSleDec ~~ PopSleInc           -0.0173780926  0.009440301
    ## 17  PopAppDec ~~ PopAppInc           -0.0087491155  0.007410719
    ## 18     PopAnh ~~    PopAnh            0.0009996403  4.123778087
    ## 19         A1 ~~        A3           -1.0000000119  0.223640520
    ## 135    PopDep ~~    PopDep            0.0239201907  2.715598418
    ## 136 PopAppDec ~~ PopAppDec            0.0337181322  0.008781122
    ## 137 PopAppInc ~~ PopAppInc            0.0627062645  0.013361339
    ## 138 PopSleDec ~~ PopSleDec            0.0235296783  0.013882201
    ## 139 PopSleInc ~~ PopSleInc            0.0296008727  0.012937359
    ## 140  PopFatig ~~  PopFatig            0.0255293723  0.017891865
    ## 141  PopGuilt ~~  PopGuilt            0.0304121608  0.012661327
    ## 142   PopConc ~~   PopConc            0.0216173081  0.017033333
    ## 143    PopSui ~~    PopSui            0.0201018964  0.009476736
    ## 196        A1 ~~        A2           -0.8823464146 21.025808147
    ## 197        A2 ~~        A3            0.8569694096 20.569677228

``` r
pop_cog_mood_neuroveg_constr.fit$modelfit
```

    ## NULL

``` r
pop_cog_mood_neuroveg_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ## NULL

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
    ## 1          A1 =~    PopDep              0.20414169 0.021805597
    ## 2          A1 =~    PopAnh              0.25222881 0.021822422
    ## 3          A1 =~  PopGuilt              0.16343162 0.021214807
    ## 4          A1 =~   PopConc              0.17576048 0.022291015
    ## 5          A1 =~    PopSui              0.11160320 0.018944335
    ## 6          A2 =~ PopAppDec              0.03170025 0.016348244
    ## 7          A2 =~ PopAppInc              0.11236623 0.027040801
    ## 8          A2 =~ PopSleDec              0.11837437 0.025122248
    ## 9          A2 =~ PopSleInc              0.10665394 0.024434677
    ## 10         A2 =~  PopFatig              0.15586184 0.030523394
    ## 13     PopDep ~~    PopAnh              0.02578165 0.009731840
    ## 109    PopDep ~~    PopDep              0.03882329 0.010197927
    ## 110    PopAnh ~~    PopAnh              0.02329887 0.011290937
    ## 111 PopAppDec ~~ PopAppDec              0.03417490 0.008844263
    ## 112 PopAppInc ~~ PopAppInc              0.06552557 0.013258839
    ## 113 PopSleDec ~~ PopSleDec              0.02799603 0.013850275
    ## 114 PopSleInc ~~ PopSleInc              0.03349031 0.013092086
    ## 115  PopFatig ~~  PopFatig              0.03123687 0.017848818
    ## 116  PopGuilt ~~  PopGuilt              0.03062639 0.009824874
    ## 117   PopConc ~~   PopConc              0.02061381 0.013897762
    ## 118    PopSui ~~    PopSui              0.01986699 0.008349319
    ## 173        A1 ~~        A2              1.13706092 0.181377450

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  16.411 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI     SRMR
    ## df 40.57811 33 0.1709599 84.57811 0.9978057 0.117905

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7357947 0.0792508060607973
    ## 2         A1 =~    PopSui    0.6213224  0.105717657956796
    ## 3         A1 =~    PopAnh    0.8736866 0.0765156558539507
    ## 4         A1 =~  PopGuilt    0.6808366 0.0890120000036697
    ## 5         A1 =~   PopConc    0.7679945 0.0986678353590292
    ## 6         A1 ~~        A1    1.0000000                   
    ## 7         A1 ~~        A2    0.9999998  0.141132336165387
    ## 8         A2 =~ PopAppDec    0.1846453 0.0939232980445869
    ## 9         A2 =~ PopAppInc    0.4363190 0.0977850774686408
    ## 10        A2 =~ PopSleDec    0.6278687  0.122381898548638
    ## 11        A2 =~ PopSleInc    0.5482416  0.117355795942911
    ## 12        A2 =~  PopFatig    0.7276584  0.126881602395355
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopSui ~~    PopSui    0.6139587  0.258567176786584
    ## 15    PopDep ~~    PopDep    0.4586060  0.131677897554377
    ## 16    PopDep ~~    PopAnh    0.2809444  0.121956881492902
    ## 17    PopAnh ~~    PopAnh    0.2366715  0.136399144725743
    ## 18 PopAppDec ~~ PopAppDec    0.9659057  0.251477285501187
    ## 19 PopAppInc ~~ PopAppInc    0.8096252   0.17290632363747
    ## 20 PopSleDec ~~ PopSleDec    0.6057697  0.330123585820198
    ## 21 PopSleInc ~~ PopSleInc    0.6994253  0.292864381021172
    ## 22  PopFatig ~~  PopFatig    0.4704967  0.323305898378957
    ## 23  PopGuilt ~~  PopGuilt    0.5364618    0.1715363882548
    ## 24   PopConc ~~   PopConc    0.4101865  0.270222449521034

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.383 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df p_chisq      AIC CFI       SRMR
    ## df 19.30104 24 0.73575 81.30104   1 0.06566416

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep   0.40288551  0.206117064013116
    ## 2         A1 =~    PopSui   0.75735768  0.370827526630537
    ## 3         A1 =~    PopAnh   0.14774125  0.171426700924138
    ## 4         A1 =~  PopGuilt   0.26332226  0.180894678184867
    ## 5         A1 =~   PopConc  -0.29735944  0.230311426119685
    ## 6         A1 ~~        A1   1.00000000                   
    ## 7         A2 =~ PopAppDec   0.11914916   0.20002746651665
    ## 8         A2 =~ PopAppInc   0.17686410  0.175194239286981
    ## 9         A2 =~ PopSleDec   0.69647515  0.358318291594383
    ## 10        A2 =~ PopSleInc  -0.63559173  0.305039379114649
    ## 11        A2 =~  PopFatig  -0.36873940   0.23212568042069
    ## 12        A2 ~~        A2   1.00000000                   
    ## 13         A =~    PopDep   0.65273490 0.0879947493287794
    ## 14         A =~    PopSui   0.48424415  0.145463090665367
    ## 15         A =~    PopAnh   0.82379047 0.0707052850702512
    ## 16         A =~ PopAppDec   0.19156961 0.0967300993241071
    ## 17         A =~ PopAppInc   0.45823800 0.0910152721602025
    ## 18         A =~ PopSleDec   0.70126734  0.110643490808361
    ## 19         A =~ PopSleInc   0.59392057  0.113075771149899
    ## 20         A =~  PopFatig   0.77320679  0.108032711745068
    ## 21         A =~  PopGuilt   0.63722296  0.101702366888651
    ## 22         A =~   PopConc   0.90658882  0.119618485817343
    ## 23         A ~~         A   1.00000000                   
    ## 24    PopSui ~~    PopSui   0.19191696  0.571442627841873
    ## 25    PopDep ~~    PopDep   0.41162024  0.181202455285498
    ## 26    PopDep ~~    PopAnh   0.32655898  0.122413606059493
    ## 27    PopAnh ~~    PopAnh   0.29954175  0.117151176702791
    ## 28 PopAppDec ~~ PopAppDec   0.94910425  0.251104080725653
    ## 29 PopAppInc ~~ PopAppInc   0.75873679  0.195977352175428
    ## 30 PopSleDec ~~ PopSleDec   0.02314588  0.602650982750261
    ## 31 PopSleInc ~~ PopSleInc   0.24328158  0.467301264901184
    ## 32  PopFatig ~~  PopFatig   0.26618140  0.359908764484926
    ## 33  PopGuilt ~~  PopGuilt   0.52460849  0.176553657162093
    ## 34   PopConc ~~   PopConc   0.08967470  0.332437480994914

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
    ## [1] "The model as initially specified failed to converge. A lower bound of 0 on residual variances has been automatically added to try and troubleshoot this."
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   4.953 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 41.47246 33 0.1478819 85.47246 0.9975468 0.1149995

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7649521 0.0877458880515884
    ## 2         A1 =~    PopSui    0.6202719  0.106120500711077
    ## 3         A1 =~    PopAnh    0.9061493 0.0904405063339581
    ## 4         A1 =~  PopGuilt    0.6746063 0.0894314921001132
    ## 5         A1 ~~        A1    1.0000000                   
    ## 6         A1 ~~        A2    0.9053709    0.1079732332199
    ## 7         A2 =~ PopAppDec    0.1909810 0.0970904966486274
    ## 8         A2 =~ PopAppInc    0.4556834 0.0946454384386007
    ## 9         A2 =~ PopSleDec    0.6552291  0.113786022748889
    ## 10        A2 =~ PopSleInc    0.5689442  0.114775012759428
    ## 11        A2 =~  PopFatig    0.7614473  0.113328458729181
    ## 12        A2 =~   PopConc    0.8091793  0.108009727085244
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopSui ~~    PopSui    0.6152586  0.258901791408545
    ## 15    PopDep ~~    PopDep    0.4148482  0.150292284007407
    ## 16    PopDep ~~    PopAnh    0.2306375  0.147433522855858
    ## 17    PopAnh ~~    PopAnh    0.1788930  0.168660707644034
    ## 18 PopAppDec ~~ PopAppDec    0.9635272  0.251249716038295
    ## 19 PopAppInc ~~ PopAppInc    0.7923540  0.174299602628215
    ## 20 PopSleDec ~~ PopSleDec    0.5706744  0.332694012772975
    ## 21 PopSleInc ~~ PopSleInc    0.6762997  0.287503353045202
    ## 22  PopFatig ~~  PopFatig    0.4201976  0.320052576512712
    ## 23  PopGuilt ~~  PopGuilt    0.5449062  0.172137848883154
    ## 24   PopConc ~~   PopConc    0.3452292  0.277230073850468

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
    ## [1] "The model as initially specified failed to converge. A lower bound of 0 on residual variances has been automatically added to try and troubleshoot this."
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2
    ##      V5 ~~ V6

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   5.133 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 16.23924 23 0.8448373 80.23924   1 0.07370012

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.4780924  0.222719468075643
    ## 2         A1 =~    PopSui    0.7683680  0.367447331323036
    ## 3         A1 =~    PopAnh    0.2018478  0.182338196770063
    ## 4         A1 =~  PopGuilt    0.3399818  0.175771995964177
    ## 5         A1 ~~        A1    1.0000000                   
    ## 6         A2 =~ PopAppDec   -0.1397071  0.204427059190325
    ## 7         A2 =~ PopAppInc    0.5471604  0.349726108401296
    ## 8         A2 =~ PopSleDec    0.2431086  0.252899292160142
    ## 9         A2 =~ PopSleInc   -0.3203907  0.231980750952747
    ## 10        A2 =~  PopFatig   -0.4029905  0.318447886359844
    ## 11        A2 =~   PopConc    0.2683276   0.23739800598392
    ## 12        A2 ~~        A2    1.0000000                   
    ## 13         A =~    PopDep    0.6044215 0.0927994316490413
    ## 14         A =~    PopSui    0.4222854  0.143066471814352
    ## 15         A =~    PopAnh    0.8108813 0.0833271984696414
    ## 16         A =~ PopAppDec    0.1984548 0.0999354759958338
    ## 17         A =~ PopAppInc    0.4775928  0.101548738001249
    ## 18         A =~ PopSleDec    0.6960472  0.115740449053474
    ## 19         A =~ PopSleInc    0.6418006  0.121180883507797
    ## 20         A =~  PopFatig    0.8235212   0.12227600122861
    ## 21         A =~  PopGuilt    0.6028376  0.105367027798162
    ## 22         A =~   PopConc    0.8399455  0.108928679966573
    ## 23         A ~~         A    1.0000000                   
    ## 24    PopSui ~~    PopSui    0.2312857  0.565284826329201
    ## 25    PopDep ~~    PopDep    0.4061025  0.201642851862343
    ## 26    PopDep ~~    PopAnh    0.3371826  0.123843133444436
    ## 27    PopAnh ~~    PopAnh    0.3017290  0.121030725047075
    ## 28 PopAppDec ~~ PopAppDec    0.9410978   0.25217770615995
    ## 29 PopAppInc ~~ PopAppInc    0.4725207  0.429667513701096
    ## 30 PopSleDec ~~ PopSleDec    0.4564148  0.343348509545384
    ## 31 PopSleDec ~~ PopSleInc   -0.3822687  0.245233877109488
    ## 32 PopSleInc ~~ PopSleInc    0.4854431  0.329095165259959
    ## 33  PopFatig ~~  PopFatig    0.1594105  0.438515563032429
    ## 34  PopGuilt ~~  PopGuilt    0.5209994  0.180172855710733
    ## 35   PopConc ~~   PopConc    0.2224929  0.319717100151645

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## Warning in lav_object_post_check(object): lavaan WARNING: the covariance matrix of the residuals of the observed
    ##                 variables (theta) is not positive definite;
    ##                 use lavInspect(fit, "theta") to investigate.

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.023 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 1         A1 =~    PopDep    0.7816756 0.0817350006188179
    ## 2         A1 =~    PopSui    0.6837606  0.116024418012115
    ## 3         A1 =~  PopGuilt    0.7371762 0.0983352342246932
    ## 4         A1 ~~        A1    1.0000000                   
    ## 5         A1 ~~        A2    0.8289466 0.0648763855973103
    ## 6         A2 =~    PopAnh    0.8801195 0.0732212439140778
    ## 7         A2 =~ PopAppDec    0.1885545  0.096620729396323
    ## 8         A2 =~ PopAppInc    0.4522554 0.0902476032712158
    ## 9         A2 =~ PopSleDec    0.6438603   0.10767119483265
    ## 10        A2 =~ PopSleInc    0.5659896  0.111444754853953
    ## 11        A2 =~  PopFatig    0.7562459  0.104205049519062
    ## 12        A2 =~   PopConc    0.7982981  0.100661116967663
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopSui ~~    PopSui    0.5324713  0.265339067835184
    ## 15    PopDep ~~    PopDep    0.3889830   0.14367558015365
    ## 16    PopDep ~~    PopAnh    0.3535095  0.108409839891775
    ## 17    PopAnh ~~    PopAnh    0.2253890  0.130825510573089
    ## 18 PopAppDec ~~ PopAppDec    0.9644476  0.251317571310416
    ## 19 PopAppInc ~~ PopAppInc    0.7954645   0.17323170565993
    ## 20 PopSleDec ~~ PopSleDec    0.5854420  0.332410748034123
    ## 21 PopSleInc ~~ PopSleInc    0.6796555  0.285924240977901
    ## 22  PopFatig ~~  PopFatig    0.4280916  0.321064641072931
    ## 23  PopGuilt ~~  PopGuilt    0.4565710  0.182906234773013
    ## 24   PopConc ~~   PopConc    0.3627193  0.269143524901783

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

    ## [1] "Calculating SRMR"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): CFI estimates below 0 should not be trusted, and
    ## indicate that the other model fit estimates should be interpreted with caution.
    ## A negative CFI estimates typically appears due to negative residual variances.

    ## elapsed 
    ##   2.255 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V1 ~~ V2

    ## Warning in lav_object_post_check(object): lavaan WARNING: the covariance matrix of the residuals of the observed
    ##                 variables (theta) is not positive definite;
    ##                 use lavInspect(fit, "theta") to investigate.

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  13.413 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0024174380515027 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975605 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 17.14824 24 0.8422711 79.14824   1 0.06815131

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep  0.451178628  0.265064088405363
    ## 2         A1 =~    PopSui  0.663758134  0.367884699599379
    ## 3         A1 =~  PopGuilt  0.283690813  0.165843803258948
    ## 4         A1 ~~        A1  1.000000000                   
    ## 5         A2 =~    PopAnh  0.210043174  0.111698145298657
    ## 6         A2 =~ PopAppDec -0.122774695  0.191515908957514
    ## 7         A2 =~ PopAppInc -0.112424220  0.175249290193913
    ## 8         A2 =~ PopSleDec -0.692902580  0.348173203445454
    ## 9         A2 =~ PopSleInc  0.642251155  0.256241222236777
    ## 10        A2 =~  PopFatig  0.378444495  0.224544929636592
    ## 11        A2 =~   PopConc -0.083853089  0.190901194020524
    ## 12        A2 ~~        A2  1.000000000                   
    ## 13         A =~    PopDep  0.627332264 0.0863700049620385
    ## 14         A =~    PopSui  0.510752133  0.108481023826643
    ## 15         A =~    PopAnh  0.856559701 0.0813453309867796
    ## 16         A =~ PopAppDec  0.203985908 0.0994410074610901
    ## 17         A =~ PopAppInc  0.467004968 0.0914896078852279
    ## 18         A =~ PopSleDec  0.764568178  0.125881683803069
    ## 19         A =~ PopSleInc  0.534946322   0.12727131402426
    ## 20         A =~  PopFatig  0.738520912  0.112053917703521
    ## 21         A =~  PopGuilt  0.641209241  0.094483291435049
    ## 22         A =~   PopConc  0.823003566  0.105970204230425
    ## 23         A ~~         A  1.000000000                   
    ## 24    PopSui ~~    PopSui  0.298556420  0.514632769007461
    ## 25    PopDep ~~    PopDep  0.402892104   0.22747197361689
    ## 26    PopDep ~~    PopAnh  0.386450482   0.11171854720895
    ## 27    PopAnh ~~    PopAnh  0.222186546  0.140481190565615
    ## 28 PopAppDec ~~ PopAppDec  0.943316530  0.250727940520103
    ## 29 PopAppInc ~~ PopAppInc  0.769266419  0.190035280654748
    ## 30 PopSleDec ~~ PopSleDec  0.001000422   0.59797908693848
    ## 31 PopSleInc ~~ PopSleInc  0.301345592  0.387420414870047
    ## 32  PopFatig ~~  PopFatig  0.311366492  0.340687466384187
    ## 33  PopGuilt ~~  PopGuilt  0.508370008  0.178850574757552
    ## 34   PopConc ~~   PopConc  0.315633818  0.276353178968614

### Model comparisons

``` r
model_fits <- 
data.frame(Model=c('1a', '1b', '1c', '1d', '1e',
                   '2a', '2a(ii)', '2b(i)', '2b(ii)', '2c(i)', '2c(ii)'),
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
              'Affect-Neuroveg (BiF)')) %>%
bind_cols(
bind_rows(
lapply(list(pop_commonfactor.fit,
            pop_commonfactor_gating.fit,
            pop_commonfactor_app.fit,
            pop_commonfactor_sle.fit,
            pop_commonfactor_app_sle.fit,
            pop_psych_soma_constr.fit,
            pop_psych_soma_bif.fit,
            pop_psych_veg.fit,
            pop_psych_veg_bif.fit,
            pop_affect_veg.fit,
            pop_affect_veg_bif_constr.fit),
       function(fit) signif(fit$modelfit))
))
rownames(model_fits) <- NULL
model_fits %>%
select(-chisq, -df) %>%
mutate(dAIC=AIC-min(AIC))
```

    ##     Model                  Name    p_chisq     AIC      CFI      SRMR    dAIC
    ## 1      1a                Common 0.00559503 99.8087 0.992817 0.1257260 25.4561
    ## 2      1b       Common (gating) 0.20295200 82.5781 0.998095 0.1179050  8.2255
    ## 3      1c          Common (App) 0.19126100 83.8628 0.998013 0.1155380  9.5102
    ## 4      1d          Common (Sle) 0.24035000 82.3327 0.998456 0.1063220  7.9801
    ## 5      1e       Common(App,Sle) 0.20976900 84.1558 0.998218 0.1036840  9.8032
    ## 6      2a         Psych-Somatic 0.17096000 84.5781 0.997806 0.1179050 10.2255
    ## 7  2a(ii)   Psych-Somatic (BiF) 0.73575000 81.3010 1.000000 0.0656642  6.9484
    ## 8   2b(i)        Psych-Neuroveg 0.14788200 85.4725 0.997547 0.1149990 11.1199
    ## 9  2b(ii)  Psych-Neuroveg (BiF) 0.84483700 80.2392 1.000000 0.0737001  5.8866
    ## 10  2c(i)       Affect-Neuroveg 0.59958600 74.3526 1.000000 0.1123390  0.0000
    ## 11 2c(ii) Affect-Neuroveg (BiF) 0.84227100 79.1482 1.000000 0.0681513  4.7956

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

    ## Warning in lav_partable_flat(FLAT, blocks = "group", meanstructure = meanstructure, : duplicated elements in model syntax have been ignored:
    ##      V9 ~~ V10
    ##      V11 ~~ V12

    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  97.508 
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
    ## df 325.4652 92 8.014343e-28 413.4652 0.9739853 0.1697943

``` r
pop_psych_veg_bif_clin_reg_a.fit$results[c(1:3,6,7,9)] %>%
filter(op == '~')
```

    ##   lhs op         rhs STD_Genotype   STD_Genotype_SE   p_value
    ## 1   A  ~  ClinAppDec   -0.1572096 0.220793965881341 0.4765777
    ## 2   A  ~  ClinAppInc    0.4392862 0.404612517012996 0.2776316
    ## 3   A  ~  ClinSleDec    0.1257157 0.260381625537067 0.6292001
    ## 4   A  ~  ClinSleInc    0.3921621 0.576637463883157 0.4964695
    ## 5   A  ~ ClinPsycInc    0.4384630  0.42194342353807 0.2987264
    ## 6   A  ~     ClinSui    1.3948227  1.22381741221557 0.2543979

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

    ##  [1]  4.440224843  1.728123652  1.091673444  0.888999868  0.752762364
    ##  [6]  0.645909858  0.293496010  0.188301343  0.006100351 -0.035591734

``` r
S_sim_pos <- S_sim[,sim_cov_keep,sim_cov_keep]

S_sim_pos_ev <- plyr::aaply(S_sim_pos, 1, function(x) eigen(cov2cor(x))$values)

summary(S_sim_pos_ev)
```

    ##        1               2               3                4         
    ##  Min.   :3.602   Min.   :1.517   Min.   :0.9986   Min.   :0.7813  
    ##  1st Qu.:4.333   1st Qu.:1.770   1st Qu.:1.1963   1st Qu.:0.9452  
    ##  Median :4.542   Median :1.910   Median :1.3094   Median :1.0063  
    ##  Mean   :4.678   Mean   :1.909   Mean   :1.3182   Mean   :1.0238  
    ##  3rd Qu.:4.880   3rd Qu.:2.032   3rd Qu.:1.4086   3rd Qu.:1.0771  
    ##  Max.   :6.580   Max.   :2.507   Max.   :1.8286   Max.   :1.4408  
    ##        5                6                7                  8            
    ##  Min.   :0.6076   Min.   :0.3094   Min.   :-0.01965   Min.   :-0.105049  
    ##  1st Qu.:0.7154   1st Qu.:0.4867   1st Qu.: 0.22877   1st Qu.: 0.005558  
    ##  Median :0.7723   Median :0.5513   Median : 0.32151   Median : 0.050654  
    ##  Mean   :0.7815   Mean   :0.5577   Mean   : 0.30368   Mean   : 0.069745  
    ##  3rd Qu.:0.8532   3rd Qu.:0.6382   3rd Qu.: 0.38630   3rd Qu.: 0.121469  
    ##  Max.   :1.0505   Max.   :0.7728   Max.   : 0.53299   Max.   : 0.349140  
    ##        9                  10          
    ##  Min.   :-0.64829   Min.   :-2.16665  
    ##  1st Qu.:-0.19294   1st Qu.:-0.61902  
    ##  Median :-0.09775   Median :-0.43412  
    ##  Mean   :-0.13334   Mean   :-0.50801  
    ##  3rd Qu.:-0.03422   3rd Qu.:-0.26716  
    ##  Max.   : 0.05460   Max.   :-0.07708

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
    ##       0.895       0.439       0.545       0.271       0.950       0.226 
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
    ## SS loadings      4.142   3.129   2.553
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
    ##       0.052       0.005       0.163       0.408       0.357       0.575 
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
