GenomicSEM of MDD symptoms sample counts
================
Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel
Nivard

Make table of counts of sypmtom presence/absence used in final analyses

``` r
library(readxl)
library(readr)
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.1.1

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(stringr)
```

File for aligning sumstats filenames to cohort/symptom

``` r
cohort_alignment <- read_tsv('meta/cohort_alignment.txt')
```

    ## Rows: 51 Columns: 3── Column specification ────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): filename, cohort, reference
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

List of `basic.xls` files from Ricopili with sample counts

``` r
basic_xls_list <- list.files('meta/distribution', pattern='basic.*.xls', full.names=TRUE)

basics <- bind_rows(lapply(basic_xls_list, read_excel)) %>%
filter(Dataset != 'SUM')
```

Merge and reformat

``` r
basics_aligned <- 
cohort_alignment %>%
filter(!reference %in% c('MDD3', 'MDD4')) %>%
mutate(Dataset=str_replace(str_remove(filename, pattern='daner_'), 'gz', 'aligned.gz')) %>%
left_join(basics, by='Dataset') 

presence_absence <- basics_aligned %>% str_glue_data("{N_cases}:{N_controls}")

basics_aligned %>%
mutate(PresenceAbsence=as.character(presence_absence)) %>%
select(cohort, reference, PresenceAbsence) %>%
pivot_wider(names_from=cohort, values_from=PresenceAbsence)
```

    ## # A tibble: 12 × 5
    ##    reference PGC        AGDS       ALSPAC    UKBB       
    ##    <chr>     <chr>      <chr>      <chr>     <chr>      
    ##  1 MDD1      11669:1152 14151:471  328:2966  66558:54647
    ##  2 MDD2      10887:1456 13601:965  610:2684  47769:73393
    ##  3 MDD3a     6520:5873  2487:8087  44:3250   24632:33816
    ##  4 MDD3b     2684:8688  4826:5748  52:3242   13813:44635
    ##  5 MDD4a     9332:3141  10036:3169 668:2626  43987:14751
    ##  6 MDD4b     3210:7163  7349:5848  464:2830  9750:48988 
    ##  7 MDD5a     5072:5032  5614:7538  NA:NA     <NA>       
    ##  8 MDD5b     5295:5911  7216:5927  NA:NA     <NA>       
    ##  9 MDD6      10913:1579 12610:600  1055:2239 50575:11159
    ## 10 MDD7      9363:2862  12270:865  386:2908  32741:31462
    ## 11 MDD8      10281:1421 12527:599  370:2924  47682:13019
    ## 12 MDD9      6721:5631  10439:2699 367:2927  33617:31108

``` r
basics_aligned %>%
mutate(sample=if_else(cohort %in% c('AGDS', 'PGC'), true='enriched', false='unselected')) %>%
group_by(sample, reference) %>%
summarise(N_cases=sum(N_cases), N_controls=sum(N_controls)) %>%
group_by(sample) %>%
summarise(minCa=min(N_cases, na.rm=T), maxCa=max(N_cases, na.rm=T), minCo=min(N_controls, na.rm=T), maxCo=max(N_controls, na.rm=T))
```

    ## `summarise()` has grouped output by 'sample'. You can override using the
    ## `.groups` argument.

    ## # A tibble: 2 × 5
    ##   sample     minCa maxCa minCo maxCo
    ##   <chr>      <dbl> <dbl> <dbl> <dbl>
    ## 1 enriched    7510 25820  1623 14436
    ## 2 unselected 10214 66886 13398 76077
