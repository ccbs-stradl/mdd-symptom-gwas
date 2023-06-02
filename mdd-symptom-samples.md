GenomicSEM of MDD symptoms sample counts
================
Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel
Nivard

Make table of counts of sypmtom presence/absence used in final analyses

``` r
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
```

File for aligning sumstats filenames to cohort/symptom

``` r
cohort_alignment <- read_tsv('meta/cohort_alignment.txt')
```

    ## Rows: 245 Columns: 3
    ## ── Column specification ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): filename, cohort, reference
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

List of `basic.xls` files from Ricopili with sample counts

``` r
basic_xls_list <- list.files('meta/distribution', pattern='basic.*.xls', full.names=TRUE, recursive=TRUE)

names(basic_xls_list) <- str_match(basename(basic_xls_list), "basic\\.([A-Za-z_]+\\.MDD[0-9A-Za-z_]+)")[,2]

basics <- bind_rows(lapply(basic_xls_list, read_excel), .id='meta')

basics_datasets <- basics %>%
    filter(Dataset != 'SUM') %>%
    select(-meta) |>
    distinct()
    
basics_meta <-
basics %>%
    filter(Dataset == 'SUM') %>%
    mutate(Dataset=meta) %>%
    select(-meta)
```

Merge and reformat

``` r
basics_aligned <- 
cohort_alignment %>%
filter(!reference %in% c('MDD3', 'MDD4')) %>%
mutate(Dataset=str_replace(str_remove(filename, pattern='daner_'), 'gz', 'align.gz')) %>%
left_join(basics_datasets, by='Dataset') |>
mutate(cohort = if_else(cohort == "janpy", true = "PGC", false=cohort)) |>
group_by(cohort, reference) |>
summarize(N_cases=sum(N_cases, na.rm=T), N_controls=sum(N_controls, na.rm=T), N_eff_half=sum(N_eff_half, na.rm=T)) |>
ungroup()
```

    ## `summarise()` has grouped output by 'cohort'. You can override using the
    ## `.groups` argument.

``` r
presence_absence <- basics_aligned %>% str_glue_data("{N_cases}:{N_controls} {round(100*N_cases/(N_cases + N_controls))}")

basics_formatted <- 
basics_aligned %>%
mutate(PresenceAbsence=as.character(presence_absence)) %>%
select(cohort, reference, PresenceAbsence) %>%
pivot_wider(names_from=cohort, values_from=PresenceAbsence)

knitr::kable(basics_formatted)
```

| reference | AGDS          | ALSPAC       | EstBB          | GenScot      | PGC          | UKBB           |
|:----------|:--------------|:-------------|:---------------|:-------------|:-------------|:---------------|
| MDD1      | 14151:471 97  | 328:2966 10  | 41070:41867 50 | 3334:159 95  | 4196:1118 79 | 66558:54647 55 |
| MDD2      | 13601:965 93  | 610:2684 19  | 32734:50090 40 | 3229:261 93  | 7902:1575 83 | 47769:73393 39 |
| MDD3a     | 2487:8087 24  | 44:3250 1    | 14777:9527 61  | 1141:462 71  | 5637:6045 48 | 24632:23720 51 |
| MDD3b     | 4826:5748 46  | 52:3242 2    | 8747:9527 48   | 431:510 46   | 2645:6909 28 | 13813:23720 37 |
| MDD4a     | 10036:3169 76 | 668:2626 20  | 28489:5182 85  | 1368:254 84  | 7513:3150 70 | 43987:12043 79 |
| MDD4b     | 7349:5848 56  | 464:2830 14  | 9911:5182 66   | 359:289 55   | 2878:4913 37 | 9750:12043 45  |
| MDD5a     | 5614:7538 43  | 113:3181 3   | NA             | 668:547 55   | 4165:4287 49 | NA             |
| MDD5b     | 7216:5927 55  | 299:2995 9   | NA             | 687:526 57   | 4798:4761 50 | NA             |
| MDD6      | 12610:600 95  | 1055:2239 32 | 33674:3338 91  | 3165:209 94  | 8166:1688 83 | 50575:11159 82 |
| MDD7      | 12270:865 93  | 386:2908 12  | 28630:9200 76  | 2847:520 85  | 6804:2503 73 | 32741:31462 51 |
| MDD8      | 12527:599 95  | 370:2924 11  | 27138:7473 78  | 3097:277 92  | 8350:1510 85 | 47682:13019 79 |
| MDD9      | 10439:2699 79 | 367:2927 11  | 13000:24850 34 | 2135:1230 63 | 5596:5680 50 | 33617:31108 52 |

``` r
basics_sum <-
basics_aligned %>%
mutate(sample=if_else(cohort %in% c('AGDS', 'PGC', 'GenScot'), true='enriched', false='unselected')) %>%
group_by(sample, reference) %>%
summarise(N_cases=sum(N_cases), N_controls=sum(N_controls)) 
```

    ## `summarise()` has grouped output by 'sample'. You can override using the
    ## `.groups` argument.

``` r
basics_sum %>%
group_by(sample) %>%
summarise(minCa=min(N_cases, na.rm=T), maxCa=max(N_cases, na.rm=T), minCo=min(N_controls, na.rm=T), maxCo=max(N_controls, na.rm=T))
```

    ## # A tibble: 2 × 5
    ##   sample     minCa  maxCa minCo  maxCo
    ##   <chr>      <dbl>  <dbl> <dbl>  <dbl>
    ## 1 enriched    7902  24732  1748  14594
    ## 2 unselected   113 107956  2995 126167

Sample prevalences

``` r
basics_sum %>% 
transmute(sample, reference, P=N_cases / (N_controls + N_cases)) %>% 
mutate(d50=abs(0.5-P)) %>%
arrange(sample, d50) %>%
print(n=Inf)
```

    ## # A tibble: 24 × 4
    ## # Groups:   sample [2]
    ##    sample     reference      P      d50
    ##    <chr>      <chr>      <dbl>    <dbl>
    ##  1 enriched   MDD4b     0.489  0.0107  
    ##  2 enriched   MDD5b     0.531  0.0311  
    ##  3 enriched   MDD5a     0.458  0.0422  
    ##  4 enriched   MDD3a     0.388  0.112   
    ##  5 enriched   MDD3b     0.375  0.125   
    ##  6 enriched   MDD9      0.654  0.154   
    ##  7 enriched   MDD4a     0.742  0.242   
    ##  8 enriched   MDD7      0.849  0.349   
    ##  9 enriched   MDD2      0.898  0.398   
    ## 10 enriched   MDD6      0.906  0.406   
    ## 11 enriched   MDD8      0.909  0.409   
    ## 12 enriched   MDD1      0.925  0.425   
    ## 13 unselected MDD4b     0.501  0.000871
    ## 14 unselected MDD3a     0.519  0.0195  
    ## 15 unselected MDD1      0.520  0.0204  
    ## 16 unselected MDD9      0.444  0.0562  
    ## 17 unselected MDD7      0.586  0.0863  
    ## 18 unselected MDD2      0.391  0.109   
    ## 19 unselected MDD3b     0.383  0.117   
    ## 20 unselected MDD8      0.763  0.263   
    ## 21 unselected MDD4a     0.787  0.287   
    ## 22 unselected MDD6      0.836  0.336   
    ## 23 unselected MDD5b     0.0908 0.409   
    ## 24 unselected MDD5a     0.0343 0.466

Datasets

``` r
meta_aligned <-
basics_meta %>%
    mutate(cohort=str_match(Dataset, "([A-Z_]+)")[,2],
           reference=str_match(Dataset, "(MDD[0-9a-b]+)")[,2])
           
write_csv(meta_aligned, "meta/meta.csv")
```
