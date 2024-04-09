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
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): filename, cohort, reference
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Symptom descriptions

``` r
mdd_symptoms <- read_tsv('dsm_mdd.tsv')
```

    ## Rows: 15 Columns: 7
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (7): Ref, Reference, h, v, abbv, Symptom, Description
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
left_join(basics_datasets, by='Dataset')

basics_summary <- basics_aligned |>
mutate(meta=case_when(cohort %in% c('PGC', 'AGDS', 'GenScot', 'janpy') ~ "Clinical",
                      cohort %in% c('ALSPAC', 'EstBB', 'UKBB') ~ "Community",
                      TRUE ~ NA_character_)) |>
group_by(meta, reference) |>
summarize(N_cases=sum(N_cases, na.rm=T), N_controls=sum(N_controls, na.rm=T), N_eff_half=sum(N_eff_half, na.rm=T)) |>
ungroup()
```

    ## `summarise()` has grouped output by 'meta'. You can override using the
    ## `.groups` argument.

``` r
presence_absence <- basics_summary |>
    str_glue_data("{N_eff_half} ({round(100*N_cases/(N_cases + N_controls))}%)")

basics_formatted <- 
basics_summary %>%
mutate(PresenceAbsence=as.character(presence_absence)) %>%
select(meta, reference, PresenceAbsence) %>%
pivot_wider(names_from=meta, values_from=PresenceAbsence) |>
left_join(mdd_symptoms, by = c('reference' = 'Reference')) |>
mutate(Symptom = str_glue("{Ref}. {Symptom}")) |>
select(Symptom, Abbr. = abbv, Clinical, Community)

knitr::kable(basics_formatted)
```

| Symptom                                    | Abbr.   | Clinical    | Community    |
|:-------------------------------------------|:--------|:------------|:-------------|
| 1\. Depressed mood                         | Dep     | 2471 (93%)  | 102071 (52%) |
| 2\. Anhedonia                              | Anh     | 4494 (90%)  | 98458 (39%)  |
| 3a. Weight loss / decrease in appetite     | AppDec  | 10119 (39%) | 35837 (52%)  |
| 3b. Weight gain / increase in appetite     | AppInc  | 9259 (38%)  | 26681 (38%)  |
| 4a. Insomnia                               | SleDec  | 9418 (74%)  | 28741 (79%)  |
| 4b. Hypersomnia                            | SleInc  | 10031 (49%) | 18377 (50%)  |
| 5a. Psychomotor agitation                  | MotoInc | 10380 (46%) | 218 (3%)     |
| 5b. Psychomotor slowing                    | MotoDec | 11130 (53%) | 543 (9%)     |
| 6\. Fatigue                                | Fatig   | 3907 (91%)  | 25790 (84%)  |
| 7\. Feelings of worthlessness / guilt      | Guilt   | 5503 (85%)  | 46694 (59%)  |
| 8\. Diminished concentration               | Conc    | 3793 (91%)  | 32827 (76%)  |
| 9\. Recurrent thoughts of death or suicide | Sui     | 10545 (65%) | 50035 (44%)  |

``` r
basics_sum <-
basics_aligned %>%
mutate(sample=if_else(cohort %in% c('AGDS', 'PGC', 'GenScot', 'janpy'), true='enriched', false='unselected')) %>%
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
    ## 1 enriched   21681  21681  1748   1748
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
    ##    sample     reference       P       d50
    ##    <chr>      <chr>       <dbl>     <dbl>
    ##  1 enriched   MDD1       0.925   0.425   
    ##  2 enriched   MDD2      NA      NA       
    ##  3 enriched   MDD3a     NA      NA       
    ##  4 enriched   MDD3b     NA      NA       
    ##  5 enriched   MDD4a     NA      NA       
    ##  6 enriched   MDD4b     NA      NA       
    ##  7 enriched   MDD5a     NA      NA       
    ##  8 enriched   MDD5b     NA      NA       
    ##  9 enriched   MDD6      NA      NA       
    ## 10 enriched   MDD7      NA      NA       
    ## 11 enriched   MDD8      NA      NA       
    ## 12 enriched   MDD9      NA      NA       
    ## 13 unselected MDD4b      0.501   0.000871
    ## 14 unselected MDD3a      0.519   0.0195  
    ## 15 unselected MDD1       0.520   0.0204  
    ## 16 unselected MDD9       0.444   0.0562  
    ## 17 unselected MDD7       0.586   0.0863  
    ## 18 unselected MDD2       0.391   0.109   
    ## 19 unselected MDD3b      0.383   0.117   
    ## 20 unselected MDD8       0.763   0.263   
    ## 21 unselected MDD4a      0.787   0.287   
    ## 22 unselected MDD6       0.836   0.336   
    ## 23 unselected MDD5b      0.0908  0.409   
    ## 24 unselected MDD5a      0.0343  0.466

Datasets

``` r
cohort_sumstats <- basics_aligned |>
    transmute(cohort, reference, dataset=str_match(Dataset, "([A-Za-z0-9_]+)")[,2],
              N_cases, N_controls, `N_eff_half`, `LAMBDA-GC`, `N-SNPs`) |>
    filter(!is.na(N_cases))

meta_sumstats <-
basics_meta %>%
    mutate(cohort=str_match(Dataset, "([A-Za-z_]+)")[,2],
           reference=str_match(Dataset, "(MDD[0-9a-b]+)")[,2]) |>
    select(cohort, reference, dataset=Dataset,
           N_cases, N_controls, `N_eff_half`,
           `LAMBDA-GC`, `N-SNPs`)
           
sumstats <- bind_rows(cohort_sumstats, meta_sumstats)
           
write_csv(sumstats, "meta/meta.csv")
```

Totals

``` r
cohort_sumstats |> filter(cohort %in% c('PGC', 'AGDS', 'GenScot', 'janpy')) |>
    mutate(N=N_cases+N_controls, cohort=if_else(cohort=='PGC', true=sapply(str_split(dataset, "_"), last), false=cohort)) |>
    group_by(cohort) |>
    summarize(N=max(N)) |>
    ungroup() |>
    summarize(N=sum(N))
```

    ## # A tibble: 1 × 1
    ##       N
    ##   <dbl>
    ## 1 30148

``` r
cohort_sumstats |> filter(cohort %in% c('ALSPAC', 'EstBB', 'UKBB')) |>
    mutate(N=N_cases+N_controls) |>
    group_by(cohort) |>
    summarize(N=max(N)) |>
    ungroup() |>
    summarize(N=sum(N))
```

    ## # A tibble: 1 × 1
    ##        N
    ##    <dbl>
    ## 1 207436
