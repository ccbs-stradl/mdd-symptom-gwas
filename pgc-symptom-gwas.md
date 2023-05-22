---
title: PGC Symptom GWAS 
author: Bradley Jermy, Mark Adams 
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

Run GWAS of MDD symptoms on the [LISA server](http://geneticcluster.org/). Because these commands need to be run on a specific cluster with specific group permissions, this notebook documents the commands that were run rather than reproducing the whole analysis internally, with `eval=FALSE`.

# Directory setup

Softlink to all MDD files. Replace `PATH` with the actual path to the MDD Wave2 directory (listed in `README.mdd2sum`):


```bash
mkdir mdd_symptoms

cd mdd_symptoms

ln -s /home/PATH/* .
```

# Install [Ricopili](https://sites.google.com/a/broadinstitute.org/ricopili/)


```bash
#Install and unpack ricopili - follow this link for custom installation https://docs.google.com/document/d/14aa-oeT5hF541I8hHsDAL_42oyvlHRC5FWR7gir4xco/edit#heading=h.y8igfs7neh22
wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.2019_Jun_25.001.tar.gz

tar -xvzf rp_bin.2019_Jun_25.001.tar.gz

#Create custom installation file in ricopili

#TEST!!!  Postimp RICOPILI code - run this home directory after all softlinks have been set-up
postimp_navi --out MDD1_case_con --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --refiex refiex.mddw2v01
```

# Symptom phenotype files

Launch an interactive session


```bash

#Go into interactive mode
srun -n 16 -t 1:00:00 --pty bash -il

#Load R

module load 2022
module load R/4.2.1-foss-2022a

R

```

Load libraries and configure input directories


```r
library(yaml)
library(readxl)
library(readr)
library(stringr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(tidyr)

# load configuration file (see config-example.yaml')
config <- yaml.load_file('config.yaml')

# input and output directories and files
mdd_v1_dir <- config$data$pgc$mdd$v1
mdd_v2_dir <- config$data$pgc$mdd$v2

output_dir = "pgc_gwas"
```

## Phenotype files 

Extract 12 MDD symptoms


```r
#Extract 12 symptoms
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
studies <- c("RADIANT", "STR", "ROTTERDAM", "QIMR", "NESDANTR", "GSK", "GENRED", "GENRED2", "CoLaus", "COFAMS", "BonnMH")
studies_v2 <- c("BiDirect")
study_noupdate <- c("STARD", "SHIP0", "MARS", "GENPOD")

dir.create(output_dir, recursive=TRUE)
```

```
## Warning in dir.create(output_dir, recursive = TRUE): 'pgc_gwas' already exists
```

```r
pheno_xls <- c(
  list.files(file.path(mdd_v1_dir, "secondary_phenotypes", "1216_update"), ".xls", full.names=TRUE),
  list.files(file.path(mdd_v2_dir, "phenotype"), ".xls", full.names=TRUE)
)

names(pheno_xls) <- basename(pheno_xls)

phenotypes <- 
bind_rows(
lapply(pheno_xls, 
  function(xls) read_excel(xls) |>
                       mutate(across(ID2:Sub.study, .fns=as.character)) |>
                       select(ID1, ID2, Study, Sub.study, MDD1:MDD9) |> 
                       mutate(across(MDD1:MDD9, as.numeric))
),
.id="filename"
) |>
separate(filename, into=c('file', 'update'))
```

```
## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion
```

```
## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion

## Warning in mask$eval_all_mutate(quo): NAs introduced by coercion
```

```
## Warning: Expecting numeric in B3615 / R3615C2: got 'E00680001'
```

```
## Warning: Expecting numeric in B3616 / R3616C2: got 'E08270150'
```

```
## Warning: Expecting numeric in B3617 / R3617C2: got 'E82170001'
```

```
## Warning: Expecting numeric in B3618 / R3618C2: got 'E82230001'
```

```
## Warning: Expecting numeric in B3619 / R3619C2: got 'E82720001'
```

```
## Warning: Expecting numeric in B3620 / R3620C2: got 'E83000001'
```

```
## Warning: Expecting numeric in B3621 / R3621C2: got 'E84110001'
```

```
## Warning: Expecting numeric in B3622 / R3622C2: got 'E88160001'
```

```
## Warning: Expecting numeric in B3623 / R3623C2: got 'M407915'
```

```
## Warning: Expecting numeric in B3624 / R3624C2: got 'M415462'
```

```
## Warning: Expecting numeric in B3510 / R3510C2: got 'GDP0202'
```

```
## Warning: Expecting numeric in B3511 / R3511C2: got 'GDP0203'
```

```
## Warning: Expecting numeric in B3512 / R3512C2: got 'GDP0205'
```

```
## Warning: Expecting numeric in B3513 / R3513C2: got 'GDP0201'
```

```
## Warning: Expecting numeric in B3514 / R3514C2: got 'GDP0207'
```

```
## Warning: Expecting numeric in B3515 / R3515C2: got 'GDP5839'
```

```
## Warning: Expecting numeric in B3516 / R3516C2: got 'GDP5840'
```

```
## Warning: Expecting numeric in B3517 / R3517C2: got 'GDP5838'
```

```
## Warning: Expecting numeric in B3518 / R3518C2: got 'GDP5860'
```

```
## Warning: Expecting numeric in B3519 / R3519C2: got 'GDP5859'
```

```
## Warning: Expecting numeric in B3520 / R3520C2: got 'GDP5858'
```

```
## Warning: Expecting numeric in B3521 / R3521C2: got 'GDP5857'
```

```
## Warning: Expecting numeric in B3522 / R3522C2: got 'GDP5841'
```

```
## Warning: Expecting numeric in B3523 / R3523C2: got 'GDP5844'
```

```
## Warning: Expecting numeric in B3524 / R3524C2: got 'GDP5845'
```

```
## Warning: Expecting numeric in B3525 / R3525C2: got 'GDP5846'
```

```
## Warning: Expecting numeric in B3526 / R3526C2: got 'GDP5848'
```

```
## Warning: Expecting numeric in B3527 / R3527C2: got 'GDP5856'
```

```
## Warning: Expecting numeric in B3528 / R3528C2: got 'GDP5855'
```

```
## Warning: Expecting numeric in B3529 / R3529C2: got 'GDP5852'
```

```
## Warning: Expecting numeric in B3530 / R3530C2: got 'GDP5851'
```

```
## Warning: Expecting numeric in B3531 / R3531C2: got 'GDP5868'
```

```
## Warning: Expecting numeric in B3532 / R3532C2: got 'GDP5870'
```

```
## Warning: Expecting numeric in B3533 / R3533C2: got 'GDP5871'
```

```
## Warning: Expecting numeric in B3534 / R3534C2: got 'GDP5872'
```

```
## Warning: Expecting numeric in B3535 / R3535C2: got 'GDP5873'
```

```
## Warning: Expecting numeric in B3536 / R3536C2: got 'GDP5888'
```

```
## Warning: Expecting numeric in B3537 / R3537C2: got 'GDP5890'
```

```
## Warning: Expecting numeric in B3538 / R3538C2: got 'GDP5891'
```

```
## Warning: Expecting numeric in B3539 / R3539C2: got 'GDP5892'
```

```
## Warning: Expecting numeric in B3540 / R3540C2: got 'GDP5894'
```

```
## Warning: Expecting numeric in B3541 / R3541C2: got 'GDP5895'
```

```
## Warning: Expecting numeric in B3542 / R3542C2: got 'GDP5896'
```

```
## Warning: Expecting numeric in B3543 / R3543C2: got 'GDP5828'
```

```
## Warning: Expecting numeric in B3544 / R3544C2: got 'GDP5829'
```

```
## Warning: Expecting numeric in B3545 / R3545C2: got 'GDP5830'
```

```
## Warning: Expecting numeric in B3546 / R3546C2: got 'GDP5831'
```

```
## Warning: Expecting numeric in B3547 / R3547C2: got 'GDP5832'
```

```
## Warning: Expecting numeric in B3548 / R3548C2: got 'GDP5834'
```

```
## Warning: Expecting numeric in B3549 / R3549C2: got 'GDP5815'
```

```
## Warning: Expecting numeric in B3550 / R3550C2: got 'GDP5816'
```

```
## Warning: Expecting numeric in B3551 / R3551C2: got 'GDP5813'
```

```
## Warning: Expecting numeric in B3552 / R3552C2: got 'GDP5817'
```

```
## Warning: Expecting numeric in B3553 / R3553C2: got 'GDP5818'
```

```
## Warning: Expecting numeric in B3554 / R3554C2: got 'GDP5819'
```

```
## Warning: Expecting numeric in B3555 / R3555C2: got 'GDP5820'
```

```
## Warning: Expecting numeric in B3556 / R3556C2: got 'GDP5821'
```

```
## Warning: Expecting numeric in B3557 / R3557C2: got 'GDP5822'
```

```
## Warning: Expecting numeric in B3558 / R3558C2: got 'GDP5823'
```

```
## Warning: Expecting numeric in B3559 / R3559C2: got 'GDP5824'
```

```
## Warning: Expecting numeric in B3560 / R3560C2: got 'GDP5811'
```

```
## Warning: Expecting numeric in B3561 / R3561C2: got 'GDP5812'
```

```
## Warning: Expecting numeric in B3562 / R3562C2: got 'GDP5814'
```

```
## Warning: Expecting numeric in B3563 / R3563C2: got 'GDP5836'
```

```
## Warning: Expecting numeric in B3564 / R3564C2: got 'GDP5854'
```

```
## Warning: Expecting numeric in B3565 / R3565C2: got 'GDP5900'
```

```
## Warning: Expecting numeric in B3566 / R3566C2: got 'GDP5853'
```

```
## Warning: Expecting numeric in B3567 / R3567C2: got 'GDP5850'
```

```
## Warning: Expecting numeric in B3568 / R3568C2: got 'GDP5842'
```

```
## Warning: Expecting numeric in B3569 / R3569C2: got 'GDP5881'
```

```
## Warning: Expecting numeric in B3570 / R3570C2: got 'GDP5882'
```

```
## Warning: Expecting numeric in B3571 / R3571C2: got 'GDP5883'
```

```
## Warning: Expecting numeric in B3572 / R3572C2: got 'GDP5885'
```

```
## Warning: Expecting numeric in B3573 / R3573C2: got 'GDP5887'
```

```
## Warning: Expecting numeric in B3574 / R3574C2: got 'GDP5886'
```

```
## Warning: Expecting numeric in B3575 / R3575C2: got 'GDP5874'
```

```
## Warning: Expecting numeric in B3576 / R3576C2: got 'GDP5875'
```

```
## Warning: Expecting numeric in B3577 / R3577C2: got 'GDP5876'
```

```
## Warning: Expecting numeric in B3578 / R3578C2: got 'GDP5867'
```

```
## Warning: Expecting numeric in B3579 / R3579C2: got 'GDP5877'
```

```
## Warning: Expecting numeric in B3580 / R3580C2: got 'GDP6408'
```

```
## Warning: Expecting numeric in B3581 / R3581C2: got 'GDP5878'
```

```
## Warning: Expecting numeric in B3582 / R3582C2: got 'GDP5879'
```

```
## Warning: Expected 2 pieces. Additional pieces discarded in 30185 rows [1, 2, 3,
## 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
```

## Code phenotypes


```r
pheno_au <- 
phenotypes |> mutate(across(MDD1:MDD9, ~ case_when(.x == 1 ~ 2L,
                                             .x == 0 ~ 1L,
                                             TRUE ~ NA_integer_))) |>
mutate(dataset = sapply(str_split(ID1, pattern = "_", n = 4), function(id) id[3]))

write_tsv(pheno_au, file.path(output_dir, "mddw2_symptoms.pheno"), na="-9")
```

## Covariates


```r
# extract covariates
v1_cov <- read_table(file.path(mdd_v1_dir, "MDD29.0515.nproj.menv.mds_cov"))
```

```
## 
## -- Column specification ---------------------------------------------------------------------------------------------------------------------------
## cols(
##   .default = col_double(),
##   FID = col_character(),
##   IID = col_character()
## )
## i Use `spec()` for the full column specifications.
```

```r
v2_cov <- read_table(file.path(mdd_v2_dir, "cobg_gw.Bidirect.menv.mds_cov"))
```

```
## 
## -- Column specification ---------------------------------------------------------------------------------------------------------------------------
## cols(
##   .default = col_double(),
##   FID = col_character(),
##   IID = col_character()
## )
## i Use `spec()` for the full column specifications.
```

```r
covar <- bind_rows(v1_cov, v2_cov)

write_tsv(covar, file.path(output_dir, "mddw2.covar"))
```

## Cases cohorts

Identify cohorts with cases that vary in symptoms.


```r
cases <- pheno_au |>
filter(str_detect(ID1, "^cas_mdd"))

cohorts_with_case_symptoms <- 
cases |>
pivot_longer(cols=MDD1:MDD9, names_to = "symptom", values_to = "status") |>
filter(status %in% 1:2) |>
distinct(dataset, status) |>
count(dataset) |>
filter(n == 2)

cases_ids <- cases |>
filter(dataset %in% pull(cohorts_with_case_symptoms, dataset)) |>
select(FID=ID1, IID=ID2)

write_tsv(cases_ids, file.path(output_dir, "cases.id"))
```

## UKB Overlap

Remove overlap with UKB


```r
ukb <- read_table("pgc_gwas/959_PGC_UKB_overlap.id", col_names=c("FID", "IID"))
```

```
## 
## -- Column specification ---------------------------------------------------------------------------------------------------------------------------
## cols(
##   FID = col_character(),
##   IID = col_character()
## )
```

## QC phenotypes

Set phenotypes to missing if there are fewer than 80 (10x number of GWAS predictors)) symptom present/absent cases.


```r
pheno_long <-
pheno_au |>
filter(ID1 %in% pull(cases_ids, FID)) |>
filter(!ID1 %in% pull(ukb, FID)) |>
pivot_longer(cols=MDD1:MDD9, names_to = "symptom", values_to = "status")

keep_symptoms <- 
pheno_long |>
filter(!is.na(status)) |>
group_by(dataset, symptom, status) |>
tally() |>
group_by(dataset, symptom) |>
summarize(samples=min(n), responses=n()) |>
filter(samples >= 80, responses == 2)
```

```
## `summarise()` has grouped output by 'dataset'. You can override using the
## `.groups` argument.
```

```r
keep_cohorts <- keep_symptoms |>
distinct(dataset)

pheno <-
pheno_long |>
inner_join(keep_symptoms, by=c("dataset", "symptom")) |>
select(-samples, responses) |>
arrange(symptom) |>
pivot_wider(names_from=symptom, values_from=status) |>
arrange(ID1) |>
select(FID=ID1, IID=ID2, starts_with("MDD"))

write_tsv(pheno, file.path(output_dir, "mddw2_symptoms_cases.pheno"), na="-9")
write_tsv(keep_cohorts, file.path(output_dir, "cohorts.keep"), col_names=F)
```

# GWAS Pipeline

Run GWAS across all symptoms for just cases


```sh

# symlink hardcall filenames

mkdir -p pgc_gwas/cobg_dir_genome_wide
while read cohort; do
  if [ -f $V1/cobg_dir_genome_wide/mdd_${cohort}_eur_sr-qc*.hg19.ch.fl.bgn.bed ]; then
    echo Linked $cohort
    ln -sf $V1/cobg_dir_genome_wide/mdd_${cohort}_*.bgn.{bed,bim,fam} pgc_gwas/cobg_dir_genome_wide/
 elif [ -f $V1/genred_case_only/mdd_${cohort}_eur_sr-qc*.hg19.ch.fl.bgn.cases.bed ]; then
 	echo Linked $cohort
 	ln -sf $V1/genred_case_only/mdd_${cohort}_*.bgn.cases.{bed,bim,fam} pgc_gwas/cobg_dir_genome_wide/
  elif [ -f $V2/cobg_dir_genome_wide/mdd_${cohort}_eur_sr-qc2.hg19.ch.fl.bgn.bed ]; then
    echo Linked $cohort
    ln -sf $V2/cobg_dir_genome_wide/mdd_${cohort}_*.bgn.{bed,bim,fam} pgc_gwas/cobg_dir_genome_wide/
  else
    echo Did not link $cohort
  fi
done<pgc_gwas/cohorts.keep

# patch in rau2 (raus)
for cohort in rau2; do
  if [ -f $V1/cobg_dir_genome_wide/mdd_${cohort}_eur_sr-qc*.hg19.ch.fl.bgn.bed ]; then
    echo Linked $cohort
    ln -sf $V1/cobg_dir_genome_wide/mdd_${cohort}_*.bgn.{bed,bim,fam} pgc_gwas/cobg_dir_genome_wide/
  else
    echo Did not link $cohort
  fi
done
```


```sh

nextflow run pgc-symptom-gwas.nf \
-config lisa.config \
-resume \
--bfile="pgc_gwas/cobg_dir_genome_wide/mdd_*.bgn*{bed,bim,fam}"

```
