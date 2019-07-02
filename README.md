# GWAS and GenomicSEM of MDD symptoms

## Prerequisites

- [R](https://r-project.org)
- [GenomicSEM](https://github.com/MichelNivard/GenomicSEM)
- [corrplot](https://cran.r-project.org/package=corrplot)
- [rmarkdown](https://rmarkdown.rstudio.com)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)

## Configuration

System configuration is specified in a [YAML](https://yaml.org/) file. Make a copy of [`config-example.yaml`](config-example.yaml) to `config.yaml` and set the options up for the system each part of the analyis is run on.

## GWAS sumstats

Sumstats are versioned in a separate internal repository using LFS. This
project is based on intermediate representations of genomic covariance
objects that can be shared and worked on as small text files.

Eventually code to generate `ldsc()` objects will be included here starting
from publicly downloadable sumstats files.

## Running GenomicSEM analyses

```r

rmarkdown::render('commonfactor.Rmd', output_format='all')


```
