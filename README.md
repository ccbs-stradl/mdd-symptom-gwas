# Genetic structure of major depression symptoms across clinical and community cohorts

- [GWAS of PGC cohorts](pgc-symptom-gwas.md)
- [Sample characteristics](mdd-symptom-samples.md)
- [Genetic covariance estimation](mdd-symptom-gsem.md)
- [Factor models](mdd-symptom-gsem-model.md)
- [Comparisons with other phenotypes](ext/mdd-symptom-gsem-ext.md)

## Prerequisites

- [R](https://r-project.org)
- [GenomicSEM](https://github.com/MichelNivard/GenomicSEM)
- [corrplot](https://cran.r-project.org/package=corrplot)
- [rmarkdown](https://rmarkdown.rstudio.com)
- [yaml](https://cran.r-project.org/web/packages/yaml/index.html)

## Configuration

System configuration is specified in a [YAML](https://yaml.org/) file. Make a copy of [`config-example.yaml`](config-example.yaml) to `config.yaml` and set the options up for the system each part of the analysis is run on (mostly used for GWAS of PGC cohorts on LISA).

## GWAS sumstats

Sumstats are versioned in a separate internal repository. This
project is based on intermediate representations of genomic covariance
objects that can be shared and worked on as small text files.

