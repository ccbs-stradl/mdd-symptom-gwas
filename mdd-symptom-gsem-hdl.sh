#! /bin/bash
#$ -N gsem_hdl
#$ -cwd
#$ -l h_vmem=8G
#$ -l h_rt=08:00:00
#$ -e logs
#$ -o logs

/exports/igmm/eddie/GenScotDepression/madams23/local/anaconda/envs/gsem/bin/Rscript --no-save --no-restore -e "rmarkdown::render('mdd-symptom-gsem-hdl.Rmd')"