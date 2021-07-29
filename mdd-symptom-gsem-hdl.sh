#! /bin/bash
#$ -N gsem_hdl
#$ -cwd
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -l h_rt=24:00:00
#$ -e logs
#$ -o logs

/exports/igmm/eddie/GenScotDepression/madams23/local/anaconda/envs/gsem/bin/Rscript --no-save --no-restore mdd-symptom-gsem-hdl.R
