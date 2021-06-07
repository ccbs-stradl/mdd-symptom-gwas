#$ gsem_sldsc
#$ -l h_vmem=8G
#$ -pe sharedmem 6
#$ -l h_rt=6:00:00
#$ -o logs
#$ -e logs

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript -e "rmarkdown::render('mdd-symptom-gsem-sldsc.Rmd')"