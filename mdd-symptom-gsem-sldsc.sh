#$ -N gsem_sldsc
#$ -l h_vmem=8G
#$ -pe sharedmem 12
#$ -l h_rt=12:00:00
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript -e "rmarkdown::render('mdd-symptom-gsem-sldsc.Rmd')"
