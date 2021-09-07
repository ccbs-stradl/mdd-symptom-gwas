#$ -N popneuroveg_gwas
#$ -t 1-24
#$ -l h_rt=10:00:00
#$ -l h_vmem=2G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript popneuroveg_gwas.R