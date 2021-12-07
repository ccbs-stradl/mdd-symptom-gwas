#$ -N clinsoma_gwas
#$ -t 1-24
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript clinsoma_gwas.R