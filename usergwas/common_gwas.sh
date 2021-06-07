#$ -N common_gwas
# -t 1
#$ -l h_rt=6:00:00
#$ -l h_vmem=5G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript common_factorgwas.R