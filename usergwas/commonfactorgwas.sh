#$ -t 1-100
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript commonfactorgwas.R $1

## submit as
# qsub -N all_common commonfactorgwas.sh all.common
# Rscript commonfactorgwas_merge.R all.common
