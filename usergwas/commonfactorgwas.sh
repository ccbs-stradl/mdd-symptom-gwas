#$ -t 1-24
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript commonfactorgwas.R $1

## submit as
# qsub -N all_affect commonfactorgwas.sh all.affect
# qsub -N all_common commonfactorgwas.sh all.common
# qsub -N all_neuro commonfactorgwas.sh all.neuroveg
