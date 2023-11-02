#$ -t 1-100
#$ -l h_rt=4:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

Rscript commonfactorgwas.R $1

## submit as
# qsub -N all_common commonfactorgwas.sh comm.ukb.common
# Rscript commonfactorgwas_merge.R comm.ukb.common
