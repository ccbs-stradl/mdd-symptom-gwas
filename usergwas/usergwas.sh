#$ -t 1-400
#$ -l h_rt=4:00:00
#$ -l h_vmem=6G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript $1

## submit as
# qsub -N mel_aty usergwas.sh usergwas_mel_aty.R
# Rscript commonfactorgwas_merge.R all.mel_aty_gate
