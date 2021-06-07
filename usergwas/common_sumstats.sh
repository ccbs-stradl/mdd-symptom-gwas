#$ -N common_sumstats
#$ -l h_rt=2:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -o logs
#$ -e logs
#$ -cwd

. /etc/profile.d/modules.sh
module load igmm/apps/R/4.1.0

Rscript common_sumstats.R