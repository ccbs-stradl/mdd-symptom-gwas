#$ -N gs_prs_all
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 8
#$ -cwd
#$ -e logs
#$ -o logs

# calculate PRS for all samples

BASE=$1
OUT=$2

# Initialise the environment modules
. /etc/profile.d/modules.sh

# load PRSice
module load igmm/apps/PRSice/2.1.11 

Rscript $(which PRSice.R) \
--prsice $(which PRSice_linux) \
--base $BASE \
--A1 A1 \
--A2 A2 \
--bp BP \
--chr CHR \
--maf-base MAF,0.05 \
--beta \
--pvalue Pval_Estimate \
--snp SNP \
--stat est \
--se se_c \
--bar-levels 0.025,0.05,0.25 \
--fastscore \
--no-regress \
--nonfounders \
--target genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA_# \
--maf 0.05 \
--type bed \
--out ${OUT}_all \
--thread 8