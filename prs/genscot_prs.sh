#$ -N gs_prs
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 8
#$ -cwd
#$ -e logs
#$ -o logs

# calculate PRS

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
--binary-target T \
--target genscot/genetics/imputed/HRC/updated_bims/GS20K_HRC_0.8_GCTA_# \
--maf 0.05 \
--pheno-file scid.pheno.txt \
--pheno-col MDD \
--cov-file genscot/genetics/genotypes/GS20K_PLINK_files/PCA_MDS_components/HM3mds.mds \
--cov-col @C[1-6] \
--prevalence 0.15 \
--keep genscot/genetics/genotypes/GS20K_PLINK_files/QCd_data/unrelated_individuals/keep_unrelated_non-italian_non-outlier_non-smr.txt \
--nonfounders \
--type bed \
--out $OUT \
--thread 8