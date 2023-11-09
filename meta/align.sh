
module load 2022
module load Java/13.0.2

ref_dir=$(head -n 1 reference_info)

nextflow run align.nf \
-config ../snellius.config \
-work-dir /scratch-shared/$USER/mdd-symptoms-gwas/align/work \
-resume \
--frq "${ref_dir}/*.EUR.frq2.gz" \
--daner "../sumstats/PGC/cohorts/daner_MDD*.gz"

nextflow run align.nf \
-config ../snellius.config \
-work-dir /scratch-shared/$USER/mdd-symptoms-gwas/align/work \
-resume \
--frq "${ref_dir}/*.EUR.frq2.gz" \
--daner "../sumstats/JRDUS/*.meta.gz"