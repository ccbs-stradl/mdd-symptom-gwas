# Munge sumstats for all cohorts symptom GWASs

for sumstats in $(ls sumstats/PGC/CasesAllCohorts/daner_MDD*.meta.gz); do

        prefix=$(basename $sumstats .gz)

        munge_sumstats.py --daner-n \
        --sumstats $sumstats \
        --merge-alleles sumstats/reference/w_hm3.snplist \
        --out sumstats/PGC/CasesAllCohorts/${prefix}.ldsc

done
