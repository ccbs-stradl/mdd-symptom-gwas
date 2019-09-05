# Munge sumstats for UKB CIDI

for sumstats in $(ls sumstats/UKB/CIDI/UKB_CIDI_MDD*.gz); do

        prefix=$(basename $sumstats .gz)

        munge_sumstats.py \
        --sumstats $sumstats \
        --N-cas-col Nca \
        --N-con-col Nco \
        --signed-sumstats OR,1 \
        --p P \
        --merge-alleles sumstats/reference/w_hm3.snplist \
        --out sumstats/UKB/CIDI/${prefix}.ldsc

done
