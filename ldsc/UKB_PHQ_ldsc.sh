# Munge sumstats for UKB PHQ (INT)

for sumstats in sumstats/UKB/PHQ9_INT/UKB_PHQ9_INT_MDD*bgenie..gz; do

        prefix=$(basename $sumstats .gz)

        munge_sumstats.py \
        --sumstats $sumstats \
        --a1 a_1 \
        --a2 a_0 \
        --info info \
        --frq af \
        --N-col N \
        --signed-sumstats beta,0 \
        --p p \
        --merge-alleles sumstats/reference/w_hm3.snplist \
        --out sumstats/UKB/PHQ9_INT/${prefix}.ldsc

done

for sumstats in sumstats/UKB/PHQ9_INT/UKB_PHQ9_INT_MDD*.ldsc.sumstats.gz; do
        OUT=$(dirname $sumstats)/$(basename $sumstats .gz).h2

        ldsc.py \
        --h2 $sumstats \
        --out $OUT \
        --ref-ld-chr sumstats/reference/eur_w_ld_chr/ \
        --w-ld-chr sumstats/reference/eur_w_ld_chr/

done

