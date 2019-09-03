# meta analyze UKB CIDI and PGC DSM GWAS
mkdir meta

for symptom in MDD1 MDD2 MDD3a MDD3b MDD4a MDD4b MDD6 MDD7 MDD8 MDD9; do
  # copy symptom file
  cp PGC/CasesAllCohorts/daner_${symptom}_*.meta.gz pgc_dsm.gz
  cp UKB/CIDI/UKB_CIDI_${symptom}.eur.bgenie.OR.tsv.gz ukb_cidi.gz

  metal ukb_pgc_meta.metal

  cat pgc_ukb_meta_1.tbl | gzip -c > meta/${symptom}.PGC_DSM.UKB_CIDI.meta.gz
  cp pgc_ukb_meta_1.tbl.info   meta/${symptom}.PGC_DSM.UKB_CIDI.meta.info

  rm pgc_ukb_meta_1.tbl pgc_ukb_meta_1.tbl.info pgc_dsm.gz ukb_cidi.gz

done

