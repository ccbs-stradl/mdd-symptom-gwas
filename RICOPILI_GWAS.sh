#bash

#Install and unpack ricopili - follow this link for custom installation https://docs.google.com/document/d/14aa-oeT5hF541I8hHsDAL_42oyvlHRC5FWR7gir4xco/edit#heading=h.y8igfs7neh22
wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.2019_Jun_25.001.tar.gz

tar -xvzf rp_bin.2019_Jun_25.001.tar.gz

#Create custom installation file in ricopili

#Softlink to all files within the testing datasets - actual path replaced by PATH 
ln -s //home/PATH/* /home/USER/ 

#TEST!!!  Postimp RICOPILI code - run this from my home directory now all softlinks have been set-up
postimp_navi --out MDD1_case_con --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --refiex refiex.mddw2v01

#Create datasets_info with information on the files you are interested in 
touch datasets_info
vim datasets_info
mdd_cof3_eur_sr-qc.hg19.ch.fl
mdd_col3_eur_sr-qc.hg19.ch.fl
mdd_rot4_eur_sr-qc.hg19.ch.fl
mdd_nes1_eur_sr-qc2.hg19.ch.fl
mdd_shp0_eur_sr-qc.hg19.ch.fl

#Run analysis
postimp_navi --out MDD1_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD1_PGCMDD2_final.txt #Done
postimp_navi --out MDD2_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD2_PGCMDD2_final.txt #Done
postimp_navi --out MDD3a_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3a_PGCMDD2_final.txt #Done
postimp_navi --out MDD3b_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3b_PGCMDD2_final.txt
postimp_navi --out MDD4a_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4a_PGCMDD2_final.txt
postimp_navi --out MDD4b_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4b_PGCMDD2_final.txt
postimp_navi --out MDD5a_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5a_PGCMDD2_final.txt
postimp_navi --out MDD5b_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5b_PGCMDD2_final.txt
postimp_navi --out MDD6_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD6_PGCMDD2_final.txt
postimp_navi --out MDD7_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD7_PGCMDD2_final.txt
postimp_navi --out MDD8_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD8_PGCMDD2_final.txt
postimp_navi --out MDD0_CaseControl --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD9_PGCMDD2_final.txt

# find names of datasets with case-only symptom data
# parse ID1s and extract name of cohort
cat MDD1_CasesOnly_PGCMDD2_final.tsv | awk '$3 != -9 {print $1}' | awk -F_ '{print $3}' | sort | uniq
cat > datasets_info_mdd1 << EOF
mdd_boma_eur_sr-qc.hg19.ch.fl
mdd_gep3_eur_sr-qc.hg19.ch.fl
mdd_gens_eur_sr-qc.hg19.ch.fl
mdd_grnd_eur_sr-qc.hg19.ch.fl
mdd_gsk2_eur_sr-qc.hg19.ch.fl
mdd_mmi2_eur_sr-qc.hg19.ch.fl
mdd_mmo4_eur_sr-qc.hg19.ch.fl
mdd_qi3c_eur_sr-qc.hg19.ch.fl
mdd_qi6c_eur_sr-qc.hg19.ch.fl
mdd_qio2_eur_sr-qc.hg19.ch.fl
mdd_rad3_eur_sr-qc.hg19.ch.fl
mdd_rage_eur_sr-qc.hg19.ch.fl
mdd_rai2_eur_sr-qc.hg19.ch.fl
mdd_rau2_eur_sr-qc.hg19.ch.fl
mdd_rde4_eur_sr-qc2.hg19.ch.fl
mdd_stm2_eur_sr-qc.hg19.ch.fl
mdd_twg2_eur_sr-qc.hg19.ch.fl
EOF

cat MDD3b_CasesOnly_PGCMDD2_final.tsv | awk '$3 != -9 {print $1}' | awk -F_ '{print $3}' | sort | uniq
cat > datasets_info_mdd3b << EOF
mdd_boma_eur_sr-qc.hg19.ch.fl
mdd_gep3_eur_sr-qc.hg19.ch.fl
mdd_gens_eur_sr-qc.hg19.ch.fl
mdd_grnd_eur_sr-qc.hg19.ch.fl
mdd_gsk2_eur_sr-qc.hg19.ch.fl
mdd_qi3c_eur_sr-qc.hg19.ch.fl
mdd_qi6c_eur_sr-qc.hg19.ch.fl
mdd_qio2_eur_sr-qc.hg19.ch.fl
mdd_rad3_eur_sr-qc.hg19.ch.fl
mdd_rage_eur_sr-qc.hg19.ch.fl
mdd_rai2_eur_sr-qc.hg19.ch.fl
mdd_rau2_eur_sr-qc.hg19.ch.fl
mdd_rde4_eur_sr-qc2.hg19.ch.fl
mdd_stm2_eur_sr-qc.hg19.ch.fl
mdd_twg2_eur_sr-qc.hg19.ch.fl
EOF

cat MDD4b_CasesOnly_PGCMDD2_final.tsv | awk '$3 != -9 {print $1}' | awk -F_ '{print $3}' | sort | uniq
cat > datasets_info_mdd4b << EOF
mdd_boma_eur_sr-qc.hg19.ch.fl
mdd_gep3_eur_sr-qc.hg19.ch.fl
mdd_gens_eur_sr-qc.hg19.ch.fl
mdd_grnd_eur_sr-qc.hg19.ch.fl
mdd_gsk2_eur_sr-qc.hg19.ch.fl
mdd_qi3c_eur_sr-qc.hg19.ch.fl
mdd_qi6c_eur_sr-qc.hg19.ch.fl
mdd_qio2_eur_sr-qc.hg19.ch.fl
mdd_rad3_eur_sr-qc.hg19.ch.fl
mdd_rage_eur_sr-qc.hg19.ch.fl
mdd_rai2_eur_sr-qc.hg19.ch.fl
mdd_rau2_eur_sr-qc.hg19.ch.fl
mdd_rde4_eur_sr-qc2.hg19.ch.fl
mdd_stm2_eur_sr-qc.hg19.ch.fl
EOF

cat MDD5a_CasesOnly_PGCMDD2_final.tsv | awk '$3 != -9 {print $1}' | awk -F_ '{print $3}' | sort | uniq
cat > datasets_info_mdd5a << EOF
mdd_boma_eur_sr-qc.hg19.ch.fl
mdd_gens_eur_sr-qc.hg19.ch.fl
mdd_gep3_eur_sr-qc.hg19.ch.fl
mdd_grnd_eur_sr-qc.hg19.ch.fl
mdd_mmi2_eur_sr-qc.hg19.ch.fl
mdd_mmo4_eur_sr-qc.hg19.ch.fl
mdd_qi3c_eur_sr-qc.hg19.ch.fl
mdd_qi6c_eur_sr-qc.hg19.ch.fl
mdd_qio2_eur_sr-qc.hg19.ch.fl
mdd_rad3_eur_sr-qc.hg19.ch.fl
mdd_rage_eur_sr-qc.hg19.ch.fl
mdd_rai2_eur_sr-qc.hg19.ch.fl
mdd_rau2_eur_sr-qc.hg19.ch.fl
mdd_rde4_eur_sr-qc2.hg19.ch.fl
mdd_stm2_eur_sr-qc.hg19.ch.fl
EOF

# run cases-only analyses
# construct postimp_navi commands
# extract "MDDn" from CasesOnly phenotype file names and use it to parameterize the postimputation pipeline
# for mdd in $(ls MDD*_CasesOnly_PGCMDD2_final.tsv | awk -F_ '{print $1}'); do echo postimp_navi --out ${mdd}_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno ${mdd}_CasesOnly_PGCMDD2_final.tsv; done
cp datasets_info_mdd1 datasets_info
postimp_navi --out MDD1_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD1_CasesOnly_PGCMDD2_final.tsv

postimp_navi --out MDD2_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD2_CasesOnly_PGCMDD2_final.tsv
postimp_navi --out MDD3a_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3a_CasesOnly_PGCMDD2_final.tsv


cp datasets_info_mdd3b datasets_info
postimp_navi --out MDD3b_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD3b_CasesOnly_PGCMDD2_final.tsv

cp datasets_info_mdd1 datasets_info
postimp_navi --out MDD4a_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4a_CasesOnly_PGCMDD2_final.tsv

cp datasets_info_mdd4b datasets_info
postimp_navi --out MDD4b_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD4b_CasesOnly_PGCMDD2_final.tsv

cp datasets_info_mdd5a datasets_info
postimp_navi --out MDD5a_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5a_CasesOnly_PGCMDD2_final.tsv

postimp_navi --out MDD5b_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD5b_CasesOnly_PGCMDD2_final.tsv

cp datasets_info_mdd1 datasets_info
postimp_navi --out MDD6_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD6_CasesOnly_PGCMDD2_final.tsv

postimp_navi --out MDD7_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD7_CasesOnly_PGCMDD2_final.tsv
postimp_navi --out MDD8_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD8_CasesOnly_PGCMDD2_final.tsv
postimp_navi --out MDD9_CasesOnly --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --popname eur --pheno MDD9_CasesOnly_PGCMDD2_final.tsv


