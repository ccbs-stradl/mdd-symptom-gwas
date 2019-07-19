#bash

#Install and unpack ricopili - follow this link for custom installation https://docs.google.com/document/d/14aa-oeT5hF541I8hHsDAL_42oyvlHRC5FWR7gir4xco/edit#heading=h.y8igfs7neh22
wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.2019_Jun_25.001.tar.gz

tar -xvzf rp_bin.2019_Jun_25.001.tar.gz

#Create custom installation file in ricopili


#Softlink to testing datasets 
ln -s //home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/* /home/jermy/


#TEST!!!  Postimp RICOPILI code - run this from my home directory now all softlinks have been set-up
postimp_navi --out MDD1_case_con --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --refiex refiex.mddw2v01

#Create datasets_info with information on the files you are interested in - not including SHIP0 yet as waiting on access. 
touch datasets_info
vim datasets_info
dasuqc1_mdd_cof3_eur_sr-qc.hg19.ch.fl
dasuqc1_mdd_col3_eur_sr-qc.hg19.ch.fl
dasuqc1_mdd_rot4_eur_sr-qc.hg19.ch.fl
dasuqc1_mdd_nes1_eur_sr-qc.hg19.ch.fl

#Actual code
postimp_navi --out MDD1_noship --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --bgscore --pheno MDD1_PGCMDD2_final.txt

 
