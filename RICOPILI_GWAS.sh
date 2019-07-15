#bash

#Install and unpack ricopili - follow this link for custom installation https://docs.google.com/document/d/14aa-oeT5hF541I8hHsDAL_42oyvlHRC5FWR7gir4xco/edit#heading=h.y8igfs7neh22
wget https://sites.google.com/a/broadinstitute.org/ricopili/download//rp_bin.2019_Jun_25.001.tar.gz

tar -xvzf rp_bin.2019_Jun_25.001.tar.gz

#Create custom installation file in ricopili

#Test to make sure it has installed correctly
#Softlink to testing datasets
ln -s /home/pgcdac/rp_bin_1018/dependencies/testing_data/PGC_cohort1.ch.fl.r4.gz
ln -s /home/pgcdac/rp_bin_1018/dependencies/testing_data/PGC_cohort2.ch.fl.r4.gz
ln -s /home/pgcdac/rp_bin_1018/dependencies/testing_data/PGC_cohort3.ch.fl.r4.gz
ln -s /home/pgcdac/rp_bin_1018/dependencies/testing_data/PGC_cohort4.ch.fl.r4.gz
ln -s /home/pgcdac/rp_bin_1018/dependencies/testing_data/PGC_meta.r4.gz

#Sorting out reference info as none available currently - softlink to file
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/reference_info

#softlink to all the imputed files within my directory - didn't have permission using a wildcard so had to do separately
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rau2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_cof3_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rde4_eur_sr-qc2.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_mmo4_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_mmi2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_qio2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_qi3c_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_qi6c_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_stm2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_twg2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_nes1_eur_sr-qc2.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_grdg_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_grnd_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rai2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_gep3_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rad3_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_gsk2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rot4_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_rage_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_i2b3_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_edi2_eur_sr-qc.hg19.ch.fl
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/dasuqc1_mdd_col3_eur_sr-qc.hg19.ch.fl

#Softlink to the pca file - NEED TO CONFIRM THIS IS THE CORRECT FILE FROM SASKIA/OLLI
ln -s /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/mdd/wave2/v1/MDD29.0515.nproj.menv.mds_cov

#Postimp RICOPILI code - run this from my home directory now all softlinks have been set-up
postimp_navi --out MDD1_case_con --mds MDD29.0515.nproj.menv.mds_cov --coco 1,2,3,4,5,6 --addout run1 --bgscore --nocon --pheno MDD1_PGCMDD2_final.tsv

 