# reformat Clin AppInc symptoms gwas into commongwas columns for 
# input into PRS calculation

# input columns

#  1	CHR
#  2	SNP
#  3	BP
#  4	A1
#  5	A2
#  6	FRQ_A_17160
#  7	FRQ_U_8330
#  8	INFO
#  9	OR
# 10	SE
# 11	P
# 12	ngt
# 13	Direction
# 14	HetISqt
# 15	HetDf
# 16	HetPVa
# 17	Nca
# 18	Nco
# 19	Neff_half

# output columns: SNP     CHR     BP      MAF     A1      A2      est     se_c         Pval_Estimate

zcat meta/distribution/AGDS_PGC.MDD3b_weightGain/daner_AGDS_PGC.MDD3b_weightGain.gz | awk '{if(NR == 1) {print "SNP", "CHR", "BP", "MAF", "A1", "A2", "est", "se_c", "Pval_Estimate"} else {print $2, $1, $3, $6, $4, $5, log($9), $10, $11}}' | gzip -c > meta/usergwas/agds_pgc.appinc.sumstats.txt.gz 