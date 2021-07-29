alspac=../sumstats/ALSPAC
daner=$alspac/daner

mkdir -p $daner

Rscript alspac_daner.R \
$alspac/RESULTS.mood.21042020.txt.gz \
$daner/daner_MDD1_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.inter.21042020.txt.gz \
$daner/daner_MDD2_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.wloss.21042020.txt.gz \
$daner/daner_MDD3a_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.wgain.21042020.txt.gz \
$daner/daner_MDD3b_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.insom.21042020.txt.gz \
$daner/daner_MDD4a_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.hypsom.21042020.txt.gz \
$daner/daner_MDD4b_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.agit.21042020.txt.gz \
$daner/daner_MDD5a_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.retar.21042020.txt.gz \
$daner/daner_MDD5b_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.fati.21042020.txt.gz \
$daner/daner_MDD6_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.guilt.21042020.txt.gz \
$daner/daner_MDD7_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.conc.21042020.txt.gz \
$daner/daner_MDD8_ALSPAC_CISR

Rscript alspac_daner.R \
$alspac/RESULTS.suic.21042020.txt.gz \
$daner/daner_MDD9_ALSPAC_CISR
