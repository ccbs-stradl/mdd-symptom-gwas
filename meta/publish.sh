# set up input files for published sumstats

mkdir inputs
# link to daner sumstats and basics spreadsheet
# rename daner_Cohort.Symptom.gz to Cohort-Symptom.gz
for distribution in distribution/*; do
    dataset=$(basename $distribution)
    input=$(echo $dataset | sed s/\\./-/ )
    ln -sf $(readlink -f $distribution/daner_${dataset}.gz) inputs/${input}.gz
    ln -sf $(readlink -f $distribution/basic.${dataset}.num.xls) inputs/${input}.xls
done

# workflow to format sumstats
nextflow run publish.nf --sumtats="inputs/*.{gz,xls}" --analyst "Mark J Adams" -resume