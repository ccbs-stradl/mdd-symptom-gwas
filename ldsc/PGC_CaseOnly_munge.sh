# save LDSC output to a CSV
echo 'Symptom,Npresent,Nabsent,Lambda GC,Intercept,h2,se' > sumstats/PGC/CaseOnly/CaseOnly.ldsc.csv

# find all daner meta sumstats
for meta in $(ls sumstats/PGC/CaseOnly/*.meta.gz); do
       # get filename prefix and construct output names
       prefix=$(basename $meta .gz)
       out_merge=$(dirname $meta)/${prefix}.ldsc.hm3
       out=$(dirname $meta)/${prefix}.ldsc
       # check if munge has already been done
       # munge and merge with HM3 
       if [ ! -f ${out_merge}.sumstats.gz ]; then
         munge_sumstats.py --daner-n --sumstats $meta --out $out_merge --merge-alleles sumstats/reference/w_hm3.snplist
      fi
      # munge w/o merge
      if [ ! -f ${out}.sumstats.gz ]; then
         munge_sumstats.py --daner-n --sumstats $meta --out $out 
      fi
      # determine number of cases/contols from daner header names
      n_pres=$(zcat $meta | head -n 1 | awk '{print $6}' | grep -o '[[:digit:]]*')
      n_abse=$(zcat $meta | head -n 1 | awk '{print $7}' | grep -o '[[:digit:]]*')
      # calculate sample prevalence from N cases and controls
      samp_prev=$(echo "scale=2; $n_pres / ($n_pres + $n_abse)" | bc)
      # approximate population prevalence as prevalence of depression * prevalence of symptom
      pop_prev=$(echo "scale=2; 0.15 * $n_pres / ($n_pres + $n_abse)" | bc)

      # calculate LDSC h2
      if [ ! -f ${out_merge}.h2.log ]; then
        ldsc.py --h2 ${out_merge}.sumstats.gz \
                --out ${out_merge}.h2 \
                --samp-prev $samp_prev \
                --pop-prev $pop_prev \
                --ref-ld-chr sumstats/reference/eur_w_ld_chr/ \
                --w-ld-chr sumstats/reference/eur_w_ld_chr/
      fi
      
      # extract information from LDSC h2 log file
      mdd=$(echo $prefix | awk -F_ '{print $2}')
      h2=$(cat ${out_merge}.h2.log | grep 'Liability' | awk '{print $5}')
      se=$( cat ${out_merge}.h2.log | grep 'Liability' | awk '{print $6}' | grep -o '[[:digit:]]*\.[[:digit:]]*')
      lambda=$(cat ${out_merge}.h2.log | grep 'Lambda' | awk '{print $3}')
      intercept=$(cat ${out_merge}.h2.log | grep 'Intercept:' | awk '{print $2}')
      # add on to CSV
      echo "${mdd},${n_pres},${n_abse},${lambda},${intercept},${h2},${se}" >> sumstats/PGC/CaseOnly/CaseOnly.ldsc.csv
      
      if [ ! -f ${out}.h2.log ]; then
        ldsc.py --h2 ${out}.sumstats.gz \
                --out ${out}.h2 \
                --samp-prev $samp_prev \
                --pop-prev $pop_prev \
                --ref-ld-chr sumstats/reference/eur_w_ld_chr/ \
                --w-ld-chr sumstats/reference/eur_w_ld_chr/
      fi
done


