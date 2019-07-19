#Bash

#Go into interactive mode
srun -n 16 -t 1:00:00 --pty bash -il

#Load R

module load R/3.4.3

Rscript pheno_pgc.R

#COMPLETE

