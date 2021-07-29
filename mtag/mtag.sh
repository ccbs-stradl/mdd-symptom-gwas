#$ -N symp_mtag
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 8
#$ -cwd
#$ -e logs
#$ -o logs


# Rscript sumstats.sh
# git clone git@github.com:JonJala/mtag.git
# conda env create --file  mtag/environment.yaml

source ~/.bashrc
conda activate mtag
mkdir clin_pop

python mtag/mtag.py \
--sumstats Clin_AppInc.txt,Clin_SleDec.txt,Clin_SleInc.txt,Pop_Dep.txt,Pop_Anh.txt,Pop_AppDec.txt,Pop_AppInc.txt,Pop_SleDec.txt,Pop_SleInc.txt,Pop_Fatig.txt,Pop_Guilt.txt,Pop_Sui.txt \
--use_beta_se \
--out clin_pop/clinpop \
--stream_stdout

#  1	Clin_AppInc.txt
#  2	Clin_SleDec.txt
#  3	Clin_SleInc.txt
#  4	Pop_Dep.txt
#  5	Pop_Anh.txt
#  6	Pop_AppDec.txt
#  7	Pop_AppInc.txt
#  8	Pop_SleDec.txt
#  9	Pop_SleInc.txt
# 10	Pop_Fatig.txt
# 11	Pop_Guilt.txt
# 12	Pop_Sui.txt

# Clin_SleDec.txt,Clin_SleInc.txt,Pop_Anh.txt,Pop_AppInc.txt,
# Pop_Dep.txt,Pop_Guilt.txt,Pop_PsycInc.txt,Pop_SleInc.txt,Pop_AppDec.txt,Pop_Conc.txt,Pop_Fatig.txt,Pop_PsycDec.txt,Pop_SleDec.txt,Pop_Sui.txt 

# Clin_AppInc.txt GC=1.059                                                      
# Clin_SleDec.txt GC=1.024                                                      
# Clin_SleInc.txt GC=1.03                                                       
# Pop_Dep.txt GC=1.118                                                          
# Pop_Anh.txt GC=1.116                                                          
# Pop_AppDec.txt GC=1.023                                                       
# Pop_AppInc.txt GC=1.054                                                       
# Pop_SleDec.txt GC=1.025                                                       
# Pop_SleInc.txt GC=1.029                                                       
# Pop_Fatig.txt GC=1.026                                                        
# Pop_Guilt.txt GC=1.055                                                        
# Pop_Sui.txt GC=1.03   