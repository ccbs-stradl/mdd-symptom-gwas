library(GenomicSEM)

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))

cat(paste("Running task", k, " of ", max_tasks, "\n"))

covstruct_prefix <- 'all.common.covstruct'
sumstats_prefix <- 'all.common.sumstats'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))

gwas_prefix <- 'all.mel_aty_gate'
gwas_rds <- here::here('meta', 'usergwas', gwas_prefix, paste(gwas_prefix, k, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)
symptoms_sumstats <- readRDS(sumstats_rds)

# get subset of SNPs to analyze in this task
nsnps <- nrow(symptoms_sumstats)
snps_per_task <- ceiling(nsnps / max_tasks)

# index of task IDs repeated for number of SNPs per ask, clip the last task
# to the actual number of snps

task_allocations <- rep(1:max_tasks, each=snps_per_task)[seq.int(nsnps)]

SNPs <- symptoms_sumstats[which(task_allocations == k),]

rm(symptoms_sumstats)

# symptom factor path model
mel_aty_gate_model <- "
MELG =~ AllDep + AllAnh + AllAppDec + AllMotoInc + AllGuilt + AllConc + AllSui 
ATYG =~ AllAppInc + AllSleInc
GATE =~ AllDep + AllAnh

MELG~~ATYG
GATE ~~ 0*MELG + 0*ATYG

MELG ~ SNP
ATYG ~ SNP
GATE ~ SNP
"

mel_aty_gate_gwas <-userGWAS(covstruc = symptoms_covstruct,
    SNPs = SNPs,
    estimation = "DWLS",
    model = mel_aty_gate_model,
    printwarn = TRUE,
    sub=c("MELG ~ SNP", "ATYG ~ SNP", "GATE ~ SNP"),
    cores = 8,
    toler = FALSE,
    SNPSE = FALSE,
    parallel = TRUE,
    GC="standard",
    MPI=FALSE,
    smooth_check=TRUE)

mel_aty_model <- "
MEL =~ AllDep + AllAnh + AllAppDec + AllMotoInc + AllGuilt + AllConc + AllSui 
ATY =~ AllAppInc + AllSleInc

MEL~~ATY

MEL ~ SNP
ATY ~ SNP
"

mel_aty_gwas <-userGWAS(covstruc = symptoms_covstruct,
    SNPs = SNPs,
    estimation = "DWLS",
    model = mel_aty_model,
    printwarn = TRUE,
    sub=c("MEL ~ SNP", "ATY ~ SNP"),
    cores = 8,
    toler = FALSE,
    SNPSE = FALSE,
    parallel = TRUE,
    GC="standard",
    MPI=FALSE,
    smooth_check=TRUE)
    
# common factor path model
mdd_model <- "
MDD =~ AllDep + AllAnh + AllAppDec + AllMotoInc + AllGuilt + AllConc + AllSui + AllAppInc + AllSleInc

MDD ~ SNP
"
    
mdd_gwas <-userGWAS(covstruc = symptoms_covstruct,
    SNPs = SNPs,
    estimation = "DWLS",
    model = mdd_model,
    printwarn = TRUE,
    sub=c("MDD ~ SNP"),
    cores = 8,
    toler = FALSE,
    SNPSE = FALSE,
    parallel = TRUE,
    GC="standard",
    MPI=FALSE,
    smooth_check=TRUE)
    
# # Calculate the chi-square for each SNP
# Q_chisq <- mdd_gwas[[1]]$chisq - mel_aty_gwas[[1]]$chisq
# 
# # Calculate the df associated with this chi-square 
# Q_df <- mdd_gwas[[1]]$chisq_df[1] - mel_aty_gwas[[1]]$chisq_df[1]
# 
# #Step 3c: Calculate the p-value associated with the Q-statistic
# Q_chisq_pval <- pchisq(Q_chisq, Q_df, lower.tail=FALSE)

gwas <- c(mel_aty_gate_gwas, mel_aty_gwas, mdd_gwas)

names(gwas) <- c('MELG', 'ATYG', 'GATE', 'MEL', 'ATY', 'MDD')


gwas_rows <- dplyr::bind_rows(gwas)

saveRDS(symptoms_gwas, gwas_rds)