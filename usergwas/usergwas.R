library(GenomicSEM)

# get task number and max task number
k <- as.numeric(Sys.getenv("SGE_TASK_ID"))
max_tasks <- as.numeric(Sys.getenv("SGE_TASK_LAST"))

cat(paste("Running task", k, " of ", max_tasks, "\n"))

covstruct_prefix <- 'agds_pgc.alspac_ukb.common.covstruct'
sumstats_prefix <- 'agds_pgc.alspac_ukb.common.sumstats'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))

gwas_prefix <- 'mdd_aty_gate'
gwas_tsv <- here::here('meta', 'usergwas', gwas_prefix, paste(gwas_prefix, k, 'rds', 'gz', sep='.'))

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

# common path model
model_common <- "
MDD =~ NA*CommAnh + UkbAnh + CommSleDec + CommAppDec + CommGuilt + CommDep + UkbDep + CommConc + CommSui 
ATY =~ NA*CommAppInc + CommSleInc + CommFatig
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

ATY ~~ 1*ATY
MDD ~~ 1*MDD
GATE ~~ 1*GATE
GATE ~~ 0*MDD + 0*ATY

ATY ~ SNP
MDD ~ SNP
GATE ~~ SNP
"

mel_aty_common <-userGWAS(covstruc = symptoms_covstruct,
    SNPs = SNPs,
    estimation = "DWLS",
    model = model_common,
    printwarn = TRUE,
    sub=c("ATY ~ SNP", "MDD ~ SNP", "GATE ~ SNP"),
    cores = 8,
    toler = FALSE,
    SNPSE = FALSE,
    parallel = TRUE,
    GC="standard",
    MPI=FALSE,
    smooth_check=TRUE)
    