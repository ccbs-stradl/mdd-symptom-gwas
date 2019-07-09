#Bash

#Go into interactive mode
srun -n 16 -t 1:00:00 --pty bash -il

#Load R

module load R/3.3.1

R

#R

###########################################################################################################################
########################################### Select phenotypes from all the PGC MDD2 files #################################
###########################################################################################################################
####### selectpheno is a basic R function to extract the phenotypes for a particular variable(s) ##########################
####### Please install the R libary gdata which requires the read.xls function implemented here  ##########################
####### variable = phenotype(s) to extract, eg "Sex" or if more than one c("Sex","Ageinterview") ##########################
####### working directory = path directory containing the xls pheno files for each cohort #################################
####### outputdirectory = directory where to save the extracted phenotype output file #####################################
###########################################################################################################################

library(yaml)
library(gdata)
library(stringr)
library(plyr)
library(dplyr)

config <- yaml.load_file('config.yaml')

selectpheno = function(variable, workingdirectory, outputdirectory) {
library(gdata)
files <- list.files(path=workingdirectory, pattern="*.xls")
res1 <- do.call(`rbind`,lapply(files, read.xls, header=T))
mddtab1 <-as.matrix(res1[,c("ID1","ID2","Study","Sub.study",paste(variable))])
filename <- paste(variable,"_PGCMDD2",sep="",".tsv")
write.table(mddtab1,file.path(outputdirectory,filename),row=F,col=T,sep="\t",quote=F)
}

##################################################################################################
#Part 1

#Extract 9 pheno files - One for each symptom - for some reason passing a vector straight to the parameter 'variable' resulted in an error so set it up in a for loop instead.
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")

dir.create(config$data$pgc$mdd$output, recursive=TRUE)

for (i in symptoms) {
selectpheno(variable=i, workingdirectory=file.path(config$data$pgc$mdd$v1, "secondary_phenotypes/1216_update"), outputdirectory=config$data$pgc$mdd$output)
}

#Rotterdam dataset to be extracted separately
rotterdam <- read.xls(file.path(config$data$pgc$mdd$output, "secondary_phenotypes/1216_update/ROTTERDAM_PGC2updated.xls", header=T)

#GSK dataset to be extracted separately
gsk <- read.xls(file.path(config$data$pgc$mdd$output, "secondary_phenotypes/1216_update/GSK_PGC2updated.xls", header=T)


#Save down these separately extracted files
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
for (i in symptoms) {
rot <- rotterdam[,c('ID1','ID2','Study','Sub.study',i)]
write.table(rot,paste(i,"rotterdam_PGC2.tsv",sep="_"),sep="\t")
gs <- gsk[,c('ID1','ID2','Study','Sub.study',i)]
write.table(gs,paste(i,"gsk_PGC2.tsv",sep="_"),sep="\t")
}

#Part 2
##Subset by study ID that is relevant to my cohorts - cases and controls
for (i in symptoms) {
rotterdam1 <- read.table(file.path(config$data$pgc$mdd$output, paste0(i,'_rotterdam_PGC2.tsv', sep="")), sep="\t", header=T)
gsk1 <- read.table(file.path(config$data$pgc$mdd$output,paste0(i,'_gsk_PGC2.tsv', sep="")), sep="\t", header=T)
mdd1 <-read.table(file.path(config$data$pgc$mdd$output, paste0('symptom_',i,'_PGCMDD2.xls', sep="")), header=T)
mdd1_norotgsk <- mdd1 %>% 
  				filter(!str_detect(ID1, '_rot4_eur_sr|_gsk2_eur_sr_'))
MDD <- rbind.fill(mdd1_norotgsk, rotterdam1, gsk1)

keepStudy <- c("COF", "Colaus", "NESDA/NTR", "ROTTERDAM", "SHIP-0", "RS")
#Change pheno file for 1=controls, 2=cases and -9=unknowns as stated on RICOPILI
MDD[,i][!MDD$Study %in% keepStudy] <- -9 #Remove cases and controls not in the correct cohort
MDD[,i][is.na(MDD[,i])] <- -9  #Recode NAs to -9
MDD[,i][MDD[,i]==1] <- 2 
MDD[,i][MDD[,i]==0] <- 1 
print(table(MDD[[i]]))
write.table(MDD, file.path(config$data$pgc$mdd$output, paste0("symptom_",i,"_PGCMDD2_recoded.tsv", sep="")), sep="\t", row.names=F)
}


#Check coding of phenos to see if any adjustments are required
for (i in symptoms) {
MDD <-read.table(file.path(config$data$pgc$mdd$output, paste0("symptom_",i,"_PGCMDD2_recoded.tsv", sep="")), sep="\t" header=T)
print(dim(MDD))
print(table(MDD[[i]]))
} #Everything seems to be working - 6th dataset needs a 9 recoded as -9

#Part 3
#Remove overlapping samples with UKB - remove according to ID1 and ID2 - UP TO HERE!
#symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
#UKB <- read.table(config$data$pgc$mdd$ukb_remove, sep="\t", header=F)
#for (i in symptoms) {
#MDD <-read.csv(file.path(config$data$pgc$mdd$output, paste0("symptom_",i,"_PGCMDD2_recoded.csv")), header=T)
#print(table(MDD[[i]]))
#MDD <- MDD[!MDD$ID1 %in% UKB$V2]

}

