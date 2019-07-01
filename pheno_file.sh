#Bash

#Softlink to the UKB Overlap dataset to remove
ln -s /home/madams/pgcdac/PGCUKBoverlap/959_PGC_UKB_overlap.txt

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
####### outpurdirectory = directory where to save the extracted phenotype output file #####################################
###########################################################################################################################

selectpheno = function(variable, workingdirectory, outputdirectory) {
setwd(paste(workingdirectory))
library(gdata)
files <- list.files(pattern="*.xls")
res1 <- do.call(`rbind`,lapply(files, read.xls, header=T))
mddtab1 <-as.matrix(res1[,c("ID1","ID2","Study","Sub.study",paste(variable))])
filename <- paste(variable,"_PGCMDD2",sep="",".xls")
write.table(mddtab1,paste(outputdirectory,filename,sep=""),row=F,col=T,sep="\t",quote=F)
    }

##################################################################################################
#Part 1

#Extract 9 pheno files - One for each symptom - for some reason passing a vector straight to the parameter 'variable' resulted in an error so set it up in a for loop instead.
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")

for (i in symptoms) {
selectpheno(variable=i, workingdirectory="/home/jermy/v1/secondary_phenotypes/1216_update", outputdirectory="/home/jermy/symptom_")
}

selectpheno(variable="MDD1", workingdirectory="/home/jermy/v1/secondary_phenotypes/1216_update", outputdirectory="/home/jermy/symptom_")

#Part 2
##Subset by study ID that is relevant to my cohorts - cases and controls
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
for (i in symptoms) {
MDD <-read.table(paste("symptom_",i,"_PGCMDD2.xls", sep=""), sep="\t", header=T)
keepStudy <- c("COF", "Colaus", "NESDA/NTR", "ROTTERDAM", "SHIP-0", "RS")
#Change pheno file for 1=controls, 2=cases and -9=unknowns as stated on RICOPILI
MDD[,i][!MDD$Study %in% keepStudy] <- -9 #Remove cases and controls not in the correct cohort
MDD[,i][is.na(MDD[,i])] <- -9  #Recode NAs to -9
MDD[,i][MDD[,i]==1] <- 2 
MDD[,i][MDD[,i]==0] <- 1 
write.csv(MDD, paste("symptom_",i,"_PGCMDD2_recoded.csv", sep=""))
}

#Check coding of phenos to see if any adjustments are required
for (i in symptoms) {
MDD <-read.csv(paste("symptom_",i,"_PGCMDD2_recoded.csv", sep=""), header=T)
print(dim(MDD))
print(table(MDD[[i]]))
} #Everything seems to be working

#Part 3
#Remove overlapping samples with UKB
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
UKB <- read.table("959_PGC_UKB_overlap.txt", sep="\t", header=F)
for (i in symptoms) {
MDD <-read.csv(paste("symptom_",i,"_PGCMDD2_recoded.csv", sep=""), header=T)
print(table(MDD[[i]]))
#MDD <- MDD[!MDD$ID1 %in% UKB$V2]
print(sum(MDD$ID2 %in% UKB$V3)) #Matches with 315
print(sum(MDD$ID1 %in% UKB$V2)) #Matches with 2287 people instead of 959 - duplicates in my dataset?
print(sum(duplicated(MDD$ID1)))  #duplicates in the datasets - ask Mark about this. 
#write.csv(MDD, paste("symptom_",i,"_PGCMDD2_recoded.csv", sep=""))
}

