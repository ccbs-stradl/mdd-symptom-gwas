#Bash

#Go into interactive mode
srun -n 16 -t 1:00:00 --pty bash -il

#Load R

module load R/3.4.3

R

#R

library(yaml)
library(gdata)
library(stringr)
library(plyr)
library(dplyr)

config <- yaml.load_file('config.yaml')


##################################################################################################
#Part 1

#Extract 9 pheno files - One for each symptom - for some reason passing a vector straight to the parameter 'variable' resulted in an error so set it up in a for loop instead.
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
studies <- c("RADIANT", "STR", "ROTTERDAM", "QIMR", "NESDANTR", "GSK", "GENRED", "GENRED2", "EDINBURGH", "CoLaus", "COFAMS", "BonnMH")
study_noupdate <- c("STARD", "SHIP0", "MARS", "Janssen", "GENPOD")

dir.create(config$data$pgc$mdd$output, recursive=TRUE)


updated <- data.frame()
for (i in studies) {
mdd <- read.xls(file.path(config$data$pgc$mdd$v1, paste0("secondary_phenotypes/1216_update/",i,"_PGC2updated.xls")), header=T)
updated <- rbind.fill(updated, mdd)
}

noupdate <- data.frame()
for (j in study_noupdate) {
mdd <- read.xls(file.path(config$data$pgc$mdd$v1, paste0("secondary_phenotypes/1216_update/",j,"_PGC2.xls")), header=T)
noupdate <- rbind.fill(noupdate, mdd)
}

for (i in symptoms){
full <- rbind.fill(updated,noupdate)
full <- full[,c("ID1","ID2","Study","Sub.study", i)]
write.table(full,file.path(config$data$pgc$mdd$output, paste0(i,"_PGCMDD2.tsv", sep="")), sep="\t", row.names=F)
}

#Part 2
##Subset by study ID that is relevant to my cohorts - cases and controls
for (i in symptoms) {
MDD <-read.table(file.path(config$data$pgc$mdd$output, paste0(i,'_PGCMDD2.tsv', sep="")), sep="\t", header=T)
keepStudy <- c("COF", "Colaus", "NESDA/NTR", "ROTTERDAM", "SHIP-0", "RS")
#Change pheno file for 1=controls, 2=cases and -9=unknowns as stated on RICOPILI
MDD[,i][!MDD$Study %in% keepStudy] <- -9 #Remove cases and controls not in the correct cohort
MDD[,i][is.na(MDD[,i])] <- -9  #Recode NAs to -9
MDD[,i][MDD[,i]==9] <- -9  #Recode 9 to -9
MDD[,i][MDD[,i]==1] <- 2 
MDD[,i][MDD[,i]==0] <- 1 
print(table(MDD[[i]]))
write.table(MDD, file.path(config$data$pgc$mdd$output, paste0(i,"_PGCMDD2_recoded.tsv", sep="")), sep="\t", row.names=F)
}


#Part 3
#Remove overlapping samples with UKB - remove according to ID1 and ID2 - UP TO HERE!

UKB <- read.table(config$data$pgc$mdd$ukb_remove, sep="\t", header=F)
UKB$ID1adj <- gsub("\\*.*","",UKB$V2) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
for (i in symptoms) {
MDD <-read.table(file.path(config$data$pgc$mdd$output, paste0(i,"_PGCMDD2_recoded.tsv", sep="")), sep="\t", header=T)
MDD$ID1adj <- gsub("\\*.*","",MDD$ID1) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
print(dim(MDD))
#filter MDD for IDs not in the UKB sample
MDD2 <- MDD %>%
filter(!MDD$ID2 %in% UKB$V3 & !MDD$ID1adj %in% UKB$V2)
print(dim(MDD2))
print(MDD %>%
filter(MDD$ID2 %in% UKB$V3 & MDD$ID1adj %in% UKB$ID1adj) %>% group_by(Study) %>% tally())
#293 removed due to overlap with UKB - 205 from RADIANT which makes sense
MDD2 <- MDD2[,-c(3,4,6)]
write.table(MDD2, file.path(config$data$pgc$mdd$output, paste0(i,"_PGCMDD2_final.tsv", sep="")), sep="\t", row.names=F)
}


#COMPLETE

