setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
library(ExtractR)
#source("Rcode/functions/Functions.R")

#Get Calcom ages and save results
if(F) {
	calcomSamps <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/CALCOMsamps.sql",db="CALCOM",uid="allan.hicks",sp="WDOW")
    save(calcomSamps,file="extractedData/CALCOMsamps_wdow.Rdat")
    file.copy("extractedData/CALCOMsamps_wdow.Rdat",paste("extractedData/CALCOMsamps_wdow_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)


    calcomAges <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/CALCOMages.sql",db="CALCOM",uid="allan.hicks",sp="WDOW")
    save(calcomAges,file="extractedData/CALCOMages_wdow.Rdat")
    file.copy("extractedData/CALCOMages_wdow.Rdat",paste("extractedData/CALCOMages_wdow_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}

setwd("C:/NOAA2015/Widow/Data")
load("extractedData/CALCOMages_wdow.Rdat")

table(substring(calcomAges$sample_no,1,4))

data.frame(best_age=sum(!is.na(calcomAges$best_age)),
		   age1=sum(!is.na(calcomAges$age1)),
		   age2=sum(!is.na(calcomAges$age2)),
		   age3=sum(!is.na(calcomAges$age3)),
		   age4=sum(!is.na(calcomAges$age4)))

x <- table(calcomAges$best_age,calcomAges$age1)
x[x==0] <- ""
x
x <- table(calcomAges$best_age,calcomAges$age2)
x[x==0] <- ""
x
x <- table(calcomAges$best_age,calcomAges$age3)
x[x==0] <- ""
x
x <- table(calcomAges$best_age,calcomAges$age4)
x[x==0] <- ""
x

#what is difference between age2 and age3
x <- table(calcomAges$age2,calcomAges$age3)
x[x==0] <- ""
x

table(!is.na(calcomAges$age2),!is.na(calcomAges$age3)
	#for every obs with age 2, there is also an age3

doubleReads <- calcomAges[!is.na(calcomAges$age2) | !is.na(calcomAges$age3),]
write.csv(doubleReads,file="Biological/Ages/CA_CALCOM/CalcomDoubleReads.csv",row.names=F)

table(substring(doubleReads$sample_no,1,4))
