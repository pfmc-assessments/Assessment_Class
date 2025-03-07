dont source
setwd("C:\\NOAA2015\\Widow\\Data\\PacFIN")


load("PacFIN.WDOW.Catch.by.Gear.05.Jun.2015.dmp")
dat <- PacFIN.WDOW.Catch.by.Gear.05.Jun.2015
rm(PacFIN.WDOW.Catch.by.Gear.05.Jun.2015)

dat$gear <- NA
dat$gear[dat$GRID%in%c("SST","SHT","PWT","DST")] <- "ShrimpTrawl"
dat$gear[dat$GRID%in%c("RLT","OTW","GFT","GFS","GFL","FTS","FFT","BTT","BMT")] <- "BottomTrawl"
dat$gear[dat$GRID%in%c("MDT")] <- "MidwaterTrawl"
dat$gear[dat$GRID%in%c("PRT","DNT")] <- "MiscTrawl"
dat$gear[dat$GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- "Pot"
dat$gear[dat$GRID%in%c("JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
dat$gear[dat$GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
dat$gear[dat$GRID%in%c("DVG","USP")] <- "Other"

dat$mt <- dat$CATCH.KG/1000
tapply(dat$mt,list(dat$YEAR,dat$SPID,dat$PCID),sum)

tapply(dat$mt,list(dat$YEAR,dat$GRID,dat$PCID),sum)
tapply(dat$mt,list(dat$YEAR,dat$gear,dat$PCID),sum)


##### By Gear Group
# load("PacFIN.WDOW.Catch.by.Gear.Group.05.Jun.2015.dmp")
# dat <- PacFIN.WDOW.Catch.by.Gear.Group.05.Jun.2015
# rm(PacFIN.WDOW.Catch.by.Gear.Group.05.Jun.2015)
# lapply(dat[,1:6],table)
### Not useful because MDT in with TWL

































##############################################################################
### OLD rougheye stuff below

#Function to extract PacFIN rougheye catch
#Needs the RODBC library
#set working directory to have subdirectories:
###sql: with the sql script called Pacfin.Catch.query
###Catches: To save the catch results in
#This assigns the fleet type and summarizes over fleet and year and period
#Writes out two files, one on catches by each fleet by year
#Possible species codes
### GRND ROCK NSLP SSLP REYE       ROUGHEYE ROCKFISH         SEBASTES ALEUTIANUS
### GRND ROCK UDW1                 SHORTRAKER+ROUGHEYE       SEBASTES SPP.         (No results using UDW1 as species code)
### GRND ROCK NSLP SSLP USLP       UNSP. SLOPE ROCKFISH      N/A


setwd("C:/NOAA2015/Widow/Data")
library(RODBC)  #need 32-bit R
source("Rcode/functions/Functions.R")

#Get PacFIN catches and save results
if(F) {
    pcatch <- queryDB(queryFilename="Pacfin.Catch",db="PACFIN",uid="hicksa",sp="REYE", querydir="sql/")
    save(pcatch,file="extractedData/PacfinCatch.REYE.Rdat")
    file.copy("extractedData/PacfinCatch.REYE.Rdat",paste("extractedData/PacfinCatch.REYE_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

    #There was no UDW1 or URK1
    pcatch.udw1 <- queryDB(queryFilename="Pacfin.Catch",db="PACFIN",uid="hicksa",sp="UDW1", querydir="sql/")  #SHortRaker/REYE???

    pcatch.urk1 <- queryDB(queryFilename="Pacfin.Catch",db="PACFIN",uid="hicksa",sp="URK1", querydir="sql/")  #SRKR+REYE+NRCK+SHRP
}

if(F) {
    pcatch3 <- queryDB(queryFilename="PARGRPexamp",db="PACFIN",uid="hicksa",sp="REYE", querydir="sql/")
    save(pcatch3,file="extractedData/PacfinCatchTribal.Rdat")
    file.copy("extractedData/PacfinCatchTribal.Rdat",paste("extractedData/PacfinCatchTribal_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}
tapply(pcatch3$MTONS,list(pcatch3$MONTH,pcatch3$PARGRP),sum)



###########################################
#  John's extraction
###########################################
load("C:\\NOAA2013\\Rougheye\\Data\\Catches\\PacFIN.REYE.Catch.by.Gear.Group.08.May.13.dmp")
pcatch <- PacFIN.REYE.Catch.by.Gear.Group.08.May.13
rm(PacFIN.REYE.Catch.by.Gear.Group.08.May.13)
lapply(pcatch,unique)
pcatch$MT <- pcatch$CATCH.KG * 0.001
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$PCID,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID),sum)
#Careful when comparing because John only provides ACA, AOR, and AWA
















load("extractedData/PacfinCatch.Rdat")

lapply(pcatch,unique)


tapply(pcatch$MT,list(pcatch$YEAR,pcatch$PCID,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$PCID),sum)

tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID),sum)




#Sum over areas to create fleets and write out results
pcatch$state <- substring(pcatch$PCID,2)
pcatch$fleet <- as.character(pcatch$GRID)
pcatch$yrFact <- factor(pcatch$YEAR,levels=min(pcatch$YEAR):max(pcatch$YEAR))
pcatch[pcatch$PCID=="ACN","fleet"] <- "USforeign"
pcatch[pcatch$PCID=="AJV","fleet"] <- "USjv"
pcatch[pcatch$PCID=="ACN","state"] <- "Foreign"
pcatch[pcatch$PCID=="AJV","state"] <- "JV"
pcatch <- pcatch[order(pcatch$fleet,pcatch$state,pcatch$YEAR,pcatch$PERIOD),]
pcatch.yr.fl.st <- tapply(pcatch$MT,list(pcatch$yrFact,pcatch$fleet,pcatch$state),sum)
#pcatch.yr.fl.st <- aggregate(pcatch$MT,list(pcatch$yrFact,pcatch$state,pcatch$fleet),sum)
out <- NULL
for(i in 1:dim(pcatch.yr.fl.st)[3]) {
    out <- rbind(out,data.frame(year=dimnames(pcatch.yr.fl.st)[[1]],state=dimnames(pcatch.yr.fl.st)[[3]][i],pcatch.yr.fl.st[,,i]))
}
out[is.na(out)] <- 0
write.csv(out,file="Catches/PacFinCatchByYearFleetState.csv",row.names = F)

par(mfrow=c(2,4),las=2)
ylims <- c(0,300)
alph <- 1
cols <- c(rgb(0,0,1,alph),rgb(0.8,0,0.8,alph),rgb(1,0.5,0,alph),rgb(0.7,0,0,alph),rgb(0,0.7,0,alph))
yrs <- as.numeric(dimnames(pcatch.yr.fl.st)[[1]])
xNames <- yrs;xNames[seq(2,length(yrs),2)]<- ""
#HKL
gear <- "HKL"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#NET
gear <- "NET"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#POT
gear <- "POT"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#TWL
gear <- "TWL"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#TWS
gear <- "TWS"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#Foreign
gear <- "USforeign"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#JV
gear <- "USjv"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))





par(mfrow=c(1,2),las=2)
ylims <- c(0,300)
alph <- 1
cols <- c(rgb(0,0,1,alph),rgb(0.8,0,0.8,alph),rgb(1,0.5,0,alph),rgb(0.7,0,0,alph),rgb(0,0.7,0,alph))
yrs <- as.numeric(dimnames(pcatch.yr.fl.st)[[1]])
xNames <- yrs;xNames[seq(2,length(yrs),2)]<- ""
#HKL
gear <- "HKL"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))
#TWL
gear <- "TWL"
gr <- pcatch.yr.fl.st[,gear,]
gr[is.na(gr)] <- 0
tmp <- apply(gr,1,sum,na.rm=T)
xx <- barplot(t(gr),beside=F,col=cols,main=gear,ylim=ylims,names.arg=xNames)
legend("topright",c("WA","OR","CA"),col=cols[c(5,4,1)],pch=15,bty="n",cex=1.1)
axis(1,at=xx[seq(1,length(yrs),2)],labels=rep("",length(xx[seq(1,length(yrs),2)])))


#This does not match the 2000 catch report produced by PacFIN in 2010. http://pacfin.psmfc.org/pacfin_pub/data_rpts_pub/pfmc_rpts_pub/r009_p00.txt








plot(yrs,tmp,type="l",lwd=3,ylim=ylims,xlab="Year",ylab="Catch (mt)",main=gear)
for(i in 1:ncol(gr)) {
    lines(yrs,gr[,i],col=cols[i])
}
legend("topleft",colnames(gr),col=cols,lty=1)



























#Get PacFIN catches with fish ticket identifier and save results
pcatch <- queryDB(queryFilename="Pacfin.Catch.fishtix",db="PACFIN",uid="hicksa",sp="PWHT", querydir="sql/")
save(pcatch,file="extractedData/PacfinCatchFishTix.Rdat")
