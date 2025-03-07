setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
library(ExtractR)
#source("Rcode/functions/Functions.R")

#Get PacFIN catches and save results
#This is a new method, which includes research catch and breaks out tribal catch
#Remove XXX fleet (foreign catch?)
if(F) {
    pcatchWDOW <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/vdrfd.sql",db="PACFIN",uid="hicksa",sp="WDOW", start=1980, end=2014)
    save(pcatchWDOW,file="extractedData/PacfinVdrfdCatch_wdow.Rdat")
    file.copy("extractedData/PacfinVdrfdCatch_wdow.Rdat",paste("extractedData/PacfinVdrfdCatch_wdow_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
    pcatchWDW1 <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/vdrfd.sql",db="PACFIN",uid="hicksa",sp="WDW1", start=1980, end=2014)
    save(pcatchWDW1,file="extractedData/PacfinVdrfdCatch_wdw1.Rdat")
    file.copy("extractedData/PacfinVdrfdCatch_wdw1.Rdat",paste("extractedData/PacfinVdrfdCatch_wdw1_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

    pcatchWDOW$SPID <- "WDOW"
    pcatchWDW1$SPID <- "WDW1"
	pcatch <- rbind(pcatchWDOW,pcatchWDW1)
	save(pcatch,file="extractedData/PacfinVdrfdCatch.Rdat")
    file.copy("extractedData/PacfinVdrfdCatch.Rdat",paste("extractedData/PacfinVdrfdCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

    hake <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/vdrfd.sql",db="PACFIN",uid="hicksa",sp="PWHT", start=1980, end=2014)
    save(hake,file="extractedData/PacfinHakeVdrfdCatch.Rdat")
    file.copy("extractedData/PacfinHakeVdrfdCatch.Rdat",paste("extractedData/PacfinHakeVdrfdCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}


#start here with extraction"
setwd("C:/NOAA2015/Widow/Data")
load("extractedData/PacfinVdrfdCatch.Rdat")
load("extractedData/PacfinHakeVdrfdCatch.Rdat")

table(pcatch$FLEET,pcatch$GRID,pcatch$IFQ_LANDING)
table(pcatch$GRID)
sort(tapply(pcatch$MT,list(pcatch$GRID),sum))
sort(tapply(pcatch$MT,list(pcatch$DAHL_SECTOR),sum))
tapply(pcatch$MT,list(pcatch$DAHL_SECTOR,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$DAHL_SECTOR),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID,pcatch$AGID),sum)
x <- tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID,pcatch$DAHL_SECTOR,pcatch$AGID),sum)
x[,c("RLT","GFT","GFS","GFL","FTS","FFT","OTW","MDT"),c("03","04"),"W"]
x[,c("GFT","MDT"),,"W"]
y <- x[,c("GFT","MDT"),c("03","04","XX"),"W"]
y[is.na(y)] <- 0
y[,,"03"]/(y[,,"03"]+y[,,"04"])
y[,"GFT","03"]/apply(x[,"GFT",,"W"],1,sum,na.rm=T)
y[,"GFT","XX"]/apply(x[,"GFT",,"W"],1,sum,na.rm=T)

#Not sure where OTW (CA only) belongs
#not sure about "OTH"
pcatch$gear <- NA
pcatch$gear[pcatch$GRID%in%c("SST","SHT","PWT","DST")] <- "ShrimpTrawl"
pcatch$gear[pcatch$GRID%in%c("RLT","OTW","GFT","GFS","GFL","FTS","FFT","BMT")] <- "BottomTrawl"
pcatch$gear[pcatch$GRID%in%c("MDT")] <- "MidwaterTrawl"
pcatch$gear[pcatch$GRID%in%c("PRT","DNT")] <- "MiscTrawl"
pcatch$gear[pcatch$GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- "Pot"
pcatch$gear[pcatch$GRID%in%c("JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
pcatch$gear[pcatch$GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
pcatch$gear[pcatch$GRID%in%c("DVG","USP")] <- "Other"


table(pcatch$gear,useNA="ifany")
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$DAHL_SECTOR,pcatch$gear),sum)
#look at how much catch from HnL fishery is open access
hnl <- tapply(pcatch$MT,list(pcatch$YEAR,pcatch$DAHL_SECTOR,pcatch$gear),sum)[,,"HnL"]
hnlTot <- tapply(pcatch$MT,list(pcatch$YEAR,pcatch$gear),sum)[,"HnL"]
oa <- apply(hnl[,c("06","08","10","11","12")],1,sum)
plot(as.numeric(names(oa)),100*oa/hnlTot,type="h",lwd=7,col="purple",ylim=c(0,100),ylab="Percentage Open Access in H&L fishery",yaxs="i")

tapply(pcatch$MT,list(pcatch$YEAR,pcatch$gear),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$FLEET),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$FLEET,pcatch$OVERAGE),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$FLEET),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$AGID),sum)
#USE RMT which is round weight (landed*conversion factor)


tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$PRMTLST),sum)

#split out shoreside hake (DAHL_SECTOR=c(3,17))
#DAHL_SECTOR only available for 1994 and onward
#He defined these basically as catch of PWHT>50% daily for gear

#DAHL_SECTOR is brought in as a character because of XX (and single digits are followed by a space)
pcatch$Dsector <- as.numeric(as.character(pcatch$DAHL_SECTOR))

tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$DAHL_SECTOR,pcatch$AGID),sum)
table(pcatch$YEAR,pcatch$gear,pcatch$DAHL_SECTOR)

tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$AGID,pcatch$Dsector%in%c(3,17)),sum)


tapply(pcatch$RMT,list(pcatch$DAHL_SECTOR),sum)

#match in hake catch
pcatch$hakeCatch <- hake$RMT[match(pcatch$FTID,hake$FTID)]

table(!is.na(pcatch$hakeCatch))
boxplot(split(pcatch$hakeCatch,pcatch$DAHL_SECTOR))

ind <- pcatch$gear%in%c("BottomTrawl","MidwaterTrawl")
par(mfrow=c(2,2))
boxplot(split(pcatch$hakeCatch,pcatch$DAHL_SECTOR))
boxplot(split(pcatch$hakeCatch,pcatch$DAHL_SECTOR),ylim=c(0,20))
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]))
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]),ylim=c(0,20))

ind <- pcatch$gear%in%c("BottomTrawl")
par(mfrow=c(2,2))
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]))
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]),ylim=c(0,20))
ind <- pcatch$gear%in%c("MidwaterTrawl")
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]))
boxplot(split(pcatch$hakeCatch[ind],pcatch$DAHL_SECTOR[ind]),ylim=c(0,20))


tapply(pcatch$RMT[ind],list(pcatch$DAHL_SECTOR[ind],pcatch$gear[ind]),sum)

pcatch$DAHL_HAKE <- FALSE
pcatch$DAHL_HAKE[pcatch$Dsector%in%c(3,17)] <- TRUE
pcatch$DAHL_HAKE[pcatch$DAHL_SECTOR=="XX"] <- NA
table(pcatch$DAHL_HAKE,useNA="ifany")

pcatch$HIX_HAKE <- FALSE
pcatch$HIX_HAKE[pcatch$gear%in%c("BottomTrawl","MidwaterTrawl") & pcatch$hakeCatch > 0.1] <- TRUE
table(pcatch$HIX_HAKE,useNA="ifany")
x <- table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE)
x

opt <- function(hakeMin) {
    pcatch$HIX_HAKE <- FALSE
    pcatch$HIX_HAKE[pcatch$gear%in%c("BottomTrawl","MidwaterTrawl") & pcatch$hakeCatch > hakeMin] <- TRUE
    pcatch$HIX_HAKE[pcatch$Month<5] <- FALSE
    pcatch$HIX_HAKE[pcatch$Month>11] <- FALSE
    x <- table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE)
    outNum <- c(x[1,2],x[2,1])
    outProp <- c(x[1,2],x[2,1])/sum(x)
    names(outNum) <- names(outProp) <- c("misClassAsHake","misClassNotHake")
    x <- tapply(pcatch$RMT,list(pcatch$DAHL_HAKE,pcatch$HIX_HAKE),sum)
    outWt <- c(x[1,2],x[2,1])
    outPrWt <- c(x[1,2],x[2,1])/sum(x)
    names(outWt) <- names(outPrWt) <- c("misClassAsHake","misClassNotHake")
    return(list(outNum=outNum,outProp=outProp,outWt=outWt,outPrWt=outPrWt))
}
hakeTons <- c(seq(0,0.08,0.02),seq(0.1,1,0.1))
misClass <- data.frame(hakeTons=hakeTons,misClassAsHakeNum=NA,misClassNotHakeNum=NA,misClassAsHakeProp=NA,misClassNotHakeProp=NA,TotalNum=NA,
                       misClassAsHakeWt=NA,misClassNotHakeWt=NA,misClassAsHakePrWt=NA,misClassNotHakePrWt=NA,TotalWt=NA)
for(i in 1:length(hakeTons)) {
    tmp <- opt(hakeTons[i])
    misClass[i,c("misClassAsHakeNum","misClassNotHakeNum","misClassAsHakeProp","misClassNotHakeProp",
                 "misClassAsHakeWt","misClassNotHakeWt","misClassAsHakePrWt","misClassNotHakePrWt")] <-
                 c(tmp$outNum,tmp$outProp,tmp$outWt,tmp$outPrWt)
    misClass$TotalNum[i] <- sum(tmp$outNum)
    misClass$TotalWt[i] <- sum(tmp$outWt)
}

par(mfrow=c(2,1),mar=c(3.5,4,1,5))
plot(misClass$hakeTons,misClass$TotalNum/1000,ylim=c(0,max(misClass$TotalNum/1000)),type="b",pch=20,ylab="# of Observations (Thousands)",las=1,xlab="")
par(new=T)
plot(misClass$hakeTons,misClass$misClassAsHakeProp,ylim=c(0,0.1),type="b",pch=17,col="orangered",yaxt="n",ylab="",xlab="",xaxt="n")
lines(misClass$hakeTons,misClass$misClassNotHakeProp,type="b",pch=15,col="blue")
axis(4,las=1,col.ticks="orangered",col="blue")
mtext("Proportion (of numbers)",side=4,line=3.2)
#legend("right",c("Misclassify as hake","Misclassify as not hake"),pch=c(17,15),col=c("orangered","blue"),lty=1)
plot(misClass$hakeTons,misClass$TotalWt/1000,ylim=c(0,max(misClass$TotalWt/1000)),type="b",pch=20,ylab="Tons of Widow Rockfish (Thousands)",las=1,xlab="")
par(new=T)
plot(misClass$hakeTons,misClass$misClassAsHakePrWt,ylim=c(0,0.04),type="b",pch=17,col="orangered",yaxt="n",ylab="",xlab="",xaxt="n")
lines(misClass$hakeTons,misClass$misClassNotHakePrWt,type="b",pch=15,col="blue")
axis(4,las=1,col.ticks="orangered",col="blue")
mtext("Proportion (of weight)",side=4,line=3.2)
legend("topright",c("Misclassify as hake","Misclassify as not hake"),pch=c(17,15),col=c("orangered","blue"),lty=1)
mtext("Hake landed weight cutoff (t)",side=1,line=2.4)





#Criteria for shoreside hake catch
### hakeCatch > 0.05:  50kg is probably about 100 hake
### Month < 5 or >11
### Year >=1991
hakeMin <- 0.05
pcatch$Month <- as.numeric(substring(pcatch$TDATE,6,7))

pcatch$HIX_HAKE <- FALSE
pcatch$HIX_HAKE[pcatch$gear%in%c("BottomTrawl","MidwaterTrawl") & pcatch$hakeCatch > hakeMin] <- TRUE  #because of looking weight of widow
#table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,pcatch$Month,useNA="ifany")
table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,useNA="ifany")

pcatch$HIX_HAKE[pcatch$Month<5] <- FALSE
pcatch$HIX_HAKE[pcatch$Month>11] <- FALSE
#table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,pcatch$Month,useNA="ifany")
table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,useNA="ifany")

pcatch$HIX_HAKE[pcatch$YEAR < 1991] <- FALSE
table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,useNA="ifany")
table(pcatch$HIX_HAKE,useNA="ifany")

pcatch$HIX_HAKE[pcatch$DAHL_HAKE & !(pcatch$HIX_HAKE)] <- TRUE   #put back in DAHL_SECTORS as the best guess
table(pcatch$HIX_HAKE,useNA="ifany")
table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,useNA="ifany")
pcatch$HIX_HAKE[pcatch$HIX_HAKE& !(pcatch$DAHL_HAKE)] <- FALSE   #put back in DAHL_SECTORS as the best guess
table(pcatch$HIX_HAKE,useNA="ifany")
table(pcatch$DAHL_HAKE,pcatch$HIX_HAKE,useNA="ifany")


tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$HIX_HAKE),sum)

pcatch$gear[pcatch$HIX_HAKE] <- "ShoresideHake"
save(pcatch,file="PacFIN/PacFinCatchWithShoresideHake.Rdat")

###Final catches from PacFIN
#remove Research catch (note that there is a space after the R)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$FLEET),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$FLEET),sum)[,"MidwaterTrawl",]
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$Month,pcatch$FLEET),sum)[,"MidwaterTrawl",,]

ind <- pcatch$FLEET != "R "   #190 obs of R
x <- tapply(pcatch$RMT[ind],list(pcatch$YEAR[ind],pcatch$gear[ind],pcatch$AGID[ind]),sum)
write.csv(x[,,"C"], file = "Catches/CA/PacFINcatches_CA.csv")
write.csv(x[,,"O"], file = "Catches/OR/PacFINcatches_OR.csv")
write.csv(x[,,"W"], file = "Catches/WA/PacFINcatches_WA.csv")


tmp <- pcatch
tmp$gear[tmp$gear%in%c("MiscTrawl","Other","Pot","ShrimpTrawl")] <- NA
xx <- tapply(tmp$RMT[ind],list(tmp$YEAR[ind],tmp$gear[ind],tmp$AGID[ind],tmp$IFQ_LANDING[ind]),sum)
xx
nNames <- unlist(lapply(dimnames(xx),length))
y <- data.frame(year  = rep(dimnames(xx)[[1]],prod(nNames[2:4])),
              fleet = rep(rep(dimnames(xx)[[2]],each=nNames[1]),prod(nNames[3:4])),
              state = rep(rep(dimnames(xx)[[3]],each=prod(nNames[1:2])),nNames[4]),
              ifq   = rep(dimnames(xx)[[4]],each=prod(nNames[1:3])),
              MT    = as.vector(xx)
)
write.csv(y,"Catches/CatchesBySectorStateIFQ.csv")


tmp <- pcatch
tmp$gear[tmp$gear%in%c("MiscTrawl","Other","Pot","ShrimpTrawl")] <- NA
xx <- tapply(tmp$RMT[ind],list(tmp$YEAR[ind],tmp$gear[ind],tmp$IFQ_LANDING[ind]),sum)
xx
nNames <- unlist(lapply(dimnames(xx),length))
y <- data.frame(year  = rep(dimnames(xx)[[1]],prod(nNames[2:3])),
              fleet = rep(rep(dimnames(xx)[[2]],each=nNames[1]),prod(nNames[3])),
              ifq   = rep(dimnames(xx)[[3]],each=prod(nNames[1:2])),
              MT    = as.vector(xx)
)
write.csv(y,"Catches/CatchesBySectorIFQ.csv")






















#did this for sablefish
xx <- tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$AGID,pcatch$IFQ_LANDING),sum)
xx
nNames <- unlist(lapply(dimnames(xx),length))
x <- data.frame(year  = rep(dimnames(xx)[[1]],prod(nNames[2:4])),
	          fleet = rep(rep(dimnames(xx)[[2]],each=nNames[1]),prod(nNames[3:4])),
	          state = rep(rep(dimnames(xx)[[3]],each=prod(nNames[1:2])),nNames[4]),
	          ifq   = rep(dimnames(xx)[[4]],each=prod(nNames[1:3])),
	          MT    = as.vector(xx)
)

write.csv(x,file="Other/CatchByIFQ.csv",row.names=F)




#for input to assessment
setwd("C:/NOAA2015/Hake/Data")
source("Rcode/functions/Functions.R")
load("extractedData/PacfinVdrfdCatch.Rdat")

pcatch <- pcatch[pcatch$FLEET != "XX",]

pcatch$Date <- as.Date(pcatch$TDATE)
pcatch$Month <- as.numeric(substr(pcatch$TDATE,6,7))
pcatch <- pcatch[order(pcatch$Date),]

pcatch$Sector <- "USshore"
pcatch$Sector[pcatch$FLEET=="R "] <- "USresearch"
#table(pcatch$FLEET,pcatch$Sector)

pcatch.yr <- tapply(pcatch$MT,list(pcatch$YEAR,pcatch$Sector),sum)
write.csv(pcatch.yr,file="Catches/USshoreCatchByYearVdrfd.csv")

#Report catch by period and year for US shore and research
pcatch.yr.per <- aggregate(pcatch$MT,list(pcatch$Month,pcatch$YEAR,pcatch$Sector),sum)
names(pcatch.yr.per) <- c("Period","Year","Sector","MT")
pcatch.yr.per <- pcatch.yr.per[,c("Sector","Period","Year","MT")]
write.csv(pcatch.yr.per,file="Catches/USshoreCatchByPeriodVdrfd.csv",row.names=F)


#Look at tribal catch in shoreside (already added in above)
tribal <- pcatch[pcatch$FLEET=="TI",]
tribal.yr <- tapply(tribal$MT,list(tribal$YEAR),sum)
#Report catch by period and year
tribal.yr.per <- aggregate(tribal$MT,list(tribal$Month,tribal$YEAR),sum)
names(tribal.yr.per) <- c("Period","Year","MT")
