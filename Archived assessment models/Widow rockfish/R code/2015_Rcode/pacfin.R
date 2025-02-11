setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
library(ExtractR)
#source("Rcode/functions/Functions.R")

#Get PacFIN catches and save results
#This is a new method, which includes research catch and breaks out tribal catch
#Remove XXX fleet (foreign catch?)
if(F) {
    pcatch <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/vdrfd.sql",db="PACFIN",uid="hicksa",sp="WDOW", start=1980, end=2014)
    save(pcatch,file="extractedData/PacfinVdrfdCatch.Rdat")
    file.copy("extractedData/PacfinVdrfdCatch.Rdat",paste("extractedData/PacfinVdrfdCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}

#start here with extraction"
setwd("C:/NOAA2015/Widow/Data")
load("extractedData/PacfinVdrfdCatch.Rdat")
table(pcatch$FLEET,pcatch$GRID,pcatch$IFQ_LANDING)
table(pcatch$GRID)
sort(tapply(pcatch$MT,list(pcatch$GRID),sum))
sort(tapply(pcatch$MT,list(pcatch$DAHL_SECTOR),sum))
tapply(pcatch$MT,list(pcatch$DAHL_SECTOR,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$DAHL_SECTOR),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID),sum)
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$GRID,pcatch$AGID),sum)


#Not sure where OTW (CA only) belongs
pcatch$gear <- NA
pcatch$gear[pcatch$GRID%in%c("SST","SHT","PWT","DST")] <- "ShrimpTrawl"
pcatch$gear[pcatch$GRID%in%c("RLT","GFT","GFS","GFL","FTS","FFT")] <- "BottomTrawl"
pcatch$gear[pcatch$GRID%in%c("OTW","MDT")] <- "MidwaterTrawl"  
pcatch$gear[pcatch$GRID%in%c("PRT","DNT","BMT")] <- "MiscTrawl"
pcatch$gear[pcatch$GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- "Pot"
pcatch$gear[pcatch$GRID%in%c("JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
pcatch$gear[pcatch$GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
pcatch$gear[pcatch$GRID%in%c("DVG","USP")] <- "Other"


table(pcatch$gear,useNA="ifany")
tapply(pcatch$MT,list(pcatch$YEAR,pcatch$DAHL_SECTOR,pcatch$gear),sum)

tapply(pcatch$MT,list(pcatch$YEAR,pcatch$gear),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear),sum)
tapply(pcatch$RMT,list(pcatch$YEAR,pcatch$gear,pcatch$AGID),sum)
#USE RMT which is round weight (landed*conversion factor)

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
