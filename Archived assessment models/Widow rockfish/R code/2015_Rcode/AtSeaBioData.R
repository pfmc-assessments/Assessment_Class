library(RODBC)

setwd("C:/NOAA2015/Widow/Data")
source("Rcode/functions/Functions.R")
#source("C:/NOAA2015/Widow/Data/Rcode/functions/workupPacFinTablesBDS.R")
source("C:/NOAA2015/Widow/Data/Rcode/functions/AtSeaComps.fn.R")
source("C:/NOAA2015/Widow/Data/Rcode/functions/Get.age.or.length.r")

#305             Widow ROCKFISH

if(F) {
    year1 <- 2008
    year2 <- 2014
    #Get NORPAC biological data and save results
    atsea.bio <- queryDB(queryFilename="atsea.detail",db="NORPAC",uid="hicksa",sp="305",start=year1,end=year2,querydir="sql/")
    save(atsea.bio,file="extractedData/atsea.bio.Rdat")
    file.copy("extractedData/atsea.bio.Rdat",paste("extractedData/atsea.bio_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

    year1 <- 1990
    year2 <- 2007
    #Get NORPAC biological data and save results for pre 2008
    atsea.bio.pre <- queryDB(queryFilename="pre.2008.atsea.length",db="NORPAC",uid="hicksa",sp="305",start=year1,end=year2,querydir="sql/")
    save(atsea.bio.pre,file="extractedData/atsea.bio.pre2008.Rdat")
    file.copy("extractedData/atsea.bio.pre2008.Rdat",paste("extractedData/atsea.bio.pre2008_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}

load(file="extractedData/atsea.bio.Rdat")
load(file="extractedData/atsea.bio.pre2008.Rdat")

table(substring(atsea.bio$RETRV_DATE_TIME,1,4))
table(substring(atsea.bio.pre$HAUL_DATE,1,4))

table(substring(atsea.bio$RETRV_DATE_TIME,1,4),atsea.bio$AGE)




#It appends to a file, so start new if necessary
Outfile.prefix = "Comps/AtSea.Length.Only/"
Report.prefix = "Comps/AtSea.Length.Only/"
Yearlist <- 2008:2014
for ( Yr in Yearlist ) {
        # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Lengths.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Lengths.report", sep="")

    tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

        Atsea.comps(tmp,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=FALSE,
              lbin.sizes=seq(14,70,2), in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=10, remove_sparse=FALSE,
              NO_LENGTH=FALSE, out.filename, report.filename)
}







load(file="extractedData/atsea.bio.pre2008.Rdat")
atsea.bio.pre$AGE <- NA

Outfile.postfix = ".for.SS.csv"
Report.postfix = ".report"
Yearlist = 2003:2007
in.season = NULL

#source("Utility.functions.r")
source("Rcode/functions/OldAtSeaComps.fn.R")

for ( Yr in Yearlist ) {
  # Setup filenames and call Old.Atsea.comps for length-comps.
  Outfile.prefix = "Comps/AtSea.Length.Only/"
  Report.prefix = "Comps/AtSea.Length.Only/"
  #Infile.prefix = "../Norpac_Annual/atsea.r.out."

  #in.filename = paste(Infile.prefix, Yr, Infile.postfix, sep="")
  out.filename = paste(Outfile.prefix, "Atsea.Lengths", Outfile.postfix, sep="")
  report.filename = paste(Report.prefix, Yr, ".Atsea.Lengths", Report.postfix, sep="")

  tmp <- atsea.bio.pre[atsea.bio.pre$YEAR==Yr,]
  tmp$RETRV_DATE_TIME <- paste(substring(tmp$HAUL_DATE,6,7),substring(tmp$HAUL_DATE,9,10),substring(tmp$HAUL_DATE,1,4),sep="/")

  Old.Atsea.comps(tmp,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=FALSE,
              lbin.sizes=seq(14,70,2), in.season=NULL, in.pctl=0.95, which="pct",
              min_Haul =0, min_T_weight=0, min_sample=15, remove_sparse=FALSE,
              NO_LENGTH=FALSE, out.filename, report.filename)

} # End for




Yearlist <- 2008:2014

Outfile.prefix = "Comps/AtSea.Age.Only/"
Report.prefix = "Comps/AtSea.Age.Only/"

for ( Yr in Yearlist ) {
    # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Ages.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Ages.report", sep="")

    tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

    Atsea.comps(tmp,BY_AGE=TRUE, BY_MONTH=FALSE, BY_GENDER=FALSE,
              lbin.sizes=NULL, in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=15, remove_sparse=FALSE,
              NO_LENGTH=TRUE, out.filename, report.filename)
}




###########################################################
## Age-at-length
library(nwfscSurveyCode)
updateSurveyCode()

setwd("C:\\NOAA2013\\Rougheye\\Data")
load("extractedData/atsea.bio.Rdat")

setwd("C:\\NOAA2013\\Rougheye\\Data\\Ages")

table(substring(atsea.bio$RETRV_DATE_TIME,1,4),!is.na(atsea.bio$AGE))
atsea.bio$Year <- substring(atsea.bio$RETRV_DATE_TIME,1,4)
reye <- atsea.bio[atsea.bio$Year%in%c(2008,2011),]
table(reye$Year,reye$AGE)
dim(reye)
reye <- reye[!is.na(reye$AGE),]
dim(reye)

x <- table(reye$Year,reye$HAUL_NUMBER)>0
apply(x,1,sum)



x <- as.data.frame(table(reye$Year,reye$LENGTH_SIZE,reye$AGE))
names(x) <- c("Year","Length","Age","AgeTallyF")
x <- x[x$AgeTallyF>0,]
x <- x[order(x$Year,x$Length,x$Age),]
x$AgeTallyM <- x$AgeTallyF
x$Length <- as.numeric(as.character(x$Length))
x$Year <- as.numeric(as.character(x$Year))
x$Age <- as.numeric(as.character(x$Age))

AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(10,80,2),ageBins=1:100,fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,"ForSS3\\AgeAtLenForSS3.AtSea.csv",row.names=F)






















































tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==2011,]
table(tmp$CRUISE,is.na(tmp$AGE))
hist(tmp[tmp$CRUISE==14420 & !is.na(tmp$AGE),"AGE"])
tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==2012,]
table(tmp$AGE,substring(tmp$RETRV_DATE_TIME,6,7))
table(tmp$AGE>0)

load("extractedData/atsea.bio.Rdat")
Yearlist <- 2008:2012

Outfile.prefix = "Comps/AtSea.Age.Only/"
Report.prefix = "Comps/AtSea.Age.Only/"

for ( Yr in Yearlist ) {
    # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Ages.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Ages.report", sep="")

    tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

    Atsea.comps(tmp,BY_AGE=TRUE, BY_MONTH=FALSE, BY_GENDER=FALSE,
              lbin.sizes=NULL, in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=15, remove_sparse=FALSE,
              NO_LENGTH=TRUE, out.filename, report.filename)
}




### By Month
Outfile.prefix = "Comps/AtSea.Age.Only.Month/"
Report.prefix = "Comps/AtSea.Age.Only.Month/"
Yearlist <- 2008:2012
for ( Yr in Yearlist ) {
        # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Ages.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Ages.report", sep="")

    tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

        Atsea.comps(tmp,BY_AGE=TRUE, BY_MONTH=TRUE, BY_GENDER=FALSE,
              lbin.sizes=NULL, in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=15, remove_sparse=FALSE,
              NO_LENGTH=TRUE, out.filename, report.filename)
}


Outfile.prefix = "Comps/AtSea.Length.Only.Month/"
Report.prefix = "Comps/AtSea.Length.Only.Month/"
Yearlist <- 2011:2012
for ( Yr in Yearlist ) {
        # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Lengths.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Lengths.report", sep="")

    tmp <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,4)==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

        Atsea.comps(tmp,BY_AGE=FALSE, BY_MONTH=TRUE, BY_GENDER=FALSE,
              lbin.sizes=NULL, in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=15, remove_sparse=FALSE,
              NO_LENGTH=TRUE, out.filename, report.filename)
}


lens <- read.csv("Comps/AtSea.Length.Only.Month/AtSeaLengthCompsMonth2011.csv")
lens <- read.csv("Comps/AtSea.Length.Only.Month/AtSeaLengthCompsMonth2012.csv")
par(mfcol=c(4,2),oma=c(1.5,1,1.5,0),mar=c(4,3,1,1))
for(i in 1:8) {
    barplot(as.numeric(lens[i,-c(1,2)]),main=lens[i,1],names=seq(20,70,2))
}
mtext("Length (cm)",side=1,outer=T)
mtext("Proportion",side=2,outer=T)
mtext("At-sea weighted length compositions",outer=T,side=3)


tmp1 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-10",]
tmp2 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-11",]
tmp3 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-12",]

par(mfrow=c(3,1),oma=c(0,0,1,0))
hist(tmp1$LENGTH_SIZE,main="October",xlim=c(0,80),xlab="")
hist(tmp2$LENGTH_SIZE,main="November",xlim=c(0,80),xlab="")
hist(tmp2$LENGTH_SIZE,main="December",xlim=c(0,80),xlab="Length (cm)")
mtext("At-sea raw length observations",outer=T,side=3)






tmp1 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-05",]
tmp2 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-06",]
tmp3 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-07",]
tmp4 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-08",]
tmp5 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-09",]
tmp6 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-10",]
tmp7 <- atsea.bio[substring(atsea.bio$RETRV_DATE_TIME,1,7)=="2011-11",]

par(mfrow=c(4,2),oma=c(0,0,1,0))
barplot(tmp1$AGE,main="May",xlim=c(0,15),xlab="")
hist(tmp2$AGE,main="June",xlim=c(0,15),xlab="")
hist(tmp3$AGE,main="July",xlim=c(0,15),xlab="")
hist(tmp4$AGE,main="August",xlim=c(0,15),xlab="")
hist(tmp5$AGE,main="September",xlim=c(0,15),xlab="")
hist(tmp6$AGE,main="October",xlim=c(0,15),xlab="")
hist(tmp7$AGE,main="November",xlim=c(0,15),xlab="")
mtext("At-sea raw age observations",outer=T,side=3)





###########################
#At-sea locations
library(nwfscMapping)
library(PBSmapping)
data(westCoastLL)
data(WCstatesInlandPBS)
source("C:/Mapping/WestCoastMapping.R")


atsea.bio$lon <- -1*(atsea.bio$RETRV_LONGITUDE_DEGREES+atsea.bio$RETRV_LONGITUDE_MINUTES/60)
atsea.bio$lat <- atsea.bio$RETRV_LATITUDE_DEGREES+atsea.bio$RETRV_LATITUDE_MINUTES/60
doPNG <- F
ht <- 6; wd<- 10
if(doPNG) {png(filename = paste(figDir,"mapREYEnBSPR.png",sep="\\"), width = wd, height = ht,units="in",res=300, pointsize = 11)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,1),mar=c(4,4,1,1)+0.1,las=1)
plotMap(westCoastLL, tck = c(-0.02), xlim=c(-130,-116.5), ylim=c(42,49.1),col="darkgreen",bg="lightblue")
addLines(westCoastLL)
addLines(WCstatesInlandPBS)
addLines(EEZ2,lwd=2)
addLines(depth1500m, col = gray(0.5))
points(atsea.bio$lon,atsea.bio$lat,pch=20,col=rgb(1,0,0,0.1))
#points(bspr$BEST_LON_DD,bspr$BEST_LAT_DD,pch=20,col=rgb(0,0,0,0.8))
if(doPNG) {dev.off()}
