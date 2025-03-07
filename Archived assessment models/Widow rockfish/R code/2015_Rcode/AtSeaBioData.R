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

atsea.bio.V <- read.csv("Biological/Widow_atsea_lengths_2000-2014.csv")

load(file="extractedData/atsea.bio.pre2008.Rdat")
atsea.bio.pre$AGE <- NA #no ages prior to 2008
table(atsea.bio.pre$YEAR)
table(atsea.bio.V$YEAR)
#I'll just use the pre-2008 because it is very close to the numbers that Vanessa supplied.
#Vanessa's extraction doesn't have sample weight, so I can't expand

table(atsea.bio.pre$LENGTH,atsea.bio.pre$SEX)
#unsexed fish are relatively few and all over, so I won't worry about them

Outfile.postfix = ".for.SS.csv"
Report.postfix = ".report"
Yearlist = 1992:2007
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

  Old.Atsea.comps(tmp,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=TRUE,
              lbin.sizes=seq(8,56,2), in.season=NULL, in.pctl=0.95, which="pct",
              min_Haul =0, min_T_weight=0, min_sample=5, remove_sparse=FALSE,
              NO_LENGTH=FALSE, out.filename, report.filename)

} # End for




atsea.bio.V <- read.csv("Biological/Widow_atsea_lengths_2000-2014.csv")
atsea.bio.V <- atsea.bio.V[atsea.bio.V$YEAR>=2008,]
load("extractedData/atsea.bio.Rdat")
atsea.bio$Year <- as.numeric(substr(atsea.bio$RETRV_DATE_TIME,1,4))

head(atsea.bio.V)
head(atsea.bio)

atsea.bio[atsea.bio$Year==2009 & atsea.bio$HAUL_NUMBER==1,]
atsea.bio[atsea.bio$Year==2008 & atsea.bio$CRUISE==11889 & atsea.bio$HAUL_NUMBER==1,]

dat <- atsea.bio[is.na(atsea.bio$BARCODE),]
table(dat$Year)
table(atsea.bio.V$YEAR)

#dat doesn't match exactly with Vanessa's, but is close

table(atsea.bio$LENGTH,atsea.bio$SEX)
#unsexed fish are few and all over, so I won't worry about them

#It appends to a file, so start new if necessary
Outfile.prefix = "Comps/AtSea.Length.Only/"
Report.prefix = "Comps/AtSea.Length.Only/"
Yearlist <- 2008:2014
for ( Yr in Yearlist ) {
        # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Lengths.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Lengths.report", sep="")

    tmp <- dat[dat$Year==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

        Atsea.comps(tmp,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=TRUE,
              lbin.sizes=seq(8,56,2), in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=5, remove_sparse=FALSE,
              NO_LENGTH=FALSE, out.filename, report.filename)
}


dat1 <- dat[!is.na(dat$AGE) & dat$Year<2011,]
sum(duplicated(dat1$SPECIMEN_NUMBER))
dat1 <- dat1[!duplicated(dat1$SPECIMEN_NUMBER),]
dat2 <- atsea.bio[!is.na(atsea.bio$AGE) & atsea.bio$Year>=2011,]
sum(duplicated(dat2$SPECIMEN_NUMBER))
dat3 <- rbind(dat1,dat2)
table(dat1$Year)
table(dat2$Year)
table(dat3$Year)

table(is.na(dat3$AGE))
table(dat3$Year,is.na(dat3$WEIGHT))

lw <- data.frame(a=1.7043e-5,b=2.9668)
calcWt <- lw$a*dat3$LENGTH^lw$b
dat3$WEIGHT[is.na(dat3$WEIGHT)] <- calcWt[is.na(dat3$WEIGHT)]
dat3$SEX <- as.character(dat3$SEX)

Yearlist <- 2008:2014  #use dat3
Outfile.prefix = "Comps/AtSea.Age.Only/"
Report.prefix = "Comps/AtSea.Age.Only/"
for ( Yr in Yearlist ) {
    # Setup filenames and call Shore.comps for age comps.
    out.filename = paste(Outfile.prefix, "Atsea.Ages.for.SS.csv", sep="")
    report.filename = paste(Report.prefix, Yr, ".Atsea.Ages.report", sep="")

    tmp <- dat3[dat3$Year==Yr,]
    tmp$RETRV_DATE_TIME <- paste(substring(tmp$RETRV_DATE_TIME,6,7),substring(tmp$RETRV_DATE_TIME,9,10),substring(tmp$RETRV_DATE_TIME,1,4),sep="/")

    Atsea.comps(tmp,BY_AGE=TRUE, BY_MONTH=FALSE, BY_GENDER=TRUE,
              lbin.sizes=NULL, in.season=NULL, which="pct",
              min_Haul=0, min_T_weight=0, in.pctl=0.95,
              min_sample=5, remove_sparse=FALSE,
              NO_LENGTH=TRUE, out.filename, report.filename)
}




###########################################################
## Age-at-length
# devtools::install_github("nwfsc-assess/nwfscSurvey")
library(nwfscSurvey)

setwd("C:/NOAA2015/Widow/Data")
load("extractedData/atsea.bio.Rdat")
atsea.bio$Year <- as.numeric(substr(atsea.bio$RETRV_DATE_TIME,1,4))
dat <- atsea.bio[is.na(atsea.bio$BARCODE),]
dat1 <- dat[!is.na(dat$AGE) & dat$Year<2011,]
dat1 <- dat1[!duplicated(dat1$SPECIMEN_NUMBER),]
dat2 <- atsea.bio[!is.na(atsea.bio$AGE) & atsea.bio$Year>=2011,]
dat3 <- rbind(dat1,dat2)
table(dat3$Year)

x <- table(dat3$Year,paste(dat3$CRUISE,dat3$HAUL_NUMBER))>0
apply(x,1,sum)

table(dat3$SEX)
dat3$SEX <- as.character(dat3$SEX)
dat3 <- dat3[dat3$SEX!="U",]
table(dat3$SEX)

x <- as.data.frame(table(dat3$Year,dat3$SEX,dat3$LENGTH_SIZE,dat3$AGE))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]

#maybe I should add in shoreside here


AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,"Comps/AgeAtLenForSS3.AtSea.Female.csv",row.names=F)
write.csv(AatL$male,"Comps/AgeAtLenForSS3.AtSea.Male.csv",row.names=F)






















































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
