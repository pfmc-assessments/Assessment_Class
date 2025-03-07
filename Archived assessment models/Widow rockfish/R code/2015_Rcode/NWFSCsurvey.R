dont source

# devtools::install_github("nwfsc-assess/nwfscMapping")
# devtools::install_github("nwfsc-assess/nwfscSurvey")
# Load package
library(nwfscMapping)
library(nwfscSurvey)
library(PBSmapping)
library(date)
data(westCoastLL)
data(WCstatesInlandPBS)
source("C:/NOAA2015/Widow/Data/Rcode/Functions/WestCoastMapping.R")  #it will error out if you don't have depth shapefiles, but some other useful things will be created
plotTheCoast <- function(xlims,ylims,plotDepth=F) {
    plotMap(westCoastLL,xlim=xlims,ylim=ylims,col="green4",bg="lightcyan2",tck=c(-0.02),cex=1.0,main="",xlab="",ylab="")
    addLines(WCstatesInlandPBS)
    #addLines(EEZ,lwd=2,col=gray(0.5))
    if(plotDepth) {
    	addLines(depth30f,col=gray(0.3))
    	addLines(depth100f,col=gray(0.3))
    }
}
f2m <- 1.8288  #fathoms to meters

wtByBiomass.fn <- function(x) {
	xx <- x[,-c(1:8)]
	xx <- t(apply(xx,1,function(y){y[length(y)] * y[-length(y)]/sum(y[-length(y)])}))
	xx <- apply(xx,2,function(y,yr){tapply(y,yr,sum)},yr=x$Year)
	cat("Check for observations in F.999 or M.999\n")
#	return(data.frame(x[!duplicated(x$Year),c("Year","Season","Fleet","gender","partition","nSamps")], xx))
	return(xx)
}
calcMeanWt <- function(x,lenWt) {
    #x is a dataframe with F and M headers for length
    fem <- x[grep("F[123456789]",names(x))]  #this doesn't select the F.999, so be careful
    mal <- x[grep("M[123456789]",names(x))]
    tmp.fn <- function(xx,ab) {
        len <- as.numeric(substring(names(xx),2))
        sum(xx * ab$a*len^ab$b)/100  # because I put them in percentages
    }  
    mnWtFem <- apply(fem,1,tmp.fn,ab=lenWt$female)
    mnWtMal <- apply(mal,1,tmp.fn,ab=lenWt$male)
    return(mnWtFem + mnWtMal)  #they are proportions that sum to one, so simply add it up
}

effN <- function(nSample,nFish,source=c("fishery","survey")) {
    if(length(nSample)!=length(nFish)) {stop("Vectors must be same length")}
    if(source=="fishery") {
        nFishCoeff <- 0.138
        nSampleCoeff <- 7.06
        breakPt <- 44
    }
    if(source=="survey") {
        nFishCoeff <- 0.0707
        nSampleCoeff <- 4.89
        breakPt <- 55
    }
    nEff <- rep(NA,length(nSample))
    ind <- nFish/nSample < breakPt
    ind[is.na(ind)] <- TRUE
    nEff[ind] <- nSample[ind] + nFishCoeff*nFish[ind]
    ind <- !ind
    nEff[ind] <- nSampleCoeff * nSample[ind]
    names(nEff) <- names(nFish)
    return(nEff)
}

setwd("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\NWFSCsurvey")





bioSS <- GetTotalBiomass.fn("BiomassAbundance.csv",headerRow=6)
plotBio.fn(bioSS,pch=16)
title(main="Design-based Biomass from NWFSC trawl survey")

lfs <- readInLengthComps.fn("LengthComps.csv",headerRow=7)
range(lfs$Length)
barplot(table(factor(lfs$Length,levels=6:55)),main="All Widow Lengths from NWFSC survey")
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,nSamps="Enter sample sizes",season=1,partition=0,NAs2zero=T)
lfsSS$F8 <- lfsSS$F8 + lfsSS$F.999
lfsSS$M8 <- lfsSS$M8 + lfsSS$M.999
lfsSS <- lfsSS[,-which(names(lfsSS)%in%c("F.999","M.999"))]

par(mfrow=c(2,1))
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

ratioF <- lfs$NumF/(lfs$NumF+lfs$NumM)
boxplot(split(ratioF,floor(lfs$Length)),ylab="Prop F",main="Sex ratios by length")
tapply(lfs$NumUnsexed,lfs$Length,sum)

doPNG <- T
wd<-6.5;ht<-3.5
if(doPNG) {png("../../../Writeup/Figures/nwfscSexRatios.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,1),mar=c(4,4,1,1))
plotSexRatio.fn(lfs,ylim=c(0,1),las=1)    #circle size is proportional to the number of observations (strata-year-length specific)
abline(h=0.5,col=rgb(0,0,0,0.3))
#abline(v=28)
if(doPNG) dev.off()
#Compare correcting for unsexed fish and adding them back in
par(mfcol=c(2,2))
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")
#correct for unsexed fish
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,NAs2zero=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

#just for giggles, look at unsexed comps
#use gender =0 to get unsexed (including male and female) proportions
par(mfcol=c(2,2))
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")
lfsSS  <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=0,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

tmpLen <- data.frame(lfs,yearAreaDepth=paste(lfs$Year,lfs$AreaName,lfs$MinStratumDepth))
lfsStrata <- SS3LFstrata.fn(tmpLen,Strata="yearAreaDepth",lgthBins=seq(8,56,2),gender=3)
lfsStrata$Strata <- substring(lfsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
lfsStrata$Strata <- strNames[match(lfsStrata$Strata,strId)]
lfsStrata <- lfsStrata[!is.na(lfsStrata$Strata),]
table(lfsStrata$Strata)
table(nwfscGlmStrata$Strata)
lfsStrata$Year <- as.numeric(substring(lfsStrata$year,1,4))

#link in biomass to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)

raw <- tapply(nwfscGlmStrata$Raw,nwfscGlmStrata$Year,sum)
raw/bioSS$Value
plotBio.fn(bioSS,pch=16)
lines(as.numeric(names(raw)),raw/1e6,col="blue",type="b")
x <- merge(lfsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))

lfsBio <- wtByBiomass.fn(x)

lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,NAs2zero=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
names(lfsBio) <- names(lfsSS)

par(mfrow=c(2,2))
plotFreqData.fn(lfsBio,ylim=c(0,55),yaxs="i",ylab="Length")
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

lfsBio$F8 <- lfsBio$F8 + lfsBio$F.999
lfsBio$M8 <- lfsBio$M8 + lfsBio$M.999
lfsBio <- lfsBio[,-which(names(lfsBio)%in%c("F.999","M.999"))]


options(digits=19)
lens <- read.csv("WidowLengths.csv",skip=8)
nSamp <- table(lens$PROJECT_CYCLE,!duplicated(lens$HAUL_IDENTIFIER))[,"TRUE"]
nFish <- table(lens$PROJECT_CYCLE)
n <- floor(effN(nSamp,nFish,"survey"))
###################################################################33
#write out lfs when ready
###Weighted by area
lfs <- readInLengthComps.fn("LengthComps.csv",headerRow=7)

lfsSS <- SS3LF.fn(lfs,lgthBins=seq(8,56,2),gender=3,NAs2zero=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28,fleet=8,partition=0,nSamps=n)
lfsSS$F8 <- lfsSS$F8 + lfsSS$F.999
lfsSS$M8 <- lfsSS$M8 + lfsSS$M.999
lfsSS <- lfsSS[,-which(names(lfsSS)%in%c("F.999","M.999"))]
lfsSS[,-c(1:6)] <- round(lfsSS[,-c(1:6)],2)
write.csv(lfsSS,"LengthCompsForSS3.NWFSC.csv",row.names=F)
###WEighted by strata biomass
tmpLen <- data.frame(lfs,yearAreaDepth=paste(lfs$Year,lfs$AreaName,lfs$MinStratumDepth))
lfsStrata <- SS3LFstrata.fn(tmpLen,Strata="yearAreaDepth",lgthBins=seq(8,56,2),gender=3)
lfsStrata$Strata <- substring(lfsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
lfsStrata$Strata <- strNames[match(lfsStrata$Strata,strId)]
lfsStrata <- lfsStrata[!is.na(lfsStrata$Strata),]
lfsStrata$Year <- as.numeric(substring(lfsStrata$year,1,4))
#link in biomass to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(lfsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
lfsBio <- wtByBiomass.fn(x)
lfsBio <- t(apply(lfsBio,1,function(x){round(100*x/sum(x),2)}))
lfsBio <- data.frame(lfsSS[c("year","Season","Fleet","gender","partition","nSamps")], lfsBio)
lfsBio$F8 <- lfsBio$F8 + lfsBio$F.999
lfsBio$M8 <- lfsBio$M8 + lfsBio$M.999
lfsBio <- lfsBio[,-which(names(lfsBio)%in%c("F.999","M.999"))]
write.csv(lfsBio,"LengthCompsForSS3.NWFSC.Bio.csv",row.names=F)


#################################
###WEighted by strata numbers (as suggested at 2015 STAR panel)
#### use properly exapnded length comp and length-weight relationship
tmpLen <- data.frame(lfs,yearAreaDepth=paste(lfs$Year,lfs$AreaName,lfs$MinStratumDepth))
lfsStrata <- SS3LFstrata.fn(tmpLen,Strata="yearAreaDepth",lgthBins=seq(8,56,2),gender=3)
lfsStrata$Strata <- substring(lfsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
lfsStrata$Strata <- strNames[match(lfsStrata$Strata,strId)]
lfsStrata <- lfsStrata[!is.na(lfsStrata$Strata),]
lfsStrata$Year <- as.numeric(substring(lfsStrata$year,1,4))
#link in biomass and mean weight to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(lfsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
#change biomass to total numbers by dividing by mean weight
lenWt <- list(female=data.frame(a=0.00001736,b=2.962),male=data.frame(a=0.00001484,b=3.005))
meanwt <- calcMeanWt(lfsStrata,lenWt=lenWt)
x$IndexMedian <- x$IndexMedian/meanwt
lfsNum <- wtByBiomass.fn(x)
lfsNum <- t(apply(lfsNum,1,function(x){round(100*x/sum(x),2)}))
lfsNum <- data.frame(lfsSS[c("year","Season","Fleet","gender","partition","nSamps")], lfsNum)
lfsNum$F8 <- lfsNum$F8 + lfsNum$F.999
lfsNum$M8 <- lfsNum$M8 + lfsNum$M.999
lfsNum <- lfsNum[,-which(names(lfsNum)%in%c("F.999","M.999"))]
write.csv(lfsNum,"LengthCompsForSS3.NWFSC.Num.csv",row.names=F)
####################################################################

options(digits=19)
lens <- read.csv("WidowLengths.csv",skip=8)
n <- table(lens$PROJECT_CYCLE,!duplicated(lens$HAUL_IDENTIFIER))[,"TRUE"]
###################################################################33
#write out lfs when ready
###    FOR BRIDGING OF 2011 ASSESSMENT
###  NEED TO MANUALLY CALCULATE THEM
###Weighted by area
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(10,64,2),gender=3,NAs2zero=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28,fleet=8,partition=0,nSamps=n)
lfsSS$F8 <- lfsSS$F8 + lfsSS$F.999
lfsSS$M8 <- lfsSS$M8 + lfsSS$M.999
lfsSS <- lfsSS[,-which(names(lfsSS)%in%c("F.999","M.999"))]
lfsSS[,-c(1:6)] <- round(lfsSS[,-c(1:6)],2)
write.csv(lfsSS,"LengthCompsForSS3.NWFSC.2011structure.csv",row.names=F)
###WEighted by strata biomass
tmpLen <- data.frame(lfs,yearAreaDepth=paste(lfs$Year,lfs$AreaName,lfs$MinStratumDepth))
lfsStrata <- SS3LFstrata.fn(tmpLen,Strata="yearAreaDepth",lgthBins=seq(8,56,2),gender=3)
lfsStrata$Strata <- substring(lfsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
lfsStrata$Strata <- strNames[match(lfsStrata$Strata,strId)]
lfsStrata <- lfsStrata[!is.na(lfsStrata$Strata),]
lfsStrata$Year <- as.numeric(substring(lfsStrata$year,1,4))
#link in biomass to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(lfsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
lfsBio <- wtByBiomass.fn(x)
lfsBio <- t(apply(lfsBio,1,function(x){round(100*x/sum(x),2)}))
lfsBio <- data.frame(lfsSS[c("year","Season","Fleet","gender","partition","nSamps")], lfsBio)
lfsBio$F8 <- lfsBio$F8 + lfsBio$F.999
lfsBio$M8 <- lfsBio$M8 + lfsBio$M.999
lfsBio <- lfsBio[,-which(names(lfsBio)%in%c("F.999","M.999"))]
write.csv(lfsBio,"LengthCompsForSS3.NWFSC.Bio.2011structure.csv",row.names=F)
####################################################################


#################################################################################################
#### AgeComps
options(digits=19)
ages <- read.csv("Widow_SexedLgthWtAge.csv",skip=8)
nSamp <- table(ages$PROJECT_CYCLE,!duplicated(ages$HAUL_IDENTIFIER))[,"TRUE"]
nFish <- table(ages$PROJECT_CYCLE)
n <- floor(effN(nSamp,nFish,"survey"))
table(ages$AGE_YRS,ages$SEX)

afs <- readInAgeComps.fn("AgeComps.csv",headerRow=7)
afs$Length <- afs$Age
afsSS  <- SS3LF.fn(afs,lgthBins=seq(0,40,1),gender=3,NAs2zero=T,
	                sexRatioUnsexed=0.5,maxSizeUnsexed=28,fleet=8,partition=0,nSamps=n)
par(mfrow=c(2,1))
plotFreqData.fn(afsSS,ylim=c(0,40),yaxs="i",ylab="Age",las=1)
plotFreqData.fn(afsSS,ylim=c(0,10),yaxs="i",ylab="Age",las=1)

afsSS <- afsSS[,-which(names(afsSS)%in%c("F.999","M.999"))]
afsSS[,-c(1:6)] <- round(afsSS[,-c(1:6)],2)
write.csv(afsSS,"AgeCompsForSS3.NWFSC.csv",row.names=F)

tmp <- afs
tmp$Length <- tmp$Age
plotSexRatio.fn(tmp,ylim=c(0,1),las=1)

#Age comps weighted by srata
tmpAge <- data.frame(afs,yearAreaDepth=paste(afs$Year,afs$AreaName,afs$MinStratumDepth))
tmpAge$Length <- tmpAge$Age
afsStrata <- SS3LFstrata.fn(tmpAge,Strata="yearAreaDepth",lgthBins=seq(0,40,1),gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
afsStrata$Strata <- substring(afsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
afsStrata$Strata <- strNames[match(afsStrata$Strata,strId)]
afsStrata <- afsStrata[!is.na(afsStrata$Strata),]
afsStrata$Year <- as.numeric(substring(afsStrata$year,1,4))
#link in biomass to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(afsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
afsBio <- wtByBiomass.fn(x)
afsBio <- t(apply(afsBio,1,function(x){round(100*x/sum(x),2)}))
afsBio <- data.frame(afsSS[c("year","Season","Fleet","gender","partition","nSamps")], afsBio)
afsBio <- afsBio[,-which(names(afsBio)%in%c("F.999","M.999"))]
write.csv(afsBio,"AgeCompsForSS3.NWFSC.Bio.csv",row.names=F)

#################################
###WEighted by strata numbers (as suggested at 2015 STAR panel)
#### use properly exapnded length comp and length-weight relationship
tmpAge <- data.frame(afs,yearAreaDepth=paste(afs$Year,afs$AreaName,afs$MinStratumDepth))
tmpAge$Length <- tmpAge$Age
afsStrata <- SS3LFstrata.fn(tmpAge,Strata="yearAreaDepth",lgthBins=seq(0,40,1),gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
afsStrata$Strata <- substring(afsStrata$year,6)
strId <- c("32To34Pt5Degs 182.88","32To34Pt5Degs 54.864","34Pt5To40Pt5Degs 182.88","34Pt5To40Pt5Degs 54.864","40Pt5To49Degs 182.88","40Pt5To49Degs 54.864")
strNames <- c(NA,NA,"Sdeep","Sshallow","Ndeep","Nshallow")
afsStrata$Strata <- strNames[match(afsStrata$Strata,strId)]
afsStrata <- afsStrata[!is.na(afsStrata$Strata),]
afsStrata$Year <- as.numeric(substring(afsStrata$year,1,4))
#link in biomass and mean weight to weight each strata
nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(afsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
#change biomass to total numbers by dividing by mean weight
lenWt <- list(female=data.frame(a=0.00001736,b=2.962),male=data.frame(a=0.00001484,b=3.005))
tmpLfs <- lfsStrata[match(afsStrata[,1],lfsStrata[,1]),]
meanwt <- calcMeanWt(tmpLfs,lenWt=lenWt)  #USE lfs
x$IndexMedian <- x$IndexMedian/meanwt
afsNum <- wtByBiomass.fn(x)
afsNum <- t(apply(afsNum,1,function(x){round(100*x/sum(x),2)}))
afsNum <- data.frame(afsSS[c("year","Season","Fleet","gender","partition","nSamps")], afsNum)
afsNum$F8 <- afsNum$F8 + afsNum$F.999
afsNum$M8 <- afsNum$M8 + afsNum$M.999
afsNum <- afsNum[,-which(names(afsNum)%in%c("F.999","M.999"))]
write.csv(afsNum,"AgeCompsForSS3.NWFSC.Num.csv",row.names=F)
####################################################################


#Age at length by strata
age <- readInAgeComps.fn("AgeComps.csv",headerRow=7,sep=",")

ageNshallow <- age[age$AreaName=="40Pt5To49Degs" & age$MinStratumDepth<55,]
AatL.NS <- SS3AgeAtLen.fn(ageNshallow,lgthBins=seq(8,56,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=F,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
AatL.NS$female$Strata <- "Nshallow"
AatL.NS$male$Strata <- "Nshallow"
ageNdeep <- age[age$AreaName=="40Pt5To49Degs" & age$MinStratumDepth>55,]
AatL.ND <- SS3AgeAtLen.fn(ageNdeep,lgthBins=seq(8,56,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=F,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
AatL.ND$female$Strata <- "Ndeep"
AatL.ND$male$Strata <- "Ndeep"
ageSshallow <- age[age$AreaName=="34Pt5To40Pt5Degs" & age$MinStratumDepth<55,]
AatL.SS <- SS3AgeAtLen.fn(ageSshallow,lgthBins=seq(8,56,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=F,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
AatL.SS$female$Strata <- "Sshallow"
AatL.SS$male$Strata <- "Sshallow"
ageSdeep <- age[age$AreaName=="34Pt5To40Pt5Degs" & age$MinStratumDepth>55,]
AatL.SD <- SS3AgeAtLen.fn(ageSdeep,lgthBins=seq(8,56,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=F,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
AatL.SD$female$Strata <- "Sdeep"
AatL.SD$male$Strata <- "Sdeep"

nwfscGlmStrata <- read.csv("NWFSCIndexByStrata.csv")
nwfscGlmStrata$Strata <- as.character(nwfscGlmStrata$Strata)
x <- merge(lfsStrata,nwfscGlmStrata[,c("Year","Strata","IndexMedian")],by=c("Year","Strata"))
#change biomass to total numbers by dividing by mean weight
lenWt <- list(female=data.frame(a=0.00001736,b=2.962),male=data.frame(a=0.00001484,b=3.005))
meanwt <- calcMeanWt(lfsStrata,lenWt=lenWt)
x$Numbers <- x$IndexMedian/meanwt
x$year <- x$Year

AatL.NS$female <- merge(AatL.NS$female,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.NS$male <- merge(AatL.NS$male,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.ND$female <- merge(AatL.ND$female,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.ND$male <- merge(AatL.ND$male,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.SS$female <- merge(AatL.SS$female,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.SS$male <- merge(AatL.SS$male,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.SD$female <- merge(AatL.SD$female,x[,c("year","Strata","Numbers")],by=c("year","Strata"))
AatL.SD$male <- merge(AatL.SD$male,x[,c("year","Strata","Numbers")],by=c("year","Strata"))

wtByBiomassAatL.fn <- function(x,ages=0:40,sex="F") {
    xx <- x[,c(paste0(sex,ages),"Numbers")]
    xx <- t(apply(xx,1,function(y){y[length(y)] * y[-length(y)]/sum(y[-length(y)])}))
    #xx <- apply(xx,2,function(y,yr){tapply(y,yr,sum)},yr=x$Year)
    cat("Check for observations in F.999 or M.999\n")
#   return(data.frame(x[!duplicated(x$Year),c("Year","Season","Fleet","gender","partition","nSamps")], xx))
    return(data.frame(x[,1:10],xx,xx))
}

AatL.NS$female <- wtByBiomassAatL.fn(AatL.NS$female)
AatL.NS$male <- wtByBiomassAatL.fn(AatL.NS$male,sex="M")
AatL.ND$female <- wtByBiomassAatL.fn(AatL.ND$female)
AatL.ND$male <- wtByBiomassAatL.fn(AatL.ND$male,sex="M")
AatL.SS$female <- wtByBiomassAatL.fn(AatL.SS$female)
AatL.SS$male <- wtByBiomassAatL.fn(AatL.SS$male,sex="M")
AatL.SD$female <- wtByBiomassAatL.fn(AatL.SD$female)
AatL.SD$male <- wtByBiomassAatL.fn(AatL.SD$male,sex="M")

ages <- 0:40
lens <- seq(8,56,2)
yrs <- 2003:2014
des <- expand.grid(yrs,lens)
names(des) <- c("year","LbinLo")
# tmp <- matrix(NA,ncol=2*length(ages),nrow=nrow(des))
# AatLfemale <- AatLmale <- data.frame(des,tmp)
# names(AatLfemale) <- names(AatLmale) <- c("year","LbinLo",paste0("F",ages),paste0("M",ages))

AatL.NS$female <- merge(AatL.NS$female,des,by=c("year","LbinLo"),all=T)
AatL.NS$female[is.na(AatL.NS$female)] <- 0
AatL.NS$male <- merge(AatL.NS$male,des,by=c("year","LbinLo"),all=T)
AatL.NS$male[is.na(AatL.NS$male)] <- 0
AatL.ND$female <- merge(AatL.ND$female,des,by=c("year","LbinLo"),all=T)
AatL.ND$female[is.na(AatL.ND$female)] <- 0
AatL.ND$male <- merge(AatL.ND$male,des,by=c("year","LbinLo"),all=T)
AatL.ND$male[is.na(AatL.ND$male)] <- 0
AatL.SS$female <- merge(AatL.SS$female,des,by=c("year","LbinLo"),all=T)
AatL.SS$female[is.na(AatL.SS$female)] <- 0
AatL.SS$male <- merge(AatL.SS$male,des,by=c("year","LbinLo"),all=T)
AatL.SS$male[is.na(AatL.SS$male)] <- 0
AatL.SD$female <- merge(AatL.SD$female,des,by=c("year","LbinLo"),all=T)
AatL.SD$female[is.na(AatL.SD$female)] <- 0
AatL.SD$male <- merge(AatL.SD$male,des,by=c("year","LbinLo"),all=T)
AatL.SD$male[is.na(AatL.SD$male)] <- 0

tmp <- list(female=NULL,male=NULL)
tmp$female <- AatL.NS$female[,paste0("F",ages)] + AatL.ND$female[,paste0("F",ages)] + 
              AatL.SS$female[,paste0("F",ages)] + AatL.SD$female[,paste0("F",ages)] 
tmp$female <- t(apply(tmp$female,1,function(x){100*x/sum(x)}))
tmp$male <- AatL.NS$male[,paste0("M",ages)] + AatL.ND$male[,paste0("M",ages)] + 
              AatL.SS$male[,paste0("M",ages)] + AatL.SD$male[,paste0("M",ages)] 
tmp$male <- t(apply(tmp$male,1,function(x){100*x/sum(x)}))

AatL$female <- data.frame(AatL.SD$female[,c("year","Season","Fleet","gender","partition","ageErr","LbinLo","LbinLo","nSamps")],
                           tmp$female,tmp$female)
AatL$female <- AatL$female[apply(AatL$female,1,function(x){!all(is.na(x[-(1:10)]))}),]
AatL$male <- data.frame(AatL.SD$male[,c("year","Season","Fleet","gender","partition","ageErr","LbinLo","LbinLo","nSamps")],
                           tmp$male,tmp$male)
AatL$male <- AatL$male[apply(AatL$male,1,function(x){!all(is.na(x[-(1:10)]))}),]

# AatL$female$Season <- AatL$male$Season <- 1
# AatL$female$Fleet <- AatL$male$Fleet <- 8
# AatL$female$gender <- 1;  AatL$male$gender <- 2
# AatL$female$partition <- AatL$male$partition <- 0

#get sample sizes from raw=T
tmp <- age[age$AreaName!="32To34Pt5Degs",]
tmp <- SS3AgeAtLen.fn(tmp,lgthBins=seq(8,56,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
tmp$female <- tmp$female[order(tmp$female$year,tmp$female$LbinLo),]
tmp$male <- tmp$male[order(tmp$male$year,tmp$male$LbinLo),]
#check that they are the same
range(tmp$female$year - AatL$female$year)
range(tmp$female$LbinLo - AatL$female$LbinLo)
range(tmp$male$year - AatL$male$year)
range(tmp$male$LbinLo - AatL$male$LbinLo)

AatL$female[,1:9] <- tmp$female[,1:9]
AatL$male[,1:9] <- tmp$male[,1:9]

write.csv(AatL$female,"AgeAtLenForSS3.Num.female.csv",row.names=F)
write.csv(AatL$male,"AgeAtLenForSS3.Num.male.csv",row.names=F)



age <- readInAgeComps.fn("AgeComps.csv",headerRow=7,sep=",")
write.csv(AatL$female,"AgeAtLenForSS3.female.2.csv",row.names=F)
write.csv(AatL$male,"AgeAtLenForSS3.male,2.csv",row.names=F)

#Age at length for 2011 structure
age <- readInAgeComps.fn("AgeComps.csv",headerRow=7,sep=",")
AatL <- SS3AgeAtLen.fn(age,lgthBins=seq(10,64,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28)
write.csv(AatL$female,"AgeAtLenForSS3.2011structure.female.csv",row.names=F)
write.csv(AatL$male,"AgeAtLenForSS3.2011structure.male.csv",row.names=F)

tmp <- age[age$Year==2003,]
SS3AgeAtLen.fn(tmp,lgthBins=seq(10,64,2),ageBins=0:40,fleet=8,season=1,partition=0,ageerr="EnterAgeErr",raw=T,sexRatioUnsexed=0.5,maxSizeUnsexed=28)

slwa <- read.csv("Widow_SexedLgthWtAge.csv",skip=8)
table(slwa$SEX)
table(slwa$AGE_YRS,slwa$SEX)


#sex ratio at length and age (to look at differential M)
slwa <- read.csv("Widow_SexedLgthWtAge.csv",skip=8)
#x <- table(slwa$LENGTH_CM,slwa$SEX)
#plot(as.numeric(rownames(x)),x[,"f"]/(x[,"f"]+x[,"m"]),type="b")

x <- table(cut(slwa$LENGTH_CM,seq(0,55,2)),slwa$SEX)
symbols(1:nrow(x),x[,"f"]/(x[,"f"]+x[,"m"]),circles=x[,"f"]+x[,"m"],inches=0.1,bg=rgb(1,0,0,0.5),xaxt="n",ylim=c(0,1),xlab="Length (cm)",ylab="Proportion female")
lines(1:nrow(x),x[,"f"]/(x[,"f"]+x[,"m"]))
axis(1,at=1:nrow(x),labels=rownames(x),las=2,cex.axis=0.7)
abline(h=0.5,lty=2)



#CV and SD of length at age
slwa <- read.csv("Widow_SexedLgthWtAge.csv",skip=8)
#ageBin currently is fixed at 1 no matter what number you enter
res <- varLengthAtAge.fn(slwa,ageBin=1,bySex=T,parStart=c(54,0.07,1),estVB=T)





#follow the 2008 and 2010 cohorts through time
setwd("C:/NOAA2015/Widow/Data")

load("ExtractedData/wdowTrawlSurveyLwsa.Rdat")
lwsa$Year <- as.numeric(substring(lwsa$PROJECT_CYCLE,7))
age <- lwsa[!is.na(lwsa$AGE_YRS),]
age$cohort <- age$Year-age$AGE_YRS
table(age$Year,age$cohort)

cohorts <- split(age,age$cohort)

dat <- cohorts[["2008"]]
table(dat$HAUL_LATITUDE_DD,dat$Year)
dat <- cohorts[["2010"]]
table(dat$HAUL_LATITUDE_DD,dat$Year)

doPNG <- T
ht <- 6;wd<-12
if(doPNG) {png(filename = "../WriteUp/Figures/cohort2008.png", width = wd, height = ht,units="in",res=300, pointsize = 11)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,5))
xlims <- c(-127.5,-116.45)
ylims <- c(32.0,48.55)
coh <- "2008"
for(yr in unique(cohorts[[coh]]$Year)) {
    dat <- cohorts[[coh]][cohorts[[coh]]$Year==yr,c("HAUL_LONGITUDE_DD","HAUL_LATITUDE_DD","AGE_YRS")]
    names(dat) <- c("X","Y","Z")
    datNoZero <- as.EventData(data.frame(EID=1:nrow(dat[dat$Z>0&!is.na(dat$Z),]),dat[dat$Z>0&!is.na(dat$Z),]))
    plotTheCoast(xlims,ylims)
    addBubbles(datNoZero,max.size=0.3,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",
                legend.pos=100,legend.cex=1.0,legend.title="",legend.breaks=c(0,1),cex=0.7,col=rgb(0,0,0,0.9))
    title(main=paste("Year",yr,"Cohort",coh,"Age",yr-as.numeric(coh)))
}
if(doPNG) {dev.off()}

doPNG <- T
ht <- 6;wd<-12
if(doPNG) {png(filename = "../WriteUp/Figures/cohort2010.png", width = wd, height = ht,units="in",res=300, pointsize = 11)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,5))
xlims <- c(-127.5,-116.45)
ylims <- c(32.0,48.55)
coh <- "2010"
plot(1,1,type="n",xlab="",xaxt="n",ylab="",yaxt="n",bty="n")
plot(1,1,type="n",xlab="",xaxt="n",ylab="",yaxt="n",bty="n")
for(yr in unique(cohorts[[coh]]$Year)) {
    dat <- cohorts[[coh]][cohorts[[coh]]$Year==yr,c("HAUL_LONGITUDE_DD","HAUL_LATITUDE_DD","AGE_YRS")]
    names(dat) <- c("X","Y","Z")
    datNoZero <- as.EventData(data.frame(EID=1:nrow(dat[dat$Z>0&!is.na(dat$Z),]),dat[dat$Z>0&!is.na(dat$Z),]))
    plotTheCoast(xlims,ylims)
    addBubbles(datNoZero,max.size=0.3,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",
                legend.pos=100,legend.cex=1.0,legend.title="",legend.breaks=c(0,1),cex=0.7,col=rgb(0,0,0,0.9))
    title(main=paste("Year",yr,"Cohort",coh,"Age",yr-as.numeric(coh)))
}
if(doPNG) {dev.off()}





#Sample sizes for Table
setwd("C:/NOAA2015/Widow/Data")
load("ExtractedData/wdowTrawlSurveyCatch.Rdat")
wdwCatch <- singleSp
rm(singleSp)
wdwCatch$Year <- as.numeric(substring(wdwCatch$PROJECT_CYCLE,7))
table(wdwCatch$Year,wdwCatch$HAUL_WT_KG>0)

load("ExtractedData/wdowTrawlSurveyLengths.Rdat")
tmp <- lengths[!duplicated(lengths$HAUL_IDENTIFIER),]
table(tmp$PROJECT_CYCLE)
table(lengths$PROJECT_CYCLE)


load("ExtractedData/wdowTrawlSurveyLwsa.Rdat")
wdwAge <- lwsa[!is.na(lwsa$AGE_YRS),]
tmp <- wdwAge[!duplicated(wdwAge$HAUL_IDENTIFIER),]
table(tmp$PROJECT_CYCLE)
table(wdwAge$PROJECT_CYCLE)


