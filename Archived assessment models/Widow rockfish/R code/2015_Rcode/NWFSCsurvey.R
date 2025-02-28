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


setwd("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\NWFSCsurvey")


bioSS <- GetTotalBiomass.fn("BiomassAbundance.csv",headerRow=6)
plotBio.fn(bioSS,pch=16)
title(main="Design-based Biomass from NWFSC trawl survey")

lfs <- readInLengthComps.fn("LengthComps.csv",headerRow=7)
range(lfs$Length)
barplot(table(factor(lfs$Length,levels=6:55)),main="All Widow Lengths from NWFSC survey")
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(6,54,2),gender=3,nSamps="Enter sample sizes",season=1,partition=0,NAs2zero=T)
par(mfrow=c(2,1))
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

ratioF <- lfs$NumF/(lfs$NumF+lfs$NumM)
boxplot(split(ratioF,floor(lfs$Length)),ylab="Prop F",main="Sex ratios by length")
tapply(lfs$NumUnsexed,lfs$Length,sum)
par(mfrow=c(1,1))
plotSexRatio.fn(lfs,ylim=c(0,1),las=1)    #circle size is proportional to the number of observations (strata-year-length specific)
abline(h=0.5)
abline(v=25)

#Compare correcting for unsexed fish and adding them back in
par(mfcol=c(2,2))
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(6,55,2),gender=3,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")
#correct for unsexed fish
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(6,55,2),gender=3,NAs2zero=T,sexRatioUnsexed=0.5,maxSizeUnsexed=25)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

#just for giggles, look at unsexed comps
#use gender =0 to get unsexed (including male and female) proportions
par(mfcol=c(2,2))
lfsSS <- SS3LF.fn(lfs,lgthBins=seq(6,55,2),gender=3,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")
lfsSS  <- SS3LF.fn(lfs,lgthBins=seq(6,55,2),gender=0,NAs2zero=T)
plotFreqData.fn(lfsSS,ylim=c(0,55),yaxs="i",ylab="Length")

###################################################################33
#write out lfs when ready
write.csv(lfsSS,"ForSS3/LengthCompsForSS3.NWFSC.csv")
####################################################################


afs <- readInAgeComps.fn("AgeComps.csv",headerRow=7)
afs$Length <- afs$Age
afsSS  <- SS3LF.fn(afs,lgthBins=seq(0,40,1),gender=3,NAs2zero=T)
par(mfrow=c(2,1))
plotFreqData.fn(afsSS,ylim=c(1,40),yaxs="i",ylab="Age",las=1)

plotFreqData.fn(afsSS,ylim=c(1,10),yaxs="i",ylab="Age",las=1)


write.csv(afsSS,"ForSS3/AgeCompsForSS3.NWFSC.csv")


AatLorig <- GetAges.fn("AgeComps.csv",headerRow=7,lgthBins=seq(6,55,2),ageBins=1:40,raw=T)

age <- readInAgeComps.fn("AgeComps.csv",headerRow=7,sep=",")
age$AgeTallyF <- age$AgeTallyF +age$AgeTallyM +age$AgeTallyU
AatL <- SS3AgeAtLen.fn(age,lgthBins=seq(10,80,2),ageBins=1:100,fleet="EnterFleet",season=1,partition=0,ageerr="EnterAgeErr",raw=T)


write.csv(AatL$female,"ForSS3\\AgeAtLenForSS3.allSex.csv",row.names=F)


#AatL.expanded <- GetAges.fn(paste(directory,"AgeComps.csv",sep="\\"),headerRow=7,lgthBins=seq(8,60,2),ageBins=1:60,raw=F)
write.csv(AatL$female,paste(directory,"ForSS3\\AgeAtLenForSS3.female.csv",sep="\\"),row.names=F)
write.csv(AatL$male,paste(directory,"ForSS3\\AgeAtLenForSS3.male.csv",sep="\\"),row.names=F)


AatL <- GetAges.fn(paste(directory,"AgeComps.csv",sep="\\"),headerRow=7,lgthBins=seq(1,14,1),ageBins=1:40,raw=T)
AatL$female[AatL$female$LbinLo<12,]

#How come the porportions of male ages seems so small?
slwa <- read.csv("Widow_SexedLgthWtAge.csv",skip=8)
table(slwa$SEX)
table(slwa$AGE_YRS,slwa$SEX)
#It seems pretty even in number aged across sex. In fact, more males are sexed.


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

