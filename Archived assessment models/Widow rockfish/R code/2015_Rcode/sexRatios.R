plotSexRatio.fn <- function (dat, circleSize = 0.1, format=c("numF","obsRow"), ...) {
	if(format=="numF") {
	    yF <- lapply(split(dat, floor(dat$Length)), function(x) {
	    				sum(x$NumF,na.rm=T)/(sum(x$NumF,na.rm=T)+sum(x$NumM,na.rm=T))})		
	}
	if(format=="obsRow") {
	    yF <- lapply(split(dat, floor(dat$Length)), function(x) {
    	            indF <- x$Sex=="f"
    	            indM <- x$Sex=="m"
    				sum(indF,na.rm=T)/(sum(indF,na.rm=T)+sum(indM,na.rm=T))})		
	}
    x <- names(split(dat, floor(dat$Length)))
    nobs <- unlist(lapply(split(dat$Length, floor(dat$Length)), length))
    plot(x, yF, type = "l", col = "red", ylab = "Fraction female", ...)
    symbols(x, yF, circles = nobs, inches = circleSize, fg = "red", 
        bg = rgb(1, 0, 0, alpha = 0.5), add = T)
    return(invisible(data.frame(X = x, fraction.female = as.numeric(yF))))
}

setwd("C:/NOAA2015/Widow/Data")
lfs <- readInLengthComps.fn("TrawlSurvey/NWFSCsurvey/LengthComps.csv",headerRow=7)
doPNG <- F
wd<-6.5;ht<-3.5
if(doPNG) {png("../Writeup/Figures/nwfscSexRatios.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,1),mar=c(4,4,1,1))
plotSexRatio.fn(lfs,ylim=c(0,1),las=1,format="numF",xlab="Length (cm)")    #circle size is proportional to the number of observations (strata-year-length specific)
abline(h=0.5,col=rgb(0,0,0,0.3))
#abline(v=28)
if(doPNG) dev.off()

bioData <- read.csv("Biological/allBioData.csv")
doPNG <- T
wd<-6.5;ht<-3.5
if(doPNG) {png("../Writeup/Figures/allSexRatios.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,1),mar=c(4,4,1,1))
plotSexRatio.fn(bioData[bioData$Length<=60,],ylim=c(0,1),las=1,format="obsRow",xlab="Length (cm)")    #circle size is proportional to the number of observations (strata-year-length specific)
abline(h=0.5,col=rgb(0,0,0,0.3))
#abline(v=28)
if(doPNG) dev.off()

bioData <- read.csv("Biological/allBioData.csv")
bioData$Length <- bioData$Age
doPNG <- F
wd<-6.5;ht<-3.5
if(doPNG) {png("../Writeup/Figures/allSexRatiosAge.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,1),mar=c(4,4,1,1))
plotSexRatio.fn(bioData[bioData$Length<=60,],ylim=c(0,1),xlim=c(0,30),las=1,format="obsRow",xlab="Age")    #circle size is proportional to the number of observations (strata-year-length specific)
abline(h=0.5,col=rgb(0,0,0,0.3))
#abline(v=28)
if(doPNG) dev.off()



