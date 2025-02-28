# devtools::install_github("nwfsc-assess/nwfscMapping")
# devtools::install_github("nwfsc-assess/nwfscSurvey")
# Load package
library(nwfscMapping)
library(PBSmapping)
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

setwd("C:/NOAA2015/Widow/Data")


load("ExtractedData/wdowTrawlSurveyCatch.Rdat")
wdwCatch <- singleSp
rm(singleSp)
wdwCatch$CPHA <- wdwCatch$HAUL_WT_KG/wdwCatch$AREA_SWEPT_HA
wdwCatch$Year <- as.numeric(substring(wdwCatch$PROJECT_CYCLE,7))
table(wdwCatch$Year,round(wdwCatch$BEST_DEPTH_M,-1))

comboC <- wdwCatch[wdwCatch$Year >= 2003,]
slopeC <- wdwCatch[wdwCatch$Year <  2003,]


dat <- comboC[,c("BEST_LON_DD","BEST_LAT_DD","CPHA")]
names(dat) <- c("X","Y","Z")
dat <- as.EventData(data.frame(EID=1:nrow(dat),dat))
datNoZero <- as.EventData(data.frame(EID=1:nrow(dat[dat$Z>0&!is.na(dat$Z),]),dat[dat$Z>0&!is.na(dat$Z),]))

ht <- 12;wd<-10
windows(height=ht,width=wd)
xlims <- c(-127.5,-116.45)
ylims <- c(32.0,48.55)
plotTheCoast(xlims,ylims)
#addBubbles(dat,max.size=0.2,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",legend.pos="topright",legend.cex=1.0,legend.title="Length (cm)",legend.breaks=c(10,100),cex=0.7,col=rgb(0,0,0,0.9))
addBubbles(datNoZero,max.size=0.3,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",
			legend.pos="topright",legend.cex=1.0,legend.title="Catch Rate (KG/HA)",legend.breaks=c(30,300),cex=0.7,col=rgb(0,0,0,0.9))

dat <- slopeC[,c("BEST_LON_DD","BEST_LAT_DD","CPHA")]
dat <- dat[!is.na(dat$BEST_LON_DD),]   #one observation in 1998 doesn't have position
names(dat) <- c("X","Y","Z")
dat <- as.EventData(data.frame(EID=1:nrow(dat),dat))
datNoZero <- as.EventData(data.frame(EID=1:nrow(dat[dat$Z>0&!is.na(dat$Z),]),dat[dat$Z>0&!is.na(dat$Z),]))

ht <- 12;wd<-10
windows(height=ht,width=wd)
xlims <- c(-127.5,-116.45)
ylims <- c(32.0,48.55)
plotTheCoast(xlims,ylims)
#addBubbles(dat,max.size=0.2,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",legend.pos="topright",legend.cex=1.0,legend.title="Length (cm)",legend.breaks=c(10,100),cex=0.7,col=rgb(0,0,0,0.9))
addBubbles(datNoZero,max.size=0.04,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero=".",
			legend.pos="topright",legend.cex=1.0,legend.title="Catch Rate (KG/HA)",legend.breaks=c(30,300),cex=0.7,col=rgb(0,0,0,0.9))


#examine lengths
load("ExtractedData/wdowTrawlSurveyLengths.Rdat")
wdwLens <- lengths
rm(lengths)
table(wdwLens$PROJECT_CYCLE)

#Length by depth and latitude
windows(height=5,width=6.5)
par(mfrow=c(2,1),mar=c(3.6,3.6,1,1))
plot(wdwLens$DEPTH_M,wdwLens$LENGTH_CM,pch=20,cex=0.8,col=rgb(0,0,0,0.5),xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0))
brks <- c(135,183)
abline(v=brks)
axis(3,at=brks,cex.axis=0.8,mgp=c(3,0.25,0))
#there are three fish smaller than 10cm, and it looks like the deeper ones could have been left from the earlier tow
#Remember to remove these or fix these
tmp <- wdwLens[wdwLens$LENGTH_CM<10,]
plot(wdwLens$HAUL_LATITUDE_DD,wdwLens$LENGTH_CM,pch=20,cex=0.8,col=rgb(0,0,0,0.5),xlab="Latitude",ylab="Length (cm)",las=1,mgp=c(2,0.7,0))
brks <- c(36,40.17)
abline(v=brks)
axis(3,at=brks,cex.axis=0.8,mgp=c(3,0.25,0))


########################
#tree regression (not extremely useful)
library(tree)
source("C:/NOAA2015/Widow/Data/Rcode/Functions/StrataFunctions.R")

#these are recommended breaks for the NWFSC survey
#put in your own breaks for your assessment
latBreaks <- c(34.5)
depthBreaks <- c(50,183,549,1300)
latLims <- c(30,49) #entire coast
depLims <- c(50,400)  #this will eliminate the two really deep small lengths which are likely from the previous tow
years <- 2003:2014  #only these years have lengths and are the combo survey

#there are three fish smaller than 10cm, and it looks like the deeper ones could have been left from the earlier tow
#Remember to remove these or fix these
tmp <- lens[lens$LENGTH_CM<10,]

lens <- wdwLens
lens$Year <- as.numeric(substring(lens$PROJECT_CYCLE,7))
lens <- lens[lens$DEPTH_M>=depLims[1] & lens$DEPTH_M<=depLims[2],]
lens <- lens[lens$HAUL_LATITUDE_DD>=latLims[1] & lens$HAUL_LATITUDE_DD<=latLims[2],]
lens <- lens[lens$Year>=min(years) &lens$Year<=max(years),]

x <- treeStrata.fn(lens[lens$LENGTH_CM>10,],"LENGTH_CM",c("DEPTH_M"))   #tree strata by only depth
   #150.45, 112.55,  85.6
x <- treeStrata.fn(lens[lens$LENGTH_CM>10,],"LENGTH_CM",c("HAUL_LATITUDE_DD"))      #tree strata by only latitude
   #36.26, 42.559, 48.0, 46.23
x <- treeStrata.fn(lens[lens$LENGTH_CM>10,],"LENGTH_CM",c("DEPTH_M","HAUL_LATITUDE_DD"))      #tree strata by depth & latitude
   #>150.45 (no additional strata), 0-85.6, 85.6-112.55, 112.55-150.45: 37.2lat

#I choose
#  135 m: somwhat stratifies small from large
#  183 m: is a sampling design break, and it breaks out more large
#  36 lat: not many large fish south of there. I don't use the samplign break of 34.5 because of the paucity of data.
#  40 10 lat: useful management break

#Length by depth and latitude
windows(height=5,width=6.5)
par(mfrow=c(2,1),mar=c(3.6,3.6,1,1))
plot(lens$DEPTH_M,lens$LENGTH_CM,pch=20,cex=0.8,col=rgb(0,0,0,0.5),xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0))
brks <- c(135,183)
abline(v=brks)
axis(3,at=brks,cex.axis=0.8,mgp=c(3,0.25,0))
#there are three fish smaller than 10cm, and it looks like the deeper ones could have been left from the earlier tow
#Remember to remove these or fix these
tmp <- lens[lens$LENGTH_CM<10,]
plot(lens$HAUL_LATITUDE_DD,lens$LENGTH_CM,pch=20,cex=0.8,col=rgb(0,0,0,0.5),xlab="Latitude",ylab="Length (cm)",las=1,mgp=c(2,0.7,0))
brks <- c(36,40.17)
abline(v=brks)
axis(3,at=brks,cex.axis=0.8,mgp=c(3,0.25,0))

#There are few data in some strata, and few observations generally south of 34.5 (25 over all years)
#so, I propose to remove that area and have the following strata
latBreaks <- c(34.5,40.17)
depthBreaks <- c(115,183)


windows(height=7.5,width=6.5)
par(mfrow=c(3,1),mar=c(3.6,3.6,1,1))
ind <- lens$HAUL_LATITUDE_DD > latBreaks[2]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(0,0,0,0.5),xlim=c(50,400),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="North of 40 10 latitude")
abline(v=depthBreaks)
axis(3,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))
ind <- lens$HAUL_LATITUDE_DD<=latBreaks[2] & lens$HAUL_LATITUDE_DD>=latBreaks[1]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(0,0,0,0.5),xlim=c(50,400),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="South of 40 10 and North of 34.5 latitude")
brks <- c(115,183)
abline(v=depthBreaks)
axis(3,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))
ind <- lens$HAUL_LATITUDE_DD < latBreaks[1]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(0,0,0,0.5),xlim=c(50,400),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="South of 34.5 latitude")
brks <- c(183)
abline(v=depthBreaks)
axis(3,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))

classify.fn <- function(x,breaks=c(36,40.5,43,47.5),labels=c("CP","MT","EK","CL","VN")) {
    #classifies areas from breaks
    x <- findInterval(x,breaks)
    return(ordered(labels[x+1],labels))
}
cbind(33:48,classify.fn(33:48,latBreaks,c("South","Central","North")))

lens$latStrat <- classify.fn(lens$HAUL_LATITUDE_DD,latBreaks,c("South","Central","North"))
lens$depthStrat <- classify.fn(lens$DEPTH_M,depthBreaks,c("shallow","midShelf","slope"))
#lens$depthStrat[lens$latStrat=="South"&lens$depthStrat=="midShelf"] <- "shallow"

table(lens$Year,lens$latStrat)
table(lens$Year,lens$depthStrat)
table(lens$Year,lens$depthStrat,lens$latStrat)

#There jsut isn't good coverage shallower than 183, so Im removing that shallow depth break.
#notice, though that there are sometimes lots of fish seen shallow and not midshlef, and vice versa
#also notice that hardly any fish were seen in 2012 in the Central area, but the most in the south area

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###Final

latBreaks <- c(34.5,40.5)
#34.5 because it is the sampling design change, and a lot of area in the SCB. REMOVE FROM ANALYSIS
#40.5 because it is a management line (and I used it instead of 40.17 so that 1 observation would be present in Central Slope in 2012)
depthBreaks <- c(183)
#183 is a smapling design break
#nothing else because of limited data within years shallower than this

windows(height=7.5,width=6.5)
par(mfrow=c(3,1),mar=c(3.6,3.6,1.2,1))
ind <- lens$HAUL_LATITUDE_DD > latBreaks[2]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(0,0,1,0.5),xlim=c(50,400),ylim=c(0,55),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="North of 40 10 latitude")
abline(v=depthBreaks)
axis(1,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))
ind <- lens$HAUL_LATITUDE_DD<=latBreaks[2] & lens$HAUL_LATITUDE_DD>=latBreaks[1]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(0,0,1,0.5),xlim=c(50,400),ylim=c(0,55),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="South of 40 10 and North of 34.5 latitude")
abline(v=depthBreaks)
axis(1,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))
ind <- lens$HAUL_LATITUDE_DD < latBreaks[1]
plot(lens$DEPTH_M[ind],lens$LENGTH_CM[ind],pch=20,cex=1.3,col=rgb(1,0,0,0.5),xlim=c(50,400),ylim=c(0,55),
		xlab="Depth (m)",ylab="Length (cm)",las=1,mgp=c(2,0.7,0),main="South of 34.5 latitude")
abline(v=depthBreaks)
axis(1,at=depthBreaks,cex.axis=0.8,mgp=c(3,0.25,0))

classify.fn <- function(x,breaks=c(36,40.5,43,47.5),labels=c("CP","MT","EK","CL","VN")) {
    #classifies areas from breaks
    x <- findInterval(x,breaks)
    return(ordered(labels[x+1],labels))
}

lens$latStrat <- classify.fn(lens$HAUL_LATITUDE_DD,latBreaks,c("South","Central","North"))
lens$depthStrat <- classify.fn(lens$DEPTH_M,depthBreaks,c("shelf","slope"))

table(lens$Year,lens$latStrat)
table(lens$Year,lens$depthStrat)
table(lens$Year,lens$depthStrat,lens$latStrat)

#look at wdwCatch
wdw <- wdwCatch[wdwCatch$HAUL_WT_KG>0,]
wdw$latStrat <- classify.fn(wdw$BEST_LAT_DD,latBreaks,c("South","Central","North"))
wdw$depthStrat <- classify.fn(wdw$BEST_DEPTH_M,depthBreaks,c("shelf","slope"))
table(wdw$Year,wdw$latStrat)
table(wdw$Year,wdw$depthStrat)
table(wdw$Year,wdw$depthStrat,wdw$latStrat)



### Triennial