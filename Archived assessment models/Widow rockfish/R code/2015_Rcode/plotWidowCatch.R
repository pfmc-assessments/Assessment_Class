#install.packages(c("maptools","PBSmapping","RgoogleMaps"))
########################################################################

library(maptools)
library(PBSmapping)
library(RgoogleMaps)

plotCellsGoogleMap <- function(lats,lons) {
    tmp1 <- LatLon2XY.centered(siteTerrain,rep(min(lats),length(lons)),lons)
    tmp2 <- LatLon2XY.centered(siteTerrain,rep(max(lats),length(lons)),lons)
    segments(tmp1$newX,tmp1$newY,tmp2$newX,tmp2$newY,lwd=3)
    tmp1 <- LatLon2XY.centered(siteTerrain,lats,rep(min(lons),length(lats)),siteTerrain$zoom)
    tmp2 <- LatLon2XY.centered(siteTerrain,lats,rep(max(lons),length(lats)),siteTerrain$zoom)
    segments(tmp1$newX,tmp1$newY,tmp2$newX,tmp2$newY,lwd=3)
}

plotCells <- function(lats,lons) {
    segments(lons,rep(min(lats),length(lons)),lons,rep(max(lats),length(lons)))
    segments(rep(min(lons),length(lats)),lats,rep(max(lons),length(lats)),lats)
}

aggCatch.fn <- function(atsea,latDiv,lonDiv,year) {
    atsea <- atsea[atsea$Year%in%year,]
    latCut <- cut(atsea$RETRV_LATITUDE,LatDiv)
    lonCut <- cut(atsea$RETRV_LONGITUDE,LonDiv)

    y <- (LatDiv[-length(LatDiv)]+LatDiv[-1])/2  #midpts of each bin
    x <- (LonDiv[-length(LonDiv)]+LonDiv[-1])/2  #midpts of each bin
    names(y) <- cut(y,LatDiv)
    names(x) <- cut(x,LonDiv)
    atsea$latAgg <- y[latCut]
    atsea$lonAgg <- x[lonCut]
    atseaTable <- apply(table(factor(atsea$latAgg,levels=y),factor(atsea$lonAgg,levels=x),atsea$VESSEL)>0, c(1,2), sum)
    print(atseaTable)
    numVessels <- data.frame(Lat = rep(rownames(atseaTable),ncol(atseaTable)),
                             Lon = rep(colnames(atseaTable),each=nrow(atseaTable)),
                             numVessels = as.vector(atseaTable))
    confidential <- numVessels[numVessels$numVessels<3 & numVessels$numVessels>0,]

    #do not aggregate by month unless checked for confidentiality
    dat <- atsea[,c("Month","lonAgg","latAgg","EXTRAPOLATED_WEIGHT")]
    dat$EXTRAPOLATED_WEIGHT <- dat$EXTRAPOLATED_WEIGHT/1000
    dat <- aggregate(dat$EXTRAPOLATED_WEIGHT,list(dat$latAgg,dat$lonAgg),sum)
    names(dat) <- c("latAgg","lonAgg","Catch.MT")
    dat <- dat[order(dat$latAgg,dat$lonAgg),]

    #remove confidential cells
    print(paste(dat$latAgg,dat$lonAgg))
    print(paste(confidential$Lat,confidential$Lon))
    rem <- paste(dat$latAgg,dat$lonAgg) %in% paste(confidential$Lat,confidential$Lon)
    print(rem)
    dat <- dat[!rem,]

    return(dat)
}



setwd("C:/NOAA2015/Widow/Data")
source("Rcode/functions/Functions.R")
source("C:/Mapping/WestCoastMapping.R")
#library(date)

#load("extractedData/atsea.bio.Type1.Rdat")
#atsea.bio$Month <- as.numeric(substr(atsea.bio$HAUL_OFFLOAD_DATE,6,7))
load("extractedData/NORPACdomesticCatchDetailed.Rdat")
wcatch <- ncatch2[ncatch2$SPECIES==305,]
wcatch <- wcatch[!is.na(wcatch$EXTRAPOLATED_WEIGHT),]
wcatch$Month <- as.numeric(substr(wcatch$HAUL_DATE,6,7))
wcatch$Year <- as.numeric(substr(wcatch$HAUL_DATE,1,4))
wcatch$RETRV_LATITUDE <- floor(wcatch$RETRV_LATITUDE/100) + 100*(wcatch$RETRV_LATITUDE/100-floor(wcatch$RETRV_LATITUDE/100))/60
wcatch$RETRV_LONGITUDE <- -1*(floor(wcatch$RETRV_LONGITUDE/100) + 100*(wcatch$RETRV_LONGITUDE/100-floor(wcatch$RETRV_LONGITUDE/100))/60)

#hotspots shapefiles
hotspots <- importShapefile("C:/NOAA2015/Widow/WidowSurvey/Shapefiles/hotspots/Widow_hotspots.shp")


#plot(atsea.bio$LONDD_START,atsea.bio$LATDD_START,pch=16)
plot(wcatch$RETRV_LONGITUDE,wcatch$RETRV_LATITUDE,pch=16)
#atsea <- wcatch[wcatch$Year==2014,]

scale <- 2
#siteTerrain <- GetMap.bbox(lonR=c(-126,-124),latR=c(41.4,48.6),size = c(320, 640), destfile = "Figures/RgoogleMap.png",
#                           MINIMUMSIZE = FALSE, RETURNIMAGE = TRUE, GRAYSCALE = FALSE,
#                           NEWMAP = TRUE, verbose = 1, SCALE = 2,maptype="hybrwcatch$RETRV_LATITUDEid")#,frame=10,zoom=6)  #,markers=mymarkers)  #zoom is another argument
bb <- qbbox(lat=c(37.5,48.5), lon=-c(127,122), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
siteTerrain <- GetMap.bbox(lonR=bb$lonR,latR=bb$latR, size = c(320, 640), destfile = "Figures/RgoogleMap.png",
                           MINIMUMSIZE = FALSE, RETURNIMAGE = TRUE, GRAYSCALE = FALSE,
                           NEWMAP = TRUE, verbose = 1, SCALE = 2,maptype="hybrid")#,frame=10,zoom=6)  #,markers=mymarkers)  #zoom is another argument



#use wcatch
maxCatch <- 20000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 5

 png("Figures/confidentialWidowAtSea.png",width=320*scale,height=640*scale)
 PlotOnStaticMap(siteTerrain,lat=wcatch$RETRV_LATITUDE,lon=wcatch$RETRV_LONGITUDE,add=FALSE,col="white",pch=16)
 tmp1 <- LatLon2XY.centered(siteTerrain,wcatch$RETRV_LATITUDE,wcatch$RETRV_LONGITUDE)
 ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
 ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
 points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
 #plotCellsGoogleMap(LatDiv,LonDiv)
 tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-128)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-127.5)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()

LatDiv <- seq(37.5,49.1,0.1)
LonDiv <- seq(-126,-123,0.1)
year <- 1991:2014
dat <- aggCatch.fn(wcatch,LatDiv,LonDiv,year)
maxCatch <- 40000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 10
png("Figures/confidentialWidowAtSea.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=0,lon=0,add=F,col=grey(0.6))
#plotCellsGoogleMap(LatDiv,LonDiv)
ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
tmp1 <- LatLon2XY.centered(siteTerrain,dat$latAgg,dat$lonAgg)
points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
#plotCellsGoogleMap(LatDiv,LonDiv)
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-128)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-127.5)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()



scale <- 2
bb <- qbbox(lat=c(43,48.5), lon=-c(127,123), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
#bb <- qbbox(lat=c(43,47), lon=-c(126,123.5), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
siteTerrain <- GetMap.bbox(lonR=bb$lonR,latR=bb$latR, size = c(320, 640), destfile = "Figures/RgoogleMapZoom.png",
                           MINIMUMSIZE = FALSE, RETURNIMAGE = TRUE, GRAYSCALE = FALSE,
                           NEWMAP = TRUE, verbose = 1, SCALE = 2,maptype="hybrid")#,frame=10,zoom=6)  #,markers=mymarkers)  #zoom is another argument

#use wcatch
maxCatch <- 20000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 10

 png("Figures/confidentialWidowAtSeaZoom.png",width=320*scale,height=640*scale)
 PlotOnStaticMap(siteTerrain,lat=wcatch$RETRV_LATITUDE,lon=wcatch$RETRV_LONGITUDE,add=FALSE,col="white",pch=16)
 PlotPolysOnStaticMap(siteTerrain, hotspots, col = rgb(0,1,0,0.1), add = TRUE)
 tmp1 <- LatLon2XY.centered(siteTerrain,wcatch$RETRV_LATITUDE,wcatch$RETRV_LONGITUDE)
 ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
 ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
 points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
 #plotCellsGoogleMap(LatDiv,LonDiv)
 tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-128)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-127.5)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()

LatDiv <- seq(37.5,49.1,0.05)
LonDiv <- seq(-126,-123,0.1)
year <- 1991:2014
dat <- aggCatch.fn(wcatch,LatDiv,LonDiv,year)
maxCatch <- 40000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 10
png("Figures/confidentialWidowAtSeaZoom.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=0,lon=0,add=F,col=grey(0.6))
#plotCellsGoogleMap(LatDiv,LonDiv)
ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
tmp1 <- LatLon2XY.centered(siteTerrain,dat$latAgg,dat$lonAgg)
points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
#plotCellsGoogleMap(LatDiv,LonDiv)
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-128)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-127.5)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()






scale <- 2
bb <- qbbox(lat=c(43.5,45), lon=-c(126,124), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
#bb <- qbbox(lat=c(43,47), lon=-c(126,123.5), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
siteTerrain <- GetMap.bbox(lonR=bb$lonR,latR=bb$latR, size = c(320, 640), destfile = "Figures/RgoogleMapZoom.png",
                           MINIMUMSIZE = FALSE, RETURNIMAGE = TRUE, GRAYSCALE = FALSE,
                           NEWMAP = TRUE, verbose = 1, SCALE = 2,maptype="hybrid")#,frame=10,zoom=6)  #,markers=mymarkers)  #zoom is another argument
#use wcatch
maxCatch <- 20000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 10
png("Figures/confidentialWidowAtSeaZoom2.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=wcatch$RETRV_LATITUDE,lon=wcatch$RETRV_LONGITUDE,add=FALSE,col="white",pch=16)
tmp1 <- LatLon2XY.centered(siteTerrain,wcatch$RETRV_LATITUDE,wcatch$RETRV_LONGITUDE)
ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
#plotCellsGoogleMap(LatDiv,LonDiv)
tmp <- LatLon2XY.centered(siteTerrain,seq(42.5,43.5,length=div),-126.5)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(42.5,43.5,length=div),-126.1)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()


scale <- 2
bb <- qbbox(lat=c(43.5,45), lon=-c(126,124), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
#bb <- qbbox(lat=c(43,47), lon=-c(126,123.5), TYPE = "all",margin=list(m=c(1,1,1,1),TYPE = c("perc", "abs")[1]))
siteTerrain <- GetMap.bbox(lonR=bb$lonR,latR=bb$latR, size = c(320, 640), destfile = "Figures/RgoogleMapZoom.png",
                           MINIMUMSIZE = FALSE, RETURNIMAGE = TRUE, GRAYSCALE = FALSE,
                           NEWMAP = TRUE, verbose = 1, SCALE = 2,maptype="hybrid")#,frame=10,zoom=6)  #,markers=mymarkers)  #zoom is another argument
#use wcatch
year <- 2014
dat <- wcatch[wcatch$Year==2014,]
maxCatch <- max(dat$EXTRAPOLATED_WEIGHT)
div <- 10
png("Figures/confidentialWidowAtSea2014.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=dat$RETRV_LATITUDE,lon=dat$RETRV_LONGITUDE,add=FALSE,col="white",pch=16)
tmp1 <- LatLon2XY.centered(siteTerrain,dat$RETRV_LATITUDE,dat$RETRV_LONGITUDE)
ind <- floor(div*dat$EXTRAPOLATED_WEIGHT/maxCatch)
ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
#plotCellsGoogleMap(LatDiv,LonDiv)
tmp <- LatLon2XY.centered(siteTerrain,seq(42.5,43.5,length=div),-126.5)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(42.5,43.5,length=div),-126.1)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
tmp <- LatLon2XY.centered(siteTerrain,42.3,-126.3)
text(tmp$newX,tmp$newY,year,col="white",cex=1.5)
dev.off()
















LatDiv <- seq(37.5,49.1,0.05)
LonDiv <- seq(-126,-123,0.1)
year <- 2014
dat <- aggCatch.fn(wcatch,LatDiv,LonDiv,year)
maxCatch <- 40000 #max(wcatch$EXTRAPOLATED_WEIGHT)
div <- 10
png("Figures/confidentialWidowAtSeaZoom.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=0,lon=0,add=F,col=grey(0.6))
#plotCellsGoogleMap(LatDiv,LonDiv)
ind <- floor(div*wcatch$EXTRAPOLATED_WEIGHT/maxCatch)
ind[ind>div] <- div
# ind[ind==0] <- 1 #this puts all the small catches into a bin to be plotted, commenting this out results in not plotting small catches
tmp1 <- LatLon2XY.centered(siteTerrain,dat$latAgg,dat$lonAgg)
points(tmp1$newX,tmp1$newY,cex=1,pch=16,col=rainbow(div,start=0.7,end=0.1)[ind])
#plotCellsGoogleMap(LatDiv,LonDiv)
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-128)
points(tmp$newX,tmp$newY,cex=5,pch=16,col=rainbow(div,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(35,39,length=div),-127.5)
text(tmp$newX,tmp$newY,paste(floor(maxCatch*(1:div)/div),"mt"),col="white",cex=0.9)
dev.off()















year <- 2014
dat <- aggCatch.fn(wcatch,LatDiv,LonDiv,year)
png("Figures/nonconfidentialAtSea2014.png",width=320*scale,height=640*scale)
PlotOnStaticMap(siteTerrain,lat=sa$Lat,lon=sa$Lon,add=F,col=grey(0.6),cex=sqrt(0.1+sa$NASC/(max(sa$NASC)/100)))
plotCellsGoogleMap(LatDiv,LonDiv)
tmp1 <- LatLon2XY.centered(siteTerrain,dat$latAgg,dat$lonAgg)
points(tmp1$newX,tmp1$newY,cex=5,pch=15,col=rainbow(10,start=0.7,end=0.1)[floor(10*dat$Catch.MT/max(dat$Catch.MT))])
tmp <- LatLon2XY.centered(siteTerrain,seq(42,45,length=10),-128)
points(tmp$newX,tmp$newY,cex=5,pch=15,col=rainbow(10,start=0.7,end=0.1))
tmp <- LatLon2XY.centered(siteTerrain,seq(42,45,length=10),-127.5)
text(tmp$newX,tmp$newY,paste(floor(max(dat$Catch.MT)*(1:5)/10),"mt"),col="white",cex=0.9)
tmp <- LatLon2XY.centered(siteTerrain,45.4,-127.7)
text(tmp$newX,tmp$newY,paste(year,"fishery"),col="white",cex=1.9)
tmp <- LatLon2XY.centered(siteTerrain,45.9,-127.7)
text(tmp$newX,tmp$newY,paste(2013,"survey"),col="white",cex=1.9)
dev.off()


