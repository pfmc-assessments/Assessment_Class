library(PBSmapping)
#lat/lon from state in the maps package
#organized by hand

CA <- matrix(c(
-124.1829, 42.00354,
-123.2375, 42.00927,
-122.3036, 42.00927,
-120.0060, 42.00927,
-119.9946, 38.98978,
-119.8685, 38.90956,
-117.8460, 37.47717,
-114.6374, 35.01918,
-114.6087, 34.88740,
-114.5515, 34.83583,
-114.5744, 34.80719,
-114.4483, 34.70978,
-114.3738, 34.46914,
-114.1389, 34.29725,
-114.1332, 34.26287,
-114.4197, 34.07953,
-114.4540, 33.99932,
-114.5400, 33.92483,
-114.5228, 33.73003,
-114.5572, 33.57533,
-114.6202, 33.52949,
-114.6087, 33.50084,
-114.7291, 33.40344,
-114.7405, 33.31176,
-114.6890, 33.26593,
-114.7004, 33.08258,
-114.5400, 33.04247,
-114.4827, 32.99664,
-114.4770, 32.93362,
-114.5515, 32.76173,
-114.7348, 32.73308,
-117.1126, 32.53827)
,ncol=2,dimnames=list(NULL,c("x","y")),byrow=T)



#tmpCA <- map('state',"california",xlim=c(-125.3,-114.2),ylim=c(32,50))
#lines(CA,col="red",lwd=2)


OR <- matrix(c(
-123.4323, 46.22623,
-123.3693, 46.14029,
-123.1172, 46.16894,
-122.9109, 46.07153,
-122.7906, 45.87673,
-122.7619, 45.64182,
-122.3208, 45.53296,
-121.7765, 45.69339,
-121.3696, 45.69339,
-121.1634, 45.60171,
-121.0545, 45.64755,
-120.9629, 45.64182,
-120.6535, 45.72776,
-120.4472, 45.68765,
-119.6966, 45.85381,
-119.5935, 45.91684,
-119.3185, 45.93402,
-119.2497, 45.92830,
-118.9575, 45.99705,
-116.9292, 45.99705,
-116.7459, 45.80798,
-116.7000, 45.83662,
-116.5453, 45.76214,
-116.4651, 45.61317,
-116.6829, 45.23502,
-116.8490, 44.99438,
-116.8318, 44.94854,
-116.9178, 44.81676,
-117.0381, 44.74800,
-117.2100, 44.43861,
-117.2042, 44.30683,
-117.1584, 44.26099,
-117.0954, 44.27245,
-117.0496, 44.23234,
-117.0037, 44.25526,
-116.9521, 44.22088,
-116.9063, 44.14640,
-116.9579, 44.08910,
-116.9350, 44.00889,
-117.0266, 43.79116,
-117.0266, 42.00927,
-119.3070, 42.00927,
-120.0060, 42.00927,
-122.3036, 42.00927,
-123.2375, 42.00927,
-124.1829, 42.00354)
,ncol=2,dimnames=list(NULL,c("x","y")),byrow=T)

#plotMap(nepacLL,xlim=c(-125.3,-123),ylim=c(45.5,46.5),col=gray(0.5))#,tck=c(-0.03),cex=1.0,mgp=c(1.9,0.7,0))#,plt=c(0.16,0.97,0.16,0.97))
#tmpOR <- map('state',"oregon",xlim=c(-125.3,-114.2),ylim=c(32,50))
#lines(OR,col="blue",lwd=1)

#par(mfrow=c(1,2))
#plotMap(nepacLL,xlim=c(-126.3,-114.2),ylim=c(41.5,49.5),col=gray(0.5))#,tck=c(-0.03),cex=1.0,mgp=c(1.9,0.7,0))#,plt=c(0.16,0.97,0.16,0.97))
#tmpWA <- map('state',"washington",xlim=c(-126.3,-114.2),ylim=c(41.5,49.5))

WA <- matrix(c(
-122.7104, 48.98789,
-122.0916, 48.99362,
-117.0209, 48.99362,
-117.0209, 48.99362,
-117.0266, 46.53563,
-117.0266, 46.42677,
-117.0324, 46.38666,
-116.9235, 46.16894,
-116.9292, 45.99705,
-118.9575, 45.99705,
-119.2497, 45.92830,
-119.3185, 45.93402,
-119.5935, 45.91684,
-119.6966, 45.85381,
-120.4472, 45.68765,
-120.6535, 45.72776,
-120.9629, 45.64182,
-121.0545, 45.64755,
-121.1634, 45.60171,
-121.3696, 45.69339,
-121.7765, 45.69339,
-122.3208, 45.53296,
-122.7619, 45.64182,
-122.7906, 45.87673,
-122.9109, 46.07153,
-123.1172, 46.16894,
-123.3693, 46.14029,
-123.4323, 46.22623)
,ncol=2,dimnames=list(NULL,c("x","y")),byrow=T)
#lines(WA,col="blue")


WCstatesInlandPBS <- as.PolySet(data.frame(PID=c(rep(1,nrow(CA)),rep(3,nrow(OR)),rep(4,nrow(WA))),
                                        SID=c(rep(1,nrow(CA)),rep(3,nrow(OR)),rep(4,nrow(WA))),
                                        POS=c(1:nrow(CA),1:nrow(OR),1:nrow(WA)),
                                        X=c(CA[,"x"],OR[,"x"],WA[,"x"]),
                                        Y=c(CA[,"y"],OR[,"y"],WA[,"y"])),projection="LL")
print("Created WCstatesInlandPBS")
#plotMap(nepacLL,xlim=c(-126.3,-113),ylim=c(31.5,50))#,col=gray(0.5))#,tck=c(-0.03),cex=1.0,mgp=c(1.9,0.7,0))#,plt=c(0.16,0.97,0.16,0.97))
#tmp <- addLines(WCstatesInlandPBS)

WCstatesInland <- list(CA=CA,OR=OR,WA=WA)       #only the inland borders, this allows different packages, such as PBSmapping to be used for coast
print("Created WCstatesInlandPBS and WCstateInland")








###################################################
#create a west coast PBS map database that goes south of 34 latitude
#uses worldLLhigh which is thinned to every 1km, whereas nepacLL is thinned every
library(PBSmapping)
data(worldLLhigh)
data(nepacLL)

#go to latitude of 34.1 to get the islands that PBSmapping cuts off
tmp <- worldLLhigh[worldLLhigh$Y>21 & worldLLhigh$Y<=34.1 & worldLLhigh$X>360-126 & worldLLhigh$X<360-103,]
#PID=1 is the west coast. there are also some islands included
tmp$X <- tmp$X-360
tmpIslands <- tmp[tmp$PID!=1,]           #this also gets the channel islands
tmp <- tmp[tmp$PID==1 & tmp$Y<=34,] #now cut off the little bit that nepacLL has
wc <- nepacLL[nepacLL$Y>=34&nepacLL$Y<65 & nepacLL$X< -112&nepacLL$X>= -145,]
wcIslands <- wc[wc$PID!=1 & wc$Y>34.1,]                   #grab islands from nepacLL, but omit the southernmost that are cutoff
tmp <- rbind(wc[wc$PID==1,],tmp)
tmp <- rbind(c(1,1,tmp[1,"X"],max(tmp$Y)),tmp)               #add in a point offscreen to connect the polygon
tmp <- rbind(c(1,1,-100,max(tmp$Y)),tmp)               #add in a point offscreen to connect the polygon
tmp <- rbind(c(1,1,-100,min(tmp$Y)),tmp)               #add in a point offscreen to connect the polygon
tmp$POS <- 1:nrow(tmp)
#just in case some of the PID numbers are the same, renumber them
tmpPID <- as.numeric(factor(tmpIslands$PID))+1            #this renumbers the PID to be between 1 and whatever, and keeps the correct groups, example: x <- c(11,12,12,13,13,13,14,14,14,14,15,15,15,15,15); as.numeric(factor(x))
wcPID <- as.numeric(factor(wcIslands$PID))+max(tmpPID)    #add to it so that other PIDs are not the same
tmpIslands$PID <- tmpPID
wcIslands$PID <- wcPID
westCoastLL <- rbind(tmp,wcIslands,tmpIslands)
#plotMap(westCoastLL)
print("Created westCoastLL")
#plotMap(westCoastLL,xlim=c(-122.3,-118.2),ylim=c(31.95,35),col=gray(0.5))

rm(tmp,wc,tmpIslands,wcIslands,tmpPID,wcPID)





wc <- westCoastLL[westCoastLL$PID==1,c("X","Y")]
wc <- wc[!(wc$Y>46.9 & wc$X> -123.3),]                                 #remove puget sound border
#PSMFC boundaries
PSMFC <- data.frame(boundary=c("Start","1A","1B","1C","2A","2B","2C","3A"),lat=c(32.5,36,40.5,42,42.933333,44.3,45.7667,47.3333))

PSMFCareasPBS <- PSMFCareaNames <- NULL
for(i in 1:(nrow(PSMFC)-1)) {
    S <- i; N <- i+1
    tmp <- wc[wc$Y>=PSMFC$lat[S] & wc$Y<=PSMFC$lat[N],c("X","Y")]
    tmp <- rbind(c(-140,PSMFC$lat[N]),c(tmp$X[1],PSMFC$lat[N]),tmp,c(tmp$X[nrow(tmp)],PSMFC$lat[S]),c(-140,PSMFC$lat[S]))
    tmp <- data.frame(PID=rep(i,nrow(tmp)),POS=1:nrow(tmp),X=tmp$X,Y=tmp$Y)
    PSMFCareasPBS <- as.PolySet(rbind(tmp,PSMFCareasPBS),projection="LL")
    PSMFCareaNames <- as.PolyData(rbind(data.frame(PID=i,label=PSMFC$boundary[N],X=-128,Y=(PSMFC$lat[S]+PSMFC$lat[N])/2),PSMFCareaNames),projection="LL")
}

PSMFC <- data.frame(boundary=c("Start","1A","1B","1C","2A","2B","2C","3A","3S/3B"),lat=c(30,36,40.5,42,42.933333,44.3,45.7667,47.3333,48.7))
xlabel <- c(-123,-127,-128,-128,-128,-128,-128,-125.2)
PSMFCareasPBS <- PSMFCareaNames <- NULL
for(i in 1:(nrow(PSMFC)-1)) {
    S <- i; N <- i+1
    tmp <- wc[wc$Y>=PSMFC$lat[S] & wc$Y<=PSMFC$lat[N],c("X","Y")]
    tmp <- rbind(c(-140,PSMFC$lat[N]),c(tmp$X[1],PSMFC$lat[N]),tmp,c(tmp$X[nrow(tmp)],PSMFC$lat[S]),c(-140,PSMFC$lat[S]))
    tmp <- data.frame(PID=rep(i,nrow(tmp)),POS=1:nrow(tmp),X=tmp$X,Y=tmp$Y)
    PSMFCareasPBS <- as.PolySet(rbind(tmp,PSMFCareasPBS),projection="LL")
    PSMFCareaNames <- as.PolyData(rbind(data.frame(PID=i,label=PSMFC$boundary[N],X=xlabel[i],Y=(PSMFC$lat[S]+PSMFC$lat[N])/2),PSMFCareaNames),projection="LL")
}
PSMFCareaNames[1,c("X","Y")] <- c(-126.2,47.7)
#PSMFCareasPBS <- joinPolys(EEZ,PSMFCareasPBS,operation="INT")
#plotMap(westCoastLL,xlim=c(-129.5,-116.2),ylim=c(30.9,49),col=gray(0.7),tck=-0.02,main="PSMFC areas")
#addLines(WCstatesInlandPBS)
#addPolys(PSMFCareasPBS,col=gray(0.9),lwd=2)
#addLabels(PSMFCareaNames)



#INPFC boundaries (no EEZ)
INPFC <- data.frame(boundary=c("Start","CP","MT","EK","CL"),lat=c(30,36,40.5,43,47.5))

INPFCareasPBS <- INPFCareaNames <- NULL  #I could probably use joinPolys
for(i in 1:(nrow(INPFC)-1)) {
    S <- i; N <- i+1
    tmp <- wc[wc$Y>=INPFC$lat[S] & wc$Y<=INPFC$lat[N],c("X","Y")]
    tmp <- rbind(c(-140,INPFC$lat[N]),c(tmp$X[1],INPFC$lat[N]),tmp,c(tmp$X[nrow(tmp)],INPFC$lat[S]),c(-140,INPFC$lat[S]))
    tmp <- data.frame(PID=rep(i,nrow(tmp)),POS=1:nrow(tmp),X=tmp$X,Y=tmp$Y)
    INPFCareasPBS <- as.PolySet(rbind(tmp,INPFCareasPBS),projection="LL")
    INPFCareaNames <- as.PolyData(rbind(data.frame(PID=i,label=INPFC$boundary[N],X=-128,Y=(INPFC$lat[S]+INPFC$lat[N])/2),INPFCareaNames),projection="LL")
}


#INPFC boundaries (with EEZ)
INPFC <- data.frame(boundary=c("Start","CP","MT","EK","CL","VN"),lat=c(30.5,36,40.5,43,47.5,48.7))

INPFCareasPBS <- INPFCareaNames <- NULL  #I could probably use joinPolys
xlabel <- c(-123,-127,-128,-128,-125.2)
for(i in 1:(nrow(INPFC)-1)) {
    S <- i; N <- i+1
    tmp <- wc[wc$Y>=INPFC$lat[S] & wc$Y<=INPFC$lat[N],c("X","Y")]
    tmp <- rbind(c(-140,INPFC$lat[N]),c(tmp$X[1],INPFC$lat[N]),tmp,c(tmp$X[nrow(tmp)],INPFC$lat[S]),c(-140,INPFC$lat[S]))
    tmp <- data.frame(PID=rep(i,nrow(tmp)),POS=1:nrow(tmp),X=tmp$X,Y=tmp$Y)
    INPFCareasPBS <- as.PolySet(rbind(tmp,INPFCareasPBS),projection="LL")
    INPFCareaNames <- as.PolyData(rbind(data.frame(PID=i,label=INPFC$boundary[N],X=xlabel[i],Y=(INPFC$lat[S]+INPFC$lat[N])/2),INPFCareaNames),projection="LL")
}
INPFCareaNames[1,c("X","Y")] <- c(-126.15,47.75)
#INPFCareasPBS <- joinPolys(EEZ,INPFCareasPBS,operation="INT")
#plotMap(westCoastLL,xlim=c(-129.5,-116.2),ylim=c(30.9,49),col=gray(0.7),tck=-0.02,main="INPFC areas")
#addLines(WCstatesInlandPBS)
#addPolys(INPFCareasPBS,col=gray(0.9),lwd=2)
#addLabels(INPFCareaNames)


print("Created PSMFCareasPBS, PSMFCareaNames, INPFCareasPBS, and INPFCareaNames")


#par(mfrow=c(1,2))
#plotMap(westCoastLL,xlim=c(-128.5,-115.9),ylim=c(31.95,49),col=gray(0.7),tck=-0.02,main="PSMFC areas")
#addLines(WCstatesInlandPBS)
#addPolys(PSMFCareasPBS,col=gray(0.9),lwd=2)
#addLabels(PSMFCareaNames)

#plotMap(westCoastLL,xlim=c(-128.5,-115.9),ylim=c(31.95,49),col=gray(0.7),tck=-0.02,main="INPFC areas")
#addLines(WCstatesInlandPBS)
#addPolys(INPFCareasPBS,col=gray(0.9),lwd=2)
#addLabels(INPFCareaNames)


#par(mfrow=c(1,2))
#plotMap(westCoastLL,xlim=c(-128.3,-122.9),ylim=c(41.5,49),col=gray(0.7),tck=-0.02,main="PSMFC areas and INPFC differences")
#addLines(WCstatesInlandPBS)
#addPolys(PSMFCareasPBS,col=gray(0.9),lwd=3)
#addLabels(PSMFCareaNames)
#tmp <- joinPolys(PSMFCareasPBS[PSMFCareasPBS$PID==7,],INPFCareasPBS[INPFCareasPBS$PID==4,],operation="XOR")
#addPolys(tmp[tmp$Y>45.8,],density=15,col=1)
#tmp <- joinPolys(PSMFCareasPBS[PSMFCareasPBS$PID==5,],INPFCareasPBS[INPFCareasPBS$PID==4,],operation="XOR")
#addPolys(tmp[tmp$Y<44,],density=15,col=1)

#plotMap(westCoastLL,xlim=c(-128.3,-122.9),ylim=c(41.5,49),col=gray(0.7),tck=-0.02,main="INPFC areas and PSMFC differences")
#addLines(WCstatesInlandPBS)
#addPolys(INPFCareasPBS,col=gray(0.9),lwd=3)
#addLabels(INPFCareaNames)
#tmp <- joinPolys(PSMFCareasPBS[PSMFCareasPBS$PID==7,],INPFCareasPBS[INPFCareasPBS$PID==4,],operation="XOR")
#addPolys(tmp[tmp$Y>45.8,],density=15,col=1)
#tmp <- joinPolys(PSMFCareasPBS[PSMFCareasPBS$PID==5,],INPFCareasPBS[INPFCareasPBS$PID==4,],operation="XOR")
#addPolys(tmp[tmp$Y<44,],density=15,col=1)




EEZ2 <- importShapefile("C:\\Mapping\\Shapefiles\\Pacific_EEZ_segmented_dd.shp",readDBF=T)
EEZ <- importShapefile("C:\\Mapping\\Shapefiles\\Pacific_EEZ_dissolve_dd.shp",readDBF=T)
tmp <- importShapefile("C:\\Mapping\\Shapefiles\\DepthContours_4Survey_dd.shp",readDBF=T)
depth30f <- tmp[tmp$PID==1,]
depth100f <- tmp[tmp$PID==2,]
depth300f <- tmp[tmp$PID==3,]
depth700f <- tmp[tmp$PID==4,]

print("Created depth contours depth30f, depth100f, depth300f, and depth700f")


INPFCareasEEZ.PBS <- joinPolys(EEZ,INPFCareasPBS,operation="INT")

tmp <- importShapefile("C:\\Mapping\\Shapefiles\\contour4allan_ETsmoothBezier5_diss.shp",readDBF=T)
depth30m <- tmp[tmp$PID==3,]
depth50m <- tmp[tmp$PID==2,]
depth1500m <- tmp[tmp$PID==1,]

print("Created depth contours depth30m, depth50m, and depth1500m")


#tmp <- importShapefile("C:\\Mapping\\Shapefiles\\INPFCbox.shp",readDBF=T)
#tmp <- read.shape("C:\\Mapping\\Shapefiles\\INPFCbox.shp")
#tmp <- importShapefile("C:\\Mapping\\Shapefiles\\Merged_WCshore_Sept2007.shp")
#tmp <- importShapefile("C:\\Mapping\\Shapefiles\\EEZ_wShoreline.shp")
#plotMap(westCoastLL,col=gray(0.8))
#addLines(EEZ2)
#plotMap(westCoastLL,col=gray(0.8))
#addLines(EEZ)
#plotMap(westCoastLL,xlim=c(-130.5,-110.9),ylim=c(30.45,49.5),col=gray(0.9))
#addLines(EEZ2)
#addLines(WCstatesInlandPBS)

#layout(matrix(c(1,2),nrow=1),widths=c(2,1))
#plotMap(westCoastLL,xlim=c(-125.5,-115.9),ylim=c(31.45,42.5),col="green3",tck=c(-0.01),cex=1.0)
#addLines(EEZ2)
#addLines(WCstatesInlandPBS)
#addLines(tmp,col=gray(0.7))
#plotMap(westCoastLL,xlim=c(-126.3,-120.9),ylim=c(41.9,49.2),col="green3",tck=c(-0.01),cex=1.0)#,mgp=c(1.9,0.7,0))
#addLines(EEZ2)
#addLines(tmp,col=gray(0.7))
#addLines(WCstatesInlandPBS)




#tmp$X <- tmp$X/100
#tmp$Y <- tmp$Y/100

#attr(tmp, "projection") <- "UTM"
#attr(tmp, "zone") <- 10
#--- convert and plot the result
#tmpUTM <- convUL(tmp)
#plotMap(nepacUTM,xlim=c(-500,1500),ylim=c(3800,5400),col=gray(0.7),tck=-0.02)
#plot(tmp$X/100,tmp$Y/100)

#plotMap(nepacLL)


tmp <- joinPolys(EEZ,INPFCareasPBS[INPFCareasPBS$PID==4,],operation="INT")
tmp <- joinPolys(EEZ,INPFCareasPBS,operation="INT")
