#investigate logbook data and create a CPUE index for widow
dont source

setwd("C://NOAA2015/Widow/Data")
logBookDir <- "C://NOAA2015/Data/Logbooks"

x <- read.csv(file.path(logBookDir,"LogBook.1982.csv"))
table(x$GRID)
table(x$SPID)
wdw <- x[x$SPID=="WDW1",]
wdw$mt <- wdw$APOUNDS*0.00045359237
table(wdw$DRVID,wdw$mt>0)
y2 <- tapply(wdw$mt,as.character(wdw$DRVID),sum,na.rm=T)


# filterOn <- function(x,vars,vals) {
	# LOOK at observer code because I do something similar
# 	if(length(vals))
# 	y <- x[]
# }


minYrCatch <- 0
yrs <- 1981:1998
yrs <- 1982:1998
logBook <- NULL
colsToSave <- c("TRIP_ID","AGID","RYEAR","RMONTH","RPORT","RPCID","DRVID","TOWNUM","ARID_PSMFC",
				"SET_LAT","SET_LONG","SET_TIME","DURATION","NET_TYPE","DEPTH1","DEPTH2","PACFIN_TARGET",
				"SPID","APOUNDS","ADJ_TOWTIME")
for(i in 1:length(yrs)) {
	cat("Year",yrs[i],"\n")
	flush.console()
	x <- read.csv(paste(logBookDir,"/LogBook.",yrs[i],".csv",sep=""),as.is=T)
	#save a smaller set of the data for analysis
	#probably want to make same criteria as vessel selection below
	#or make sure to filter when filtering out vessels later
	#probably better to filter as soon as possible so it doesn't get too big
	logBook <- rbind(logBook,x[,colsToSave])

	#work with only widow catch and any other criteria
	wdw <- x[x$SPID=="WDW1",]
wdw <- wdw[wdw$SET_LAT>=42 & wdw$SET_LAT<=46,]
wdw <- wdw[wdw$GRID != "MDT",]
	wdw$mt <- wdw$APOUNDS*0.00045359237

	y <- tapply(wdw$mt,as.character(wdw$DRVID),sum,na.rm=T)
	y <- y[y>minYrCatch]
	yy <- data.frame(vessel=names(y),mt=y)
	names(yy)[length(names(yy))] <- paste("Yr",yrs[i],sep="")
	if(i > 1) {
		vessels <- merge(vessels,yy,by="vessel",all=T)
	}else {
		vessels <- yy
	}
}
origVessels <- vessels

#vessel criteria
numYrs <- 5
totCatchMin <- 50
yrCatchMin <- 5
numYrsWithCatchMin <- 7

numYrs <- 3
totCatchMin <- 20
yrCatchMin <- 1
numYrsWithCatchMin <- 5

vessels <- origVessels
#first make vessle ID the rowname
rownames(vessels) <- vessels[,1]
vessels <- vessels[,-1]

#numYrs
vesselNumYrs <- apply(!is.na(vessels),1,sum)
vessels <- vessels[vesselNumYrs>=numYrs,]
#totCatchMin
vesselTotCatch <- apply(vessels,1,sum,na.rm=T)
vessels <- vessels[vesselTotCatch>=totCatchMin,]
#YrsWithCatchMin
vesselCatchMin <- apply(vessels>=yrCatchMin,1,sum,na.rm=T)
vessels <- vessels[vesselCatchMin>=numYrsWithCatchMin,]

round(vessels)
#make a plot of vessel coverage
windows(height=5,width=6.5)
yrs <- as.numeric(substring(colnames(vessels),3))
legX <- rep(min(yrs),3)
legY <- seq(nrow(vessels),nrow(vessels)-4,-2)
legVal <- round(max(vessels,na.rm=T)/c(1,2,5),-1)
symbols(c(legX,rep(yrs,each=nrow(vessels))),
	    c(legY,rep(1:nrow(vessels),ncol(vessels))),
	    circle=sqrt(c(legVal,unlist(vessels))),
	    inch=0.04,bg="blue",fg=NA,
	    xlab="Year",ylab="Vessel",xlim=range(yrs))
text(legX,legY,paste(legVal,"t"),cex=0.8,pos=2,offset=-2)


#CPUE analysis
origLogBook <- logBook
#filter out same as for vessels
#first choose only selected vessels
logBook <- logBook[as.character(logBook$DRVID)%in%rownames(vessels),]
#other filters
logBook <- logBook[logBook$SET_LAT>=42 & logBook$SET_LAT<=46,]
logBook <- logBook[!is.na(logBook$TRIP_ID),]
logBook$Trip_Tow <- paste(logBook$TRIP_ID,logBook$TOWNUM,sep="_")

#positive catches of widow
pos <- logBook[logBook$SPID=="WDW1",]
pos$CPUE <- pos$APOUNDS/pos$DURATION   #may want to try catch per tow rather than distance
boxplot(split(log(pos$CPUE),pos$RYEAR),ylab="log(LBS/HR)",xlab="Year")




#proportion postive
tmp <- split(logBook$SPID,logBook$Trip_Tow)
tmp <- unlist(lapply(tmp,function(x){any(x=="WDW1")}))
pos <- logBook[!duplicated(logBook$Trip_Tow),]
pos <- cbind(pos[match(names(tmp),pos$Trip_Tow),],tmp)
tmp <- table(pos$RYEAR,pos$tmp)
propPos <- tmp[,2]/(tmp[,1]+tmp[,2])
plot(as.numeric(names(propPos)), propPos, type="b", pch=16, ylab="Proportion positive", xlab="Year")
