dont source

setwd("C:/NOAA2015/Widow/Data/Catches/WA")

early <- read.csv("WashFishTix1935-1969.csv")
early$mt <- round(early$Species_RoundWeight*0.00045359237,2)
early$FTL_SPID <- as.character(early$FTL_SPID)
early$FTL_SPID <- gsub("\\s","",early$FTL_SPID)

early$PacFIN_GRID <- as.character(early$PacFIN_GRID)
early$PacFIN_GRID <- gsub("\\s","",early$PacFIN_GRID)
early$gear <- NA
early$gear[early$PacFIN_GRID%in%c("SST","SHT","PWT","DST")] <- "ShrimpTrawl"
early$gear[early$PacFIN_GRID%in%c("RLT","GFT","GFS","GFL","FTS","FFT")] <- "BottomTrawl"
early$gear[early$PacFIN_GRID%in%c("OTW","MDT")] <- "MidwaterTrawl"
early$gear[early$PacFIN_GRID%in%c("PRT","DNT","BMT")] <- "BottomTrawl"   #"MiscTrawl"
early$gear[early$PacFIN_GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- "Pot"
early$gear[early$PacFIN_GRID%in%c("JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
early$gear[early$PacFIN_GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
early$gear[early$PacFIN_GRID%in%c("DVG","USP","")] <- "Other"

table(early$gear,useNA="ifany")
tapply(early$mt,early$gear,sum,na.rm=T)
head(early)

tapply(early$mt,list(early$YearLanded_Orig,early$PacFIN_SPID),sum)
tapply(early$mt,list(early$YearLanded_Orig,early$PacFIN_SPID,early$gear),sum)

unique(early$AreaRegion)
sort(tapply(early$mt,early$AreaRegion,sum,na.rm=T))
early$area <- NA
early$area[early$AreaRegion%in%c("Tillamook Head/Cape Elizabeth","Cape Elizabeth/Cape Shoalwater",
									 "Cape Ellizabeth","Cape Johnson/Cape Elizabeth","Cape Shoalwater/Cape Johnson",
									 "Cape Shoalwater/Point Grenville","Cape Shoalwater/Tillamook Head",
									 "Grays Harbor - Outside","La Push - Outside","Westport - Outer/Outside",
									 "Willapa Bay/Cape Johnson",
									 "Oregon - Waters","California /Coast/Waters",
									 "Tillamook Head - South")]   <- "OuterUS"
early$area[early$AreaRegion%in%c("Sail Rock - Tatosh","Suez River - (Sooes River)","Willapa Bay",
									 "Quinault - General","Makah Bay")]  <- "Inside"
early$area[early$AreaRegion%in%c("Cape Flattery","Outside State Territorial Waters",
									 "Tillamook Head/Barkley Sound")]  <- "UsCan"

table(early$area,useNA="ifany")
tapply(early$mt,early$area,sum,na.rm=T)
tapply(early$mt,list(early$YearLanded_Orig,early$PacFIN_SPID,early$area,early$gear),sum,na.rm=T)

early$gear2 <- factor(early$gear,levels=c("BottomTrawl","HnL"))
early$area2 <- factor(early$area,levels=c("OuterUS","UsCan"))

tapply(early$mt,early$gear2,sum,na.rm=T)
tapply(early$mt,early$area2,sum,na.rm=T)
tapply(early$mt,list(early$YearLanded_Orig,early$FTL_SPID,early$area2,early$gear2),sum)

x <- aggregate(early$mt,list(early$YearLanded_Orig,early$gear2,early$area2,early$FTL_SPID),sum)
names(x) <- c("Year","gear2","area2","SPID","mt")
x <- x[order(x$Year,x$gear2,x$area2,x$SPID),]
write.csv(x,"WaCatchesFromTeresaFishTix.csv",row.names=F)




###############################################
### Look at logbooks for MDT
logBook <- NULL
colsToSave <- c("TRIP_ID","AGID","RYEAR","RMONTH","RPORT","RPCID","DRVID","TOWNUM","ARID_PSMFC",
				"SET_LAT","SET_LONG","SET_TIME","DURATION","NET_TYPE","GRID","DEPTH1","DEPTH2","PACFIN_TARGET",
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
}
tapply(logBook$APOUNDS,list(logBook$RYEAR,logBook$GRID),sum,na.rm=T)
x <- logBook[logBook$SPID=="WDW1" & logBook$AGID=="W",]
y <- tapply(x$APOUNDS,list(x$RYEAR,x$GRID),sum,na.rm=T)
cbind(y[,"MDT"]/apply(y,1,sum,na.rm=T))