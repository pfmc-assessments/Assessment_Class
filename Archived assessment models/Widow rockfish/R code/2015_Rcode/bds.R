dont source
setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
#library(ExtractR)
source("Rcode/functions/Functions.R")

year1 <- 1970
year2 <- 2014

if(F) {
    #Get PacFIN bds data and save results
    bds.age <- queryDB(queryFilename="bds.age",db="PACFIN",uid="hicksa",sp="WDOW", start=year1,end=year2)
    save(bds.age,file="extractedData/bds.age.Rdat")
    file.copy("extractedData/bds.age.Rdat",paste("extractedData/bds.age_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
    bds.allsp.cluster <- queryDB(queryFilename="bds.allsp.cluster",db="PACFIN",uid="hicksa",sp="WDOW",start=year1,end=year2)
    save(bds.allsp.cluster,file="extractedData/bds.allsp.cluster.Rdat")
    file.copy("extractedData/bds.allsp.cluster.Rdat",paste("extractedData/bds.allsp.cluster_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
    bds.fish <- queryDB(queryFilename="bds.fish",db="PACFIN",uid="hicksa",sp="WDOW",start=year1,end=year2,asis=T) #c(rep(F,11),T,rep(F,32)))   #I read in the FISH_LENGTH_TYPE as.is=T because if it is all "F" it is changed to logical
	bds.fish$FISH_LENGTH <- as.numeric(bds.fish$FISH_LENGTH)
    save(bds.fish,file="extractedData/bds.fish.Rdat")
    file.copy("extractedData/bds.fish.Rdat",paste("extractedData/bds.fish_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
    bds.sp.cluster <- queryDB(queryFilename="bds.sp.cluster",db="PACFIN",uid="hicksa",sp="WDOW",start=year1,end=year2,querydir="sql/")
    save(bds.sp.cluster,file="extractedData/bds.sp.cluster.Rdat")
    file.copy("extractedData/bds.sp.cluster.Rdat",paste("extractedData/bds.sp.cluster_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

	#Load the database queries
	setwd("C:/NOAA2015/Widow/Data")
	source("Rcode/functions/Functions.R")
	source("Rcode/functions/workupPacFinTablesBDS.R")

	load("extractedData/bds.age.Rdat")
	load("extractedData/bds.allsp.cluster.Rdat")
	load("extractedData/bds.fish.Rdat")
	load("extractedData/bds.sp.cluster.Rdat")

	bds.fish.worked <- workupPacFinTablesBDS(bds_fish=bds.fish,age_temp=bds.age,sp_cluster=bds.sp.cluster,all_cluster=bds.allsp.cluster)
	bds.fish.worked$SEX <- factor(bds.fish.worked$SEX)
	save(bds.fish.worked,file="PacFIN/PacFIN.WDOW.bds.Rdat")
    file.copy("PacFIN/PacFIN.WDOW.bds.Rdat",paste("PacFIN/PacFIN.WDOW.bds_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}

###########
#   START HERE if you have already done the extractions and worked up data
###########

setwd("C:/NOAA2015/Widow/Data")
load("PacFIN/PacFIN.WDOW.bds.Rdat")
dat <- bds.fish.worked
rm(bds.fish.worked)

table(dat$GRID)

dat$gear <- NA
dat$gear[dat$GRID%in%c("SST","SHT","PWT","DST","DSG")] <- "ShrimpTrawl"
dat$gear[dat$GRID%in%c("RLT","GFT","GFS","GFL","FTS","FFT","BTT","BMT","56")] <- "BottomTrawl"
dat$gear[dat$GRID%in%c("OTW","MDT","54")] <- "MidwaterTrawl"
dat$gear[dat$GRID%in%c("PRT","DNT")] <- "MiscTrawl"
dat$gear[dat$GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- "Pot"
dat$gear[dat$GRID%in%c("HKL","JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
dat$gear[dat$GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
dat$gear[dat$GRID%in%c("DVG","USP","MPT","UNK","XXX")] <- "Other"   #MPT is CP midwater

table(dat$gear,useNA="ifany")




ind <- !is.na(dat$FISH_LENGTH)
table(dat$SAMPLE_YEAR[ind],dat$SAMPLE_AGENCY[ind])

ind <- !is.na(dat$FISH_AGE_YEARS_FINAL)
table(dat$SAMPLE_YEAR[ind],dat$SAMPLE_AGENCY[ind])
table(dat$DATE_AGED,dat$AGED_BY,useNA="ifany")
table(substring(dat$DATE_AGED,1,4),dat$AGED_BY,useNA="ifany")

table(dat$SOURCE_AGID,dat$SAMPLE_AGENCY)  #PW in sample_agency is Pacific Whiting
table(dat$SAMPLE_YEAR,dat$SAMPLE_AGENCY)
table(dat$SAMPLE_YEAR,dat$gear,dat$SAMPLE_AGENCY)

#investigate state specific length comps to see how different they are
lenBoxplots <- function(x,...) {
	x.yr <- split(x,x$SAMPLE_YEAR)
	x.yr.state <- lapply(x.yr,function(xx){split(xx$FISH_LENGTH,factor(xx$SAMPLE_AGENCY,levels=c("CA","OR","PW","W")))})
	lapply(x.yr.state,names)
	length(x.yr.state)
	makePlot <- function(yy,...) {
		boxplot(yy,...)
		text(1:4,650,unlist(lapply(yy,length)),cex=0.8)
	}
	y <- lapply(x.yr.state,makePlot,...)
	invisible(y)
}
windows(height=10,width=15)
par(mfrow=c(5,8),mar=c(3,3,1,1))
lenBoxplots(dat[dat$gear=="BottomTrawl",],ylim=range(dat$FISH_LENGTH,na.rm=T))

lenBoxplots(dat[dat$gear=="MidwaterTrawl",],ylim=range(dat$FISH_LENGTH,na.rm=T))
lenBoxplots(dat[dat$gear=="HnL",],ylim=range(dat$FISH_LENGTH,na.rm=T))
lenBoxplots(dat[dat$gear=="Net",],ylim=range(dat$FISH_LENGTH,na.rm=T))


table(!is.na(dat$age1)&dat$age1>0,!is.na(dat$age2)&dat$age2>0,!is.na(dat$age3)&dat$age3>0)
doubleReads <- dat[!is.na(dat$age1)&dat$age1>0 & !is.na(dat$age2)&dat$age2>0,]
head(doubleReads)
table(doubleReads$SOURCE_AGID)
table(doubleReads$SAMPLE_YEAR,doubleReads$SOURCE_AGID)

table(dat$AGE_METHOD)



#Length at age
plot(dat$FISH_AGE_YEARS_FINAL,dat$FISH_LENGTH,pch=20)

par(mfrow=c(1,1))
plot(dat$FISH_AGE_YEARS_FINAL[dat$SEX=="F"],dat$FISH_LENGTH[dat$SEX=="F"],pch=20,col=rgb(1,0,0,0.6))
points(dat$FISH_AGE_YEARS_FINAL[dat$SEX=="M"],dat$FISH_LENGTH[dat$SEX=="M"],pch=20,col=rgb(0,0,1,0.6))


