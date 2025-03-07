dont source
#see also Widow.R in the ObserverDatabase folder

setwd("C:/NOAA2015/Widow/Data")

#Pikitch data
load("Discards\\Pikitch\\Pikitch.et.al.WDOW.Discard.Rates.x.over.y.with.MDT 17 Jun 2015.dmp")
dat <- Pikitch.et.al.WDOW.Discard.Rates.x.over.y
rm(Pikitch.et.al.WDOW.Discard.Rates.x.over.y)

table(dat$Year)
table(dat$Year,dat$Study)
table(dat$Year,dat$Gear.Type)
table(dat$Year,dat$Areas)

#I believe that "2B 2C 3A" is all areas combined
gft <- dat[dat$Gear.Type=="GFT" & dat$Areas=="2B 2C 3A",]

land <- read.csv("Catches/LandingsByFleet.csv")
land <- land[land$Year%in%gft$Year,c("Year","BottomTrawl")]
dis <- land$BottomTrawl * gft$DiscardRate.Sp.Wt.Wgting
dis.se <- sqrt(land$BottomTrawl^2 * gft$SD.DiscardRate.Sp.Wt.Wgting^2/gft$Num.Study.Tows)
dis.seLog <- sqrt(log((dis.se/dis)^2+1))

#ask Wallace if my interpretation of SD is correct
#Think about correct clauclation of sd, se, and cv for a total


Pikitch <- data.frame(Year=gft$Year, Season=1, Fleet=1,
			DiscardLog=log(dis)-dis.seLog^2/2,
			SElog=dis.seLog,
			PtEst=dis,
			Median=dis*exp(-dis.seLog^2/2),
			CV=dis.se/dis,
			Rate=dis/(dis+land$BottomTrawl),
			Ratio=dis/(land$BottomTrawl))

mdt <- dat[dat$Gear.Type=="MDT" & dat$Areas=="2B 2C 3A",]
mdt <- mdt[!is.na(mdt$DiscardRate.Sp.Wt.Wgting),]

land <- read.csv("Catches/LandingsByFleet.csv")
land <- land[land$Year%in%mdt$Year,c("Year","MidwaterTrawl")]
dis <- land$MidwaterTrawl * mdt$DiscardRate.Sp.Wt.Wgting
dis.se <- sqrt(land$MidwaterTrawl^2 * mdt$SD.DiscardRate.Sp.Wt.Wgting^2/mdt$Num.Study.Tows)
dis.seLog <- sqrt(log((dis.se/dis)^2+1))

#ask Wallace if my interpretation of SD is correct
#Think about correct clauclation of sd, se, and cv for a total


Pikitch <- rbind(Pikitch,data.frame(Year=mdt$Year, Season=1, Fleet=2,
			DiscardLog=log(dis)-dis.seLog^2/2,
			SElog=dis.seLog,
			PtEst=dis,
			Median=dis*exp(-dis.seLog^2/2),
			CV=dis.se/dis,
			Rate=dis/(dis+land$MidwaterTrawl),
			Ratio=dis/(land$MidwaterTrawl)))

write.csv(Pikitch,file="Discards/Pikitch.csv",row.names=F)



#EDCP
load("Discards\\EDCP\\EDCP DATA - R\\R - John Wallace\\EDCP; obs_tckt & obs_ttow.dmp")
head(obs_tckt)
head(obs_ttow)

dat <- obs_tckt[,c("GEAR","LAND_DATE","NTOWS","NDIS_REC","NMKT","PMFCAREA","WDOW_DIS","WDOW_TCKT")]
dat<- dat[(dat$WDOW_DIS+dat$WDOW_TCKT)>0,]
dat$Year <- as.numeric(substring(dat$LAND_DATE,nchar(dat$LAND_DATE)-3,nchar(dat$LAND_DATE)))
dat$Rate <- dat$WDOW_DIS/(dat$WDOW_DIS+dat$WDOW_TCKT)

table(dat$GEAR)
dat$fleet <- NA
dat$fleet[dat$GEAR==360] <- "Midwater Trawl"
dat$fleet[dat$GEAR==390] <- "Bottom Trawl"
#GEAR==380 is shrimp double

dat <- dat[!is.na(dat$fleet),]

table(dat$Year,dat$fleet) #number of trips observed?
tapply(dat$NTOWS,list(dat$Year,dat$fleet),sum)  #number of tows observed
hist(dat$Rate)

land <- read.csv("Catches/LandingsByFleet.csv")
land <- land[land$Year%in%dat$Year,c("Year","BottomTrawl","MidwaterTrawl")]
names(land) <- c("Year","Bottom Trawl","Midwater Trawl")

#aggregated by trip data, not tow specific data
x <- split(dat,paste(dat$Year,dat$fleet))
xx <- data.frame(Year=NULL,Fleet=NULL,Ntrips=NULL,Ntows=NULL,NtowsObs=NULL,Nbio=NULL,dis=NULL,ret=NULL,rate=NULL,ratio=NULL)
for(i in 1:length(x)) {
	#individual ratios to get se of ratio
	dis <- x[[i]]$WDOW_DIS
	ret <- x[[i]]$WDOW_TCKT
	ratio <- dis/ret
	ratio[ratio>10] <- 10   #max Expansion of 10
	se.ratio <- sd(ratio)/sqrt(length(ratio))
	cv.dis <- (sd(dis)/sqrt(length(dis)))/mean(dis)

	ret <- land[land$Year==x[[i]]$Year[1],x[[i]]$fleet[1]] #ttoal retention (landings)
	ratio <- sum(dis)/ret    #an annual ratio
	ratio[ratio>10] <- 10   #max Expansion of 10
	dis <- ret * ratio
	dis.seLog <- sqrt(log((cv.dis)^2+1))
	xx <- rbind(xx, data.frame(
		Year = x[[i]]$Year[1],
		Season = 1,
		Fleet = x[[i]]$fleet[1],
		DiscardLog = log(dis)-dis.seLog^2/2,
		SElog=dis.seLog,
		PtEst=dis,
		Median=dis*exp(-dis.seLog^2/2),
		CV=cv.dis,
		Rate = dis/(dis+ret),
		Ratio = dis/ret,
		Ntrips = nrow(x[[i]]),
		Ntows = sum(x[[i]]$NTOWS),
		NtowsObs = sum(x[[i]]$NDIS_REC),
		Nbio = sum(x[[i]]$NMKT)))
}
xx <- xx[order(xx$Fleet),]

write.csv(xx,file="Discards/EDCP.csv",row.names=F)


#WCGOP (also see widow.r in ObserverDatabase folder)
setwd("C:/NOAA2015/ObserverDatabase/Widow")
load("RdatFiles/combinedDiscardsBoot.Rdat")

wcgop.fn <- function(res,yrs,stratName,Fleet=NA,fleetName=NULL,Season=1,Yaxis=0) {
	forSS <- NULL
	for(i in 1:length(yrs)) {
		ptEst <- res[[as.character(yrs[i])]][[stratName]][["PtEst"]]
		dis <- res[[as.character(yrs[i])]][[stratName]][["dis"]]
		hist(dis,ylab="",xlab="",main=yrs[i],yaxt="n")
		if(i == Yaxis) {mtext(fleetName,side=2,outer=F,line=0.3)}
		forSS <- rbind(forSS, data.frame(Year=yrs[i],Season=1,Fleet=Fleet,
					DiscardLog=log(median(dis)),
					SElog=sqrt(log((sd(dis)/mean(dis))^2+1)),
					PtEst = ptEst["discard"],
					Median = median(dis),
					CV = sd(dis)/mean(dis),
					Rate = ptEst["discard"]/(ptEst["discard"]+ptEst["retained"]),
					Ratio = ptEst["discard"]/ptEst["retained"]
				))
	}
#	mtext(paste("Total discards for",fleetName,"(t)"),side=1,outer=T)
	return(forSS)
}



windows(width=6.5,height=9)
par(mfrow=c(7,3),mar=c(3,1,3,0),oma=c(1.5,3,0,0))
#Bottom Trawl
forSS <- wcgop.fn(res,2002:2010,"Bottom Trawl.CA:OR:WA.FALSE",Fleet=1,fleetName="Bottom Trawl",Season=1,Yaxis=4)
forSS <- rbind(forSS, data.frame(Year=2011:2013,Season=1,Fleet=1,
			DiscardLog=NA, SElog=NA, PtEst=NA, Median=NA, CV=NA, Rate=NA, Ratio=NA))
#Using a maxExpansion of 10 resulted in the exact same log(median discard), but smaller CV for some years.
#Midwater Trawl
forSS <- rbind(forSS,wcgop.fn(res,2002,"Midwater Trawl.CA:OR:WA.FALSE",Fleet=2,fleetName="Midwater\nTrawl",Season=1,Yaxis=1))
mtext("Midwater\nTrawl",side=2)
forSS <- rbind(forSS, data.frame(Year=2012:2013,Season=1,Fleet=2,
			DiscardLog=NA, SElog=NA, PtEst=NA, Median=NA, CV=NA, Rate=NA, Ratio=NA))
plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
#Hook and Line
forSS <- rbind(forSS,wcgop.fn(res,c(2004:2008,2010:2013),"HookAndLine.CA:OR:WA.FALSE",Fleet=5,fleetName="Hook & Line",Season=1,Yaxis=4))
abline(h=2.85e4,xpd=NA)
abline(h=3.9e4,xpd=NA)
mtext("Total Discards (t)",side=1,outer=T)

write.csv(forSS,file="C:/NOAA2015/Widow/Data/Discards/WCGOP.csv",row.names=F)



############################################################
### LENGTHS
setwd("C:/NOAA2015/Widow/Data/Discards")

#Pikitch
load("Pikitch/Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm 19 Jun 2015.dmp")
dat <- Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm
rm(Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm)

load("Pikitch/Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex 22 Jun 2015.dmp")
datNoSex <- Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex
rm(Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex)

load("Pikitch/Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex.2C3A  22 Jun 2015.dmp")
datNoSex2c3a <- Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex.2C3A
rm(Pikitch.et.al.Widow.Lengths.wt.PacFIN.assm.No.Sex.2C3A)

#Use only female and male discards
dis <- dat[dat$Disposition=="Discarded",]
out <- matrix(NA,nrow=0,ncol=56,
	         dimnames=list(NULL,c("Year","Season","Fleet","Gender","Partition","nSamps",
	         	                    paste0("F",seq(8,56,2)),paste0("F",seq(8,56,2),".1"))))
out <- as.data.frame(out)
yrs <- sort(unique(dis$Year))
for(i in 1:length(yrs)) {
	tmpF <- dis[dis$Year==yrs[i] & dis$Sex=="Female",]
	Ls <- data.frame(F8=sum(tmpF[,c("L.6","L.8")]), tmpF[,paste("L",seq(10,54,2),sep=".")],
	      F56=sum(tmpF[,paste("L",seq(56,64,2),sep=".")]))
	names(Ls) <- gsub("L.","F",names(Ls))
	tmp <- data.frame(Year=tmpF$Year, Season=1, Fleet=1, Gender=1, Partition=1,
		               nSamps=tmpF$Num.Study.Lengths, Ls, Ls)
	out<-rbind(out,tmp)
	tmpM <- dis[dis$Year==yrs[i] & dis$Sex=="Male",]
	Ls <- data.frame(F8=sum(tmpM[,c("L.6","L.8")]), tmpM[,paste("L",seq(10,54,2),sep=".")],
	      F56=sum(tmpM[,paste("L",seq(56,64,2),sep=".")]))
	names(Ls) <- gsub("L.","F",names(Ls))
	tmp <- data.frame(Year=tmpM$Year, Season=1, Fleet=1, Gender=2, Partition=1,
		               nSamps=tmpM$Num.Study.Lengths, Ls, Ls)
	out<-rbind(out,tmp)
}
write.csv(out,file="PikitchDiscardLengthComps.csv",row.names=F)

## WCGOP discard length comps
## weighted by state discard totals
dis <- read.csv("Comps/WDOW.lengthComps.csv")

tmp <- tapply(dis$Weighted,list(dis$Year,dis$Lenbin,dis$Gear),sum)
minSize=8; maxSize=56
gears <- c("Bottom Trawl","H&L")
out <- NULL
for(i in gears) {
    tmp[,as.character(minSize),i] <- tmp[,as.character(minSize),i] + apply(tmp[,which(as.numeric(colnames(tmp[,,i]))<minSize),i],1,sum)
    tmp[,as.character(maxSize),i] <- tmp[,as.character(maxSize),i] + apply(tmp[,which(as.numeric(colnames(tmp[,,i]))>maxSize),i],1,sum)
    tmpF <- tmp[,which(as.numeric(colnames(tmp[,,i]))>=minSize & as.numeric(colnames(tmp[,,i]))<=maxSize),i]
    tmpF <- t(apply(tmpF,1,function(x){x/sum(x)}))
    tmpM <- tmpF
    colnames(tmpF) <- paste("F",colnames(tmpF),sep="")
    colnames(tmpM) <- paste("M",colnames(tmpM),sep="")
    out <- rbind(out,data.frame(year=rownames(tmpF),season=1,fleet=i,gender=0,partition=1,nSamps="nSamps",tmpF,tmpM))
}
write.csv(out,"Comps\\discardLFs.csv",row.names=F)


dis <- read.csv("Comps\\discardLFs.csv")
dis[dis==0] <- NA
dis <- split(dis,dis$fleet)
lens <- seq(8,56,2)
years <- 2004:2013
gears <- c("Bottom Trawl","H&L")
doPNG <- T
wd<-6.5;ht<-4.5
if(doPNG) {png("../../Writeup/Figures/Observer/discardLengthComps.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(2,1),mar=c(2,2,1.5,1),oma=c(1.5,1.5,0,0),las=1)
for(i in 1:2) {
    symbols(rep(dis[[i]]$year,length(lens)),rep(lens,each=length(dis[[i]]$year)),circles=unlist(dis[[i]][,paste("F",lens,sep="")]),inches=0.1,xlab="",ylab="",main=gears[i],xlim=c(min(years)-0.5,max(years+0.5)))
}
mtext("Year",outer=T,side=1,line=0)
mtext("Length (cm)",outer=T,side=2,line=0.5,las=0)
if(doPNG) dev.off()
