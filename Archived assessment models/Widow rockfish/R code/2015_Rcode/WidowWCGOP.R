dont source
setwd("C:/NOAA2015/ObserverDatabase")
source("../ExtractR/Rcode/functions/Functions.R")
files <- dir("Rcode/functions")
for(i in 1:length(files)) {
	source(file.path("Rcode/functions", files[i]))
}


# gearLevels <- list(c("Bottom Trawl","Midwater Trawl","Shrimp Trawl"),"Hook & Line","Pot")
# gearNames <- c("Trawl","HookAndLine","Pot")
# stateLevels <- c("OR")  #r_state
# stateNames <- c("OR")
# latLevels <- c(0,43.9,100)  #r_state
# latNames <- c("OrS","OrN")


# #######################################################################################################################
# ### check for confidentiality using all observations
# out <- NULL
# load("extractedData/ob.pre2011.Rdat") #object named ob
# out$'pre2011' <- checkConfidentiality(ob,colnms=c("gear2","r_state"),
# 				  colLevs=list(gearLevels,stateLevels),
# 				  strNms=list(gearNames,stateNames))

# load("extractedData/ob.2011.Rdat") #object named ob
# out$'2011' <- checkConfidentiality(ob,colnms=c("gear2","r_state"),
# 				  colLevs=list(gearLevels,stateLevels),
# 				  strNms=list(gearNames,stateNames))

# load("extractedData/ob.2012.Rdat") #object named ob
# out$'2012' <- checkConfidentiality(ob,colnms=c("gear2","r_state"),
# 				  colLevs=list(gearLevels,stateLevels),
# 				  strNms=list(gearNames,stateNames))

# load("extractedData/ob.2013.Rdat") #object named ob
# out$'2013' <- checkConfidentiality(ob,colnms=c("gear2","r_state"),
# 				  colLevs=list(gearLevels,stateLevels),
# 				  strNms=list(gearNames,stateNames))

# save(out,file="KelpGreenling/RdatFiles/allForConfidentiality_AllOR.Rdat")




###########################################################################
#Point estimates
setwd("C:/NOAA2015/ObserverDatabase")
source("Rcode/functions/observRfunctions.R")
source("../ExtractR/Rcode/functions/Functions.R")
source("Rcode/functions/checkConfidentiality.R")
source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\functions\\Boot_summary_fxn_jej.r")
source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\functions\\bootstrapDiscardBiomass_ach_V3.r")
source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\functions\\bootstrapDiscardRatio_ach_V3.r")
source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\BootstrapDiscardData.R")

source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\discardsNonCatchShares.R")
source("C:\\NOAA2015\\ObserverDatabase\\Rcode\\discardsCatchShares.R")

dte=Sys.Date()

load("extractedData/ob.pre2011.Rdat") #object named ob
ob.all <- ob
load("extractedData/ob.2011.Rdat") #object named ob
ob.all <- rbind(ob.all,ob)
load("extractedData/ob.2012.Rdat") #object named ob
ob.all <- rbind(ob.all,ob)
load("extractedData/ob.2013.Rdat") #object named ob
ob.all <- rbind(ob.all,ob)

ob <- ob.all
rm(ob.all)


unique(ob$species[grep("Widow",ob$species)])

#strata definitions
gearLevels <- list("Bottom Trawl","Midwater Trawl","Shrimp Trawl","Hook & Line","Pot")
gearNames <- c("Bottom Trawl","Midwater Trawl","Shrimp Trawl","HookAndLine","Pot")
sectorLevels <- list("LE CA Halibut","LE Fixed Gear DTL","Limited Entry Sablefish","Limited Entry Trawl","Nearshore","OA CA Halibut","OA Fixed Gear","Pink Shrimp","Catch Shares","Shoreside Hake")
sectorNames <- c("LE CA Halibut","LE Fixed Gear DTL","Limited Entry Sablefish","Limited Entry Trawl","Nearshore","OA CA Halibut","OA Fixed Gear","Pink Shrimp","Catch Shares","Shoreside Hake")
stateLevels <- c("WA","OR","CA")  #r_state
stateNames <- c("WA","OR","CA")
#latLevels <- c(0,43.9,100)  #r_state
#latNames <- c("OrS","OrN")

outBoot <- bootstrapDiscardData(ob,sp="Widow Rockfish",B=1,
	            colnms=c("gear2","r_state","sector"),
	            colnms.new=c("gear3","State","sector2"),
				colLevs=list(gearLevels,stateLevels,sectorLevels),
				stratNms=list(gearNames,stateNames,sectorNames),
				bootFile='Widow/RdatFiles/noboot.out',
				resultsFile='Widow/RdatFiles/noboot_dat.out')

write.csv(outBoot$ncs[,c(1:13,29)],paste('Widow/Widow_OB_DisRatios_noboot_ncs_',dte,'.csv',sep=''),row.names=F)
write.csv(outBoot$cs,paste('Widow/Widow_OB_DisRatios_noboot_cs_',dte,'.csv',sep=''),row.names=F)


wdow <- ob[ob$species=="Widow Rockfish",]
table(wdow$ryear,wdow$gear2)
table(wdow$gear2,wdow$sector)
hist(wdow$dis/(wdow$dis+wdow$ret))
#widow is either all discarded or all retained NO NO NO. 
head(wdow[duplicated(wdow$haul_id),])
wdow[wdow$haul_id==81221,]
#there is a single line for retained and a single line for discarded

#find tow specific discard rates
dat <- split(wdow[,c("ryear","dis","ret")],wdow$haul_id)
dat[['78873']]
wdow[c("200921","200923","200937","201044"),]
write.csv(ob[ob$haul_id=="78873",],file="ExampleTowWithMultipleLines.csv")

dr <- unlist(lapply(dat,function(x){sum(x$dis)/sum(x$dis+x$ret)}))
hist(dr)

towSpecDr <- function(x) {
	for(yr in sort(unique(x$ryear))) {
		tmp <- x[x$ryear==yr,]
		dat <- split(tmp[,c("ryear","dis","ret")],tmp$haul_id)
		dr <- unlist(lapply(dat,function(x){sum(x$dis)/sum(x$dis+x$ret)}))
		levs <- cut(dr,c(-1,seq(0,0.9,0.1),0.999999,1))
		barplot(table(levs)/length(levs),width=c(0.1,rep(1,10),0.1),
			    names=c("0",rep(NA,10),"1"),col="blue",ylim=c(0,1),
			    main=yr,las=1)
		#hist(dr,freq=FALSE,main=yr,xlab="")
	}
	mtext("Discard Rate",outer=T,line=-0.1,side=1)
	mtext("Proportion",outer=T,line=0.2,side=2)
}

wd<-6.5;ht<-4.5
png("C:/NOAA2015/Widow/Data/Discards/Figures/HnLrates.png",width=wd,height=ht,units="in",res=300)
par(mfrow=c(3,4),mar=c(3,3,1,1),oma=c(1.5,2,2,0))
towSpecDr(wdow[wdow$gear2=="Hook & Line",])
mtext("Hook& Line",outer=T,side=3,line=0.5)
dev.off()

png("C:/NOAA2015/Widow/Data/Discards/Figures/BottomTrawlRates.png",width=wd,height=ht,units="in",res=300)
par(mfrow=c(3,4),mar=c(3,3,1,1),oma=c(1.5,2,2,0))
towSpecDr(wdow[wdow$gear2=="Bottom Trawl",])
mtext("Bottom Trawl",outer=T,side=3,line=0.5)
dev.off()

png("C:/NOAA2015/Widow/Data/Discards/Figures/MidwaterRates.png",width=wd,height=ht,units="in",res=300)
par(mfrow=c(3,4),mar=c(3,3,1,1),oma=c(1.5,2,2,0))
towSpecDr(wdow[wdow$gear2=="Midwater Trawl",])
mtext("Midwater Trawl",outer=T,side=3,line=0.5)
dev.off()



ind <- wdow$ryear==2013&wdow$sector=="Catch Shares"&wdow$gear2=="Bottom Trawl"&wdow$r_state=="CA"
hist(wdow$dis[ind])
unique(round(wdow$dis[ind],-1)
sum(wdow$dis[ind])/sum(wdow$dis[ind]+wdow$ret[ind])
cbind(wdow$dis[ind],wdow$ret[ind])
hist(wdow$dis[ind]/(wdow$dis[ind]+wdow$ret[ind]))

#there is one haul with over 2500lbs of discards that drives the high rate
#It is either all discarded or all retained
wdow[wdow$dis>2000,]


