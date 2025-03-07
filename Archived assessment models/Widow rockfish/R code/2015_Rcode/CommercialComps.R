setwd("C:\\NOAA2015\\Widow\\Data")
source("Rcode\\Functions\\BDS_filterData.R")
source("Rcode\\Functions\\BDS_comps.R")

effN <- function(nSample,nFish,source=c("fishery","survey")) {
	if(length(nSample)!=length(nFish)) {stop("Vectors must be same length")}
	if(source=="fishery") {
		nFishCoeff <- 0.138
		nSampleCoeff <- 7.06
		breakPt <- 44
	}
	if(source=="survey") {
		nFishCoeff <- 0.0707
		nSampleCoeff <- 4.89
		breakPt <- 55
	}
	nEff <- rep(NA,length(nSample))
	ind <- nFish/nSample < breakPt
	ind[is.na(ind)] <- TRUE
	nEff[ind] <- nSample[ind] + nFishCoeff*nFish[ind]
	ind <- !ind
	nEff[ind] <- nSampleCoeff * nSample[ind]
	names(nEff) <- names(nFish)
	return(nEff)
}


load("PacFIN/PacFIN.WDOW.bds.05.Jun.2015.dmp")
BDS <- PacFIN.WDOW.bds.05.Jun.2015
rm(PacFIN.WDOW.bds.05.Jun.2015)
load("PacFIN/PacFIN.WDOW.bds.Rdat")
BDS <- bds.fish.worked
rm(bds.fish.worked)

load("extractedData/PacfinVdrfdCatch.Rdat")
load("PacFIN/PacFinCatchWithShoresideHake.Rdat") #loads pcatch (created in file pacfin.R)
hakeFTID <- as.character(unique(pcatch$FTID[pcatch$HIX_HAKE]))

load("ExtractedData/CALCOMsamps_wdow.Rdat") #loads pcatch (created in file pacfin.R)

#figure out what CALCOM data is already in BDS and what needs to be added
CA <- BDS[BDS$SAMPLE_AGENCY=="CA",]

if(F) {
	calcomSamps[calcomSamps$pink_ticket=="K131959" & !is.na(calcomSamps$pink_ticket),]
	BDS[BDS$FTID=="K131959",]

	apply(calcomSamps,2,function(x){sum(is.na(x))})
	apply(CA,2,function(x){sum(is.na(x))})

	table(substring(calcomSamps$sample_no,1,4),is.na(calcomSamps$pink_ticket))

	x <- table(factor(BDS$SAMPLE_YEAR,levels=1971:2014),BDS$SAMPLE_AGENCY)[,"CA"]
	y <- table(factor(substring(calcomSamps$sample_no,1,4),levels=1971:2014))
	data.frame(BDS=x,CALCOM=y,BDS_CALCOM=x-y)
	#there are always more in CALCOM for some years, so add them to BDS

	x <- table(BDS$FTID,BDS$SAMPLE_AGENCY)
	y <- table(calcomSamps$pink_ticket)

	BDS[BDS$SAMPLE_NO=="1982606161",]
	calcomSamps[substring(calcomSamps$sample_no,1,4)==1982 & !is.na(calcomSamps$pink_ticket) & calcomSamps$pink_ticket=="N964803",]
	table(nchar(BDS$SAMPLE_NO[BDS$SAMPLE_AGENCY=="CA"]))
	table(nchar(calcomSamps$sample_no))

	calcomSamps[calcomSamps$sample_no==197800313,]
	CA[CA$SAMPLE_NO==19785501313,]

	table(CA$SAMPLE_YEAR,CA$SPECIES_WGT!=CA$CLUSTER_WGT)
}
calcomSamps$year <- as.numeric(substring(calcomSamps$sample_no,1,4))
calcomSamps$month <- match(months(calcomSamps$sample_date),month.name)
calcomSamps$day <- unlist(lapply(strsplit(as.character(calcomSamps$sample_date),"-"),function(x){as.numeric(x[3])}))


bdsString <- CA$String <- paste0(CA$SAMPLE_YEAR,
					CA$SAMPLE_MONTH,
					CA$SAMPLE_DAY,
	                CA$PORT,
	                CA$GRID,
	                CA$CLUSTER_NO,
	                CA$FISH_NO,
	                CA$FISH_LENGTH)

calcomString <- calcomSamps$String <- paste0(calcomSamps$year,
	                   calcomSamps$month,
	                   calcomSamps$day,
	                   calcomSamps$cal_port,
	                   calcomSamps$gear,
	                   calcomSamps$clust_no,
	                   calcomSamps$fish_no,
	                   calcomSamps$flength)
length(calcomString) - length(bdsString)
inds <- match(calcomString,bdsString)
sum(is.na(inds))
tmp <- calcomSamps[is.na(inds),-c(3,5,9,10,11,12,13,14,15,16,17,19,20,21,22,23,24,25,30,31,32,34,38,40)]
tmp <- data.frame(SPID=tmp$species,SAMPLE_NO=tmp$sample_no,SAMPLE_YEAR=tmp$year,SOURCE_AGID="C",
	              SAMPLE_AGENCY="CA",CLUSTER_NO=tmp$clust_no,FISH_AGE_YEARS_FINAL=tmp$age,
	              FISH_AGE_CODE_FINAL=NA,FISH_NO=tmp$fish_no,FREQ=1,FISH_LENGTH=tmp$flength,
	              FISH_LENGTH_TYPE=NA,FORK_LENGTH_ESTIMATED=NA,FORK_LENGTH=NA,MATURITY=tmp$maturity,
	              MATURITY_AGCODE=NA,FISH_WEIGHT=NA,SEX=tmp$sex,DATA_TYPE="C",DEPTH_AVG=NA,DEPTH_MIN=NA,DEPTH_MAX=NA,
	              DRVID=tmp$boat_no,GEAR=tmp$gear_grp,GRID=tmp$gear,INPFC_AREA="CA",PSMFC_AREA=NA,PSMFC_ARID=NA,
	              SAMPLE_AGID="C",SAMPLE_MONTH=tmp$month,SAMPLE_DAY=tmp$day,SAMPLE_METHOD="R",SAMPLE_TYPE=NA,
	              MALES_WGT=NA,MALES_NUM=NA,FEMALES_NUM=NA,FEMALES_WGT=NA,TOTAL_WGT=tmp$total_wgt,
	              EXP_WT=NA,PCID=NA,PORT=tmp$cal_port,FTID=tmp$pink_ticket,COND=NA,COND_AGCODE=NA,
	              AGE_STRUCT_AGCODE=NA,AGE_METHOD=NA,AGE_READABILITY=NA,AGED_BY=NA,DATE_AGED=NA,age1=tmp$age,
	              age2=NA,age3=NA,all_cluster_sum=NA,SPECIES_WGT=tmp$weight,CLUSTER_WGT=NA,
	              FRAME_CLWT=NA,ADJ_CLWT=NA,String=tmp$String,CALCOM=TRUE)
tmp$SEX <- factor(tmp$SEX, levels=c(1,2,9), labels=c("M","F","U"))
CA$CALCOM <- FALSE
#there are some samples that appear to be in there twice, but with different sample numbers
if(F) {
	calcomSamps[calcomSamps$sample_no=="197707721"  & !is.na(calcomSamps$pink_ticket),]
	calcomSamps[calcomSamps$sample_no=="197777427"  & !is.na(calcomSamps$pink_ticket),]
	CA[CA$FTID=="141060" & !is.na(CA$FTID),]
}
tmp.split <- split(tmp,tmp$String)
x <- unlist(lapply(tmp.split,nrow))
#table(x) #29 doubles
y <- names(tmp.split)[x>1]
#tmp.split[[y[1]]]
duplSampNo <- matrix(NA,nrow=length(y),ncol=max(x))
for(i in 1:length(y)) {
	xx <- unique(tmp.split[[y[i]]]$SAMPLE_NO)
	duplSampNo[i,1:length(xx)] <- xx
}
duplSampNo <- duplSampNo[!duplicated(duplSampNo[,1]),]
deleteSampNos <- duplSampNo[,2]   #these are sample numbers to delete because they are in there twice (all 1977)
cat("Deleteing the following sample numbers because they appear to be duplicated with different sample numbers:\n",deleteSampNos,"\n")
tmp <- tmp[-which(tmp$SAMPLE_NO%in%deleteSampNos),]  #but this selects 30 rows

combined <- rbind(CA,tmp)
combined <- combined[,-which(names(combined)=="String")]

if(F) {
	x <- table(tmp$SAMPLE_YEAR,tmp$SAMPLE_NO,tmp$GEAR)>0
	apply(x,c(1,3),sum)  #number of trips sampled

	table(tmp$SAMPLE_YEAR,tmp$GEAR)


}

#add in the allClustSpWt
# combinedString <- paste0(combined$SAMPLE_YEAR,
# 					combined$SAMPLE_MONTH,
# 					combined$SAMPLE_DAY,
# 	                combined$PORT,
# 	                combined$GRID,
# 	                combined$DRVID)

# tmp.split <- split(combined,combinedString)
# if(F) {
# 	table(unlist(lapply(tmp.split,nrow))) #num of samples per sample
# 	x <- unlist(lapply(tmp.split,function(x){length(unique(x$SAMPLE_NO))}))
# 	table(x)
# 	y <- names(tmp.split)[x>1]
# 	tmp.split[[y[1]]]
# 	tmp.split[[y[10]]] #some of these are the same fish ticket, but multiple sample numbers (RWT is a sum of totalwts)
# 	length(tmp.split)
# }

#just do it by sample_no
tmp.split <- split(combined,combined$SAMPLE_NO)
if(F) {
	tmp.split[[250]]; cat("\n"); pcatch[pcatch$FTID=="T050327",]  #I'm not sure why totalwt is less than rwt (rwt matches fish ticket)
	tmp.split[[2500]]; cat("\n"); pcatch[pcatch$FTID=="U705292",]
	tmp.split[[2507]]; cat("\n"); pcatch[pcatch$FTID=="T177953",]
	#looks like use species wt and RWT_LBS (from fishtickets) when possible, and fill in otherwise. This will be some issue pre-1981 because RWT is not available
}
x.fn <- function(x) {
	sum(x$SPECIES_WGT[!duplicated(x$CLUSTER_NO)])
}
allSpSum <- unlist(lapply(tmp.split,x.fn))
combined$allSpSum <- allSpSum[match(combined$SAMPLE_NO,names(allSpSum))]

if(F) {
	#check that it is correct
	tmp.split <- split(combined,combined$SAMPLE_NO)
	tmp.split[[250]]
	tmp.split[[2500]]
	tmp.split[[2508]]
}

if(F) {
	table(BDS$FISH_LENGTH)
	plot(BDS$FISH_AGE_YEARS_FINAL,BDS$FISH_LENGTH,pch=20)
	BDS[BDS$FISH_LENGTH>600 & !is.na(BDS$FISH_LENGTH),]
	#FIVE FISH ABOVE 600MM LOOK SUSPICIOUS

	#100 fish with total length from 1991, MDT, no ages
	BDS[BDS$FISH_LENGTH_TYPE & !is.na(BDS$FISH_LENGTH_TYPE),1:20]
	sort(BDS[BDS$FISH_LENGTH_TYPE & !is.na(BDS$FISH_LENGTH_TYPE),"FISH_LENGTH"])
	BDS[BDS$FISH_LENGTH_TYPE & !is.na(BDS$FISH_LENGTH_TYPE),"GRID"]

	table(BDS$SAMPLE_AGENCY,useNA="ifany")
	head(BDS[BDS$SAMPLE_AGENCY=="PW",])   #PW is Pacific Whiting
	table(BDS$SAMPLE_AGENCY,BDS$SOURCE_AGID,useNA="ifany") #all PW is in OR
}

allBDS <- BDS[BDS$SAMPLE_AGENCY!="CA",]
allBDS$CALCOM <- FALSE
allBDS$allSpSum <- NA
allBDS <- rbind(allBDS,combined)
if(F) {
	write.csv(allBDS,file="Biological/allBDS.csv",row.names=F)
}


wdow <- SetUpWidowBDS.fn(allBDS, verbose=T, max.mmLength=600,
	         dataTypes=c("C"), sampleMethods=c("R"), sampleTypes=c(NA,"C","M"),
	         states=c("CA","OR","WA","PW"))

#add in RWT from fish tickets
rwtLbs <- tapply(pcatch$RMT,pcatch$FTID,sum)/0.00045359237
wdow$RWT_LBS <- rwtLbs[match(wdow$FTID,names(rwtLbs))]
#find which fish tickets do not have an entry
x <- table(wdow$SAMPLE_YEAR[is.na(wdow$RWT_LBS)],wdow$FTID[is.na(wdow$RWT_LBS)],wdow$state[is.na(wdow$RWT_LBS)])
apply(x,c(1,3),sum)
unique(wdow$FTID[is.na(wdow$RWT_LBS)])

table(wdow$FTID[is.na(wdow$RWT_LBS)],wdow$CALCOM[is.na(wdow$RWT_LBS)],useNA='ifany')
x <- table(wdow$FTID[is.na(wdow$RWT_LBS)],wdow$SAMPLE_YEAR[is.na(wdow$RWT_LBS)],useNA='ifany')
x
apply(x>0,c(2),sum)  #Number of missing fish tickets

if(F) {
	tmp <- names(rwtLbs)[match(wdow$FTID,names(rwtLbs))]
	table(is.na(wdow$allSpSum) | is.na(wdow$RWT_LBS),wdow$SAMPLE_AGENCY)
	table(is.na(wdow$allSpSum) | is.na(wdow$RWT_LBS),is.na(wdow$totalWt),wdow$SAMPLE_AGENCY,useNA='ifany')
}

if(F) {
	table(wdow[wdow$FTID %in% hakeFTID,"gear"])
	table(wdow$gear)
}
wdow[wdow$FTID %in% hakeFTID,"gear"] <- "ShoresideHake"
#table(wdow$gear,wdow$state)  #a few misclassifications
wdow[wdow$SAMPLE_AGENCY=="PW","gear"] <- "ShoresideHake"
#table(wdow$gear,wdow$state)  #a few misclassifications

if(F) {
	table(wdow$gear)
	table(wdow$GRID[wdow$gear=="ShoresideHake"])
	table(wdow$FISH_LENGTH)

	table(is.na(wdow$totalWt))
	table(is.na(wdow$totalWt),wdow$state)
	table(wdow$SAMPLE_YEAR,is.na(wdow$totalWt))
}

tmp <- unique(wdow$FTID[is.na(wdow$totalWt)])  #these are observations where totalWt was not available
tmp <- pcatch[pcatch$FTID%in%tmp,]    #these are fish tickets that could be used to fill those in
tmp <- tapply(tmp$LBS,factor(tmp$FTID),sum)
for(i in 1:length(tmp)) {
	wdow$totalWt[wdow$FTID==names(tmp[i]) & !is.na(wdow$FTID)] <- tmp[i]
}

if(F) {
	table(is.na(wdow$totalWt))
	table(is.na(wdow$totalWt),wdow$state)
	table(wdow$SAMPLE_YEAR,is.na(wdow$totalWt))
	wdow[wdow$SAMPLE_YEAR==2013 & is.na(wdow$totalWt),]

	x <- table(wdow$state,wdow$FTID)["WA",]
	x[x>0]
	pcatch[pcatch$FTID=="3927571",]
	head(wdow[wdow$FTID=="3927571",])
	pcatch[pcatch$FTID=="EA002016",]
	head(wdow[wdow$FTID=="EA002016",],2)

}

x <- table(wdow$SAMPLE_YEAR,wdow$SAMPLE_NO,wdow$gear)>0
apply(x,c(1,3),sum)  #number of trips sampled

x <- table(wdow$SAMPLE_YEAR,wdow$SAMPLE_NO,wdow$state,wdow$gear)>0
apply(x,c(1,3,4),sum)

table(wdow$SAMPLE_YEAR,wdow$state,wdow$gear)


if(F) {
	#sex ratios
	sexRatio.fn <- function(dat,...) {
	    tmp.split <- split(dat$SEX,dat$length.cm)
	    x.fn <- function(x) {
	    	xx <- table(x)
	    	out <- xx["F"]/(xx["F"]+xx["M"])
	    	names(out) <- NULL
	    	return(out)
	    }
	    propF <- unlist(lapply(tmp.split,x.fn))
	    nobs <- unlist(lapply(tmp.split,function(x){sum(x=="M" | x=="F")}))
	    x <- as.numeric(names(propF))
	    plot(x,propF,type="l",col="red",xlab="Length (cm)",ylab="Fraction female",ylim=c(0,1),...)
	    symbols(x,propF,circles=nobs,inches=0.1,fg="red",bg=rgb(1,0,0,alpha=0.5),add=T)
	   	return(list(propF=propF,nobs=nobs))
	}
	windows(height=7,width=10)
	par(mfrow=c(2,3))
	sexRatio.fn(wdow,main="All Years")$propF
	sexRatio.fn(wdow[substring(wdow$SAMPLE_YEAR,1,3)=="197",],main="1970's")
	sexRatio.fn(wdow[substring(wdow$SAMPLE_YEAR,1,3)=="198",],main="1980's")
	sexRatio.fn(wdow[substring(wdow$SAMPLE_YEAR,1,3)=="199",],main="1990's")
	sexRatio.fn(wdow[substring(wdow$SAMPLE_YEAR,1,3)=="200",],main="2000's")
	sexRatio.fn(wdow[substring(wdow$SAMPLE_YEAR,1,3)=="201",],main="2010's")
	#there is no obvious trend, so use all years

	par(mfrow=c(1,1))
	out <- sexRatio.fn(wdow,main="All Years")
	propF <- out$propF
	nobs <- out$nobs
	sexRatio <- rep(NA,diff(range(wdow$length.cm))+1)
    names(sexRatio) <- min(wdow$length.cm):max(wdow$length.cm)
    sexRatio[names(propF)] <- propF
    lens <- as.numeric(names(propF))
    propF2 <- propF
    propF2[lens<=28] <- 0.5
    propF2[lens>55] <- 1
    nobs2 <- nobs
    nobs2[lens<=28] <- max(nobs2)
    lines(lens,propF2,col="green")
    lo <- loess(propF2~lens,weights=nobs2,span=0.3)
    lo.pred <- predict(lo,newdata=data.frame(lens=min(wdow$length.cm):max(wdow$length.cm)))
    names(lo.pred) <- names(sexRatio)
    lo.pred[lo.pred>1] <- 1
    lines(as.numeric(names(lo.pred)),lo.pred,col="blue")
    wdow$sexRatio <- lo.pred[match(wdow$length.cm,as.numeric(names(lo.pred)))]
}

write.csv(wdow,file="Biological/wdow.csv",row.names=F)

##################################################################################################
### LENGTHS
#for table
x <- table(wdow$SAMPLE_YEAR,wdow$SAMPLE_NO,wdow$SAMPLE_AGENCY,wdow$gear)
nSamp <- apply(x,c(1,3,4),function(x){sum(x>0)})
nFish <- table(wdow$SAMPLE_YEAR,wdow$SAMPLE_AGENCY,wdow$gear)

x <- table(wdow$SAMPLE_YEAR,wdow$SAMPLE_NO,wdow$gear)
nSamp <- apply(x,c(1,3),function(x){sum(x>0)})
nFish <- table(wdow$SAMPLE_YEAR,wdow$gear)

lw <- data.frame(WA=c(1.7043e-5,2.9668),OR=c(1.7043e-5,2.9668),CA=c(1.7043e-5,2.9668))

#Bottom Trawl
LFs <- commLFs.fn(wdow,lw,gear="BottomTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowBottomTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=1,gender=3,partition=2)
n <- floor(effN(nSamp[,"BottomTrawl"],nFish[,"BottomTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
apply(tmp[,-(1:6)],1,sum)
btComps <- tmp[order(tmp[,"year"]),]
write.csv(btComps,file="Comps/bottomTrawlLengthComps.csv",row.names=F)

#Midwater Trawl
LFs <- commLFs.fn(wdow,lw,gear="MidwaterTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowMidwaterTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=2,gender=3,partition=2)
n <- floor(effN(nSamp[,"MidwaterTrawl"],nFish[,"MidwaterTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
mwComps <- tmp[order(tmp[,"year"]),]
write.csv(mwComps,file="Comps/midwaterTrawlLengthComps.csv",row.names=F)

#Hook and Line
LFs <- commLFs.fn(wdow,lw,gear="HnL",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowHnL.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=5,gender=3,partition=2)
n <- floor(effN(nSamp[,"HnL"],nFish[,"HnL"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
hnlComps <- tmp[order(tmp[,"year"]),]
write.csv(hnlComps,file="Comps/HnLLengthComps.csv",row.names=F)

#Net
graphics.off()
LFs <- commLFs.fn(wdow,lw,gear="Net",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowNet.csv",
 	            maxExpansion=100,verbose=T,loessSpan=0.9)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=4,gender=3,partition=2)
n <- floor(effN(nSamp[,"Net"],nFish[,"Net"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
netComps <- tmp[order(tmp[,"year"]),]
write.csv(netComps,file="Comps/netLengthComps.csv",row.names=F)

#Shoreside Hake
wdow$state[wdow$gear=="ShoresideHake"] <- "PW"
LFs <- commLFs.fn(wdow,lw,gear="ShoresideHake",state="PW",
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowShoresideHake.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=3,gender=3,partition=2)
n <- floor(effN(nSamp[,"ShoresideHake"],nFish[,"ShoresideHake"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
sshComps <- tmp[order(tmp[,"year"]),]
write.csv(sshComps,file="Comps/shoresideHakeLengthComps.csv",row.names=F)


if(F) {
	table(is.na(LFs$sampleWgt))
	table(is.na(LFs$usetot_wgt))
	table(is.na(LFs$sampleWgt),LFs$state)
	table(is.na(LFs$usetot_wgt),LFs$state)
}

if(F) {
	#look at 1983 net
	net <- wdow[wdow$gear=="Net",]
	net[net$SAMPLE_YEAR==1983,]   #it is double counting the CALCOM samples (there are two sets of same samples)
	net[net$SAMPLE_YEAR==1984,]   #again, doubel samples
	bt <- wdow[wdow$gear=="BottomTrawl",]
    bt[bt$SAMPLE_YEAR==1981 & bt$totalWt==39270,]   #it is double counting the CALCOM samples (there are two sets of same samples)

}

#No CALCOM (may be double counting)
wdow2 <- wdow[!(wdow$CALCOM),]
x <- table(wdow2$SAMPLE_YEAR,wdow2$SAMPLE_NO,wdow2$gear)
nSamp <- apply(x,c(1,3),function(x){sum(x>0)})
nFish <- table(wdow2$SAMPLE_YEAR,wdow2$gear)

#Bottom Trawl
LFs <- commLFs.fn(wdow2,lw,gear="BottomTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowBottomTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow2$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=1,gender=3,partition=2)
n <- floor(effN(nSamp[,"BottomTrawl"],nFish[,"BottomTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
apply(tmp[,-(1:6)],1,sum)
btComps <- tmp[order(tmp[,"year"]),]
write.csv(btComps,file="Comps/bottomTrawlNoCalcomLengthComps.csv",row.names=F)

#Midwater Trawl
LFs <- commLFs.fn(wdow2,lw,gear="MidwaterTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowMidwaterTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow2$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=2,gender=3,partition=2)
n <- floor(effN(nSamp[,"MidwaterTrawl"],nFish[,"MidwaterTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
mwComps <- tmp[order(tmp[,"year"]),]
write.csv(mwComps,file="Comps/midwaterTrawlNoCalcomLengthComps.csv",row.names=F)

#Hook and Line
LFs <- commLFs.fn(wdow2,lw,gear="HnL",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowHnL.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow2$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=5,gender=3,partition=2)
n <- floor(effN(nSamp[,"HnL"],nFish[,"HnL"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
hnlComps <- tmp[order(tmp[,"year"]),]
write.csv(hnlComps,file="Comps/HnLNoCalcomLengthComps.csv",row.names=F)

#Net
LFs <- commLFs.fn(wdow2,lw,gear="Net",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowNet.csv",
 	            maxExpansion=100,verbose=T,loessSpan=0.9)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow2$SAMPLE_YEAR),lens=seq(8,56,2),season=1,fleet=4,gender=3,partition=2)
n <- floor(effN(nSamp[,"Net"],nFish[,"Net"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
netComps <- tmp[order(tmp[,"year"]),]
write.csv(netComps,file="Comps/netNoCalcomLengthComps.csv",row.names=F)



###############################################################################
# AGES
#Run above to get allBDS and wdow (just before lengths)
#trick it to think that ages are lengths
allAges <- wdow[!is.na(wdow$FISH_AGE_YEARS_FINAL),]
allAges$FISH_AGE_YEARS_FINAL <- as.numeric(as.character(allAges$FISH_AGE_YEARS_FINAL))
x <- table(allAges$SAMPLE_YEAR,allAges$SAMPLE_NO,allAges$gear)
nSamp <- apply(x,c(1,3),function(x){sum(x>0)})
nFish <- table(allAges$SAMPLE_YEAR,allAges$gear)

#For tables
x <- table(allAges$SAMPLE_YEAR,allAges$SAMPLE_NO,allAges$SAMPLE_AGENCY,allAges$gear)
nSamp <- apply(x,c(1,3,4),function(x){sum(x>0)})
nFish <- table(allAges$SAMPLE_YEAR,allAges$SAMPLE_AGENCY,allAges$gear)



lw <- data.frame(WA=c(1.7043e-5,2.9668),OR=c(1.7043e-5,2.9668),CA=c(1.7043e-5,2.9668))

#Bottom Trawl
LFs <- commLFs.fn(allAges,lw,gear="BottomTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowBottomTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=0:40,season=1,fleet=1,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"BottomTrawl"],nFish[,"BottomTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/bottomTrawlAgeComps.csv",row.names=F)

#Midwater Trawl
LFs <- commLFs.fn(allAges,lw,gear="MidwaterTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowMidwaterTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=0:40,season=1,fleet=2,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"MidwaterTrawl"],nFish[,"MidwaterTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/midwaterTrawlAgeComps.csv",row.names=F)

#Hook and Line
LFs <- commLFs.fn(allAges,lw,gear="HnL",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowHnL.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.8,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=0:40,season=1,fleet=5,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"HnL"],nFish[,"HnL"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/hookAndLineAgeComps.csv",row.names=F)

#Net
LFs <- commLFs.fn(allAges,lw,gear="Net",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowNet.csv",
 	            maxExpansion=100,verbose=T,loessSpan=0.8,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=0:40,season=1,fleet=4,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"Net"],nFish[,"Net"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/netAgeComps.csv",row.names=F)

#Shoreside Hake
allAges$state[allAges$gear=="ShoresideHake"] <- "PW"
LFs <- commLFs.fn(allAges,lw,gear="ShoresideHake",state="PW",
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowShoresideHake.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(wdow$SAMPLE_YEAR),lens=0:40,season=1,fleet=3,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"ShoresideHake"],nFish[,"ShoresideHake"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/shoresideHakeAgeComps.csv",row.names=F)


#################################################################
### No CALCOM
allAges2 <- wdow[!is.na(wdow$FISH_AGE_YEARS_FINAL) & !(wdow$CALCOM),]
allAges2$FISH_AGE_YEARS_FINAL <- as.numeric(as.character(allAges2$FISH_AGE_YEARS_FINAL))
x <- table(allAges2$SAMPLE_YEAR,allAges2$SAMPLE_NO,allAges2$gear)
nSamp <- apply(x,c(1,3),function(x){sum(x>0)})
nFish <- table(allAges2$SAMPLE_YEAR,allAges2$gear)

#Bottom Trawl
LFs <- commLFs.fn(allAges2,lw,gear="BottomTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowBottomTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(allAges2$SAMPLE_YEAR),lens=0:40,season=1,fleet=1,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"BottomTrawl"],nFish[,"BottomTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/bottomTrawlNoCalcomAgeComps.csv",row.names=F)

#Midwater Trawl
LFs <- commLFs.fn(allAges2,lw,gear="MidwaterTrawl",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowMidwaterTrawl.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.3,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(allAges2$SAMPLE_YEAR),lens=0:40,season=1,fleet=2,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"MidwaterTrawl"],nFish[,"MidwaterTrawl"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/midwaterTrawlNoCalcomAgeComps.csv",row.names=F)

#Hook and Line
LFs <- commLFs.fn(allAges2,lw,gear="HnL",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowHnL.csv",
 	            maxExpansion=300,verbose=T,loessSpan=0.8,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(allAges2$SAMPLE_YEAR),lens=0:40,season=1,fleet=5,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"HnL"],nFish[,"HnL"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/hookAndLineNoCalcomAgeComps.csv",row.names=F)

#Net
LFs <- commLFs.fn(allAges2,lw,gear="Net",state=NULL,
 	            catchFile="C:/NOAA2015/Widow/Data/Catches/ForComps/WidowNet.csv",
 	            maxExpansion=100,verbose=T,loessSpan=0.8,ageComp=TRUE)
tmp <- lfsForSS3_gender3.fn(LFs$both,years=unique(allAges2$SAMPLE_YEAR),lens=0:40,season=1,fleet=4,gender=3,partition=2,ageComp=TRUE)
tmp <- tmp[apply(tmp[,-(1:6)],1,sum)>0,]
n <- floor(effN(nSamp[,"Net"],nFish[,"Net"],"fishery"))
tmp[,"nSamps"] <- n[as.character(tmp[,"year"])]
tmp <- tmp[tmp[,"nSamps"]>0,]
comps <- tmp[order(tmp[,"year"]),]
comps <- data.frame(comps[,1:5],AgeError="ENTER",Lbin_lo=-1,Lbin_hi=-1,
	                   comps[,6:ncol(comps)])
write.csv(comps,file="Comps/netNoCalcomAgeComps.csv",row.names=F)



########################################################
#Age-at-length
#Run above to get allBDS and wdow (just before lengths)

# devtools::install_github("nwfsc-assess/nwfscSurvey")
library(nwfscSurvey)

allAges <- wdow[!is.na(wdow$FISH_AGE_YEARS_FINAL),]
allAges$FISH_AGE_YEARS_FINAL <- as.numeric(as.character(allAges$FISH_AGE_YEARS_FINAL))

table(allAges$SEX)
table(allAges$FISH_AGE_YEARS_FINAL,allAges$SEX)
#the unsexed arefew and not in places where it would matter, so remove
allAges$SEX <- as.character(allAges$SEX)
allAges <- allAges[allAges$SEX!="U",]
table(allAges$SEX)

table(allAges$gear)

gr <- "BottomTrawl"
fl <- 1
dat3 <- allAges[allAges$gear==gr,]
x <- as.data.frame(table(dat3$SAMPLE_YEAR,dat3$SEX,dat3$length.cm,dat3$FISH_AGE_YEARS_FINAL))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]
AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet=fl,season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,paste("Comps/AgeAtLenForSS3",gr,"Female.csv",sep="."),row.names=F)
write.csv(AatL$male,paste("Comps/AgeAtLenForSS3",gr,"Male.csv",sep="."),row.names=F)


gr <- "MidwaterTrawl"
fl <- 2
dat3 <- allAges[allAges$gear==gr,]
x <- as.data.frame(table(dat3$SAMPLE_YEAR,dat3$SEX,dat3$length.cm,dat3$FISH_AGE_YEARS_FINAL))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]
AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet=fl,season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,paste("Comps/AgeAtLenForSS3",gr,"Female.csv",sep="."),row.names=F)
write.csv(AatL$male,paste("Comps/AgeAtLenForSS3",gr,"Male.csv",sep="."),row.names=F)

gr <- "HnL"
fl <- 5
dat3 <- allAges[allAges$gear==gr,]
x <- as.data.frame(table(dat3$SAMPLE_YEAR,dat3$SEX,dat3$length.cm,dat3$FISH_AGE_YEARS_FINAL))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]
AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet=fl,season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,paste("Comps/AgeAtLenForSS3",gr,"Female.csv",sep="."),row.names=F)
write.csv(AatL$male,paste("Comps/AgeAtLenForSS3",gr,"Male.csv",sep="."),row.names=F)

gr <- "Net"
fl <- 4
dat3 <- allAges[allAges$gear==gr,]
x <- as.data.frame(table(dat3$SAMPLE_YEAR,dat3$SEX,dat3$length.cm,dat3$FISH_AGE_YEARS_FINAL))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]
AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet=fl,season=1,partition=0,ageerr="EnterAgeErr",raw=T)
write.csv(AatL$female,paste("Comps/AgeAtLenForSS3",gr,"Female.csv",sep="."),row.names=F)
write.csv(AatL$male,paste("Comps/AgeAtLenForSS3",gr,"Male.csv",sep="."),row.names=F)


#I'll bring in at-sea observations and do them together
gr <- "ShoresideHake"
fl <- 3
datS <- allAges[allAges$gear==gr,]
datS$state[datS$gear=="ShoresideHake"] <- "PW"
datS <- datS[,c("SAMPLE_YEAR","SEX","length.cm","FISH_AGE_YEARS_FINAL")]
load("extractedData/atsea.bio.Rdat")
atsea.bio$Year <- as.numeric(substr(atsea.bio$RETRV_DATE_TIME,1,4))
dat <- atsea.bio[is.na(atsea.bio$BARCODE),]
dat1 <- dat[!is.na(dat$AGE) & dat$Year<2011,]
dat1 <- dat1[!duplicated(dat1$SPECIMEN_NUMBER),]
dat2 <- atsea.bio[!is.na(atsea.bio$AGE) & atsea.bio$Year>=2011,]
dat3 <- rbind(dat1,dat2)
dat3$SEX_CODE <- as.character(dat3$SEX_CODE)
dat3 <- dat3[dat3$SEX!="U",]
datA <- data.frame(SAMPLE_YEAR=dat3$Year,SEX=dat3$SEX_CODE,length.cm=dat3$LENGTH_SIZE,FISH_AGE_YEARS_FINAL=dat3$AGE)

dat3 <- rbind(datS,datA)

x <- as.data.frame(table(dat3$SAMPLE_YEAR,dat3$SEX,dat3$length.cm,dat3$FISH_AGE_YEARS_FINAL))
x <- cbind(x[x$Var2=="F",],x[x$Var2=="M",])
names(x) <- c("Year","Sex","Length","Age","AgeTallyF","YearM","SexM","LengthM","AgeM","AgeTallyM")
sum(as.numeric(as.character(x$Length)) - as.numeric(as.character(x$LengthM)))  #check that they match
x <- x[,c("Year","Length","Age","AgeTallyF","AgeTallyM")]
x <- as.data.frame(apply(x,2,function(x){as.numeric(as.character(x))}))
x <- x[x$AgeTallyF>0 | x$AgeTallyM>0,]
x <- x[order(x$Year,x$Length,x$Age),]
AatL <- SS3AgeAtLen.fn(x,lgthBins=seq(8,56,2),ageBins=0:40,fleet=fl,season=1,partition=0,ageerr="EnterAgeErr",raw=T)
gr <- "Hake"
write.csv(AatL$female,paste("Comps/AgeAtLenForSS3",gr,"Female.csv",sep="."),row.names=F)
write.csv(AatL$male,paste("Comps/AgeAtLenForSS3",gr,"Male.csv",sep="."),row.names=F)





