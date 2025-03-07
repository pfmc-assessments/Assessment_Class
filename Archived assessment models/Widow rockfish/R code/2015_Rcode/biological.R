dont source
#an investigation of biological relationships
#l-w, age-length, etc
#from all dta sources (fishery, survey, port sampling, etc)

setwd("C:\\NOAA2015\\Widow\\Data")

#nwfscCombo survey
nwfscBio <- read.csv("TrawlSurvey/NWFSCsurvey/Widow_SexedLgthWtAge.csv",skip=8)
table(is.na(nwfscBio$AGE_YRS),is.na(nwfscBio$WEIGHT_KG),is.na(nwfscBio$LENGTH_CM))  #all lengths are there, 1 observation without weight, 17 without age, none without age or weight
nwfscBioSm <- data.frame(Year=as.numeric(substring(nwfscBio$PROJECT_CYCLE,7)),
						 Lat = nwfscBio$HAUL_LATITUDE_DD,
						 Lon = nwfscBio$HAUL_LONGITUDE_DD,
						 Depth = nwfscBio$DEPTH_M,
						 Sex = nwfscBio$SEX,
						 Length = nwfscBio$LENGTH_CM,
						 Weight = nwfscBio$WEIGHT_KG,
						 Age = nwfscBio$AGE_YRS,
						 Maturity = nwfscBio$MATURE_YN,
						 State=NA,
						 Sample=nwfscBio$HAUL_IDENTIFIER,
						 Source = "nwfscCombo")
plot(nwfscBio$LENGTH_CM,nwfscBio$WEIGHT_KG,pch=16)

#Triennial survey
load("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\Triennial\\AK.Surveys.Bio.widow.04.Jun.2015.dmp")
triBio <- AK.Surveys.Bio.widow.04.Jun.2015$Ages
rm(AK.Surveys.Bio.widow.04.Jun.2015)
table(triBio$SEX)
triBio$SEX[triBio$SEX==2] <- "f"
triBio$SEX[triBio$SEX==1] <- "m"
table(is.na(triBio$AGE),is.na(triBio$INDV_WGHT),is.na(triBio$LENGTH))  ##all ahve a length, 22 ages, 55 weights; no fish has all three
triBio <- triBio[!(is.na(triBio$AGE)&is.na(triBio$INDV_WGHT)),]
table(is.na(triBio$AGE),is.na(triBio$INDV_WGHT),is.na(triBio$LENGTH))  ##all ahve a length, 22 ages, 55 weights; no fish has all three
triBioSm <- data.frame(Year=triBio$YEAR,
						 Lat = triBio$START_LATITUDE,
						 Lon = triBio$START_LONGITUDE,
						 Depth = triBio$BOTTOM_DEPTH,
						 Sex = triBio$SEX,
						 Length = as.numeric(triBio$LENGTH)/10,
						 Weight = as.numeric(triBio$INDV_WGHT_G)/1000,
						 Age = as.numeric(triBio$AGE),
						 Maturity = NA,
						 State = NA,
						 Sample = triBio$HAULJOIN,
						 Source = "triennial")


#PacFIN BDS
allBDS <- read.csv("Biological/allBDS.csv")
table(is.na(allBDS$FISH_AGE_YEARS_FINAL),is.na(allBDS$FISH_WEIGHT),is.na(allBDS$FISH_LENGTH))  ##46 do not have a length, 9134 have all three
table(is.na(allBDS$FISH_LENGTH),is.na(allBDS$MATURITY))
table(allBDS$MATURITY,useNA="ifany")
allBDS <- allBDS[!(is.na(allBDS$FISH_LENGTH)&is.na(allBDS$FISH_AGE_YEARS_FINAL)&is.na(allBDS$FISH_WEIGHT)),] #removed only 32
allBDS <- allBDS[!(is.na(allBDS$FISH_AGE_YEARS_FINAL)&is.na(allBDS$FISH_WEIGHT)),] #removed 94659 observations

allBDS$SEX <- as.character(allBDS$SEX)
table(allBDS$SEX)
allBDS$SEX[allBDS$SEX=="F"] <- "f"
allBDS$SEX[allBDS$SEX=="M"] <- "m"
allBDS$SEX[allBDS$SEX=="U"] <- "u"

table(allBDS$FISH_WEIGHT)
plot(allBDS$FISH_LENGTH,allBDS$FISH_WEIGHT,pch=16)
ind <- allBDS$FISH_WEIGHT>8 & !is.na(allBDS$FISH_WEIGHT)
allBDS$FISH_WEIGHT[ind] <- allBDS$FISH_WEIGHT[ind]/1000
allBDS$FISH_WEIGHT[!ind] <- allBDS$FISH_WEIGHT[!ind] * 0.453592
plot(allBDS$FISH_LENGTH,allBDS$FISH_WEIGHT,pch=16)

allBdsSm <- data.frame(Year=allBDS$SAMPLE_YEAR,
						 Lat = NA,
						 Lon = NA,
						 Depth = NA,
						 Sex = allBDS$SEX,
						 Length = allBDS$FISH_LENGTH/10,
						 Weight = allBDS$FISH_WEIGHT,
						 Age = allBDS$FISH_AGE_YEARS_FINAL,
						 Maturity = allBDS$MATURITY,
						 State = allBDS$SOURCE_AGID,
						 Sample = allBDS$SAMPLE_NO,
						 Source = "bds")



#at-sea hake
atseaBio <- read.csv("Biological/Widow_atsea_age.csv")
table(is.na(atseaBio$AGE),is.na(atseaBio$WEIGHT),is.na(atseaBio$LENGTH))  ##all have a length, 4256 have all three
atseaBio <- atseaBio[!(is.na(atseaBio$AGE)&is.na(atseaBio$WEIGHT)),] #removed only 8
atseaBio$SEX <- as.character(atseaBio$SEX)
table(atseaBio$SEX)
atseaBio$SEX[atseaBio$SEX=="F"] <- "f"
atseaBio$SEX[atseaBio$SEX=="M"] <- "m"
atseaBio$SEX[atseaBio$SEX=="U"] <- "u"

atseaBioSm <- data.frame(Year=atseaBio$YEAR,
						 Lat = atseaBio$LATDD_START,
						 Lon = atseaBio$LONDD_START,
						 Depth = atseaBio$BOTTOM_DEPTH_FATHOMS * 1.8288,
						 Sex = atseaBio$SEX,
						 Length = atseaBio$LENGTH,
						 Weight = atseaBio$WEIGHT,
						 Age = atseaBio$AGE,
						 Maturity = NA,
						 State = NA,
						 Sample = atseaBio$HAUL_JOIN,
						 Source = "atsea")

#WCGOP
wcgopBio <- read.csv("Discards/WcgopBioData_semicolonSep.csv",sep=";")
table(is.na(wcgopBio$age),is.na(wcgopBio$specimen_weight),is.na(wcgopBio$length))  ##all have a length, no ages, 81 with weight
wcgopBio <- wcgopBio[!(is.na(wcgopBio$age)&is.na(wcgopBio$specimen_weight)),] #kept only 81
wcgopBio$sex <- as.character(wcgopBio$sex)
table(wcgopBio$sex)
wcgopBio$sex[wcgopBio$sex=="F"] <- "f"
wcgopBio$sex[wcgopBio$sex=="M"] <- "m"

wcgopBioSm <- data.frame(Year=wcgopBio$ryear,
						 Lat = wcgopBio$set_lat,
						 Lon = wcgopBio$set_long,
						 Depth = wcgopBio$set_depth * 1.8288,
						 Sex = wcgopBio$sex,
						 Length = wcgopBio$length,
						 Weight = wcgopBio$specimen_weight * 0.453592,
						 Age = wcgopBio$age,
						 Maturity = NA,
						 State = NA,
						 Sample = wcgopBio$haul_id,
						 Source = "wcgop")



bioData <- rbind(nwfscBioSm, triBioSm, allBdsSm, atseaBioSm, wcgopBioSm)
write.csv(bioData, file="Biological/allBioData.csv",row.names=F)

#######################################################################################
EcheverraiLenarz 1984: conversions from fork tot toal length in mm
a<- 6.954; b<-0.931  #TL to FL
FL <- a + b*seq(100,400,10)
a<- 6.845; b<-1.072  #FL to TL
TL <- a + b*seq(100,400,10)


setwd("C:\\NOAA2015\\Widow\\Data")
bioData <- read.csv("Biological/allBioData.csv")

table(bioData$Sex,useNA='ifany')
table(bioData$Year,useNA='ifany')


plot(bioData$Length,bioData$Weight,pch=20)
plot(jitter(bioData$Length,factor=7),bioData$Weight,pch=20)

doPNG <- T
wd<-6.5;ht<-4.5
if(doPNG) {png("../Writeup/Figures/weightAtLengthBySource.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(4,3,1,1))
sources <- unique(bioData$Source)
cols <- c("black","blue","green3","red","purple")
cols <- adjustcolor(cols,alpha.f<- 0.6)
plot(bioData$Length,bioData$Weight,pch=20,type="n",xlab="Length (cm)",ylab="Weight (kg)",las=1)
for(i in length(sources):1) {
	tmp <- bioData[bioData$Source==sources[i],]
	points(jitter(tmp$Length,factor=7),tmp$Weight,col=cols[i],pch=16,cex=0.7)

}
legend("topleft",c("NWFSC","Triennial","BDS","At-Sea","WCGOP"),col=cols,pch=20)
if(doPNG) dev.off()
#all outliers are from fisheries data
#WCGOP data are biased low
#the >60cm fish seems way out there
table(round(bioData$Length[bioData$Weight<0.2],0))
#a group of fish between 30 and 45 that are less than 0.2

ind <- bioData$Length>30 & bioData$Length<45 & bioData$Weight<0.2 & !is.na(bioData$Weight)   #bioData[ind,]
bioDat2 <- bioData[bioData$Weight>0 & bioData$Source!="wcgop" & bioData$Length<=60 & !ind,]
bioDat2 <- bioDat2[!is.na(bioDat2$Weight),]
sources <- unique(bioDat2$Source)
cols <- c("black","blue","green3","red","purple")
cols <- adjustcolor(cols,alpha.f<- 0.6)
plot(bioDat2$Length,bioDat2$Weight,pch=20,type="n")
for(i in length(sources):1) {
	tmp <- bioDat2[bioDat2$Source==sources[i],]
	points(jitter(tmp$Length,factor=7),tmp$Weight,col=cols[i],pch=16)
}
legend("topleft",c("NWFSC","Triennial","BDS","At-Sea"),col=cols,pch=20)



logLength <- log(bioDat2$Length)
lw.lm <- lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit)
lw.lm.nwfsc <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="nwfscCombo")
lw.lm.tri <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="triennial")
lw.lm.bds <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="bds")
lw.lm.atsea <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="atsea")
#lw.lm.wcgop <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="wcgop")
allSex <- rbind(coefficients(lw.lm),coefficients(lw.lm.nwfsc),coefficients(lw.lm.tri),coefficients(lw.lm.bds),coefficients(lw.lm.atsea))# ,coefficients(lw.lm.wcgop))
rownames(allSex) <- c("all","nwfsc","tri","bds","atsea")# ,"wcgop")
c(paste("a =",format(exp(coefficients(lw.lm)[1]),digits=5)),paste("b =",round(coefficients(lw.lm)[2],4)))



doPNG <- T
wd<-6.5;ht<-3
if(doPNG) {png("../Writeup/Figures/weightAtLengthPred.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(1,2),mar=c(1.5,1.5,2,0.5),oma=c(2,2,0,0))
lens <- 1:70
logLength <- log(bioDat2$Length)

lw.lm <- lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Sex=="f")
lw.lm.nwfsc <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="nwfscCombo"&bioDat2$Sex=="f")
lw.lm.tri <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="triennial"&bioDat2$Sex=="f")
lw.lm.bds <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="bds"&bioDat2$Sex=="f")
lw.lm.atsea <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="atsea"&bioDat2$Sex=="f")
#lw.lm.wcgop <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="wcgop"&bioDat2$Sex=="f")
fem <- rbind(coefficients(lw.lm),coefficients(lw.lm.nwfsc),coefficients(lw.lm.tri),coefficients(lw.lm.bds),coefficients(lw.lm.atsea))# ,coefficients(lw.lm.wcgop))
rownames(fem) <- c("all","nwfsc","tri","bds","atsea")# ,"wcgop")

plot(lens,exp(predict(lw.lm,newdata=data.frame(logLength=log(lens)))), type="l", lwd=7,col="grey",xlab="",ylab="",las=1,main="Female",xlim=c(0,70),ylim=c(0,5.1))
lines(lens,exp(predict(lw.lm.nwfsc,newdata=data.frame(logLength=log(lens)))),col=cols[1],lwd=3)
lines(lens,exp(predict(lw.lm.tri,newdata=data.frame(logLength=log(lens)))),col=cols[2],lwd=3)
lines(lens,exp(predict(lw.lm.bds,newdata=data.frame(logLength=log(lens)))),col=cols[3],lwd=3)
lines(lens,exp(predict(lw.lm.atsea,newdata=data.frame(logLength=log(lens)))),col=cols[4],lwd=3)
lines(lens,0.00000545*lens^3.28781,col="black",lty=2,lwd=3)
legend("topleft",c("All","NWFSC","Triennial","BDS","At-Sea","2011 Assessment"),col=c("grey",cols[-5],"black"),lwd=c(5,3,3,3,3,3),lty=c(1,1,1,1,1,2),cex=0.7,seg.len=3)
legend("bottomright",c(paste("a =",format(exp(coefficients(lw.lm)[1]),digits=5)),paste("b =",round(coefficients(lw.lm)[2],4))),pch=c(NA,NA),cex=0.8,bty="n")

lw.lm <- lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Sex=="m")
lw.lm.nwfsc <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="nwfscCombo"&bioDat2$Sex=="m")
lw.lm.tri <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="triennial"&bioDat2$Sex=="m")
lw.lm.bds <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="bds"&bioDat2$Sex=="m")
lw.lm.atsea <-lm(log(bioDat2$Weight) ~ logLength,na.action=na.omit,subset=bioDat2$Source=="atsea"&bioDat2$Sex=="m")
#logLength <- log(bioData$Length)
#lw.lm.wcgop <-lm(log(bioData$Weight) ~ logLength,na.action=na.omit,subset=bioData$Source=="wcgop"&bioData$Sex=="m")
#logLength <- log(bioDat2$Length)
male <- rbind(coefficients(lw.lm),coefficients(lw.lm.nwfsc),coefficients(lw.lm.tri),coefficients(lw.lm.bds),coefficients(lw.lm.atsea))# ,coefficients(lw.lm.wcgop))
rownames(male) <- c("all","nwfsc","tri","bds","atsea") #,"wcgop")

lens <- 1:70
plot(lens,exp(predict(lw.lm,newdata=data.frame(logLength=log(lens)))), type="l", lwd=7,col="grey",xlab="",ylab="",las=1,main="Male",ylim=c(0,5.1),xlim=c(0,70))
lines(lens,exp(predict(lw.lm.nwfsc,newdata=data.frame(logLength=log(lens)))),col=cols[1],lwd=3)
lines(lens,exp(predict(lw.lm.tri,newdata=data.frame(logLength=log(lens)))),col=cols[2],lwd=3)
lines(lens,exp(predict(lw.lm.bds,newdata=data.frame(logLength=log(lens)))),col=cols[3],lwd=3)
lines(lens,exp(predict(lw.lm.atsea,newdata=data.frame(logLength=log(lens)))),col=cols[4],lwd=3)
lines(lens,0.00001188*lens^3.06631,col="black",lty=2,lwd=3)
legend("topleft",c("All","NWFSC","Triennial","BDS","At-Sea","2011 Assessment"),col=c("grey",cols[-5],"black"),lwd=c(5,3,3,3,3,3),lty=c(1,1,1,1,1,2),cex=0.7,seg.len=3)
legend("bottomright",c(paste("a =",format(exp(coefficients(lw.lm)[1]),digits=5)),paste("b =",round(coefficients(lw.lm)[2],4))),pch=c(NA,NA),cex=0.8,bty="n")
mtext("Length (cm)",outer=T,side=1,line=0.8)
mtext("Weight (kg)",outer=T,side=2,line=0.8)

if(doPNG) dev.off()



#############################################
# Length at age
setwd("C:\\NOAA2015\\Widow\\Data")
library(nwfscSurvey)
VB.fn <- function(age,Linf,k,t0) {
    out <- Linf*(1-exp(-k*(age-t0)))
    return(out)
}
VBopt.fn <- function(x,age,lengths) {sum((lengths-VB.fn(age,Linf=x[1],k=x[2],t0=x[3]))^2)}


bioData <- read.csv("Biological/allBioData.csv")
table(bioData$Age,useNA="ifany")
table(!is.na(bioData$Length),!is.na(bioData$Age))

ageData <- bioData[!is.na(bioData$Age) & !is.na(bioData$Length),]

plot(bioData$Age,bioData$Length,pch=20,xlim=c(0,60))
ageData[ageData$Age<2,]
#I'm not sure I trust any of these observations. The age fish is 41 cm. The age 0 fish came from very deep water (which could have been hold overs from previous although none seen in previous)
#I may want to remove these from the age-at-length

sources <- unique(ageData$Source)
sourceNames <- c("NWFSC survey","Triennial","BDS","At-Sea")
cols <- c("black","blue","green3","red")
cols <- adjustcolor(cols,alpha.f<- 0.5)

doPNG <- T
wd<-6.5;ht<-8.5
if(doPNG) {png("../Writeup/Figures/LengthAgeAll.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfrow=c(2,1),mar=c(3,3,3,1),oma=c(1,1,0,0))
a.split <- split(ageData,ageData$Sex)
plot(a.split[["f"]]$Age,a.split[["f"]]$Length,pch=20,type="n",xlim=c(0,60),ylim=c(0,60),las=1,xlab="",ylab="",main="Female")
for(i in length(sources):1) {
	tmp <- a.split[["f"]][a.split[["f"]]$Source==sources[i],]
	points(jitter(tmp$Age,factor=0.8),tmp$Length,col=cols[i],pch=16)
}
xpar <- optim(c(48,0.18,0),VBopt.fn,age=a.split[["f"]]$Age,lengths=a.split[["f"]]$Length)$par
lines(0:60,VB.fn(0:60,Linf=xpar[1],k=xpar[2],t0=xpar[3]),col="yellow2",lwd=3)
legend("bottomright",c("NWFSC","Triennial","BDS","At-Sea"),col=cols,pch=20)
legend("bottomleft",c(paste("Linf=",round(xpar[1],2)),
	                  paste("k=",round(xpar[2],2)),
	                  paste("t0=",round(xpar[3],2))),pch=NA,bty="n",cex=0.8)
plot(a.split[["m"]]$Age,a.split[["m"]]$Length,pch=20,type="n",xlim=c(0,60),ylim=c(0,60),las=1,xlab="",ylab="",main="Male")
for(i in length(sources):1) {
	tmp <- a.split[["m"]][a.split[["m"]]$Source==sources[i],]
	points(jitter(tmp$Age,factor=0.8),tmp$Length,col=cols[i],pch=16)
}
xpar <- optim(c(48,0.18,0),VBopt.fn,age=a.split[["m"]]$Age,lengths=a.split[["m"]]$Length)$par
lines(0:60,VB.fn(0:60,Linf=xpar[1],k=xpar[2],t0=xpar[3]),col="yellow2",lwd=3)
legend("bottomright",c("NWFSC","Triennial","BDS","At-Sea"),col=cols,pch=20)
legend("bottomleft",c(paste("Linf=",round(xpar[1],2)),
	                  paste("k=",round(xpar[2],2)),
	                  paste("t0=",round(xpar[3],2))),pch=NA,bty="n",cex=0.8)
if(doPNG) dev.off()


doPNG <- T
wd<-6.5;ht<-8.5
if(doPNG) {png("../Writeup/Figures/LengthAgeEach.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mfcol=c(4,2),mar=c(2,4,2,1),oma=c(2,2,0,0))
for(i in length(sources):1) {
	tmp <- ageData[ageData$Source==sources[i] & ageData$Sex=="f",]
	plot(jitter(tmp$Age,factor=0.8),tmp$Length,col=cols[i],pch=16,xlim=c(0,60),ylim=c(0,60),ylab="",xlab="")
	if(i==4) title(main="Female")
	title(ylab=sourceNames[i],cex.lab=1.2)
	xpar <- optim(c(48,0.18,0),VBopt.fn,age=tmp$Age,lengths=tmp$Length)$par
    lines(0:60,VB.fn(0:60,Linf=xpar[1],k=xpar[2],t0=xpar[3]),col="yellow2",lwd=3)
    legend("bottomright",c(paste("Linf=",round(xpar[1],2)),
	                  paste("k=",round(xpar[2],2)),
	                  paste("t0=",round(xpar[3],2))),pch=NA,bty="n",cex=0.8)
}
for(i in length(sources):1) {
	tmp <- ageData[ageData$Source==sources[i] & ageData$Sex=="m",]
	plot(jitter(tmp$Age,factor=0.8),tmp$Length,col=cols[i],pch=16,xlim=c(0,60),ylim=c(0,60),ylab="",xlab="")
	if(i==4) title(main="Male")
	xpar <- optim(c(48,0.18,0),VBopt.fn,age=tmp$Age,lengths=tmp$Length)$par
    lines(0:60,VB.fn(0:60,Linf=xpar[1],k=xpar[2],t0=xpar[3]),col="yellow2",lwd=3)
    legend("bottomright",c(paste("Linf=",round(xpar[1],2)),
	                  paste("k=",round(xpar[2],2)),
	                  paste("t0=",round(xpar[3],2))),pch=NA,bty="n",cex=0.8)
}
mtext("Age",side=1,outer=T,line=0.5)
mtext("Length (cm)",side=2,outer=T,line=0.5)
if(doPNG) dev.off()


doPNG <- T
wd<-6.5;ht<-4.5
if(doPNG) {png("../Writeup/Figures/VarLengthAge.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(4,4,2,4))
tmp <- ageData[ageData$Age<60,]
res <- varLengthAtAge.fn(data.frame(AGE_YRS=tmp$Age,LENGTH_CM=tmp$Length,SEX=tmp$Sex),
	                     ageBin=1,bySex=T,parStart=c(48,0.18,0),estVB=T,legX="topright")
if(doPNG) {dev.off()}


######################################################################
#Maturity
library(gtools)
setwd("C:\\NOAA2015\\Widow\\Data")
source("C:\\NOAA2015\\Widow\\Data\\Rcode\\Functions\\maturityFunctions.R")
mat.fn <- function(x,len) {
    omega3 <- x[1]
    omega4 <- x[2]
    den <- 1+exp(omega3*(len-omega4))
    return(1/den)
}

mat <- read.csv("Biological/Maturity/BarssEchevarriaMaturityData.csv")
matIndiv <- data.frame(Length=NULL,I.M=NULL,Source=NULL)
for(i in 1:nrow(mat)) {
	if(!is.na(mat[i,"ORimmature"])&mat[i,"ORimmature"]>0)	matIndiv <- rbind(matIndiv,data.frame(Length=mat[i,"Length"], I.M=rep("I",mat[i,"ORimmature"]), Source="OR"))
	if(!is.na(mat[i,"ORmature"])&mat[i,"ORmature"]>0)	matIndiv <- rbind(matIndiv,data.frame(Length=mat[i,"Length"], I.M=rep("M",mat[i,"ORmature"]), Source="OR"))
	if(!is.na(mat[i,"CAimmature"])&mat[i,"CAimmature"]>0)	matIndiv <- rbind(matIndiv,data.frame(Length=mat[i,"Length"], I.M=rep("I",mat[i,"CAimmature"]), Source="CA"))
	if(!is.na(mat[i,"CAmature"])&mat[i,"CAmature"]>0)	matIndiv <- rbind(matIndiv,data.frame(Length=mat[i,"Length"], I.M=rep("M",mat[i,"CAmature"]), Source="CA"))
}
matIndiv$Mature <- NA
matIndiv$Mature[matIndiv$I.M=="I"] <- 0
matIndiv$Mature[matIndiv$I.M=="M"] <- 1

doPNG <- F
wd<-6.5;ht<-4.5
if(doPNG) {png("Biological/Maturity/maturityByState.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(5,4,3,3))
out <- plotMaturity.fn(matIndiv, theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, labHadj=0.65, legYinter=1.5, legTitle="Source")
if(doPNG) dev.off()

out <- plotPredMaturityLength(lens=23:54,matIndiv[matIndiv$Source=="CA",],start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F)
outPredCA <- predLogitLen(out$par,lengths=23:54)
cbind(23:54,outPred)
omega3 <- -1*out$par[2]
omega4 <- -1*out$par[1]/out$par[2]
c(omega3,omega4)
out <- plotPredMaturityLength(lens=23:54,matIndiv[matIndiv$Source=="OR",],start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F)
outPredOR <- predLogitLen(out$par,lengths=23:54)
cbind(23:54,outPred)
omega3 <- -1*out$par[2]
omega4 <- -1*out$par[1]/out$par[2]
c(omega3,omega4)
av <- apply(cbind(outPredCA,outPredOR),1,mean)
cbind(23:54,outPredCA,outPredOR,av)
lines(23:54,av)

#set under 30cm to immature (age 2 and less fish)
tmp <- matIndiv
tmp$I.M[tmp$Length<30] <- "I"
tmp$Mature[tmp$Length<30] <- 0
doPNG <- F
wd<-6.5;ht<-5
if(doPNG) {png("Biological/Maturity/maturityByStateUnder30Im.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(5,4,3,3))
out <- plotMaturity.fn(tmp, theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, labHadj=0.65, legYinter=1.5, legTitle="Source")
if(doPNG) dev.off()

#####################################################################################################
### Equal weighting to each state
tmp <- mat
tmp$ORimmature <- tmp$ORimmature*(sum(mat$CAimmature,na.rm=T)+sum(mat$CAmature,na.rm=T))
tmp$ORmature <- tmp$ORmature*(sum(mat$CAimmature,na.rm=T)+sum(mat$CAmature,na.rm=T))
tmp$CAimmature <- tmp$CAimmature*(sum(mat$ORimmature,na.rm=T)+sum(mat$ORmature,na.rm=T))
tmp$CAmature <- tmp$CAmature*(sum(mat$ORimmature,na.rm=T)+sum(mat$ORmature,na.rm=T))
tmpIndiv <- data.frame(Length=NULL,I.M=NULL,Source=NULL)
for(i in 1:nrow(tmp)) {
	if(!is.na(tmp[i,"ORimmature"])&tmp[i,"ORimmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Length=tmp[i,"Length"], I.M=rep("I",tmp[i,"ORimmature"]), Source="OR"))
	if(!is.na(tmp[i,"ORmature"])&tmp[i,"ORmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Length=tmp[i,"Length"], I.M=rep("M",tmp[i,"ORmature"]), Source="OR"))
	if(!is.na(tmp[i,"CAimmature"])&tmp[i,"CAimmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Length=tmp[i,"Length"], I.M=rep("I",tmp[i,"CAimmature"]), Source="CA"))
	if(!is.na(tmp[i,"CAmature"])&tmp[i,"CAmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Length=tmp[i,"Length"], I.M=rep("M",tmp[i,"CAmature"]), Source="CA"))
}
tmpIndiv$Mature <- NA
tmpIndiv$Mature[tmpIndiv$I.M=="I"] <- 0
tmpIndiv$Mature[tmpIndiv$I.M=="M"] <- 1
#tmpIndiv$I.M[tmpIndiv$Length<30] <- "I"
#tmpIndiv$Mature[tmpIndiv$Length<30] <- 0

doPNG <- T
wd<-6.5;ht<-5
if(doPNG) {png("Biological/Maturity/maturityByStateEqualWeight.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(5,4,3,3))
out <- plotMaturity.fn(matIndiv, theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, labHadj=0.65, legYinter=1.5, legTitle="Source")
if(doPNG) dev.off()

out <- plotPredMaturityLength(lens=23:54,tmpIndiv,start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F)
outPred <- predLogitLen(out$par,lengths=23:54)
cbind(23:54,outPred)
omega3 <- -1*out$par[2]
omega4 <- -1*out$par[1]/out$par[2]
c(omega3,omega4)




matAge <- read.csv("Biological/Maturity/BarssEchevarriaAgeMaturityData.csv")
matAgeIndiv <- data.frame(Age=NULL,I.M=NULL,Source=NULL)
for(i in 1:nrow(matAge)) {
	if(!is.na(matAge[i,"OrImmature"])&matAge[i,"OrImmature"]>0)	matAgeIndiv <- rbind(matAgeIndiv,data.frame(Age=matAge[i,"Age"], I.M=rep("I",matAge[i,"OrImmature"]), Source="OR"))
	if(!is.na(matAge[i,"OrMature"])&matAge[i,"OrMature"]>0)	matAgeIndiv <- rbind(matAgeIndiv,data.frame(Age=matAge[i,"Age"], I.M=rep("M",matAge[i,"OrMature"]), Source="OR"))
	if(!is.na(matAge[i,"CaImmature"])&matAge[i,"CaImmature"]>0)	matAgeIndiv <- rbind(matAgeIndiv,data.frame(Age=matAge[i,"Age"], I.M=rep("I",matAge[i,"CaImmature"]), Source="CA"))
	if(!is.na(matAge[i,"CaMature"])&matAge[i,"CaMature"]>0)	matAgeIndiv <- rbind(matAgeIndiv,data.frame(Age=matAge[i,"Age"], I.M=rep("M",matAge[i,"CaMature"]), Source="CA"))
}
matAgeIndiv$Mature <- NA
matAgeIndiv$Mature[matAgeIndiv$I.M=="I"] <- 0
matAgeIndiv$Mature[matAgeIndiv$I.M=="M"] <- 1

matAgeIndiv$Length <- matAgeIndiv$Age
n <- 1000
tmp2 <- data.frame(Age=rep(2,n), I.M=rep("I",n), Source="CA", Mature=0,Length=rep(2,n))
doPNG <- F
wd<-6.5;ht<-5
if(doPNG) {png("Biological/Maturity/maturityAgeByState.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(5,4,3,3))
out <- plotMaturity.fn(matAgeIndiv, theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, labHadj=0.65, legYinter=1.5, legTitle="Source", xlabel="Age")
out <- plotPredMaturityLength(lens=2:23,rbind(tmp,tmp2),start=c(-6,1.2,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="blue",lwd=3)
lines(matAge$Age,matAge$X2011Assess,col="red",lwd=2)
if(doPNG) dev.off()
outPred <- predLogitLen(out$par,lengths=2:23)
cbind(matAge$X2011Assess,outPred)


outPred <- predLogitLen(c(-6,1.2,1),lengths=2:22)
lines(2:22,outPred)
out <- plotPredMaturityLength(lens=2:23,matAgeIndiv[matAgeIndiv$Source=="CA",],start=c(-6,1.2,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=T)
tmp <- matAgeIndiv[matAgeIndiv$Source=="CA",]
n <- 1000
tmp2 <- data.frame(Age=rep(2,n), I.M=rep("I",n), Source="CA", Mature=0,Length=rep(2,n))
out <- plotPredMaturityLength(lens=2:23,rbind(tmp,tmp2),start=c(-6,1.2,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=T)
out <- plotMaturity.fn(rbind(matAgeIndiv,tmp2), theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, labHadj=0.65, legYinter=1.5, legTitle="Source", xlabel="Age")

#equal weighting of CA and OR
tmp <- matAge
tmp$OrImmature <- tmp$OrImmature*(sum(matAge$CaImmature,na.rm=T)+sum(matAge$CaMature,na.rm=T))
tmp$OrMature <- tmp$OrMature*(sum(matAge$CaImmature,na.rm=T)+sum(matAge$CaMature,na.rm=T))
tmp$CaImmature <- tmp$CaImmature*(sum(matAge$OrImmature,na.rm=T)+sum(matAge$OrMature,na.rm=T))
tmp$CaMature <- tmp$CaMature*(sum(matAge$OrImature,na.rm=T)+sum(matAge$OrMature,na.rm=T))
tmpIndiv <- data.frame(Age=NULL,I.M=NULL,Source=NULL)
for(i in 1:nrow(tmp)) {
	if(!is.na(tmp[i,"OrImmature"])&tmp[i,"OrImmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Age=tmp[i,"Age"], I.M=rep("I",tmp[i,"OrImmature"]), Source="OR"))
	if(!is.na(tmp[i,"OrMature"])&tmp[i,"OrMature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Age=tmp[i,"Age"], I.M=rep("M",tmp[i,"OrMature"]), Source="OR"))
	if(!is.na(tmp[i,"CaImmature"])&tmp[i,"CaImmature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Age=tmp[i,"Age"], I.M=rep("I",tmp[i,"CaImmature"]), Source="CA"))
	if(!is.na(tmp[i,"CaMature"])&tmp[i,"CaMature"]>0)	tmpIndiv <- rbind(tmpIndiv,data.frame(Age=tmp[i,"Age"], I.M=rep("M",tmp[i,"CaMature"]), Source="CA"))
}
tmpIndiv$Mature <- NA
tmpIndiv$Mature[tmpIndiv$I.M=="I"] <- 0
tmpIndiv$Mature[tmpIndiv$I.M=="M"] <- 1
#tmpIndiv$I.M[tmpIndiv$Length<30] <- "I"
#tmpIndiv$Mature[tmpIndiv$Length<30] <- 0
tmpIndiv$Length <- tmpIndiv$Age
n <- 1000
tmp2 <- data.frame(Age=rep(2,n), I.M=rep("I",n), Source="CA", Mature=0,Length=rep(2,n))


doPNG <- T
wd<-6.5;ht<-4.5
if(doPNG) {png("Biological/Maturity/maturityAgeByStateEqual.png",height=ht,width=wd,pointsize=10,units="in",res=300)}
if(!doPNG) {windows(height=ht,width=wd)}
par(mar=c(5,4,3,3))
outState <- plotMaturity.fn(rbind(matAgeIndiv,data.frame(Age=0,I.M=NA,Source=NA,Mature=NA,Length=0)), theFactor="Source", cols=c("orange","black"),
                factorLab=c("CA","OR"),
                alpha=0.6, inch=0.1, ht=0.1, legx=100, legy=100, labHadj=0.65,,xlabel="Age")  #, legYinter=1.5, legTitle="Source")
outAll <- plotPredMaturityLength(lens=0:23,tmpIndiv,start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="purple",lwd=3)
legend("right",c("CA","OR","All"),col=c("orange","black","purple"),lty=1,lwd=3,y.intersp = 1.5,cex=0.8,title="Source")
if(doPNG) dev.off()

outPred <- predLogitLen(outAll$par,lengths=2:23)
cbind(2:23,outPred)
omega3 <- -1*outAll$par[2]
omega4 <- -1*outAll$par[1]/outAll$par[2]
c(omega3,omega4)

outCA <- plotPredMaturityLength(lens=0:23,matAgeIndiv[matAgeIndiv$Source=="CA",],start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="orange",lwd=3)
outPredCA <- predLogitLen(outCA$par,lengths=0:23)
c(-1*outCA$par[2], -1*outCA$par[1]/outCA$par[2])
outOR <- plotPredMaturityLength(lens=0:23,matAgeIndiv[matAgeIndiv$Source=="OR",],start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="black",lwd=3)
outPredOR <- predLogitLen(outOR$par,lengths=0:23)
c(-1*outOR$par[2], -1*outOR$par[1]/outOR$par[2])
round(cbind(0:23,outPredCA,outPredOR),4)



#with Immature at zero added in
#out <- plotPredMaturityLength(lens=2:23,rbind(matAgeIndiv[matAgeIndiv$Source=="CA",],tmp2),start=c(-6,1.2,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="blue",lwd=3)
out <- plotPredMaturityLength(lens=2:23,rbind(tmpIndiv,tmp2),start=c(-6.2,0.179,1),low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,col="purple")
outPred <- predLogitLen(out$par,lengths=2:23)
cbind(2:23,outPred)
omega3 <- -1*out$par[2]
omega4 <- -1*out$par[1]/out$par[2]
c(omega3,omega4)



or.out = glm(cbind(ORimmature,ORmature) ~ Length, family=binomial(logit), data=mat)
omega3 <- -1*or.out$coeff["Length"]
omega4 <- -1*or.out$coeff["(Intercept)"]/or.out$coeff["Length"]
c(omega3,omega4)

ca.out = glm(cbind(CAimmature,CAmature) ~ Length, family=binomial(logit), data=mat)
omega3 <- ca.out$coeff["Length"]
omega4 <- -1*ca.out$coeff["(Intercept)"]/ca.out$coeff["Length"]
c(omega3,omega4)

all.out <- glm(cbind(CAimmature,CAmature) ~ Length, family=binomial(logit), data=mat)
omega3 <- -1*all.out$coeff["Length"]
omega4 <- -1*all.out$coeff["(Intercept)"]/all.out$coeff["Length"]
c(omega3,omega4)

#use OR samples only for maturity curve because it agrees more with later sources (that didn't have data)



#Fecundity
#From Dick (2009)
#Intercept = 286.2 * 1000 = 286200 kg
#Slope = 0.0026 * 1000000 = 2600.0 kg

ag <- 286.2
bg <- 0.0026
wtg <- seq(0,4000,length=100)
plot(wtg,wtg*(ag+bg*wtg),ylab="Eggs",xlab="Weight (g)",type="l")

a <- ag*1000
b <- bg*1e6
wt <- wtg/1000
#wt <- seq(0,4,length=100)  #kg
plot(wt,wt*(a+b*wt),ylab="Eggs",xlab="Weight (kg)",type="l")

eggs <- wt*(a+b*wt)
summary(lm(eggs~wt))

#fromBoehlert et al (1982) reported in barss & Echeverria
#F = 605.71*W - 261830  (wt in grams) (based on 64 observations from Dec1980-Jan1981)
eggs <- 605.71*wtg - 261830

ag <- 286.2
bg <- 0.0026
wtg <- seq(0,4000,length=100)
plot(wtg,wtg*(ag+bg*wtg),ylab="Eggs",xlab="Weight (g)",type="l",ylim=range(eggs))
lines(wtg,eggs,col="blue")

plot(wtg,wtg*(ag+bg*wtg)/wtg,ylab="Eggs/gram",xlab="Weight (g)",type="l",ylim=c(0,540))
lines(wtg,eggs/wtg,col="blue")




































#eliminate those that are beyond 1 and 99% quanitles at length
bioDat2 <- bioData[!is.na(bioData$Weight),c("Length","Weight")]
tmp.split <- split(bioDat2,round(bioDat2$Length,0))
tmp <- lapply(tmp.split,function(x){
							   if(length(x)>30) {
							   	 cat(sum(x<=quantile(x,0.01,na.rm=T) & x>=quantile(x,0.99,na.rm=T)),"observations removed for length",round(x$Length[1],0),"\n")
							   	 return(x[x>quantile(x,0.01,na.rm=T) & x<quantile(x,0.99,na.rm=T)])
							   } else{ return(x) }
})
plot(bioDat2$Length,bioDat2$Weight,pch=16,type="n")
for(i in 1:length(sources)) {
	tmp <- bioDat2[bioDat2$Source==sources[i],]
	points(jitter(tmp$Length,factor=7),tmp$Weight,col=cols[i],pch=16)

}
