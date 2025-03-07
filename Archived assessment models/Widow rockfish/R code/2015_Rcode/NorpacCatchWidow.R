#Extract NORPAC domestic and foreign rougheye catch from at-sea fishery
#Needs the RODBC library
#Start R 2.14.0 in 32-bit mode for drivers to work
#set working directory to have subdirectories:
###extractedData: folder to hold raw extractions
###sql: with the sql script called Pacfin.Catch.query
###Catches: To save the catch results in
#This assigns the fleet type and summarizes over fleet and year

setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
library(ExtractR)
#source("Rcode/functions/Functions.R")
source("Rcode/functions/processNorpacCatchSpecies.R")
#options(digits=19)  #set this so that the full haul join number is displayed

if(F) {
convert to extractR
    ncatch <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/NORPAC.domestic.catch.query",db="NORPAC",uid="hicksa",start=1966,end=2014)
    save(ncatch,file="extractedData/NORPACdomesticCatch.Rdat")
    file.copy("extractedData/NORPACdomesticCatch.Rdat",paste("extractedData/NORPACdomesticCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}

if(F) {
    ncatch2 <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/NORPAC.domestic.catch.detailed.query",db="NORPAC",uid="hicksa",start=1990,end=2014)
    save(ncatch2,file="extractedData/NORPACdomesticCatchDetailed.Rdat")
    file.copy("extractedData/NORPACdomesticCatchDetailed.Rdat",paste("extractedData/NORPACdomesticCatchDetailed_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)
}






#Get NORPAC catches
setwd("C:/NOAA2015/Widow/Data")
#source("Rcode/functions/Functions.R")
source("Rcode/functions/processNorpacCatchSpecies.R")

load(file="extractedData/NORPACdomesticCatch.Rdat")
load(file="extractedData/NORPACdomesticCatchDetailed.Rdat")

ncatchAll <- ncatch
#ncatch <- ncatch[substring(ncatch$HAUL_DATE,1,4)%in%2003:2011,]


#Get species list
#nspecies <- queryDB(queryFilename="NORPACspecies",db="NORPAC",uid="hicksa",start=1966,end=2011,querydir="sql/")
#nspecies[grep("ROCK",nspecies$SPECIES_NAME),]
#354  SHORTRAKER/ROUGHEYE ROCKFISH
#307             ROUGHEYE ROCKFISH
#357         BLACKSPOTTED ROCKFISH
#308           RED BANDED ROCKFISH
#310           SILVERGRAY ROCKFISH
#326           SHORTRAKER ROCKFISH
#334               AURORA ROCKFISH
#305                WIDOW ROCKFISH
#511            Humboldt Squid


#process results, calculate bycatch rate, and fill in unsampled hauls
out <- processNorpacCatch(ncatch,species=305)
widowHauls <- NorpacUnsampHauls(ncatch,species=305)
widowHauls$percUnsampled <- round(100*widowHauls$Unsampled/widowHauls$Hauls,2)
save(out,file="extractedData/NorpacAtSeaCatch.Rdat")
table(ncatch$SPECIES) #no blackspotted identified
table(ncatch$SPECIES,substring(ncatch$HAUL_DATE,1,4))

load(file="extractedData/NorpacAtSeaCatch.Rdat")
out.yr <- aggregate(out$Catch,list(out$Sector,out$Year),sum)
names(out.yr) <- c("Sector","Year","Catch.MT")
widow <- out.yr
barplot(widow$Catch.MT,beside=T,names=1990:2014,ylab="Catch (mt)",xlab="Year",main="Widow catch in the hake at-sea fishery")

#the extrapolated to haul catch but not extrapolated to unsampled hauls
wd <- ncatchAll[ncatchAll$SPECIES==305,]
tapply(wd$EXTRAPOLATED,substring(wd$HAUL_DATE,1,4),sum)
table(substring(wd$HAUL_DATE,1,4))

#Look at extrapolations
load(file="extractedData/NORPACdomesticCatchDetailed.Rdat") #ncatch2
wd <- ncatch2[ncatch2$SPECIES==305 & !is.na(ncatch2$SPECIES),]
tapply(wd$EXTRAPOLATED,substring(wd$HAUL_DATE,1,4),sum)
widow.df <- data.frame(widow,sampledOnly=tapply(wd$EXTRAPOLATED,substring(wd$HAUL_DATE,1,4),sum)/1000)
widow.df$percWtNotSamp <- (widow.df$Catch.MT-widow.df$sampledOnly)/widow.df$Catch.MT

wd$expand.species <- round(wd$EXTRAPOLATED_WEIGHT/wd$SAMPLE_WEIGHT,2)
wd$expand.tow <- round(wd$OFFICIAL_TOTAL_CATCH/(wd$SAMPLE_SIZE/1000),2)
range(wd$expand.species-wd$expand.tow,na.rm=T)
boxplot(wd$expand.species-wd$expand.tow)

boxplot(wd$expand.species)
boxplot(wd$expand.species,log="y")
range(wd$expand.species,na.rm=T)

# windows(height=5,width=6.5)
# boxplot(split(wd$expand.species,factor(wd$YEAR,levels=2002:2012)),log="y",xlab="Year",ylab="Within tow expansion factor",main="Expansion factors for rougheye rockfish")

widow.df$medWithinTowExp <- round(tapply(wd$expand.species,wd$YEAR,median),2)
tapply(wd$expand.tow,wd$YEAR,median)














#do this by sector (CP and MS)
#From http://www.afsc.noaa.gov/FMA/domestic_haul_v.htm
#VESSEL_TYPE
#A one-digit numeric code that indicates whether the vessel processes fish or delivers it to a processing plant where:
#1 = a catcher processor vessel,
#2 = a mothership or a ship that receives unsorted codends from other vessels,
#3 = a catcher only vessel that delivers unprocessed fish to a shoreside or floating plant or vessel,
#4 = a mothership that receives sorted codends,
#5 = a vessel that sells the majority of their catch over the side to other fishing vessels who will utilize the fish for bait,
#6 = vessels that discard all catch from a haul; would be used for codend dumping of an entire haul (added January 2004).

wd.sec <- split(wd,wd$VESSEL_TYPE)
boxplot(split(wd$expand.species,wd$VESSEL_TYPE),log="y")
x <- tapply(wd$EXTRAPOLATED_WEIGHT,list(wd$YEAR,wd$VESSEL_TYPE),sum)
x
png("../WriteUp/Figures/CPandMScatch.png",width=6.5,height=4,units="in",res=600)
par(mar=c(4,4,1,1))
plot(as.numeric(row.names(x))-0.1, x[,1]/1000,type="h",lwd=3,ylab="Catch of Widow Rockfish (thousand t)",xlab="Year",las=1)
lines(as.numeric(row.names(x))+0.1, x[,2]/1000,type="h",lwd=3,col="green4")
legend("topright",c("Catcher-Processor","Mothership"),lwd=5,lty=1,col=c("black","green4"))
dev.off()

windows(height=5,width=6.5)
boxplot(split(wd.sec[[1]]$expand.species,factor(wd.sec[[1]]$YEAR,levels=2002:2012)),log="y",xlab="Year",ylab="Within tow expansion factor",main="Expansion factors for CP fleet")
windows(height=5,width=6.5)
boxplot(split(wd.sec[[2]]$expand.species,factor(wd.sec[[2]]$YEAR,levels=2002:2012)),log="y",xlab="Year",ylab="Within tow expansion factor",main="Expansion factors for MS fleet")

wd$exSpCat <- cut(wd$expand.species,c(0,1,2,3,4,5,10,20,100,Inf))
table(wd$YEAR,wd$exSpCat)

cbind(tapply(wd$SAMPLE_WEIGHT,wd$exSpCat,sum))
tapply(wd$EXTRAPOLATED_WEIGHT,list(wd$YEAR,wd$exSpCat),sum)

tapply(wd$EXTRAPOLATED_WEIGHT,list(wd$YEAR,wd$VESSEL_TYPE),sum)





















pAtSeaCatch$CP <- pAtSeaCatch$PROCESSOR == pAtSeaCatch$CATCHER
#ms <- pAtSeaCatch[pAtSeaCatch$PROCESSOR != pAtSeaCatch$CATCHER,]
tapply(pAtSeaCatch$TOTAL_WEIGHT,list(pAtSeaCatch$YEAR,pAtSeaCatch$CP),sum)



wd$SPECIFICHAUL <- paste(format(re[,1],digits=19),re$HAUL,sep="_")
hauls.fn <- function(x) {
    x <- x[(!is.na(x$OFFICIAL_TOTAL_CATCH)) & x$OFFICIAL_TOTAL_CATCH != 0,]  #only observations with an official catch
    Nhauls <- length(unique(x$SPECIFICHAUL))
    unsampled <- x[is.na(x$HAUL_SAMPLED_BY) | x$HAUL_SAMPLED_BY == 0,]
    Un.hauls <- length(unique(unsampled$SPECIFICHAUL))
    return(data.frame(hauls=Nhauls,unsampled=Un.hauls))
}
hauls <- lapply(split(re,re$YEAR),hauls.fn)


hauls.fn(re[re$Year==2012,])

table(re$YEAR)
table(pAtSeaCatch$YEAR)


table(pAtSeaCatch$YEAR,pAtSeaCatch$PROCESSOR != pAtSeaCatch$CATCHER)
t(table(re$VESSEL_TYPE,re$YEAR))


tmpP <- pAtSeaCatch[pAtSeaCatch$YEAR==1991,]
tmpR <- re[re$YEAR==1991 & !is.na(re$SPECIES),]
# I wonder if these are the days that unobserved tows are aggregated to? I don't understand why there is a difference in number of rows

split(tmpP,tmpP$PROCESSOR)
split(tmpR,factor(tmpR$VESSEL))

split(tmpP,tmpP$DAY)
split(tmpR,substring(tmpR$HAUL_DATE,9))

lapply(split(tmpP$TOTAL_WEIGHT,tmpP$PROCESSOR),sum)
lapply(split(tmpR$EXTRAPOLATED_WEIGHT,factor(tmpR$VESSEL)),sum)


tmpP <- pAtSeaCatch[pAtSeaCatch$YEAR==2009 & pAtSeaCatch$MONTH==6 & pAtSeaCatch$DAY==1,]
tmpR <- re[as.character(re$HAUL_DATE)=="2009-06-01" & !is.na(re$SPECIES),]




#NPAC4900_SPCOMP
load("extractedData/PacfinAtSeaCatch.Rdat")
data.frame(rougheye,sampledOnly=tapply(re$EXTRAPOLATED,substring(re$HAUL_DATE,1,4),sum)/1000,NPAC4900=tapply(pAtSeaCatch$TOTAL_WEIGHT,pAtSeaCatch$YEAR,sum))



table(re$PERCENT_RETAINED,re$YEAR)


load(file="extractedData/NORPACdomesticCatchDetailed.Rdat")

#first find unique tows and widow catch and tows with zero widow
dim(ncatch2)
ncatch2 <- ncatch2[!is.na(ncatch2$EXTRAPOLATED_WEIGHT),]
dim(ncatch2)
ncatch2$yr <- as.numeric(substring(ncatch2$HAUL_DATE,1,4))  #some records of YEAR are NA when there is a date, but those are removed above now
ncatch2$uniqueHaul <- paste(ncatch2$yr,ncatch2$CRUISE,ncatch2$VESSEL,ncatch2$VESSEL_TYPE,ncatch2$HAUL_DATE,ncatch2$HAUL,sep=".")


tows <- ncatch2
tows$widow <- 0
tows$widow[tows$SPECIES==305] <- 1
tows <- tows[rev(order(tows$widow)),]
tows <- tows[!duplicated(tows$uniqueHaul),]
table(tows$yr,tows$widow,useNA="ifany")
tows$EXTRAPOLATED_WEIGHT[tows$widow==0] <- 0
table(tows$yr,tows$EXTRAPOLATED_WEIGHT>0)

propPos <- table(tows$yr,tows$widow==1)
propPos <- apply(propPos,1,function(x){x[2]/sum(x)})
yr <- as.numeric(names(propPos))
plot(yr,propPos,pch=16)

widowTows <- tows[tows$widow==1,]
boxplot(split(widowTows$EXTRAPOLATED_WEIGHT,widowTows$yr),log="y")
