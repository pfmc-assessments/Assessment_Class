#install.packages("devtools")
devtools::install_github("nwfsc-data/ExtractR")  #,ref="Adding")

setwd("C:/NOAA2015/Widow/Data/TrawlSurvey")
library(RODBC)
library(ExtractR)

#survDbDir <- "C:/NOAA2015/SurveyDatabase"
#source(paste(survDbDir,"Rcode/functions/Functions.R",sep="/"))  #eventually, this will be a package

#SPECIES_CODE          SCIENTIFIC_NAME          SPECIES ORGANISM_TYPE FMP_YN
#661         30220                  Sebastes entomelas                             widow rockfish          Fish    FMP

sp <- "Sebastes entomelas"
start <- 1998
end <- 2014

#Haul catch/weight for widow
singleSp <- queryDB(queryFilename="SingleSpecies.sql",
	                  db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
table(singleSp$PROJECT_CYCLE)
save(singleSp,file="../extractedData/wdowTrawlSurveyCatch.Rdat")
file.copy("../extractedData/wdowTrawlSurveyCatch.Rdat",
	       paste("../extractedData/wdowTrawlSurveyCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)

singleSpUnSat <- queryDB(queryFilename="SingleSpecies_UnSat.sql",
	                  db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
head(singleSpUnSat)
save(singleSp,file="../extractedData/wdowUnsatTrawlSurveyCatch.Rdat")
file.copy("../extractedData/wdowUnsatTrawlSurveyCatch.Rdat",
	       paste("../extractedData/wdowUnsatTrawlSurveyCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)


#Lengths
lengths <- queryDB(queryFilename="Lengths.sql",
	                  db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
save(lengths,file="../extractedData/wdowTrawlSurveyLengths.Rdat")
file.copy("../extractedData/wdowTrawlSurveyLengths.Rdat",
	       paste("../extractedData/wdowTrawlSurveyLengths_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)

lengthsUnSat <- queryDB(queryFilename="Lengths_UnSat.sql",
	                  db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
save(lengths,file="../extractedData/wdowTrawlSurveyLengthsUnSat.Rdat")
file.copy("../extractedData/wdowTrawlSurveyLengthsUnSat.Rdat",
	       paste("../extractedData/wdowTrawlSurveyLengthsUnSat_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)


#Length,Sex,Wt,Age
lwsa <- queryDB(queryFilename="LenAgeWt.sql",
                      db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
save(lwsa,file="../extractedData/wdowTrawlSurveyLwsa.Rdat")
file.copy("../extractedData/wdowTrawlSurveyLwsa.Rdat",
	       paste("../extractedData/wdowTrawlSurveyLwsa_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)

lwsaUnSat <- queryDB(queryFilename="LenAgeWt_UnSat.sql",
	                  db = "FRAM",
	                  uid="hicksal",
	                  sp=sp,start=start,end=end)
save(lwsaUnSat,file="../extractedData/wdowTrawlSurveyLwsaUnSat.Rdat")
file.copy("../extractedData/wdowTrawlSurveyLwsaUnSat.Rdat",
	       paste("../extractedData/wdowTrawlSurveyLwsaUnSat_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
	       overwrite=T)


#Length,Sex,Wt,Age, Maturity
lwsam <- queryDB(queryFilename="LenAgeWtMat.sql",
                      db = "FRAM",
	                  uid="hicksal",
                      sp=sp,start=start,end=end)
table(lwsam$MATURE_YN,useNA="ifany")
# save(lwsam,file="../extractedData/wdowTrawlSurveyLwsam.Rdat")
# file.copy("../extractedData/wdowTrawlSurveyLwsam.Rdat",
# 	       paste("../extractedData/wdowTrawlSurveyLwsam_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
# 	       overwrite=T)

lwsamUnSat <- queryDB(queryFilename="LenAgeWtMat_UnSat.sql",
	                  db = "FRAM",
	                  uid="hicksal",
                  sp=sp,start=start,end=end)
table(lwsamUnSat$MATURE_YN,useNA="ifany")
# save(lwsamUnSat,file="../extractedData/wdowTrawlSurveyLwsamUnSat.Rdat")
# file.copy("../extractedData/wdowTrawlSurveyLwsamUnSat.Rdat",
# 	       paste("../extractedData/wdowTrawlSurveyLwsamUnSat_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),
# 	       overwrite=T)






#Is there an issue with two 2013 tows?
devtools::install_github("nwfsc-assess/nwfscMapping")
# Load package
library(nwfscMapping)
library(PBSmapping)
data(westCoastLL)
data(WCstatesInlandPBS)
source("C:/Mapping/WestCoastMapping.R")


setwd("C:/NOAA2015/Widow/Data")
load("ExtractedData/wdowTrawlSurveyLengths.Rdat")
wdwLens <- lengths
rm(lengths)
#Length by depth and latitude
plot(wdwLens$DEPTH_M,wdwLens$LENGTH_CM,pch=16)

tmp <- wdwLens[wdwLens$LENGTH_CM<10,]
distance(tmp$HAUL_LONGITUDE_DD[1],tmp$HAUL_LATITUDE_DD[1],tmp$HAUL_LONGITUDE_DD[2],tmp$HAUL_LATITUDE_DD[2])



plotTheCoast <- function(xlims,ylims,plotDepth=F) {
    plotMap(westCoastLL,xlim=xlims,ylim=ylims,col="green4",bg="lightcyan2",tck=c(-0.02),cex=1.0,main="",xlab="",ylab="")
    addLines(WCstatesInlandPBS)
    #addLines(EEZ,lwd=2,col=gray(0.5))
    if(plotDepth) {
    	addLines(depth50m,col=gray(0.3))
    	addLines(depth1500m,col=gray(0.3))
    }
}

dat <- wdwLens[wdwLens$PROJECT_CYCLE=="Cycle 2013",c("HAUL_LONGITUDE_DD","HAUL_LATITUDE_DD","LENGTH_CM")]
names(dat) <- c("X","Y","Z")
dat <- as.EventData(data.frame(EID=1:nrow(dat),dat))
#This is dumb, probably best to plot mena length
ht <- 12;wd<-10
windows(height=ht,width=wd)
xlims <- c(-127.5,-116.45)
ylims <- c(32.0,48.5)
plotTheCoast(xlims,ylims)
addBubbles(dat,max.size=0.2,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero="x",legend.pos="topright",legend.cex=1.0,legend.title="Length (cm)",legend.breaks=c(15,30),cex=0.7,col=rgb(0,0,0,0.9))


#but, look at the issue of the smallest fish
dat <- tmp[tmp$PROJECT_CYCLE=="Cycle 2013",c("HAUL_LONGITUDE_DD","HAUL_LATITUDE_DD","LENGTH_CM")]
names(dat) <- c("X","Y","Z")
dat <- as.EventData(data.frame(EID=1:nrow(dat),dat))

ht <- 12;wd<-10
windows(height=ht,width=wd)
xlims <- c(-125.5,-122.45)
ylims <- c(38.0,40.5)
plotTheCoast(xlims,ylims)
addLines(depth300f)
addLines(depth100f)
addBubbles(dat,max.size=0.2,symbol.bg=rgb(0,0,0,0.5),symbol.fg="black",symbol.zero="x",legend.pos="topright",legend.cex=1.0,legend.title="Length (cm)",legend.breaks=c(15,30),cex=0.7,col=rgb(0,0,0,0.9))
