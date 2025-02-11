###################################
########### TRIENNIAL #############
###################################
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM", ref="1.0.0")
library(nwfscDeltaGLM)

library(stats)
library(runjags)
library(R2jags)
library(coda)
library(superdiag)
library(pscl)
load.module("glm") # this is loading a specific library in JAGS / BUGS that implements a conditional sampler
runif(1)

Dir<-"C:/Assessments/Widow2015/Data/TrawlSurvey/Triennial"
load(paste0(Dir, "/Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015.dmp"))

Data = Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015
colnames(Data)[colnames(Data)=="HAULJOIN"]<-"OP_CODE"
colnames(Data)[colnames(Data)=="START_LATITUDE"]<-"BEST_LAT_DD"
colnames(Data)[colnames(Data)=="GEAR_DEPTH"]<-"BEST_DEPTH_M"
Data$BEST_DEPTH_M[5703]<-Data$BOTTOM_DEPTH[5703]
Data$TIME.GMT <- Data$DISTANCE_FISHED *1000* Data$NET_WIDTH
colnames(Data)[colnames(Data)=="TIME.GMT"]<-"AREA_SWEPT_MSQ"

#Separate shelf and slope surveys
triShelf <- subset(Data,SURVEY=="Tri.Shelf")
triSlope <- subset(Data,SURVEY=="AFSC.Slope")

Which = which(triShelf[,'WEIGHT']>500 )
triShelf[Which,c('YEAR','WEIGHT', 'NUMBER_FISH')]
#Which = which(triShelf[,'HAUL_WT_KG']>500 )
#triShelf[Which,c('PROJECT_CYCLE','YEAR','HAUL_WT_KG')]
#library(maps)
#library(mapdata)
#map("state", c("Oregon","Washington","California"), xlim=c(-126,-117), ylim=c(32,49), col="grey90", fill=TRUE, main="", mar=c(0,0,2,0),interior=TRUE)#, main=Year_Set[t])
#points( x=triShelf[Which,"BEST_LON_DD"], y=triShelf[Which,"BEST_LAT_DD"])

#Select years
triShelf <- subset(triShelf, YEAR > 1977)
tri.early<- subset(triShelf, YEAR < 1995)
tri.late <- subset(triShelf, YEAR > 1992)

temp.triSlope <- subset(triSlope, YEAR > 1996)
triSlope 	  <- temp.triSlope[!is.na(temp.triSlope$AREA_SWEPT),]


setwd("C:/Assessments/Widow2015/Data/Index standardization/chain_GLMM_shelf")
masterDat = triShelf
strata.limits.tri <- readIn(ncol=5,nlines=3)
  STRATA   NLat SLat MinDepth MaxDepth
  shallow 49.0 34.5  55       183
  deep    49.0 34.5  183      400

species = "widow_rockfish"
names(masterDat)[which(names(masterDat)=="WEIGHT")] = species
strata.limits<-strata.limits.tri
DataList = processData()

#mcmc.control = list(chains=3, thin=1e1, burnin=1e4, iterToSave=1e4)
mcmc.control = list(chains=3, thin=100, burnin=1e5, iterToSave=1e5)
Parallel = TRUE   # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="randomExpanded", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="randomExpanded", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure2 = list("StrataYear.positiveTows"="random", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="random", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure3 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
mods = list()
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[2]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[3]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[4]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[5]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[6]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[7]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[8]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[9]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[10]]= fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

# Process MCMC output
data(SA3)
#cut out the first two model structure since they are not running
new.mods = list(mods[[3]], mods[[4]], mods[[5]], mods[[6]], mods[[7]], mods[[8]], mods[[9]], mods[[10]])
#new.mods = list(mods[[7]], mods[[8]], mods[[9]], mods[[10]])
doMCMCDiags(my.wd, new.mods)

# DIC
DIC = data.frame("DIC"=sapply( new.mods, FUN=function(List){ List$BUGSoutput$DIC }), "WAIC"=NA, 
                "Likelihood"=sapply( new.mods, FUN=function(List){ List$likelihood }), 
                "StrataYear"=sapply( new.mods, FUN=function(List){ List$modelStructure$StrataYear.zeroTows }) )
for( StrucI in 1:2){
for( LikeI in 1:2){
  NumI = 4*(StrucI-1) + LikeI
  File = paste0( "widow rockfish_FinalDiagnostics/Model=",NumI,"/")
  WAIC = read.csv( paste0(File,"WAIC.csv") )
  DIC[NumI,'WAIC'] = WAIC[1,'WAIC']  
}}
write.csv( DIC, file="DIC_table.csv")


###########################################################
###			  	Triennial Early
###########################################################
setwd("C:/Assessments/Widow2015/Data/Index standardization/chain_GLMM_shelf_early")                                    

masterDat = tri.early
strata.limits <- readIn(ncol=5,nlines=3)
  STRATA   NLat SLat MinDepth MaxDepth
  shallow 49.0 34.5  55       183
  deep    49.0 34.5  183      400

species = "widow_rockfish"
names(masterDat)[which(names(masterDat)=="WEIGHT")] = species
DataList = processData()

#mcmc.control = list(chains=3, thin=1e1, burnin=1e4, iterToSave=1e4)
Parallel = TRUE   # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="randomExpanded", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="randomExpanded", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure2 = list("StrataYear.positiveTows"="random", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="random", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure3 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
mods = list()
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[2]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[3]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[4]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[5]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[6]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[7]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[8]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[9]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[10]] =fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)

# Process MCMC output
data(SA3)
new.mods = list(mods[[3]], mods[[4]], mods[[5]], mods[[6]], mods[[7]], mods[[8]], mods[[9]], mods[[10]])
#new.mods = list(mods[[7]], mods[[8]], mods[[9]], mods[[10]])
doMCMCDiags(my.wd, mods)

# DIC
DIC = data.frame("DIC"=sapply( mods, FUN=function(List){ List$BUGSoutput$DIC }), "WAIC"=NA, 
                "Likelihood"=sapply( mods, FUN=function(List){ List$likelihood }), 
                "StrataYear"=sapply( mods, FUN=function(List){ List$modelStructure$StrataYear.zeroTows }) )
for( StrucI in 1:2){
for( LikeI in 1:3){
  NumI = 4*(StrucI-1) + LikeI
  File = paste0( "widow_rockfish_FinalDiagnostics/Model=",NumI,"/")
  WAIC = read.csv( paste0(File,"WAIC.csv") )
  DIC[NumI,'WAIC'] = WAIC[1,'WAIC']  
}}
write.csv( DIC, file="DIC_table.csv")

###########################################################
###			  	Triennial Late
###########################################################
setwd("C:/Assessments/Widow2015/Data/Index standardization/chain_GLMM_shelf_late")
masterDat = tri.late
strata.limits <- readIn(ncol=5,nlines=3)
  STRATA   NLat SLat MinDepth MaxDepth
  shallow 49.0 34.5  55       183
  deep    49.0 34.5  183      400

species = "widow_rockfish"
names(masterDat)[which(names(masterDat)=="WEIGHT")] = species
DataList = processData()

#mcmc.control = list(chains=3, thin=1e1, burnin=1e4, iterToSave=1e4)
Parallel = TRUE # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure2 = list("StrataYear.positiveTows"="random", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="random", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure3 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
mods = list()
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[2]] = fitDeltaGLM(modelStructure=modelStructure1,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[3]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[4]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[5]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[6]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[7]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[8]] = fitDeltaGLM(modelStructure=modelStructure2,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[9]] = fitDeltaGLM(modelStructure=modelStructure3,likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[10]] =fitDeltaGLM(modelStructure=modelStructure3,likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)


# Process MCMC output
data(SA3)
new.mods = list(mods[[3]], mods[[4]], mods[[5]], mods[[6]], mods[[7]], mods[[8]], mods[[9]], mods[[10]])
#new.mods = list(mods[[7]], mods[[8]], mods[[9]], mods[[10]])
doMCMCDiags(my.wd, new.mods)

# DIC
DIC = data.frame("DIC"=sapply( new.mods, FUN=function(List){ List$BUGSoutput$DIC }), "WAIC"=NA, 
                "Likelihood"=sapply( new.mods, FUN=function(List){ List$likelihood }), 
                "StrataYear"=sapply( new.mods, FUN=function(List){ List$modelStructure$StrataYear.zeroTows }) )
for( StrucI in 1:2){
for( LikeI in 1:2){
  NumI = 4*(StrucI-1) + LikeI
  File = paste0( "widow_rockfish_FinalDiagnostics/Model=",NumI,"/")
  WAIC = read.csv( paste0(File,"WAIC.csv") )
  DIC[NumI,'WAIC'] = WAIC[1,'WAIC']  
}}
write.csv( DIC, file="DIC_table.csv")




#Deviation Calculations
#Load up the Save.RData file for each run
mods = Save$mods
dev <- cbind(1:length(mods),unlist(lapply(mods,function(x){mean(x$BUGSoutput$sims.list$deviance)})))
dev[order(dev[,2]),]