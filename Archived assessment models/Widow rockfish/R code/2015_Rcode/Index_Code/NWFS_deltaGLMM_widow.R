
rm(list=ls())                                    
library(stats)
library(runjags)
library(R2jags)
library(coda)
library(superdiag)
library(pscl)
load.module("glm") # this is loading a specific library in JAGS / BUGS that implements a conditional sampler
runif(1)

# Install package
library("devtools")
install_github("nwfsc-assess/nwfscDeltaGLM", ref="1.0.0")
library(nwfscDeltaGLM)


#####################
### Load/mod data ###
#####################
#THis file had all rows without weights already removed
setwd("C:/Assessments/Widow2015/Data/TrawlSurvey/NWFSCsurvey")
species = "widow rockfish"
SA3 <- read.csv("SA3.csv")
nwfsc.widow<-read.csv("Widow2015_HaulCatchWt&Effort.csv")
nwfsc.widow<-nwfsc.widow[!is.na(nwfsc.widow$AREA_SWEPT_HA),]
nwfsc.widow = cbind(nwfsc.widow, "YEAR"=as.numeric(substr(nwfsc.widow[,'PROJECT_CYCLE'],7,10)))
nwfsc.widow = cbind(nwfsc.widow, "AREA_SWEPT_MSQ"=nwfsc.widow[,'AREA_SWEPT_HA']*1e4)
colnames(nwfsc.widow)[colnames(nwfsc.widow)=="SURVEY_PASS"]<-"PASS"
nwfsc.widow$PASS = nwfsc.widow$PASS #-1.5

# ECEs
#Which = which(nwfsc.widow[,'HAUL_WT_KG']>500 )
#nwfsc.widow[Which,c('PROJECT_CYCLE','YEAR','HAUL_WT_KG')]
#library(maps)
#library(mapdata)
#map("state", c("Oregon","Washington","California"), xlim=c(-126,-117), ylim=c(32,49), col="grey90", fill=TRUE, main="", mar=c(0,0,2,0),interior=TRUE)#, main=Year_Set[t])
#points( x=nwfsc.widow[Which,"BEST_LON_DD"], y=nwfsc.widow[Which,"BEST_LAT_DD"])

#####################################
#####################################
#####################################

# Load data and strata
masterDat = nwfsc.widow
strata.limits.nwfsc <- readIn(ncol=5,nlines=5)
  STRATA    NLat SLat MinDepth MaxDepth
  Sshallow 40.5 34.5  55       183
  Sdeep    40.5 34.5  183      400
  Nshallow 49.0 40.5  55       183
  Ndeep    49.0 40.5  183      400

#strata.limits.nwfsc <- readIn(ncol=5,nlines=3)
#  STRATA    NLat SLat MinDepth MaxDepth
#  shallow 49.0 34.5  55       183
#  deep    49.0 34.5  183      400


# Modify data slightly
#This directory will store the model results
setwd("C:/Assessments/Widow2015/Data//Index standardization/chain_4_strata_nwfsc_no_covariates/")
names(masterDat)[which(names(masterDat)=="HAUL_WT_KG")] = species
strata.limits<-strata.limits.nwfsc
if( TRUE ){
  X.bin = X.pos = as.matrix(masterDat[,'PASS',drop=FALSE]) - 1.5 # pass_1=(-0.5); pass_2=(0.5)
  nX.pos = nX.binomial = 1
  Covariates = list(positive=TRUE, binomial=TRUE)
}else{
  Covariates = list(positive=FALSE, binomial=FALSE)
}

# Preliminary data processing
DataList = processData()

# Define settings
#mcmc.control = list(chains=3, thin=1e1, burnin=1e4, iterToSave=1e4)
mcmc.control = list(chains=3, thin=100, burnin=1e5, iterToSave=1e5)
Parallel = TRUE   # If having trouble, try turning off parallel
modelStructure1 = list("StrataYear.positiveTows"="fixed", "VesselYear.positiveTows"="randomExpanded", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="fixed", "VesselYear.zeroTows"="randomExpanded", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure2 = list("StrataYear.positiveTows"="random", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", 
                        "StrataYear.zeroTows"="random", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero",
                         "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")
modelStructure3 = list("StrataYear.positiveTows"="zero", "VesselYear.positiveTows"="random", "Vessel.positiveTows"="zero", "StrataYear.zeroTows"="zero", "VesselYear.zeroTows"="random", "Vessel.zeroTows"="zero", "Catchability.positiveTows"="one", "Catchability.zeroTows"="zero", "year.deviations"="fixed", "strata.deviations"="fixed")

# Define models
mods = list()
set.seed(7594)
mods[[1]] = fitDeltaGLM(modelStructure=modelStructure1, #covariates=Covariates,
            likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[2]] = fitDeltaGLM(modelStructure=modelStructure1, #covariates=Covariates,
            likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str1.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[3]] = fitDeltaGLM(modelStructure=modelStructure2, #covariates=Covariates,
            likelihood = "gamma",  model.name = "deltaGLM.AFT.gamma.str2.txt",  mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[4]] = fitDeltaGLM(modelStructure=modelStructure2, #covariates=Covariates,
            likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[5]] = fitDeltaGLM(modelStructure=modelStructure3, #covariates=Covariates,
            likelihood = "gamma", model.name = "deltaGLM.AFT.gamma.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[6]] = fitDeltaGLM(modelStructure=modelStructure3, #covariates=Covariates,
            likelihood = "lognormal", model.name = "deltaGLM.AFT.log.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[7]] = fitDeltaGLM(modelStructure=modelStructure2, #covariates=Covariates,
            likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[8]] = fitDeltaGLM(modelStructure=modelStructure2, #covariates=Covariates, 
            likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str2.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[9]] = fitDeltaGLM(modelStructure=modelStructure3, #covariates=Covariates,
            likelihood = "gammaECE", model.name = "deltaGLM.AFT.gammaECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)
mods[[10]] =fitDeltaGLM(modelStructure=modelStructure3, #covariates=Covariates,
            likelihood = "lognormalECE", model.name = "deltaGLM.AFT.logECE.str3.txt", mcmc.control=mcmc.control, Parallel=Parallel, Species=species)


# Process MCMC output
data(SA3)
#cut out the first two model structure since they are not running
new.mods = list(mods[[3]], mods[[4]], mods[[5]], mods[[6]], mods[[7]], mods[[8]], mods[[9]], mods[[10]])
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

