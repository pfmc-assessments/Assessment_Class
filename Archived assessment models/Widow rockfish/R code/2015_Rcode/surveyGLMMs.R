#library(stats)
#library(runjags)
library(coda)
#library(superdiag)
#library(R2jags)
#library(pscl)
#library(statmod)
#load.module("glm")


setwd("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey")
load("Index standardization\\stratified shelf-slope 2003-2014\\widow rockfish_FinalDiagnostics\\Save.RData")

modelNames <- c("Gamma - random","Lognormal - random",
				"Gamma - no strata year","Lognormal - no strata year",
				"Gamma - random ECE","Lognormal- random ECE",
                "Gamma - no strata year ECE","Lognormal - no strata year ECE")
shortNames <- c("G","LN",
                "G\nnoSY","LN\nnoSY",
                "G\n\nECE","LN\n\nECE",
                "G\nnoSY\nECE","LN\nnoSY\nECE")

dev <- data.frame(num=1:length(Save$mods),name=modelNames,
                  deviance=unlist(lapply(Save$mods,function(x){mean(x$BUGSoutput$sims.list$deviance)})))
dev[order(dev[,"deviance"]),]
  num                           name deviance
5   5             Gamma - random ECE 3179.680
6   6          Lognormal- random ECE 3217.455
8   8 Lognormal - no strata year ECE 3262.704
7   7     Gamma - no strata year ECE 3349.712
2   2             Lognormal - random 3473.686
4   4     Lognormal - no strata year 3490.303
1   1                 Gamma - random 3522.781
3   3         Gamma - no strata year 3574.934


tmp <- vector(mode="list",length=length(modelNames))
for(i in 1:length(tmp)) {
	tmp[[i]] <- Save$mods[[i]]$BUGSoutput$sims.list$deviance
}
names(tmp) <- modelNames

boxplot(tmp, names=rep("",length(tmp)),ylab="Deviance")
text(1:8,3048,shortNames,adj=c(0.5,1),xpd=NA)




my.wd <- "C:/NOAA2013/Rougheye/Data/Survey/GLMM/"
source(paste(my.wd,"bayesGLM v2.15.R",sep=""))
Letters = apply(MARGIN=1,FUN=paste,collapse="",expand.grid(letters,letters))






attach(Save$Data)
doMCMCDiags(my.wd,Save$mods)

#Triennial
load("Index standardization\\GLMM_shelf\\widow_rockfish_FinalDiagnostics\\Save.RData")
modelNames <- c("GammaFixed","LognormalFIxed",
                "GammaRandom","LognormalRandom",
                "GammaZero","LognormalZero")
shortNames <- c("GF","LNF",
                "GR","LNR",
                "GZ","LNZ")

dev <- data.frame(num=1:length(Save$mods),name=modelNames,
                  deviance=unlist(lapply(Save$mods,function(x){mean(x$BUGSoutput$sims.list$deviance)})))
dev[order(dev[,"deviance"]),]
