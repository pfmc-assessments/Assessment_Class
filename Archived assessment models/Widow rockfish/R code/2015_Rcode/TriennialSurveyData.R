dont source
# devtools::install_github("nwfsc-assess/nwfscSurvey")
library(nwfscSurvey)

setwd("C:/NOAA2015/Widow/Data/TrawlSurvey/Triennial")
load("../SA3.dmp")

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

if(F) {
    load("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\Triennial\\Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015.dmp")
    dat <- Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015
    rm(Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015)
    write.csv(dat[dat$SURVEY=="Tri.Shelf",],file="TriennialCatch.csv")

    load("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\Triennial\\AK.Surveys.Bio.widow.04.Jun.2015.dmp")
    dat <- AK.Surveys.Bio.widow.04.Jun.2015
    rm(AK.Surveys.Bio.widow.04.Jun.2015)
    write.csv(dat$Lengths[dat$Lengths$SURVEY=="Tri.Shelf",],file="TriennialLengths.csv")
    write.csv(dat$Ages[dat$Ages$SURVEY=="Tri.Shelf",],file="TriennialAges.csv")
}

#triennial only
tri <- ReadInBiomass.EWC.fn("TriennialCatch.csv",directory=getwd(),species=c(30220,NA),verbose=T)  #all tows

triennial <- ReadInBiomass.EWC.fn(dataFile="TriennialCatch.csv",directory=getwd(), species=30220, verbose=T)  #positive tows only
triennial$catchPerArea <- triennial$kgPerKm2
triennial$year <- triennial$YEAR
triennial <- triennial[triennial$year>1977,]

triennial$latCuts <- cut(triennial$START_LATITUDE,c(0,50))     #c(0,40.5,50)
triennial$depCuts <- cut(triennial$BOTTOM_DEPTH,c(0,183,2000))
table(triennial$latCuts,triennial$depCuts,triennial$YEAR)

table(triennial$SPECIES,useNA="ifany")

hist(triennial$BOTTOM_DEPTH)
triennial[triennial$BOTTOM_DEPTH>400,]

#only 2 fish beyond 400 m
triennial <- triennial[triennial$BOTTOM_DEPTH<=400,]

strat <- read.table("WidowStrataTri.txt",header=T)
names(strat) <- c("name","START_LATITUDE.2","START_LATITUDE.1","BOTTOM_DEPTH.1","BOTTOM_DEPTH.2")
strata.tri <- data.frame(name=strat[,1],area=NA,strat[,c(4,5,2,3)])
strata.tri <- StrataAreas.fn(strata.tri,SA3)

biomass <- DesignBasedEstBiomass.EWC.fn(triennial,strat.vars=c("BOTTOM_DEPTH","START_LATITUDE"),strat.df=strata.tri)
tri.bio <- data.frame(Year=biomass$LNtons$year,Value=biomass$Total$Bhat,seLogB=biomass$LNtons$SElogBhat)
plotBio.fn(tri.bio,pch=16)
#early tri 1980:1992    late tri 1995:2004
abline(v=mean(c(1992,1995)),lwd=3)


####  Lengths
xB <- ReadInBiomass.EWC.fn(dataFile="TriennialCatch.csv",directory=getwd(), species=30220, verbose=T)  #positive tows only
xB <- xB[xB$YEAR%in%1980:2004,]
xB$year <- xB$YEAR
#calculate the density (kg/km^2) using net width and distance fished
xB$catchPerArea <- xB$kgPerKm2

table(xB$YEAR)

strat <- read.table("WidowStrataTri.txt",header=T)
names(strat) <- c("name","START_LATITUDE.2","START_LATITUDE.1","BOTTOM_DEPTH.1","BOTTOM_DEPTH.2")
strata.tri <- data.frame(name=strat[,1],area=NA,strat[,c(4,5,2,3)])
strata.tri <- StrataAreas.fn(strata.tri,SA3)

len <- ReadInLengths.EWC.fn(dataFile="TriennialLengths.csv",directory=getwd(), species=30220, verbose=T)
xL <- len[len$SURVEY=="Tri.Shelf" & len$YEAR%in%1980:2004,]
xL <- xL[!is.na(xL$LENGTH),]
xL$LENGTH <- xL$LENGTH/10
xL$WEIGHT <- xL$SP_TOW_WGHT_KG
xL$NUMBER_FISH <- xL$SP_TOW_NUM
xL$year <- xL$YEAR

range(xL$LENGTH)
table(xL$YEAR) #numbe of samples per year
table(xL$YEAR,!duplicated(paste(xL$CRUISEJOIN,xL$HAULJOIN,xL$HAUL,xL$START_LATITUDE,xL$BOTTOM_DEPTH)))  #The TRUE column is the number of tows sampled
table(xL$YEAR,xL$SEX)
table(xL$LENGTH,xL$SEX)

nSamp <- table(xL$YEAR,!duplicated(paste(xL$CRUISEJOIN,xL$HAULJOIN,xL$HAUL,xL$START_LATITUDE,xL$BOTTOM_DEPTH)))[,"TRUE"]  #The TRUE column is the number of tows sampled
nFish <- table(xL$YEAR)
n <- floor(effN(nSamp,nFish,"survey"))


#the length comps ready for SS3, by sex (use gender=0 for unsexed, gender3 for sexed)
LFs <- SurveyLFs.EWC.fn(xL, xB, strat.vars=c("BOTTOM_DEPTH","START_LATITUDE"), strat.df=strata.tri,
          femaleMale=c(2,1),lgthBins=seq(8,56,2),SS3out=T,gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28,
          fleet=7, ,nSamps=n, season=1)
LFs$F8 <- LFs$F8 + LFs$F.999
LFs$M8 <- LFs$M8 + LFs$M.999
LFs <- LFs[,-which(names(LFs)%in%c("F.999","M.999"))]
write.csv(LFs,"LengthCompsForSS3.Triennial.csv",row.names=F)

LFs <- read.csv("LengthCompsForSS3.Triennial.csv")
par(mfrow=c(2,1))
plotFreqData.fn(LFs,ylim=c(0,65),yaxs="i",ylab="Length",main="Triennial survey")

lfsStrata <- SurveyLFs.EWC.fn(xL, xB, strat.vars=c("BOTTOM_DEPTH","START_LATITUDE"), strat.df=strata.tri,
          femaleMale=c(2,1),lgthBins=seq(8,56,2),SS3out=F,gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28,
          fleet=7, ,nSamps=n, season=1)$L.year.str
triGlmStrata <- read.csv("TriennialIndexByStrata.csv")
triGlmStrata$Strata <- as.character(triGlmStrata$Strata)
#change biomass to total numbers by dividing by mean weight
lenWt <- list(female=data.frame(a=0.00001736,b=2.962),male=data.frame(a=0.00001484,b=3.005))

wtByBiomass.EWC.fn <- function(x,bio,Lengths,lw=NULL) {
    theLs.yr <- as.numeric(as.character(x$LENGTH))
    TotalLjhF <- as.numeric(as.character(x$TotalLjhF))
    TotalLjhM <- as.numeric(as.character(x$TotalLjhM))

    #meanWt
    if(!is.null(lw)) {
        wf <- lw$female$a*theLs.yr^lw$female$b
        wm <- lw$male$a*theLs.yr^lw$male$b
        pf <- TotalLjhF/(sum(TotalLjhF)+sum(TotalLjhM))
        pm <- TotalLjhM/(sum(TotalLjhF)+sum(TotalLjhM))
        mnWt <- sum(wf*pf + wm*pm)
        bio <- bio/mnWt        
    }

    allLs <- Lengths[findInterval(theLs.yr,Lengths,all.inside=T)]    #finds the interval that the length falls in and floors it (so 23.2 would be in 23 if 23 was a level in Lengths, all.inside puts maximum age group into N-1 group, thus I padded with Inf.)
    TotalLjF <- tapply(TotalLjhF,allLs,sum,na.rm=T)
    TotalLjM <- tapply(TotalLjhM,allLs,sum,na.rm=T)
    out <- data.frame(Length=Lengths,TotalLjF=rep(NA,length(Lengths)),TotalLjM=rep(NA,length(Lengths)))
    row.names(out) <- out$Length
    out[names(TotalLjF),"TotalLjF"] <- 100*TotalLjF/(sum(TotalLjF,na.rm=T)+sum(TotalLjM,na.rm=T))
    out[names(TotalLjM),"TotalLjM"] <- 100*TotalLjM/(sum(TotalLjF,na.rm=T)+sum(TotalLjM,na.rm=T))
    out <- out[-nrow(out),]   #remove last row because Inf and always NA due to inside.all=T

    lgths <- as.character(out$Length)
    Ls <- bio*c(out$TotalLjF,out$TotalLjM)
    Ls[is.na(Ls)] <- 0
    Ls <- matrix(Ls,nrow=1,byrow=T,
        dimnames=list(NULL,paste(c(rep("F",length(lgths)),rep("M",length(lgths))),lgths,sep="")))
    # out <- data.frame(year=as.numeric(names(out)),season=season,fleet=fleet,gender=rep(3,length(out)),
    #     partition=partition,Nsamp=nSamps,Ls)
    return(Ls)
}


lenBins <- seq(8,56,2)
lfsBioStr <- data.frame(matrix(nrow=0,ncol=2*length(lenBins)+7,
                dimnames=list(NULL,c("Stratum","year","Season","Fleet","gender","partition","nSamps",
                    paste0("F",lenBins),paste0("M",lenBins)))))
for(yr in names(lfsStrata)) {
    for(str in names(lfsStrata[[yr]])) {

        tmp <- wtByBiomass.EWC.fn(lfsStrata[[yr]][[str]],
                                 triGlmStrata[triGlmStrata$Year==yr & triGlmStrata$Strata==str,"IndexMedian"],
                                 Lengths=c(-999,lenBins,Inf), lw=lenWt)
        lfsBioStr <- rbind(lfsBioStr,data.frame(Stratum=str, year=yr, Season=1, Fleet = 7,
                   gender=3, partition=0, nSamps="ENTER", tmp))
    }
}
xx <- apply(lfsBioStr[,-(1:7)],2,function(y,yr){tapply(y,yr,sum)},yr=lfsBioStr$year)
lfsBio <- data.frame(LFs[c("year","season","fleet","gender","partition","Nsamp")], xx)
lfsBio$F8 <- lfsBio$F8 + lfsBio$F.999
lfsBio$M8 <- lfsBio$M8 + lfsBio$M.999
lfsBio <- lfsBio[,-which(names(lfsBio)%in%c("F.999","M.999"))]
#write.csv(lfsBio,"LengthCompsForSS3.Tri.Bio.csv",row.names=F)  #put lw=NULL above for bio
write.csv(lfsBio,"LengthCompsForSS3.Tri.Num.csv",row.names=F)


#################################################################################
## For 2011 structure
####  Lengths
xB <- ReadInBiomass.EWC.fn(dataFile="TriennialCatch.csv",directory=getwd(), species=30220, verbose=T)  #positive tows only
xB <- xB[xB$YEAR%in%1980:2004,]
xB$year <- xB$YEAR
#calculate the density (kg/km^2) using net width and distance fished
xB$catchPerArea <- xB$kgPerKm2

strat <- read.table("WidowStrataTri.txt",header=T)
names(strat) <- c("name","START_LATITUDE.2","START_LATITUDE.1","BOTTOM_DEPTH.1","BOTTOM_DEPTH.2")
strata.tri <- data.frame(name=strat[,1],area=NA,strat[,c(4,5,2,3)])
strata.tri <- StrataAreas.fn(strata.tri,SA3)

len <- ReadInLengths.EWC.fn(dataFile="TriennialLengths.csv",directory=getwd(), species=30220, verbose=T)
xL <- len[len$SURVEY=="Tri.Shelf" & len$YEAR%in%1980:2004,]
xL <- xL[!is.na(xL$LENGTH),]
xL$LENGTH <- xL$LENGTH/10
xL$WEIGHT <- xL$SP_TOW_WGHT_KG
xL$NUMBER_FISH <- xL$SP_TOW_NUM
xL$year <- xL$YEAR

n <- table(xL$YEAR,!duplicated(paste(xL$CRUISEJOIN,xL$HAULJOIN,xL$HAUL,xL$START_LATITUDE,xL$BOTTOM_DEPTH)))[,"TRUE"]  #The TRUE column is the number of tows sampled


#the length comps ready for SS3, by sex (use gender=0 for unsexed, gender3 for sexed)
LFs <- SurveyLFs.EWC.fn(xL, xB, strat.vars=c("BOTTOM_DEPTH","START_LATITUDE"), strat.df=strata.tri,
          femaleMale=c(2,1),lgthBins=seq(10,64,2),SS3out=T,gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28,
          fleet=7, ,nSamps=n, season=1)
LFs$F10 <- LFs$F10 + LFs$F.999
LFs$M10 <- LFs$M10 + LFs$M.999
LFs <- LFs[,-which(names(LFs)%in%c("F.999","M.999"))]
write.csv(LFs,"LengthCompsForSS3.Triennial.2011structure.csv",row.names=F)

lfsStrata <- SurveyLFs.EWC.fn(xL, xB, strat.vars=c("BOTTOM_DEPTH","START_LATITUDE"), strat.df=strata.tri,
          femaleMale=c(2,1),lgthBins=seq(10,64,2),SS3out=F,gender=3,sexRatioUnsexed=0.5,maxSizeUnsexed=28,
          fleet=7, ,nSamps=n, season=1)$L.year.str
triGlmStrata <- read.csv("TriennialIndexByStrata.csv")
triGlmStrata$Strata <- as.character(triGlmStrata$Strata)


wtByBiomass.EWC.fn <- function(x,bio,Lengths) {
    theLs.yr <- as.numeric(as.character(x$LENGTH))
    TotalLjhF <- as.numeric(as.character(x$TotalLjhF))
    TotalLjhM <- as.numeric(as.character(x$TotalLjhM))
    allLs <- Lengths[findInterval(theLs.yr,Lengths,all.inside=T)]    #finds the interval that the length falls in and floors it (so 23.2 would be in 23 if 23 was a level in Lengths, all.inside puts maximum age group into N-1 group, thus I padded with Inf.)
    TotalLjF <- tapply(TotalLjhF,allLs,sum,na.rm=T)
    TotalLjM <- tapply(TotalLjhM,allLs,sum,na.rm=T)
    out <- data.frame(Length=Lengths,TotalLjF=rep(NA,length(Lengths)),TotalLjM=rep(NA,length(Lengths)))
    row.names(out) <- out$Length
    out[names(TotalLjF),"TotalLjF"] <- 100*TotalLjF/(sum(TotalLjF,na.rm=T)+sum(TotalLjM,na.rm=T))
    out[names(TotalLjM),"TotalLjM"] <- 100*TotalLjM/(sum(TotalLjF,na.rm=T)+sum(TotalLjM,na.rm=T))
    out <- out[-nrow(out),]   #remove last row because Inf and always NA due to inside.all=T

    lgths <- as.character(out$Length)
    Ls <- bio*c(out$TotalLjF,out$TotalLjM)
    Ls[is.na(Ls)] <- 0
    Ls <- matrix(Ls,nrow=1,byrow=T,
        dimnames=list(NULL,paste(c(rep("F",length(lgths)),rep("M",length(lgths))),lgths,sep="")))
    # out <- data.frame(year=as.numeric(names(out)),season=season,fleet=fleet,gender=rep(3,length(out)),
    #     partition=partition,Nsamp=nSamps,Ls)
    return(Ls)
}

lenBins <- seq(10,64,2)
lfsBioStr <- data.frame(matrix(nrow=0,ncol=2*length(lenBins)+7,
                dimnames=list(NULL,c("Stratum","year","Season","Fleet","gender","partition","nSamps",
                    paste0("F",lenBins),paste0("M",lenBins)))))
for(yr in names(lfsStrata)) {
    for(str in names(lfsStrata[[yr]])) {

        tmp <- wtByBiomass.EWC.fn(lfsStrata[[yr]][[str]],
                                 triGlmStrata[triGlmStrata$Year==yr & triGlmStrata$Strata==str,"IndexMedian"],
                                 Lengths=c(-999,lenBins,Inf))
        lfsBioStr <- rbind(lfsBioStr,data.frame(Stratum=str, year=yr, Season=1, Fleet = 7,
                   gender=3, partition=0, nSamps="ENTER", tmp))
    }
}
xx <- apply(lfsBioStr[,-(1:7)],2,function(y,yr){tapply(y,yr,sum)},yr=lfsBioStr$year)
lfsBio <- data.frame(LFs[c("year","season","fleet","gender","partition","Nsamp")], xx)
write.csv(lfsBio,"LengthCompsForSS3.Tri.2011structure.Bio.csv",row.names=F)














####  Ages
age <- ReadInAges.EWC.fn(dataFile="TriennialAges.csv",directory=getwd(), species=30220, verbose=T)
xA <- age[age$SURVEY=="Tri.Shelf" & age$YEAR%in%1980:2004,]
xA <- xA[!is.na(xA$LENGTH),]
xA$LENGTH <- xA$LENGTH/10
xA$WEIGHT <- xA$SP_TOW_WGHT_KG
xA$NUMBER_FISH <- xA$SP_TOW_NUM
xA$year <- xA$YEAR
### 22 ages from a single tow in 1980







DesignBasedEstBiomass.EWC.fn <- function(dat,strat.vars,strat.df)  {
#Calculates design based estimates from survey data
#dat is a dataframe of the data. It must have catchrate called catchPerArea and columns corresponding to strata variable names. It must also have a column called year.
#strat.vars is a vector of the strata variable names (i.e., c("BEST_LATITUDE","BEST_DEPTH"))
#strat.df is a dataframe with the first column the name of the stratum, the second column the area of the stratum, and the remaining columns are the high and low variables defining the strata
#The variables defining the strata must begin with the name in strat.vars and end with ".1" or ".2" (i.e., BEST_DEPTH.1)
#the strata are assumed to be continuous variables, thus have a lower and upper value defining them. The lower value does not necessarily have to be the same as the previous upper value.
#the stat.df dataframe is difficult to build up with more than one variable becuase it turns into a design where you have to define all areas, thus repeat the variables for one (like a design)
#An example strat.df
### data.frame(name=c("shallowS","shallowN","deepS","deepN"),area=c(1200,1400,2000,2500),BEST_DEPTH.1=c(rep(55,2),rep(100,2)),BEST_DEPTH.2=c(rep(100,2),rep(183,2)),INPFC=c(32,42,32,42),BEST_LATITUDE.2=c(42,49,42,49))
### I should think about adding in a routine that automatically puts in the latitude if stratifying by INPFC (or simply finds the strata by INPFC)
#The code below splits out the data by stratum and calculates the average density in each stratum. It then expands up by area to give an estimated weight in each stratum.
#The stratum weights are added together to get a total index for that year
#I calculate the variance given stratified sampling theory
#I work in normal space, then calculate the statistics if B is lognormal
#This is the Mean Ratio Estimate

    if(is.null(dat$catchPerArea)) stop("There must be a column called catchPerArea in the dataframe")
    row.names(strat.df) <- strat.df[,1]     #put in rownmaes to make easier to index later
    numStrata <- nrow(strat.df)

    #first create strata factors
    dat <- data.frame(dat,stratum=StrataFactors.fn(dat,strat.vars,strat.df))        #create a new column for the stratum factor
    dat.yr <- split(dat,dat$year)

    yr.fn <- function(x) {
        x <- split(x,x$stratum)
        namesStrat <- names(x)
        nobs <- unlist(lapply(x,function(x){nrow(x)}))
        if(any(nobs<=1)) {
            cat("*****\nWARNING: At least one stratum in year",x[[1]][1,"year"],"has fewer than one observation.\n*****\n")
        }
        meanCatchRateInStrata <- unlist(lapply(x,function(x){mean(x$catchPerArea)}))
        varCatchRateInStrata <- unlist(lapply(x,function(x){var(x$catchPerArea)}))
        stratStats <- data.frame(name=namesStrat,area=strat.df[namesStrat,2],ntows=nobs,meanCatchRate=meanCatchRateInStrata,varCatchRate=varCatchRateInStrata)
        stratStats$Bhat <- stratStats$area*stratStats$meanCatchRate
        stratStats$varBhat <- stratStats$varCatchRate*(stratStats$area*stratStats$area)/stratStats$ntows
        stratStats
    }
    yearlyStrataEsts <- lapply(dat.yr,yr.fn)
    names(yearlyStrataEsts) <- paste("Year",names(yearlyStrataEsts),sep="")

    yrTotal.fn <- function(x) {
        data.frame(Bhat=sum(x$Bhat),seBhat=sqrt(sum(x$varBhat)),cv=sqrt(sum(x$varBhat))/sum(x$Bhat))
    }
    ests <- as.data.frame(t(as.data.frame(lapply(lapply(yearlyStrataEsts,yrTotal.fn),t))))             #some crazy stuff to put into a dataframe with years as rows
    logVar <- log(ests$cv^2+1)
    ln <- data.frame(year=substring(row.names(ests),5),meanBhat=ests$Bhat/1000,medianBhat=ests$Bhat*exp(-0.5*logVar)/1000,SElogBhat=sqrt(logVar))

    list(Strata=yearlyStrataEsts,Total=ests,LNtons=ln)
}

plotBio.fn <- function(bio,CI=0.95,scalar=1e6,gap=0.03,ylab="Biomass ('000 mt)",xlab="Year",ylim=NULL,segCol="black",...) {
    #Plots the biomass with confidence intervals
    #uses data in the format of SS3, so you can use the GetTotalBiomass.fn or your own data.frame
    #scalar is simply the divisor for the biomass
    #gap is a value that introduces a slight gap between the point estimate and the start of the line for the CI
    # careful because a gap too large will invert the CI, making it look huge. You should know when this happens
    y <- as.numeric(as.character(bio$Value/scalar))
    x <- as.numeric(as.character(bio$Year))
    se <- as.numeric(as.character(bio$seLogB))
    logB <- log(bio$Value)
    ci <- exp(rbind(c(logB+qnorm(1-(1-CI)/2)*se),c(logB-qnorm(1-(1-CI)/2)*se)))/scalar
    if(is.null(ylim)) {
        ylim <- c(0,1.05*max(ci))
    }

    gap <- gap*max(y)
    plot(x,y,,ylab=ylab,xlab=xlab,ylim=ylim,...)
    segments(x,y+gap,x,ci[1,],col=segCol)
    segments(x,y-gap,x,ci[2,],col=segCol)
}
