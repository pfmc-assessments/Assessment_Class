commLFs.fn <- function(bds,lw,gear="TWL",state=NULL,catchFile=NULL,
                  maxExpansion=300,verbose=T,normalize=T,
                  doSexRatio=TRUE,sexRatioUnsexed=NA,maxSizeUnsexed=NA,loessSpan=0.3,
                  ageComp=FALSE)  {
    #####################################################################################
    # lw is a data.frame of length-weight params with columns as parameters for different sexes (and unsexed) and rows as WA, OR, CA (1, 2, 3)
    #      # create a predicted fish weight based on state and length
    #   Use unsexed because of uncertainty
    #   use cm and kg for length and weight for fitted regression coefficients!)
    #   lw=data.frame(WA=c(5.6992e-6,3.2523),OR=c(6.8960e-6,3.1873),CA=c(9.4427e-6,3.0855)),
    # sexRatioUnsexed is the ratio to assume for unsexed fish less than maxSizeUnsexed
    # maxExpansion of 300 is about where 20% ofthe samples would be capped, it is also where there is a kink in cumsum
    ######################################################################################

    if(!is.null(gear)) {
        bds <- bds[!is.na(bds$gear)&bds$gear==gear,]             #use only a specific gear code (see BDS_filterData.R for gear code assignments)
    }
    if(!is.null(state)) {
        bds <- bds[bds$state==state,]             #use only a specific state (see BDS_filterData.R for state assignments)
    }

    # set up the expansion based on landed weight to cluster weight
    #if that isn't available, check for species landed weight (RWT_LBS) to species weight
    #if that isn't available, guess
    #use totalWt that was created in BDS_filterData.R, unless missing
    #bds$TOTAL_WGT[is.na(bds$TOTAL_WGT)] <- 0
    bds$usetot_wgt <- NA
    bds$sampleWgt <- NA

    # create a predicted fish weight based on sex and length
    # uses cm and kg ans length and weight units
    # these will be summed to give the sample weight, however
    #don't need to use this because expand cluster weight to total weight because they are cluster sampling a mixture of species
    bds$predwt <- NA
    bds$predwt[bds$state %in% "WA"] <- lw$WA[1]*((bds$FISH_LENGTH[bds$state %in% "WA"]/10)^lw$WA[2]) * 2.20462  #convert kg to pounds
    bds$predwt[bds$state %in% "OR"] <- lw$OR[1]*((bds$FISH_LENGTH[bds$state %in% "OR"]/10)^lw$OR[2]) * 2.20462  #convert kg to pounds
    bds$predwt[bds$state %in% "PW"] <- lw$OR[1]*((bds$FISH_LENGTH[bds$state %in% "PW"]/10)^lw$OR[2]) * 2.20462  #convert kg to pounds
    bds$predwt[bds$state %in% "CA"] <- lw$CA[1]*((bds$FISH_LENGTH[bds$state %in% "CA"]/10)^lw$CA[2]) * 2.20462  #convert kg to pounds

    tmp.split <- split(bds$predwt,bds$SAMPLE_NO)
    predWtSum <- unlist(lapply(tmp.split,sum))
    bds$predWtSum <- predWtSum[match(bds$SAMPLE_NO,names(predWtSum))]
    if(any(is.na(bds$predWtSum))) {
        stop("There are some predicted weights that are NA. This means that there are NA's in lengths, the parameters, or somewhere else.")
    }

    #WA does not have a sample weight
    #I will calcualte sample weight from l-w parameters
    ### This is likely incorrect, but may not cause too much bias. The totalWt is probably the landing of mixed species, and the sample was taken from the mixed species
    ### Theoretically, the best weighting would be totalWt:SampleWt, not totalWT:SpeciesWt.
    ### So, expansion factors are probably bigger than actual when landed wt includes more than jsut Widow (early years)
    #missing fields will assume an expansion of 1
    ind <- bds$SOURCE_AGID %in% c("W")
    bds$usetot_wgt[ind] <- bds$totalWt[ind]

    bds$sampleWgt[ind] <- bds$all_cluster_sum[ind]
    ind <- ind & is.na(bds$sampleWgt)
    bds$sampleWgt[ind] <- bds$CLUSTER_WGT[ind]
    ind <- ind & is.na(bds$sampleWgt)
    bds$sampleWgt[ind] <- bds$predWtSum[ind]

    ind <- (bds$SOURCE_AGID %in% c("W")) & is.na(bds$usetot_wgt)
    bds$usetot_wgt[ind] <- bds$sampleWgt[ind]


    #OR uses expanded weight, lbs
    ind <- !is.na(bds$EXP_WT) & bds$SOURCE_AGID=="O"
    bds$usetot_wgt[ind] <- bds$EXP_WT[ind]

    ind <- bds$SOURCE_AGID=="O"
    bds$sampleWgt[ind] <- bds$all_cluster_sum[ind]
    ind <- ind & is.na(bds$sampleWgt)
    bds$sampleWgt[ind] <- bds$CLUSTER_WGT[ind]
    ind <- ind & is.na(bds$sampleWgt)
    bds$sampleWgt[ind] <- bds$predWtSum[ind]


    #CA uses total weight, lbs and all_cluster_wgt because multiple clusters
    #all cluster wt is the weight of the sample, which might not be all one species.
    #So, use specieswt
    #TOTAL_WGT is the weight of the landing for what was sampled in cluster
    #it looks like the ratio should be either RWT_LBS/allSPECIES_WGT or TOTAL_WGT/all_cluster_sum
    #the difference is when the cluster sample is not specifically widow, but a mix
    ind <- bds$SOURCE_AGID == "C" & !is.na(bds$TOTAL_WGT)
    #first try to fill in TOTAL_WGT/all_cluster_sum
    bds$usetot_wgt[ind] <- bds$TOTAL_WGT[ind]
    bds$sampleWgt[ind] <- bds$all_cluster_sum[ind]
    #then go for RWT_LBS/allSPECIES_WGT for ones that don't have both of above

    ind2 <- ind & is.na(bds$usetot_wgt) | is.na(bds$sampleWgt) & !is.na(bds$RWT_LBS)
    bds$usetot_wgt[ind2] <- bds$RWT_LBS[ind2]
    ind2 <- ind & (is.na(bds$usetot_wgt) | is.na(bds$sampleWgt)) &!is.na(bds$allSpSum)
    bds$sampleWgt[ind] <- bds$allSpSum[ind]

    ind <- is.na(bds$sampleWgt)
    bds$sampleWgt[ind] <- bds$predWtSum[ind]

    #using predWt here with total landing is no different than other states, even though expansion may be larger than it should
    #bds$usetot_wgt[ind & is.na(bds$usetot_wgt)] <- bds$all_cluster_sum[ind & is.na(bds$usetot_wgt)]*0.999*0.45359237


    # possibly use median from the right agency where total landed weight missing (if weights are missing)
    #fill in missing landings weights with predicted weights (could do medians or something also, before expansion) I will set expansion to 0.999 because
    #only 2 samples from OR and 63 from CA (mostly early years) has missing values (WA likely has a lot of missing species weights, but I used predWt)
    #so I'm just going to use 0.999 as the expansion (filled in later)
    # ind <- bds$SOURCE_AGID == "O"
    # bds$usetot_wgt[ind & is.na(bds$usetot_wgt)] <- bds$all_cluster_sum[ind & is.na(bds$usetot_wgt)]*0.999*0.45359237

    bds$usetot_wgt <- as.numeric(bds$usetot_wgt)
    bds$sampleWgt <-  as.numeric(bds$sampleWgt)

    #calculate expansion factors
    # compute an expansion factor for each length observation
    #Owen's recommended power function to account for non-homogeneity within a trip
    bds$expand <- (bds$usetot_wgt/bds$sampleWgt)^0.9
    bds$effN <- bds$sampleWgt^0.8

    bds$expand[is.na(bds$expand)] <- 0.999  #arbitrary, but I can see which ones

    samps <- bds[!duplicated(bds$SAMPLE_NO),]
    row.names(samps) <- samps$SAMPLE_NO
    if(verbose) {
        cat("Table of number of samples with missing weights to expand with by source agency ID. (TRUE or NA means missing)\n")
        print(table(is.na(samps$usetot_wgt),samps$SOURCE_AGID,useNA="ifany"))
        print(table(is.na(samps$usetot_wgt),samps$SAMPLE_YEAR,samps$SOURCE_AGID,useNA="ifany"))
        cat("If there are missing weights, you may want to decide how to fill them in.\n\n")
        print(samps[is.na(samps$usetot_wgt),])
    }

    if(verbose) {
        cat("There are",sum(bds$expand < 0.999),"expansion factors less than 0.999.\n")
        cat("There are",sum(bds$expand > maxExpansion),"expansion factors greater than",maxExpansion,"\n")
        par(mfrow=c(2,1))
        #####expansion factors by state
        boxplot(split(bds$expand,bds$state),xlab="State",ylab="expansion to landing factor")
        boxplot(split(bds$expand,bds$state),xlab="State",ylab="expansion to landing factor",ylim=c(0,maxExpansion))
        windows(height=5,width=6.5)
        plot(cumsum(table(round(samps$expand,0),useNA='ifany'))/nrow(samps))
    }

    bds$expand[bds$expand < 0.999] <- 0.998
    bds$expand[bds$expand > maxExpansion] <- maxExpansion

    #####expansion factors by state
    windows()
    boxplot(split(bds$expand,bds$state),xlab="State",ylab="expansion to landing factor")

    #USE sex ratio to assign unsexed fish to M or F
    if(doSexRatio) {
        tmp.split <- split(bds$SEX,bds$length.cm)
        x.fn <- function(x) {
            xx <- table(x)
            out <- xx["F"]/(xx["F"]+xx["M"])
            names(out) <- NULL
            return(out)
        }
        propF <- unlist(lapply(tmp.split,x.fn))
        nobs <- unlist(lapply(tmp.split,function(x){sum(x=="M" | x=="F")}))
        lens <- as.numeric(names(propF))
        if(verbose) {
            windows(height=5,width=6.5)
            plot(lens,propF,type="l",col="red",xlab="Length (cm)",ylab="Fraction female",ylim=c(0,1),main="Sex Ratio")
            symbols(lens,propF,circles=nobs,inches=0.1,fg="red",bg=rgb(1,0,0,alpha=0.5),add=T)
        }
        propF[lens<=28] <- 0.5
        propF[lens>55] <- 1
        nobs[lens<=28] <- max(nobs)
        lo <- loess(propF~lens,weights=nobs,span=loessSpan)  #determined by eye
        sexRatio <- predict(lo,newdata=data.frame(lens=min(bds$length.cm):max(bds$length.cm)))
        names(sexRatio) <- as.numeric(min(bds$length.cm):max(bds$length.cm))
        sexRatio[sexRatio>1] <- 1
        if(verbose) {
            lines(lens,propF,col="darkgreen")
            lines(as.numeric(names(sexRatio)),sexRatio,col="blue")
            legend("topleft",c("Data","Data w/ Assumed","LOESS"),col=c("red","darkgreen","blue"),lty=1)
        }
        bds$sexRatio <- sexRatio[match(bds$length.cm,as.numeric(names(sexRatio)))]
        bds$FREQ <- as.numeric(bds$FREQ)

        bdsSexRatio <- bds
        tmpF <- tmpM <- bdsSexRatio[bdsSexRatio$SEX=="U",]
        tmpF$SEX <- "F"
        tmpM$SEX <- "M"
        tmpF$FREQ <- tmpF$sexRatio
        tmpM$FREQ <- 1-tmpM$sexRatio
        bdsSexRatio <- rbind(bdsSexRatio,tmpM,tmpF)
        # for(i in 1:nrow(bdsSexRatio)) {
        #     if(bdsSexRatio$SEX[i] == "U") {
        #         tmpM <- tmpF <- bdsSexRatio[i,]
        #         tmpM$SEX <- "M"
        #         tmpF$SEX <- "F"
        #         tmpM$FREQ <- 1-tmpM$sexRatio
        #         tmpF$FREQ <- tmpF$sexRatio
        #         bdsSexRatio <- rbind(bdsSexRatio,tmpM,tmpF)
        #     }
        # }
    }else {
        bdsSexRatio <- bds
    }

    if(ageComp) {
        #switch age into length.cm to trick the code
        #Do this now so that sex ratio is determined using lengths
        origLength <- bds$length.cm
        bds$length.cm <- bds$FISH_AGE_YEARS_FINAL
        origLengthSexRatio <- bdsSexRatio$length.cm
        bdsSexRatio$length.cm <- bdsSexRatio$FISH_AGE_YEARS_FINAL
    }

    # get the number of observed lengths for each sample
    tmp <- table(bds$SAMPLE_NO,!is.na(bds$length.cm))[,"TRUE"]   #numer of lengths in each sample
    samps$obsSampNum <- NA
    samps[names(tmp),"obsSampNum"] <- tmp

    bdsSampsInd <- match(bds$SAMPLE_NO,samps$SAMPLE_NO)

    bds$numLens <- NA
    bds$numLens <- samps[bdsSampsInd,"obsSampNum"]

    samps <- bds[!duplicated(bds$SAMPLE_NO),]
    row.names(samps) <- samps$SAMPLE_NO


    #THINK ABOUT COMBINING SAMPLES WITHIN PORTS OR SIMILAR TO GET LRGER NUMBER OF FISH PER SAMPLE (LOTS OF 1'S)


    # the effective number of fish in each length is the expansion factor (b/c only one per length)
    # except if fish are put in from sex ratio. Then it is expand*FREQ
    # aggregate the effective number to build the composition data for each sex
    tmp <- bdsSexRatio[bdsSexRatio$SEX == "M",]
    if(nrow(tmp)>0) {
        maleLenComps <- aggregate(tmp$expand*tmp$FREQ,
                        list(state=tmp$state, year=tmp$SAMPLE_YEAR, sex=tmp$SEX, length=tmp$length.cm),
                        sum)
    }
    else {
        maleLenComps <- data.frame(state=NA,year=NA,sex=NA,length=NA,x=NA)
    }

    tmp <- bdsSexRatio[bdsSexRatio$SEX == "F",]
    if(nrow(tmp)>0) {
    femaleLenComps <- aggregate(tmp$expand*tmp$FREQ,
                        list(state=tmp$state, year=tmp$SAMPLE_YEAR, sex=tmp$SEX, length=tmp$length.cm),
                        sum)
    }
    else {
        maleLenComps <- data.frame(state=NA,year=NA,sex=NA,length=NA,x=NA)
    }

    tmp <- bds[bds$SEX == "U",]   #FREQ should be all be 1
    if(nrow(tmp)>0) {
    unsexedLenComps <- aggregate(tmp$expand,
                        list(state=tmp$state, year=tmp$SAMPLE_YEAR, sex=tmp$SEX, length=tmp$length.cm),
                        sum)
    }
    else {
        unsexedLenComps <- data.frame(state=NA,year=NA,sex=NA,length=NA,x=NA)
    }

    tmp <- bds
    if(nrow(tmp)>0) {  #use only the originally since M and F added on for sex ratio stuff
    allSexLenComps <- aggregate(tmp$expand,
                        list(state = tmp$state, year = tmp$SAMPLE_YEAR, length = tmp$length.cm),
                        sum)
    allSexLenComps$sex <- "FMU"
    }
    else {
        allSexLenComps <- data.frame(state=NA,year=NA,sex=NA,length=NA,x=NA)
    }

    tmpF <- femaleLenComps
    tmpF$sex <- "F"
    tmpM <- maleLenComps
    tmpM$sex <- "M"
    bothSexLenComps <- rbind(tmpF,tmpM)

    # summary of the number of fish of each sex sampled
    maleLens <- table(bds$SAMPLE_YEAR[bds$SEX == "M"],bds$state[bds$SEX == "M"])
    femaleLens <- table(bds$SAMPLE_YEAR[bds$SEX == "F"],bds$state[bds$SEX == "F"])
    unsexedLens <- table(bds$SAMPLE_YEAR[bds$SEX == "U"],bds$state[bds$SEX == "U"])
    allSexLens <- table(bds$SAMPLE_YEAR,bds$state)

    # calculate number of samples by year and state
    table(samps$SAMPLE_YEAR,samps$state)

    #the total weights are of the entire species mix, sometimes
    #so, instead of calcualting an expansion factor, I'll renormailize the LF and multiply by the catch
    #this is only important if combining states
    maleLenComps <- split(maleLenComps,maleLenComps$state)
    maleLenComps <- lapply(maleLenComps,function(x){split(x,x$year)})
    femaleLenComps <- split(femaleLenComps,femaleLenComps$state)
    femaleLenComps <- lapply(femaleLenComps,function(x){split(x,x$year)})
    unsexedLenComps <- split(unsexedLenComps,unsexedLenComps$state)
    unsexedLenComps <- lapply(unsexedLenComps,function(x){split(x,x$year)})
    allSexLenComps <- split(allSexLenComps,allSexLenComps$state)
    allSexLenComps <- lapply(allSexLenComps,function(x){split(x,x$year)})
    bothSexLenComps <- split(bothSexLenComps,bothSexLenComps$state)
    bothSexLenComps <- lapply(bothSexLenComps,function(x){split(x,x$year)})

    if(is.null(state)) { #only expand up to state if all states are included
        #You want to expand by the TotalCatch/TotalWtSampledFrom
        #But the catch is in the particular species, whereas the samples come from possibly a mix of species
        #So, instead, of calculating an expansion factor, renormailize the LF and multiply by the catch
        #Then add all states together and renormalize to get a weighted proportion by catch in each state
        #Remember to normalize males and females together to retain sex ratio

        #read in total catches by year, area, season for final expansion
        Catch <- read.csv(catchFile,header=T)
        tmp <- unlist(Catch[,-1])
        Catch <- data.frame(year=Catch[,1],state=substring(names(tmp),1,2),catch=tmp)

        normalizeLF.fn <- function(xx,catch) {
            lapply(xx,function(x){
              prop <- x$x/sum(x$x)
              prop <- catch[catch$state==x$state[1] & catch$year==x$year[1],"catch"]*prop
              ret <- data.frame(state=x$state,year=x$year,length=x$length,sex=x$sex,lf=prop)
              if(ageComp) names(ret)[3] <- "age"
              return(ret)
            })
        }
        #Only bothSexLenComps has both sexes normalized together.
        #maleLenComps and femaleLenComps set the sex ratio to 50:50 is used together
        maleLenComps <- lapply(maleLenComps,normalizeLF.fn,catch=Catch)
        femaleLenComps <- lapply(femaleLenComps,normalizeLF.fn,catch=Catch)
        unsexedLenComps <- lapply(unsexedLenComps,normalizeLF.fn,catch=Catch)
        allSexLenComps <- lapply(allSexLenComps,normalizeLF.fn,catch=Catch)
        bothSexLenComps <- lapply(bothSexLenComps,normalizeLF.fn,catch=Catch)

        cat("Comps weighted by state specific catches.\nThey can simply be added together for a coastwide comp\n")
    }

    if(!is.null(state)) { #only expand up to state if all states are included, but format all lists here for a single state
        if(normalize) {
            normalizeLF.fn <- function(xx) {
              lapply(xx,function(x){
                prop <- x$x/sum(x$x)
                ret <- data.frame(state=x$state,year=x$year,sex=x$sex,length=x$length,lf=prop)
                if(ageComp) names(ret)[4] <- "age"
                return(ret)
              })
            }
        }

        if(!normalize) {
            normalizeLF.fn <- function(xx) {
              lapply(xx,function(x){
                prop <- x$x
                ret <- data.frame(state=x$state,year=x$year,sex=x$sex,length=x$length,lf=prop)
                if(ageComp) names(ret)[4] <- "age"
                return(ret)
              })
            }
        }
        maleLenComps <- lapply(maleLenComps,normalizeLF.fn)
        femaleLenComps <- lapply(femaleLenComps,normalizeLF.fn)
        unsexedLenComps <- lapply(unsexedLenComps,normalizeLF.fn)
        allSexLenComps <- lapply(allSexLenComps,normalizeLF.fn)
        bothSexLenComps <- lapply(bothSexLenComps,normalizeLF.fn)
    }

    return(list(female=femaleLenComps,male=maleLenComps,both=bothSexLenComps,unsexed=unsexedLenComps,all=allSexLenComps,bds=bds[,c("SAMPLE_NO","SAMPLE_YEAR","state","expand","effN")]))
}


lfsForSS3_gender3.fn <- function(LFs,years,lens,season=1,fleet="Enter",
                         gender=3,partition=2,normalize=T,ageComp=FALSE) {
    #LFs is a list of years and from above would be specific to a state
    #This is for both sexes where the sex ratio is retained (gender=3)

    #it looks for all states and adds them together for a specific year
    #There are two levels to the list. This first is state, the second in year
    out <- matrix(0,nrow=length(years),ncol=2*length(lens)+6,dimnames=list(years,c("year","Season","Fleet","gender","partition","nSamps",paste("F",lens,sep=""),paste("M",lens,sep=""))))
    out[,"year"] <- as.numeric(years)
    out[,"Season"] <- season
    out[,"Fleet"] <- fleet
    out[,"gender"] <- gender
    out[,"partition"] <- partition
    out[,"nSamps"] <- NA
    for(yr in as.character(years)) {
        #cat("Year:",yr,"\n")
        for(state in names(LFs)) {
            #cat("State:",state,"\n")
            tmp <- LFs[[state]][[yr]]
            if(!is.null(tmp)) {
                if(ageComp) names(tmp)[which(names(tmp)=="age")] <- "length"
                tmp <- split(tmp,tmp$sex)
                sex <- "F"
                if(!is.null(tmp[[sex]])) {
                    tmp[[sex]]$length[tmp[[sex]]$length > max(lens)] <- max(lens)
                    tmp[[sex]]$length[tmp[[sex]]$length < min(lens)] <- min(lens)
                    lgths <- lens[findInterval(tmp[[sex]]$length,lens)]
                    lgths <- paste(sex,lgths,sep="")
                    x <- tapply(tmp[[sex]]$lf,lgths,sum)
                    out[yr,names(x)] <- out[yr,names(x)] + x
                }

                sex <- "M"
                if(!is.null(tmp[[sex]])) {
                    tmp[[sex]]$length[tmp[[sex]]$length > max(lens)] <- max(lens)
                    tmp[[sex]]$length[tmp[[sex]]$length < min(lens)] <- min(lens)
                    lgths <- lens[findInterval(tmp[[sex]]$length,lens)]
                    lgths <- paste(sex,lgths,sep="")
                    x <- tapply(tmp[[sex]]$lf,lgths,sum)
                    out[yr,names(x)] <- out[yr,names(x)] + x
                }
            }
        }
    }
    if(normalize) {
        out[,-(1:6)] <- t(round(100*apply(out[,-(1:6)],1,function(x){x/sum(x)}),2))
    }
    out[,-(1:6)][is.na(out[,-(1:6)])] <- 0
    return(out)
}




































lfsForSS3_combinedsex.fn <- function(LFs,years,season=1,fleet,gender=0,partition=2,lens,combineStates=T) {
    #combineStates means to add states together (expanded by catch)
    #this is not generalized and assumes you want to combine the three states
    #if(combineStates) {
    #    yrs <- sort(unique(lapply(LFs,function(x){as.numeric(names(x))})))
    #    comb <- list()
    #    for(i in 1:length(yrs)) {
    #        tmp <- merge(LFs$CA,LFs$OR,by="length")
    #        comb[[i]] <-
    #    }

    #LFs is a list of years and from above would be specific to a state
    #years <- sort(unique(c(names(LFs))))
    out <- matrix(0,nrow=length(years),ncol=2*length(lens)+6,dimnames=list(years,c("year","Season","Fleet","gender","partition","nSamps",paste("F",lens,sep=""),paste("M",lens,sep=""))))
    out[,"year"] <- as.numeric(years)
    out[,"Season"] <- season
    out[,"Fleet"] <- fleet
    out[,"gender"] <- gender
    out[,"partition"] <- partition
    out[,"nSamps"] <- NA
    for(i in names(LFs)) {
        tmp <- LFs[[i]]
        lfsplus <- sum(tmp$lf[tmp$length>=max(lens)])
        lfsminus <- sum(tmp$lf[tmp$length<=min(lens)])
        tmp[tmp$length==max(lens),"lf"] <- lfsplus
        tmp[tmp$length==min(lens),"lf"] <- lfsminus
        tmp <- tmp[tmp$length>=min(lens)&tmp$length<=max(lens),]
        lgths <- lens[findInterval(tmp$length,lens)]
        lgths <- paste("F",lgths,sep="")
        x <- tapply(tmp$lf,lgths,sum)
        lgths <- names(x)
        out[i,lgths] <- x
        lgths <- lens[findInterval(tmp$length,lens)]
        lgths <- paste("M",lgths,sep="")
        x <- tapply(tmp$lf,lgths,sum)
        lgths <- names(x)
        out[i,lgths] <- x
    }
    return(out)
}














lfsForSS3.fn <- function(femLFs,maleLFs,season,fleet,gender,partition,lens) {
    #LFs is a list of years and from above would be specific to a state
    #This is for individual sexes where the sex ratio is 1:1
    years <- sort(unique(c(names(femLFs),names(maleLFs))))
    out <- matrix(0,nrow=length(years),ncol=2*length(lens)+6,dimnames=list(years,c("year","Season","Fleet","gender","partition","nSamps",paste("F",lens,sep=""),paste("M",lens,sep=""))))
    out[,"year"] <- as.numeric(years)
    out[,"Season"] <- season
    out[,"Fleet"] <- fleet
    out[,"gender"] <- gender
    out[,"partition"] <- partition
    out[,"nSamps"] <- NA
    for(i in names(femLFs)) {
        tmp <- femLFs[[i]]
        lfsplus <- sum(tmp$lf[tmp$length>=max(lens)])
        lfsminus <- sum(tmp$lf[tmp$length<=min(lens)])
        tmp[tmp$length==max(lens),"lf"] <- lfsplus
        tmp[tmp$length==min(lens),"lf"] <- lfsminus
        tmp <- tmp[tmp$length>=min(lens)&tmp$length<=max(lens),]
        lgths <- paste("F",tmp$length,sep="")
        out[i,lgths] <- tmp$lf
    }
    for(i in names(maleLFs)) {
        tmp <- maleLFs[[i]]
        lfsplus <- sum(tmp$lf[tmp$length>=max(lens)])
        lfsminus <- sum(tmp$lf[tmp$length<=min(lens)])
        tmp[tmp$length==max(lens),"lf"] <- lfsplus
        tmp[tmp$length==min(lens),"lf"] <- lfsminus
        tmp <- tmp[tmp$length>=min(lens)&tmp$length<=max(lens),]
        lgths <- paste("M",tmp$length,sep="")
        out[i,lgths] <- tmp$lf
    }
    return(out)
}
lfsForSS3_combinedsex.fn <- function(LFs,years,season=1,fleet,gender=0,partition=2,lens,combineStates=T) {
    #combineStates means to add states together (expanded by catch)
    #this is not generalized and assumes you want to combine the three states
    #if(combineStates) {
    #    yrs <- sort(unique(lapply(LFs,function(x){as.numeric(names(x))})))
    #    comb <- list()
    #    for(i in 1:length(yrs)) {
    #        tmp <- merge(LFs$CA,LFs$OR,by="length")
    #        comb[[i]] <-
    #    }

    #LFs is a list of years and from above would be specific to a state
    #years <- sort(unique(c(names(LFs))))
    out <- matrix(0,nrow=length(years),ncol=2*length(lens)+6,dimnames=list(years,c("year","Season","Fleet","gender","partition","nSamps",paste("F",lens,sep=""),paste("M",lens,sep=""))))
    out[,"year"] <- as.numeric(years)
    out[,"Season"] <- season
    out[,"Fleet"] <- fleet
    out[,"gender"] <- gender
    out[,"partition"] <- partition
    out[,"nSamps"] <- NA
    for(i in names(LFs)) {
        tmp <- LFs[[i]]
        lfsplus <- sum(tmp$lf[tmp$length>=max(lens)])
        lfsminus <- sum(tmp$lf[tmp$length<=min(lens)])
        tmp[tmp$length==max(lens),"lf"] <- lfsplus
        tmp[tmp$length==min(lens),"lf"] <- lfsminus
        tmp <- tmp[tmp$length>=min(lens)&tmp$length<=max(lens),]
        lgths <- lens[findInterval(tmp$length,lens)]
        lgths <- paste("F",lgths,sep="")
        x <- tapply(tmp$lf,lgths,sum)
        lgths <- names(x)
        out[i,lgths] <- x
        lgths <- lens[findInterval(tmp$length,lens)]
        lgths <- paste("M",lgths,sep="")
        x <- tapply(tmp$lf,lgths,sum)
        lgths <- names(x)
        out[i,lgths] <- x
    }
    return(out)
}
