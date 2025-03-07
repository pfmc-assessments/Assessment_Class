propMat <- function(dat) {
	x <- table(round(dat$Length,0),dat$Mature)
	x <- cbind(as.numeric(rownames(x)),x[,"1"]/(x[,"0"]+x[,"1"]))
	colnames(x) <- c("Length","PrMat")
	return(x)
}

maturity.fn <- function(dat,lens=dat$Length,matur=dat$I.M) {
    if(is.factor(lens)) {
        tmp <- strsplit(as.character(lens),",")
        tmp1 <- unlist(lapply(tmp,function(x){as.numeric(substring(x[1],2))}))
        tmp2 <- unlist(lapply(tmp,function(x){as.numeric(substring(x[2],1,nchar(x[2])-1))}))
        lens <- (tmp1+tmp2)/2
    }
    mats <- table(lens,matur)
    tots <- mats[,1]+mats[,2]
    prop <- mats[,2]/tots
    len <- as.numeric(names(prop))
    invisible(data.frame(len,tots,prop))
}



fitDataLogit <- function(pars,dat) {
    #fits to lengths and ages and does not asymtote to 1
    #inspired by code written by Jim Thorson
    tmp <- (pars[1] + pars[2]*dat$Length + pars[3]*dat$Age)  #in Real space
    pred <- inv.logit(tmp) * pars[4]
    loglike <- sum(dat$Mature*log(pred) + (1-dat$Mature)*log(1-pred),na.rm=T)   #Bernoulli log-likelihood
    return(-2*loglike)
}

fitDataLogitLen <- function(pars,dat) {
    #fits to lengths and ages and does not asymtote to 1
    #inspired by code written by Jim Thorson
    tmp <- (pars[1] + pars[2]*dat$Length)  #in Real space
    pred <- inv.logit(tmp) * pars[3]
    loglike <- sum(dat$Mature*log(pred) + (1-dat$Mature)*log(1-pred),na.rm=T)   #Bernoulli log-likelihood
    return(-2*loglike)
}

predLogit <- function(pars,lengths,ages) {
    pred <- matrix(NA,nrow=length(lengths),ncol=length(ages),dimnames=list(as.character(lengths),as.character(ages)))
    for(j in ages) {
        tmp <- (pars[1] + pars[2]*lengths + pars[3]*j)  #in Real space
        pred[,j] <- inv.logit(tmp) * pars[4]
    }
    return(pred)
}

predLogitLen <- function(pars,lengths) {
    pred <- rep(NA,length(lengths))
    pred <- (pars[1] + pars[2]*lengths)  #in Real space
    pred <- inv.logit(pred) * pars[3]
    return(pred)
}

plotPredMaturityLength <- function(lens,dat,start,low=c(-Inf,-Inf,1),upp=c(Inf,Inf,1),fitFnc=fitDataLogitLen,doL50=F,...) {
	if(length(start)>3) stop("This is a three parameter model\n")
	if(start[3]>1 | low[3]>1 | upp[3]>1) stop("The asymptote (3rd parameter) should be 1 or less\n")
	if(low[3] > start[3]) {
		cat("WARNING: Three parameters are supplied but the asymptote is not allowed to be less than 1 in the low argument\n")
	}
	out <- nlminb(fitFnc, start=start,dat=dat,lower=low,upper=upp)
	outPred <- predLogitLen(out$par,lengths=lens)
	propMatureAtLen <- outPred   #proportion mature at alength with asymptte estimated
	lines(lens,propMatureAtLen,...)

    if(doL50) {
        #find l50 and l95
        prop <- 0.5
        tmp <- lens[order(abs(outPred-prop))][1]
        lens2 <- seq(tmp-1,tmp+1,0.01)
        outPred2 <- predLogitLen(out$par,lengths=lens2)
        l50 <- lens2[order(abs(outPred2-prop))][1]
        prop <- 0.80
        tmp <- lens[order(abs(outPred-prop))][1]
        lens2 <- seq(tmp-1,tmp+1,0.01)
        outPred2 <- predLogitLen(out$par,lengths=lens2)
        l80 <- lens2[order(abs(outPred2-prop))][1]
        legend("bottomright",
               c(paste("l50",l50,sep="="),paste("l80",l80,sep="="),paste("Asymptote",round(out$par[3],2),sep="=")),
               pch=NA, bty="n")
    }

	invisible(out)
}


plotMaturity.fn <- function(x, theFactor, cols, factorLab=sort(as.character(unique(x[,theFactor]))), lens=sort(unique(x$Length)),
                            alpha=0.6, inch=0.1, ht=0.1, labHadj=0.3, legYinter=1, legx="right", legy=NULL, legTitle=NULL, xlabel="Length (cm)") {
    #x is the raw data, a single data frame
    #theFactor is the column name of the factor to split on (plot separately)
    #cols is a vector of colors that is length of unique values in
    #alpha is the transperancy
    #inch scales the circle size relative to number of samples
    #ht is the interval for the rug plots

    dat <- split(x,x[,theFactor])
    PrMat <- lapply(dat,maturity.fn)
    theMax <- max(unlist(lapply(PrMat,function(x){x$tots})))
    numMats <- length(PrMat)  #number of different items to plot

    cat("Check that the labels match the factors\n")
    print(rbind(names(dat),factorLab,cols))

    plot(1,1,xlim=range(x$Length,na.rm=T),ylim=c(-ht*numMats,1.05),type="n",xlab=xlabel,ylab="",yaxt="n",yaxs="i")
    axis(2,at=seq(0.1,1,0.1),las=1)
    mtext("Proportion Mature", side=2, at=0.5, adj=0.5, line=3, outer=F)
    axis(2,at=seq(-1*ht,-numMats*ht,-1*ht)+0.5*ht,label=factorLab,las=1,cex.axis=0.7,tick=F,hadj=labHadj)
    axis(3)
    axis(4,at=c(-numMats*ht,-(numMats-1)*ht),label=c(0,theMax),cex.axis=0.7,las=1,hadj=0.5)
    mtext("# of samples", side=4, line=1.3, at=mean(seq(-1*ht,-numMats*ht,-1*ht)+0.5*ht), adj=0.5, cex=0.7)
    #axis(4,at=mean(seq(-1*ht,-numMats*ht,-1*ht)+0.5*ht),label="# of samples",las=0,cex.axis=0.8,tick=F,hadj=0
    abline(h=0,lwd=3)
    abline(h=seq(-1*ht,-numMats*ht,-1*ht))

    for(i in 1:numMats) {
        out <- PrMat[[i]]
        theCol <- col2rgb(cols[i])[,1]/255
        symbols(c(-1,out$len), c(-1,out$prop), circles=sqrt(c(theMax,out$tots)),
                 inches=inch,fg="black",bg=rgb(theCol["red"],theCol["green"],theCol["blue"],alpha),add=T)
#        segments(out$len,-i*ht,out$len,out$tots/(max(out$tots)/ht)-i*ht,col=cols[i],lwd=3)  #this will plot each factor to have its own max (and at least 1 segment reaching the top)
        segments(out$len,-i*ht,out$len,out$tots/(theMax/ht)-i*ht,col=cols[i],lwd=3)
        tmp <- plotPredMaturityLength(lens, dat[[i]],
                           start=c(-6.2,0.179,1), low=c(-Inf,-Inf,1), upp=c(Inf,Inf,1),
                           col=cols[i], lwd=3)
    }
    legend(x=legx,y=legy,factorLab,col=cols,lty=1,lwd=3,y.intersp = legYinter,cex=0.8,title=legTitle)
    invisible(PrMat)
}



