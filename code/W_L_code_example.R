# estimate weight-length parameters

# read in data (this example uses data in CSV file)
wl <- read.csv("enter location of your csv file here")

# subset based on sex
wlF <- wl[wl$SEX=="f",]
wlM <- wl[wl$SEX=="m",]
wlB <- wl[wl$SEX%in%c("f","m"),]

# apply linear models in log-log space
lmF <- lm(log(WEIGHT_KG)~log(LENGTH_CM),data=wlF)
lmM <- lm(log(WEIGHT_KG)~log(LENGTH_CM),data=wlM)
lmB <- lm(log(WEIGHT_KG)~log(LENGTH_CM),data=wlB)

# get output, including correction for median vs. mean in lognormal distribution
getline <- function(model){
  Amed <- exp(model$coefficients[1])
  B <- model$coefficients[2]
  sdres <- sd(model$residuals)
  Amean <- Amed*exp(0.5*sdres^2)
  return(as.numeric(c(Amed,sdres,Amean,B)))
}

# get parameter values for each model
resultsF <- getline(lmF)
resultsM <- getline(lmM)
resultsB <- getline(lmB)

# make table of all models
results <- data.frame(rbind(resultsF,resultsM,resultsB))
names(results) <- c("median_intercept","SD_resids","A","B")
print(results)

##          median_intercept SD_resids            A        B
## resultsF     2.289612e-06 0.1212388 2.306502e-06 3.152602
## resultsM     3.468706e-06 0.1133655 3.491068e-06 3.034950
## resultsB     2.761746e-06 0.1214178 2.782179e-06 3.098613

# make empty plot
plot(wlB$LENGTH_CM, wlB$WEIGHT_KG,type='n',las=1,
     xlab="Length (cm)", ylab="Weight (kg)",ylim=c(0,1.1*max(wlF$WEIGHT_KG)))
# add line at 0
abline(h=0,col='grey',lty=3)
# add points
points(wlB$WEIGHT_KG~wlB$LENGTH_CM,pch=16,
       col=ifelse(wlB$SEX=="f",rgb(1,0,0,.1),rgb(0,0,1,.1)))

# get sequences of x-values spanning range of data for females and males
xF <- seq(min(wlF$LENGTH_CM),max(wlF$LENGTH_CM),length=200)
xM <- seq(min(wlM$LENGTH_CM),max(wlM$LENGTH_CM),length=200)
# add lines to plot
lines(xF,results$A[1]*xF^results$B[1],col=rgb(1,0,0,0.5),lwd=5)
lines(xM,results$A[2]*xM^results$B[2],col=rgb(0,0,1,0.5),lwd=5)

