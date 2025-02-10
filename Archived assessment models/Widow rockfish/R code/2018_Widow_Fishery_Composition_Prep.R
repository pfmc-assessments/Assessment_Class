# Navigate to directory where data is
old_Wd <- getwd()
setwd("..")
setwd("2019 update/Data")

# Read data
length_data <- read.csv("WDOW.PacFIN.LengthComps2015-2018.csv", header = TRUE, skip = 1)
age_data <- read.csv("WDOW.PacFIN.AgeComps2015-2018.csv", skip = 1)

# Remove length 0
length_data <- length_data[,-which(colnames(length_data) %in% c("L0", "L0.1"))]

# Calculate samples size
n_samp_fish <- function(nfish, nsample){
  neff = ifelse(nfish/nsample < 44, 
                nsample + 0.138 * nfish,
                7.06 * nsample)
  return(neff)
}

# Set up length data for saving
neff <- n_samp_fish(nfish = length_data$Nsamps, nsample = length_data$Ntrips)
ss_length_data <- data.frame(yr = length_data$fishyr,  month = length_data$season, fleet = length_data$fleet, sex = length_data$gender, part = length_data$partition, Nsamp = round(neff, 0))
ss_length_data <- cbind(ss_length_data, round(length_data[,9:ncol(length_data)], 2))
write.table(ss_length_data, "WIDOW_PACFIN_LENGTH_COMPS_2015_2018_FORMATTED.csv", sep = " ")

# Set up age data for saving
neff <- n_samp_fish(nfish = age_data$Nsamps, nsample = age_data$Ntrips)
ss_age <- data.frame(yr = age_data$fishyr,  month = age_data$season, fleet = age_data$fleet, sex = age_data$gender, part = age_data$partition, Nsamp = round(neff, 0))
ss_age <- cbind(ss_age, round(age_data[,9:ncol(age_data)], 2))
write.table(ss_age, "WIDOW_PACFIN_AGE_COMPS_2015_2018_FORMATTED.csv", sep = " ")

# setwd to old
setwd(old_Wd)


# May 2019 - Grant Adams
library(r4ss)
library(data.table)

####################################################
## GET BASE CASE VALUES FOR ASSESSMENT
####################################################
# directories where models were run need to be defined
oldwd <- getwd()
setwd("../") # Move back one
dir = 'Model Runs/2019_Base_case' # Updated reweighting of the model


mod1 <- SS_output(dir = dir)
data1 <- SS_readdat_3.30(file = paste0(dir, "/2019widow.dat"))

# Plot age comp
for(i in 1:length(unique(data1$agecomp$FltSvy))){
  fltsrv <- sort(unique(data1$agecomp$FltSvy))[i]
  if(fltsrv > 0 ){
    dat_sub <- data1$agecomp[which(data1$agecomp$FltSvy == fltsrv),]
    dat_sub <- dat_sub[which(dat_sub$Yr > 0),]
    dat_sub <- dat_sub[which(dat_sub$Part %in% c(0, 2)),]
    
    ylab="Ages"; xlab="Year"
    
    nages <- 40
    x <- as.numeric(as.character(dat_sub$Yr))
    dat <- dat_sub[,10:ncol(dat_sub)]
    inch <- 0.15
    y <- as.numeric(substring(names(dat),2))
    y <- y[1:nages]
    xlim <- range(x)
    
    par(mfrow=c(2,1))
    
    name <- "Female"
    plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = c(0, nages), main = paste(name, fltsrv))
    for(j in 1:nrow(dat_sub)){
      if(dat_sub$Yr[j] > 0){
        if(dat_sub$Lbin_lo[j] < 0){
          if(dat_sub$Gender[j] %in% c(1, 3)){
            symbols(rep(dat_sub$Yr[j], length(0:nages)),0:nages,circles=dat[j,grep("f", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
    }
    
    name <- "Male"
    plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = c(0, nages), main = paste(name, fltsrv))
    for(j in 1:nrow(dat_sub)){
      if(dat_sub$Yr[j] > 0){
        if(dat_sub$Lbin_lo[j] < 0){
          if(dat_sub$Gender[j] >= 2){
            symbols(rep(dat_sub$Yr[j], length(0:nages)),0:nages,circles=dat[j,grep("m", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
    }
  }
}

# Total age historgrams
dat_sub <- data1$agecomp[which(data1$agecomp$FltSvy > 0),]
dat_sub <- dat_sub[which(dat_sub$Yr > 0),]
dat <- dat_sub[,10:ncol(dat_sub)]
fdat <- dat[,grep("f", colnames(dat))]
mdat <- dat[,grep("m", colnames(dat))]

par(mfrow = c(2, 1))
# females
plot( y = colSums(fdat), x =  0:40, type = "l", xlab = "Age", ylab = "Number observed", main = "Females")
# males
plot( y = colSums(mdat), x =  0:40, type = "l", xlab = "Age", ylab = "Number observed", main = "Males")

# Plot length comp
for(i in 1:length(unique(data1$lencomp$FltSvy))){
  fltsrv <- sort(unique(data1$lencomp$FltSvy))[i]
  if(fltsrv > 0 ){
    dat_sub <- data1$lencomp[which(data1$lencomp$FltSvy == fltsrv),]
    dat_sub <- dat_sub[which(dat_sub$Yr > 0),]
    dat_sub <- dat_sub[which(dat_sub$Part %in% c(0, 2)),]
    ylab="Length (cm)"; xlab="Year"
    if(nrow(dat_sub) > 0){
      
      nages <- 25
      x <- as.numeric(as.character(dat_sub$Yr))
      dat <- dat_sub[,7:ncol(dat_sub)]
      inch <- 0.1
      y <- as.numeric(substring(names(dat),2))
      y <- y[1:nages]
      xlim <- range(x)
      
      par(mfrow=c(2,1))
      
      name <- "Female"
      plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = range(y), main = paste(name, "catch", fltsrv))
      for(j in 1:nrow(dat_sub)){
        if(dat_sub$Yr[j] > 0){
          if(dat_sub$Gender[j] %in% c(1, 3)){
            symbols(rep(dat_sub$Yr[j], nages),y,circles=dat[j,grep("f", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
      
      name <- "Male"
      plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = range(y), main = paste(name, "catch", fltsrv))
      for(j in 1:nrow(dat_sub)){
        if(dat_sub$Yr[j] > 0){
          if(dat_sub$Gender[j] == 2){
            symbols(rep(dat_sub$Yr[j], nages),y,circles=dat[j,grep("m", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
    }
  }
}




# Plot discard length comp
for(i in 1:length(unique(data1$lencomp$FltSvy))){
  fltsrv <- sort(unique(data1$lencomp$FltSvy))[i]
  if(fltsrv > 0 ){
    dat_sub <- data1$lencomp[which(data1$lencomp$FltSvy == fltsrv),]
    dat_sub <- dat_sub[which(dat_sub$Yr > 0),]
    dat_sub <- dat_sub[which(dat_sub$Part %in% c(1)),]
    ylab="Length (cm)"; xlab="Year"
    if(nrow(dat_sub) > 0){
      
      nages <- 25
      x <- as.numeric(as.character(dat_sub$Yr))
      dat <- dat_sub[,7:ncol(dat_sub)]
      inch <- 0.1
      y <- as.numeric(substring(names(dat),2))
      y <- y[1:nages]
      xlim <- range(x)
      
      par(mfrow=c(2,1))
      
      name <- "Female"
      plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = range(y), main = paste(name, "discard", fltsrv))
      for(j in 1:nrow(dat_sub)){
        if(dat_sub$Yr[j] > 0){
          if(dat_sub$Gender[j] %in% c(1, 3)){
            symbols(rep(dat_sub$Yr[j], nages),y,circles=dat[j,grep("f", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
      
      name <- "Male"
      plot(NA, NA,xlab=xlab,ylab=ylab,xlim=xlim, ylim = range(y), main = paste(name, "discard", fltsrv))
      for(j in 1:nrow(dat_sub)){
        if(dat_sub$Yr[j] > 0){
          if(dat_sub$Gender[j] == 2){
            symbols(rep(dat_sub$Yr[j], nages),y,circles=dat[j,grep("m", colnames(dat))],inches=inch, add = TRUE)
          }
        }
      }
    }
  }
}
