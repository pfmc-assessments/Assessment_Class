
setwd("C:\\NOAA2013\\Rougheye\\Data\\AtSea")

foreign <- read.csv("foreignAtSeaCatch.csv")
domestic <- read.csv("doemsticAtSeaCatch.csv")

foreign$catch <- foreign$EXPANDED_WEIGHT_SECTOR_LEVEL..kg.
domestic$catch <- domestic$EXPANDED_WEIGHT_SECTOR_LEVEL..kg.

par(mfrow=c(2,1))
foreign <- foreign[foreign$YEAR>1975,]
boxplot(split(foreign$catch,factor(foreign$YEAR,levels=1976:2012)),log="y",main="Foreign",ylim=c(0.1,5e3))

domestic <- domestic[domestic$YEAR>1990,]
boxplot(split(domestic$catch,factor(domestic$YEAR,levels=1976:2012)),log="y",main="Domestic",ylim=c(0.1,5e3))
