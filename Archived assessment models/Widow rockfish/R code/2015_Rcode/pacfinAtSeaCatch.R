#Function to extract PacFIN at-sea catches from NPAC4900_SPCOMP
#Needs the RODBC library
#Start R in 32-bit mode for drivers to work
#set working directory to have subdirectories:
###sql: with the sql script called Pacfin.AtSeaCatch.query
###Catches: To save the catch results in
#This assigns the fleet type and summarizes over fleet and year and period
#Writes out two files, one on catches by each fleet by year
#The second file is USshoreside only summarized by period and year

setwd("C:/NOAA2015/Widow/Data")
library(RODBC)
library(ExtractR)
#source("Rcode/functions/Functions.R")

#Get PacFIN catches and save results
pAtSeaCatch <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/Pacfin.AtSeaCatch.sql",db="PACFIN",uid="hicksa",sp="WDOW", start="1980", end="2014")
save(pAtSeaCatch,file="extractedData/PacfinAtSeaCatch.Rdat")
    file.copy("extractedData/PacfinAtSeaCatch.Rdat",paste("extractedData/PacfinAtSeaCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

# NPAC4900SpComp <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/NPAC4900_SPCOMP.sql",db="PACFIN",uid="hicksa",sp="WDOW", start="1980", end="2014")
# save(pAtSeaCatch,file="extractedData/PacfinAtSeaCatch.Rdat")
#     file.copy("extractedData/PacfinAtSeaCatch.Rdat",paste("extractedData/PacfinAtSeaCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)

# NPAC4900V <- queryDB(queryFilename="C:/NOAA2015/Widow/Data/SQL/NPAC4900_V.sql",db="PACFIN",uid="hicksa",sp="WDOW", start="1980", end="2014")
# save(pAtSeaCatch,file="extractedData/PacfinAtSeaCatch.Rdat")
#     file.copy("extractedData/PacfinAtSeaCatch.Rdat",paste("extractedData/PacfinAtSeaCatch_",format(Sys.time(),"%Y.%m.%d"),".Rdat",sep=""),overwrite=T)



setwd("C:/NOAA2015/Widow/Data")
load("extractedData/PacfinAtSeaCatch.Rdat")

table(pAtSeaCatch$YEAR)

wdowTons <- tapply(pAtSeaCatch$TOTAL_WEIGHT,pAtSeaCatch$YEAR,sum)
cbind(wdowTons)
tapply(pAtSeaCatch$TOTAL_WEIGHT,list(pAtSeaCatch$MONTH,pAtSeaCatch$YEAR),sum)


pAtSeaCatch <- queryDB(queryFilename="Pacfin.AtSeaCatch",db="PACFIN",uid="hicksa",sp="PWHT", start="1980", end="2012", querydir="sql/")
hakeAtSeaTons <- c(4388.776, 204687.636, 151210.546,  93510.597, 166248.170, 100383.341, 125944.612, 145840.998, 144546.383, 140913.215, 120712.587,
                100329.690, 84750.168, 86611.965,116972.202, 151067.740, 139786.631, 126240.612, 180681.188,  72349.572, 106313.917, 128073.417, 93736.273)
names(hakeAtSeaTons) <- as.character(1990:2012)
#From 2013 assessment
#c(4.54,205.82,154.74,98.04,179.87,102.31,128.11,146.05,145.16,141.02,120.92,100.53,84.75,86.61,117.07,151.07,139.79,126.24,180.64,72.35,106.31,128.07,93.78)
plot(hakeAtSeaTons,reyeTons)
