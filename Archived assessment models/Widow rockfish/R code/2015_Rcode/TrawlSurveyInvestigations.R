setwd("C:/NOAA2015/Widow/Data")

#NWFSC slope and combo
load("extractedData/wdowTrawlSurveyCatch.Rdat")
table(singleSp$PROJECT_CYCLE)
table(singleSp$PROJECT_CYCLE,singleSp$HAUL_WT_KG>0)
tapply(singleSp$HAUL_WT_KG,singleSp$PROJECT_CYCLE,sum)

load("extractedData/wdowTrawlSurveyLengths.Rdat")
table(lengths$PROJECT_CYCLE)

#AK Slope and triennial
load("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\Triennial\\Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015.dmp")
dat <- Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015
rm(Tri.Shelf.and.AFSC.Slope.widow.Catch.04.Jun.2015)

table(dat$YEAR,dat$SPECIES_CODE,dat$SURVEY)



load("C:\\NOAA2015\\Widow\\Data\\TrawlSurvey\\Triennial\\AK.Surveys.Bio.widow.04.Jun.2015.dmp")
dat <- AK.Surveys.Bio.widow.04.Jun.2015
rm(AK.Surveys.Bio.widow.04.Jun.2015)

table(dat$Lengths$YEAR,dat$Lengths$SURVEY)
table(dat$Lengths$YEAR,dat$Lengths$HAUL,dat$Lengths$SURVEY)
table(dat$Ages$YEAR,dat$Ages$SURVEY)
