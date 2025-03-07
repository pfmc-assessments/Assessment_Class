library(nwfscSurvey)

setwd("C:/NOAA2015/Widow/Data")

#triennial
x <- ReadInBiomass.EWC.fn(dataFile="TriennialCatch.csv",directory="TrawlSurvey/Triennial", species=30220, verbose=T)

DesignBasedEstBiomass.EWC.fn(x, strat.vars= , strat.df=)