
###Installation

install.packages("remotes")
remotes::install_github("pfmc-assessments/nwfscSurvey")


#####List of functions

ls("package:nwfscSurvey")


#######Pull the data

catch = PullCatch.fn(Name = "rex sole", 
                     SurveyName = "NWFSC.Combo")
                     
bio   = PullBio.fn(Name = "rex sole", 
                   SurveyName = "NWFSC.Combo")
 
                   
head(catch)
head(bio)


####Save the data to explore the data
write.csv(catch,"catch.csv")
write.csv(bio,"bio.csv")

####Visualize the data

plot_cpue(
  dir = getwd(), 
  catch = catch)

plot_bio_patterns(
  dir = getwd(), 
  bio = bio, 
  col_name = "Length_cm")

wh_plot_proportion(
  data_catch = catch,
  data_bio = bio
)

######Set spatial strata

  
  strata = CreateStrataDF.fn(names=c("shallow_CA", "middle_CA", "deepCA","shallow_OR", "middle_OR", "deep_OR", "shallow_WA", "middle_WA",  "deep_WA"), 
                           depths.shallow = c(55,    183,  549,  55,   183,  549,		55, 	183,   549),
                           depths.deep    = c(183,   549,  1280, 183,  549,  1280,	183,	549,  1280),
                           lats.south     = c(32,    32,   32,   42,   42,   42,   	46, 	46,    46),
                           lats.north     = c(42,    42,   42,   46,   46,   46,   	49, 	49,    49))
                           
                           
###calculate design-based index

biomass = Biomass.fn(dir = getwd(), 
                     dat = catch,  
                     strat.df = strata)
                     
PlotBio.fn(
  dir = getwd(), 
  dat = biomass)
  
PlotBioStrata.fn(
  dir = getwd(), 
  dat = biomass)
                       
#####Length comps                  

###########Changed from GetN to GetN.fn
n <- GetN.fn(dir = getwd(), 
         dat = bio, 
         type = "length", 
         
         species = "flatfish")
         
len_bins <- seq(2, 60, 2)

Length_Freq <- SurveyLFs.fn(dir = getwd(), 
                        datL =  bio, 
                        datTows = catch,
                       strat.df = strata,
                        lgthBins = len_bins)


PlotFreqData.fn(dir = getwd(), 
                dat = Length_Freq)

plot_comps(dir = getwd(), data = Length_Freq)




### Age comps
n <- GetN.fn(dir = getwd(), 
         dat = bio,
         type = "age", 
         species = "flatfish")
         
age_bins <- 1:22

Ages <- SurveyAFs.fn(dir = getwd(), 
                     datA = bio,
                     datTows = catch, 
                     strat.df = strata,
                     ageBins = age_bins,
                     nSamps = n)
                     
  PlotFreqData.fn(
  dir = getwd(), 
  dat = Ages)

plot_comps(dir = getwd(), 
  data = Ages)
  
#CAAL

caal <- SurveyAgeAtLen.fn(dir = getwd(), 
                          datAL = bio, 
                          datTows = catch,
                          strat.df = strata,
                          lgthBins = len_bins, 
                          ageBins = age_bins)
                          
                                  
##MAP                               
PlotMap.fn(dir = getwd(),
  dat = catch)
  
  
#Additional plots  
PlotVarLengthAtAge.fn(dat=bio,(dir = getwd()))


PlotSexRatio.fn(dat=bio,(dir = getwd()))

