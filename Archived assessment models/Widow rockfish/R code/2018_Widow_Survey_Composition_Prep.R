#===============================================================================
# NWFSC Survey Data Package
#===============================================================================

# https://github.com/nwfsc-assess/nwfscSurvey
devtools::install_github("nwfsc-assess/nwfscSurvey", build_vignettes = TRUE)

# Load the packaged
library(nwfscSurvey)
# Look at the vignette
vignette("nwfscSurvey")
# Look at all the functions in the package
ls("package:nwfscSurvey")
?PullCatch.fn

#===============================================================================
#=============          NWFSC Combo          ===================================
#===============================================================================
# Change wd
old_Wd <- getwd()
setwd("..")
setwd("2019 update/Data")

catch = PullCatch.fn(Name = "widow rockfish", SurveyName = "NWFSC.Combo", SaveFile = TRUE, Dir = getwd()) 

bio   = PullBio.fn(Name = "widow rockfish", SurveyName = "NWFSC.Combo", SaveFile = TRUE, Dir = getwd())

head(catch)
head(bio)

# Can pull data based on the general name (Name) of the scientific name(SciName). The default year range (YearRange)
# is set to cover all potential years.  The SurveyName options are: Triennial, AFSC.Slope, NWFSC.Slope, NWFSC.Shelf
# NWFSC.Combo, NWFSC.Hypoxia, NWFSC.Santa.Barb.Basin, or NWFSC.Video. These data pulls can also be saved to a specified 
# directory using SaveFile = TRUE and Dir = "directory to save file". 
# load("Catch_2018-08-06__NWFSC.Combo_2018-08-06.rda")


# Create Stratafication:
# The stratafication areas are calculated from the SA3 file which is attached to the package. This is consistent with the 2015 widow assessment
strata = CreateStrataDF.fn(names=c("shallow_s", "deep_s", "shallow_n", "deep_n"), 
                           depths.shallow = c( 55,   183,    55,  183),
                           depths.deep    = c(183,   400,   183,  400),
                           lats.south     = c(34.5,  34.5, 40.5, 40.5),
                           lats.north     = c(40.5,  40.5, 49.0, 49.0))

strata

#============================================================================================
#Length Biological Data 
#============================================================================================
len = bio
len.bins = seq(8, 56, 2)

# Calculate the input sample size
# Calculate samples size
n_samp_srv <- function(ntows, nsample){
  neff = ifelse(nsample/ntows < 55, 
                ntows + 0.0707 * nsample,
                4.89 * ntows)
  return(neff)
}

n <- c()
ntows <- c()
ncatch <- c()
nages = c()
nlengths <- c()
for(i in 1:length(unique(bio$Year))){
  yr <- sort(unique(bio$Year))[i]
  bio_sub <- bio[which(bio$Year == yr),]
  catch_sub <- catch[which(catch$Year == yr),]
  ntows[i] <- length(unique(bio_sub$Trawl_id))
  nsample = nrow(bio_sub)
  nages[i] <- length(which(!is.na(bio_sub$Age)))
  nlengths[i] <- length(which(!is.na(bio_sub$Length_cm)))
  n[i] <- n_samp_srv(ntows[i], nsample)
  # 2003 should be 33
  ncatch[i] <- length(which(catch_sub$total_catch_numbers != 0))
}

# Expand and format length composition data for SS
LFs <- SurveyLFs.fn(dir = getwd(), datL = len, datTows = catch,  
                    strat.df = strata, lgthBins = len.bins, gender = 3, 
                    sexRatioStage = 2, sexRatioUnsexed = 0.5, maxSizeUnsexed = 28, 
                    nSamps = round(n, 0), month = 8.8, fleet = 8)

# The code offers two options for applying the sex ratio based on expansion stage. The sex ratio will be
# applied based on a tow basis first if sexRatioStage = 1. The other option applies the sex ratio to the
# expanded numbers of fish across a whole strata (sexRatioStage = 2, this was the option applied to the
# NWFSC combo survey data in the past).




PlotFreqData.fn(dir = getwd(), dat = LFs, survey = "NWFSC Shelf-Slope Survey", ylim=c(0, max(len.bins) + 4), yaxs="i", ylab="Length (cm)", dopng = TRUE)
PlotSexRatio.fn(dir = getwd(), dat = len, data.type = "length", dopng = TRUE, survey = "NWFSC Shelf-Slope Survey", main = "NWFSC Shelf-Slope Survey")

#============================================================================================
# Marginal Age Data 
#============================================================================================
age = bio
age.bins = 0:40

n = GetN.fn(dir = getwd(), dat = age, type = "age", species = "shelfrock", printfolder = "forSS")

# Exand and format the marginal age composition data for SS
Ages <- SurveyAFs.fn(dir = getwd(), datA = age, datTows = catch,  
                     strat.df = strata, ageBins = age.bins, 
                     sexRatioStage = 2, sexRatioUnsexed = 0.50, maxSizeUnsexed = 5, 
                     gender = 3, nSamps = n, month = 7, fleet = 8)



PlotFreqData.fn(dir = getwd(), dat = Ages, survey = "NWFSC Shelf-Slope Survey", ylim=c(0, max(age.bins) + 2), yaxs="i", ylab="Age (yr)", dopng=TRUE)
PlotVarLengthAtAge.fn(dir = getwd(), dat = age, survey = "NWFSC Shelf-Slope Survey", dopng = TRUE) 
PlotSexRatio.fn(dir = getwd(), dat = age, data.type = "age", survey = "NWFSC Shelf-Slope Survey", dopng = TRUE, main = "NWFSC Shelf-Slope Survey",)

#============================================================================================
# Conditional Age-at-Length Data
#============================================================================================
Ages <- SurveyAgeAtLen.fn (dir = getwd(), datAL = age, datTows = catch, 
                           strat.df = strata, lgthBins = len.bins, ageBins = age.bins, partition = 0,
                           month = 8.8, fleet = 8)

par(mfrow = c(2, 8))

# Plot it
for(i in 1:length(sort(unique(Ages$female$year)))){
  yr <- sort(unique(Ages$female$year))[i]
  female_dat <- Ages$female[which(Ages$female$year == yr),]
  male_dat <- Ages$male[which(Ages$male$year == yr),]
  
  plot(NA, NA, ylab = "Length (cm)", xlab = "Age", ylim = c(0, 60), xlim = c(0, 40))
  legend("topleft", legend = yr, bty = "n")
  
  # Plot females
  for(j in 1:nrow(female_dat)){
    symbols(y = rep(female_dat$LbinHi[j], 41), x = 0:40, circles = female_dat[j,10:50]/sum(female_dat[j,10:50]), inches = 0.1/length(which(female_dat[j,10:50] != 0)), fg = 2, add = TRUE)
  }

  symbols(y = rep(55, 3), x = c(10, 20, 30), circles = c(0.001, 0.5, 1), inches = 0.1, add = TRUE)
  text(y = rep(55, 3), x = c(15, 25, 35), labels = c(0.001, 0.5, 1))
}


par(mfrow = c(2, 8))

# Plot males
for(i in 1:length(sort(unique(Ages$male$year)))){
  yr <- sort(unique(Ages$male$year))[i]
  female_dat <- Ages$male[which(Ages$male$year == yr),]
  
  plot(NA, NA, ylab = "Length (cm)", xlab = "Age", ylim = c(0, 60), xlim = c(0, 40))
  legend("topleft", legend = yr, bty = "n")
  
  # Plot females
  for(j in 1:nrow(female_dat)){
    symbols(y = rep(female_dat$LbinHi[j], 41), x = 0:40, circles = female_dat[j,10:50]/sum(female_dat[j,10:50]), inches = 0.1/length(which(female_dat[j,10:50] != 0)), fg = 3, add = TRUE)
  }
  
  symbols(y = rep(55, 3), x = c(10, 20, 30), circles = c(0.001, 0.5, 1), inches = 0.1, add = TRUE)
  text(y = rep(55, 3), x = c(15, 25, 35), labels = c(0.001, 0.5, 1))
}
