## ----eval=FALSE--------------------------------------------------------------------------------------------------------
# pak::pkg_install("pfmc-assessments/nwfscSurvey")


## ----plot_cpue, fig.show='hide'----------------------------------------------------------------------------------------
library(nwfscSurvey)
# Pull survey data from NWFSC data warehouse
catch <- pull_catch(common_name = "rex sole", # all lowercase 
                    survey = "NWFSC.Combo",
                    dir = NULL) # this is the default
plot_cpue(dir = NULL, catch)




## ----strata, results='hide'--------------------------------------------------------------------------------------------
strata <- CreateStrataDF.fn(
  names = c("shallow_s", "mid_s", "deep_s", "shallow_n", "mid_n", "deep_n"), 
  depths.shallow = c( 55,   200, 300,    55, 200, 300),
  depths.deep    = c(200,   300, 400,   200, 300, 400),
  lats.south     = c( 32,    32,  32,    42,  42,  42),
  lats.north     = c( 42,    42,  42,    49,  49,  49))
strata




## ----index_plot, fig.show='hide'---------------------------------------------------------------------------------------
biomass <- get_design_based(data = catch,  
                            strata = strata)

plot_index(data = biomass, plot = 1)

