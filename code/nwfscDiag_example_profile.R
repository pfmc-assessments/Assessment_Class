library(nwfscDiag)
library(here)

profile.settings <- get_settings_profile(parameters = rep('SR_LN(R0)', 2),
                                         # as discovered, the range here is too wide
                                         # suggest ~9.75-11 for Rex
                                         low = rep(-2,2), high = rep(2,2),
                                         step_size = rep(0.25,2),
                                         param_space = rep('relative',2),
                                         use_prior_like = c(0,1)) 
settings <- get_settings(settings = list(base_name = 'default_settings',
                             run = 'profile',
                             profile_details = profile.settings))
  
settings$exe <- 'ss_win'
run_diagnostics(mydir = here('models/Sensitivity_Anal/model8a two fleet Francis'), 
                model_settings = settings)
