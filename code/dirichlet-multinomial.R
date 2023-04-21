library(here)
library(r4ss)

# Example
# 
# l.skate <- SS_read('Q:/Assessments/Archives/LongnoseSkate/LongnoseSkate_2019/2_base_model')
# 
# l.skate$dat$len_info |> View()
# unique(l.skate$dat$lencomp$FltSvy) # zeros do not have length comp data
# 
# l.skate$dat$age_info
# unique(l.skate$dat$agecomp$FltSvy)
# 
# l.skate$ctl$DoVar_adjust
# l.skate$ctl$dirichlet_parms
# l.skate$ctl$Variance_adjustment_list


copy_SS_inputs(dir.old = here('models/sensitivity_anal/5.7_Update_Selex_and_Biology/23.fix.M.males'),
                 #here('models/Sensitivity_Anal/5.7_Update_Selex_and_Biology/23.fix.M.males'),
               dir.new = here('models/Sensitivity_Anal/dirichlet-multinomial'), 
               overwrite = TRUE)

rex <- SS_read(here('models/Sensitivity_Anal/dirichlet-multinomial'))

# The function r4ss::tune_comps() can adjust control and data files for D-M data weighting
# See https://r4ss.github.io/r4ss/articles/r4ss-intro-vignette.html#tuning-composition-data
# However, it cannot handle 1) models without age data
#                           2) separate retained and discard sample sizes
# So, we will do it manually

# First update data file
rex$dat$use_lencomp
rex$dat$use_lencomp <- 2 # allow different data weights by partition

rex$dat$len_info

# need to add new columns and row
rex$dat$len_info$Fleet <- (1:4)
rex$dat$len_info$Partition <- 0
rex$dat$len_info <- rex$dat$len_info[c('FISHERY', 'FISHERY', 'SURVEY1', 
                                       'SURVEY2', 'SURVEY3'),
                                     c(8,9,1:7)]

rex$dat$len_info$Partition[1:2] <- 1:2

rex$dat$len_info$CompError <- 1
rex$dat$len_info$ParmSelect <- 1:5

# add terminator row
rex$dat$len_info[6,] <- c(-9999, rep(0,8))
rownames(rex$dat$len_info)[6] <- 'TERMINATOR'

rex$dat$len_info

# Now update control file
# make D-M parameter table, structured just like other parameter tables
n_dm_pars <- max(rex$dat$len_info$ParmSelect)
dm_table <- matrix(0, nrow = n_dm_pars,
                   ncol = ncol(rex$ctl$size_selex_parms),
                   dimnames = list(paste0('ln(DM_theta)_', 1:n_dm_pars),
                                   colnames(rex$ctl$size_selex_parms))) |>
  as.data.frame()
dm_table

dm_table$PHASE <- r4ss:::get_last_phase(rex$ctl) + 1

# increase last phase in starter file if necessary
if(rex$start$last_estimation_phase < max(dm_table$PHASE)) {
  rex$start$last_estimation_phase <- max(dm_table$PHASE)
}

dm_table$LO <- -5
dm_table$HI <- 5
dm_table$INIT <- 0
# normal(0, 1.813) prior to counteract parameter transformation 
# i.e., ensure theta ~ uniform
dm_table$PR_type <- 6
dm_table$PRIOR <- 0
dm_table$PR_SD <- 1.813

dm_table

# Add to control file list
rex$ctl$dirichlet_parms <- dm_table

# Get rid of Francis data weights
# (In this case, they are not there, so this is not necessary)
rex$ctl$DoVar_adjust <- 0
rex$ctl$Variance_adjustment_list <- NULL 

# other playing around
rex$ctl$SR_parms['SR_LN(R0)', 'INIT'] <- 15
rex$ctl$SR_parms['SR_sigmaR', 'INIT'] <- 0.4

SS_write(rex, dir = here('models/Sensitivity_Anal/dirichlet-multinomial'),
         overwrite = TRUE)

run(dir = here('models/Sensitivity_Anal/dirichlet-multinomial'),
    exe = here('Executables/SS_V3_30_21/ss_win.exe'), 
    show_in_console = TRUE,  
    extras = '-nohess',
    skipfinished = FALSE)

rex.out <- SS_output(here('models/Sensitivity_Anal/dirichlet-multinomial'), 
          verbose = FALSE)
rex.out$Dirichlet_Multinomial_pars

SS_plots(rex.out, verbose = FALSE)