### R script that Ian walked through during 2025-02-21 class

### optional commands to install r4ss and load library
# # install r4ss
# pak::pak("r4ss/r4ss")
# # or
# remotes::install_github("r4ss/r4ss")

# library(r4ss) # loading library avoids requiring `r4ss::` in front of commands below

# define outer directory where all your models are stored
model_dir <- "Archived assessment models/Widow rockfish/"
# read model output
widow2019 <- r4ss::SS_output(dir = file.path(model_dir, "2019 base model/Base_45/"))
# read model output with less info printed to R console
widow2019 <- r4ss::SS_output(dir = "", verbose = FALSE, printstats = FALSE)

# make all plots
r4ss::SS_plots(widow2019) # optionally add verbose = FALSE

# make specific plot group
r4ss::SS_plots(widow2019, plot = 1) # only make plot group 1 (biology)

### comparing two models
# create summary of model results for multiple models 
# (in this case the same model compared to itself just as an example)
widow_summary <- r4ss::SSsummarize(list(widow2019, widow2019))
# also see r4ss::SSgetoutput() if you have vector of directories

# make table comparing models 
# (note that elements of the table can be changed via the `names` and `likenames` inputs)
r4ss::SStableComparisons(widow_summary)
# make plots comparing models
r4ss::SSplotComparisons(widow_summary)
