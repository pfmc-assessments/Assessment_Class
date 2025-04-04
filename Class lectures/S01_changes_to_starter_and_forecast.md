# Changes to starter and forecast files

## starter file suggestions
* don't use .par file unless absolutely necessary (better to use updated initial values from a previous model run in control.ss_new)
* set `Include prior_like for non-estimated parameters (0,1)` to 1 (makes it easier to compare likelihood in profiles--to be covered in a later class) 
* set `Number of datafiles to produce` to 1 or greater to produce `.ss_new` files

## forecast file suggestions
* update benchmark years (convert to negative value representing years before the ending year of the model)
* consider using the -12345 code and the new format for expanded fcast year controls
* set `-1 # Buffer:  enter Control rule target as fraction of Flimit (e.g. 0.75), negative value invokes list of [year, scalar] with filling from year to YrMax` 
* update buffer values using `PEPtools::get_buffer()` https://github.com/pfmc-assessments/PEPtools/blob/main/R/buffer_fxn.R
  * install via `remotes::install_github("pfmc-assessments/PEPtools")`
  * command is something like `PEPtools::get_buffer(years = 2025:2036, sigma = 
     0.5, pstar = 0.45)`
* ignore "FirstYear for caps and allocations" because we have no caps or allocations
* change "stddev of log(realized catch/target catch) in forecast" to 0
* update fixed forcast catches at the bottom for assumed catches in 2025 and 2026 (values will likely be provided by GMT)
