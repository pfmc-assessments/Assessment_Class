############################################################################
#
# Atsea.comps.r
#
# Build the size- or age-compositions for the hake at-sea fishery.
#
# Andi Stephens, August 2010
#
# The output data are (optionally) copied from the inputs, so that
# "original.data" can always be examined for comparison (and reality checks)
# with intermediate and final products.
#
# Input files:  
#
#    atsea.query.out.<year>.csv
#
# Arguments:
#
#    BY_AGE        Logical.  Do age- rather than the default length-comps.
#    BY_MONTH      Logical.  Do monthly summaries.
#    BY_GENDER     Logical.  Summarize by gender
#    NO_LENGTH     Logical.  Summarize by age-only.
#    lbin.sizes    Vector, sorted smallest to largest.  Not required;
#                  default is 2-cm bins.
#    in.season     List whose elements describe the months to include
#                  in a SEASON for seasonal summaries.  Not required.
#                  Number and size of seasons is arbitrary.  See example
#                  in DEBUG section below.
#    in.pctl       Threshhold quantile for expansion factors.  Default:  0.95.
#    min_Haul      Numeric.  Minimum haul weight.
#    min_T_weight  Numeric.  Minimum total weight.
#    min_sample    Numeric.  Minimum fish in a sample.
#    remove_sparse Logical.  Remove samples without the minimum # fish.
#    which         Output mode.  "pct", "num", or "both".  Default:  "pct".
#    in.filename   Data file.
#    out.filename  Results file.
#    rpt.filename  Report file.
#
#
# Process:
#
#    Evaluate input arguments
#    Extract the input columns of interest
#    Delete records containing missing or bad values
#    Develop the appropriate expansion factor to normalize sampling effort
#        (i.e. Some hauls are more extensively sampled than others).
#
#    haul_weight
#    Use the OBSVR_EST_CATCH where possible, and the VESSEL_EST_CATCH where
#    the former is not available as the haul_weight.
#
#    haul_exp_factor
#    The ratio of the haul_weight to the TOTAL_SAMPLE_WEIGHT per haul is the
#    haul_exp_factor.
#
#    species_est_haul_weight
#    The species_est_haul_weight is the sum of the per-haul species weight
#    multiplied by the haul_exp_factor.
#
#    Age_exp_factor
#    The Age_exp_factor is the per-haul ratio of the species_est_haul_weight
#    to the aged_sample_weight, which is the per-haul WEIGHT of aged samples.
#
#    For lengths, we use the TOTAL_SAMPLE_WEIGHT/SAMPLE_NUMBER for the median
#    per-haul individual weight.
#
#    Extrap_Median_Indiv_Wt
#    Median individual weights per haul are divided by the number of individuals
#    weighed per haul to provide Weighted_Medians.  These are aggregated per-
#    trip to provide By_Trip_Median_Indiv_Wt.
#
#    Extrap_Median_Indiv_Wt is the Weighted Medians where available, and the
#    By_Trip_Median_Indiv_Wt where it is not.  When neither value is available,
#    a search is instituted to find a By_Trip value within the three days
#    prior or following, searching outward from the day in question.  If no
#    value is found, the daily median individual weight for that day is used,
#    and finally, any remaining empty values are filled with the Annual median
#    of the Weighted_Medians.
#
#    Eff_sample_weight
#    Eff_sample_weight is the product of Extrap_Mean_Indiv_Wt and the number of
#    fish lengthed per haul.
#
#    Len_exp_factor
#    Finally, the expansion factor for the lengths is the per-haul ratio of the
#    species_est_haul_weight and the Eff_sample_weight.
#
#    The expansion factor (either age or length) is the normalized catch.
#    Tabulate these according to the summary-control arguments to the function.
#    Report summary statistics and write out results.
#
##############################################################################

Atsea.comps = function(x,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=FALSE,
                       lbin.sizes=NULL, in.season=NULL, in.pctl=0.95,
                       which="pct", min_Haul=0, min_T_weight=0,
                       min_sample=15, remove_sparse=FALSE, NO_LENGTH=FALSE,
                       out.filename, report.filename,ages=0:40) {

# Stop turning my character vectors into factors!
 #   x$Month = substring(x$RETRV_DATE_TIME,6,7)
 #   x$YEAR = substring(x$RETRV_DATE_TIME,1,4)

options(stringsAsFactors = FALSE)

# DEBUG Controls whether we hang onto output data, and whether we run as a
# function with inputs or as a script with defaults

if (!exists("DEBUG") ) {

  DEBUG = FALSE

} # End if

if (DEBUG) {

  # Start clean
  rm(list=ls())

  cat("DEBUG is ON\n"); flush.console();

  # Hah!  Just deleted!
  DEBUG = TRUE

  which = "both"
  min_Haul      = 0
  min_T_weight  = 0
  min_sample    = 15
  remove_sparse = FALSE

  in.filename     = "../Norpac_Annual/atsea.r.out.2008.csv"
  out.filename    = "../results/Atsea.Length.2008.forSS.csv"
  report.filename = "../results/Atsea.2008.report"

  # Note input examples here

  in.pctl    = 0.95
  lbin.sizes = NULL
  #in.season  = NULL
  in.season = list(c(1,2,3),4,5,6,7,8,9,10,c(11,12))

  BY_AGE     = FALSE
  BY_MONTH   = FALSE
  BY_GENDER  = TRUE
  NO_LENGTH  = FALSE

} # End if DEBUG

##############################################################################
#
# Read in utility functions unless they're already there (running in a loop).
#
##############################################################################

if (!exists("try.detach")) {

  source("Utility.functions.r")
  source("Get.age.or.length.r")

}

# Read data and binsizes if given, otherwise use defaults.

input.data = x#read.csv(in.filename)

if ( ! is.null(lbin.sizes) ) { 

  Lbins = as.integer(lbin.sizes)

print("USER ENTERED LBIN.SIZES")
} # End if

BY_SEASON = (! is.null(in.season))

if (BY_SEASON && BY_MONTH) {

  stop("Choose either monthly or seasonal summaries, not both!\n")

} # End if

TEMPORAL = FALSE

if (BY_SEASON || BY_MONTH) {

  TEMPORAL = TRUE

} # End if

##############################################################################
#
# Perform preliminary processing.
#
# Keep only those columns used in the analysis.
#
# Remove bad records and record statistics.
#
##############################################################################

# Before attaching things, make sure they're not already there.

try.detach(atsea.data)
try.detach(input.data)

# 'Attach' exports all the column names into the environment
# so that we don't have to refer to each 'x' as input.data$x.
# Note that you still assign to yyy$x, and after doing so,
# you need to refer to yyy$x.  Do not ask me why this is so
# damned screwy!  Actually, it is because "attach" shows you
# a copy of the columns rather than the true ones.

attach(input.data)

# Use package "date" to convert SAS dates. 

require(date)

Date = as.date(RETRV_DATE_TIME)

Month = date.mdy(Date)[[1]]
YEAR = date.mdy(Date)[[3]]

YEAR[YEAR < 1945] = YEAR + 100


if (length(sort(unique(YEAR))) > 1 ) {

  stop("\n\n *** Fatal error: called Atsea.comps with multi-year dataset.\n\n")

} # End if

# Create trip, haul identifiers

Trip.ID = paste(CRUISE, CRUISE_VESSEL_SEQ, TRIP_SEQ, sep=".")
Haul.ID = paste(CRUISE, CRUISE_VESSEL_SEQ, TRIP_SEQ, HAUL_SEQ, sep=".")

if (DEBUG) {

  original.data = data.frame(CRUISE_VESSEL_SEQ, TRIP_SEQ, TOTAL_SAMPLE_WEIGHT, 
                             OBSVR_EST_CATCH, VESSEL_EST_CATCH, HAUL_SEQ, CRUISE,
                             SPECIES_CODE, COMMON_NAME, SEX_CODE, LENGTH_SIZE, 
                             FREQUENCY, WEIGHT, AGE, RETRV_LATITUDE_DEGREES,
                             RETRV_DATE_TIME, SAMPLE_NUMBER, SPECIES_WEIGHT)

} # End if DEBUG

atsea.data = data.frame(TOTAL_SAMPLE_WEIGHT, OBSVR_EST_CATCH, VESSEL_EST_CATCH, 
                        LENGTH_SIZE, FREQUENCY, WEIGHT, AGE, SEX_CODE,
                        SAMPLE_NUMBER, SPECIES_WEIGHT, Trip.ID, Haul.ID, Date, Month)


detach(input.data)

rm(input.data, Trip.ID, Haul.ID, Date, Month)

attach(atsea.data)

Original.N.records = nrow(atsea.data)

# Collect initial summary data

pre_length.summary = summary(atsea.data$LENGTH_SIZE)
pre_obs_est.summary = summary(atsea.data$OBSVR_EST_CATCH)
pre_ves_est.summary = summary(atsea.data$VESSEL_EST_CATCH)
pre_total.wgt.summary = summary(atsea.data$TOTAL_SAMPLE_WEIGHT)
pre_age.summary = summary(atsea.data$AGE)

# Set up seasons

if (TEMPORAL) {

  # Convert Month to SEASON

  if (BY_SEASON) {

    SEASON = as.character(atsea.data$Month)

    for ( i in 1:length(in.season) ) {

      SEASON[SEASON %in% in.season[[i]]] = paste("Season.", i, sep="")

    } # End for

  } else {

    SEASON = as.numeric(atsea.data$Month)

  } # End if-else

  atsea.data = cbind(atsea.data, SEASON)
  rm(SEASON)

} # End if

# Flag bad records by setting their LENGTH_SIZE to 0 for deletion

# Data without length

atsea.data$LENGTH_SIZE[is.na(LENGTH_SIZE)] = 0
No.length = sum(atsea.data$LENGTH_SIZE == 0, na.rm=T)

# Data without gender

if (BY_GENDER) {

  atsea.data$SEX_CODE[SEX_CODE == 0] = 0
  atsea.data$SEX_CODE[SEX_CODE == "U"] = 0
  atsea.data$SEX_CODE[is.na(SEX_CODE)] = 0

  No.gender = sum(atsea.data$SEX_CODE == 0, na.rm=T)
  atsea.data$LENGTH_SIZE[atsea.data$SEX_CODE == 0] = 0

} else {

  #############################################################################%
  #
  # This is a terrible hack, but it helps when summing the sexes.
  # Note that if we are going by gender, then those not sexed have already
  # been deleted.
  #
  #############################################################################%

  atsea.data$SEX_CODE[atsea.data$SEX_CODE != "F"] = "M"

} # End if

# Data without Date

if (TEMPORAL) {

  atsea.data$Date[is.na(Date)] = 0
  No.month = sum(atsea.data$Date == 0, na.rm=T)
  atsea.data$LENGTH_SIZE[atsea.data$Date == 0] = 0

} # End if TEMPORAL

# Unaged fish

if (BY_AGE) {

  atsea.data$AGE[is.na(AGE)] = 0
  atsea.data$AGE[AGE == -9] = 0
  No.age = sum(atsea.data$AGE == 0, na.rm=T)
  atsea.data$LENGTH_SIZE[atsea.data$AGE == 0] = 0

} # End if BY_AGE


##############################################################################
#
# Now delete the flagged data.  This includes all unlengthed records.
#
# Detach while changing the dataset.  Add the columns we need to save,
# and remove them from the environment.  This ensures that as we remove
# bad records the column lengths stay the same.
#
##############################################################################

detach(atsea.data)

atsea.data = atsea.data[atsea.data$LENGTH_SIZE != 0 , ]

# Find too-sparse samples

Too_few_fish = NULL

if (min_sample > 0) {

  samples = sort(unique(atsea.data$SAMPLE_NUMBER))
  sample_count = rep(0, length(samples))

  for ( i in 1:length(samples) ) {

    sample_count[i] = sum(atsea.data$SAMPLE_NUMBER == samples[i], na.rm=T)

  } # End for

  Too_few_fish = cbind(samples[sample_count < min_sample], sample_count[sample_count < min_sample])
  colnames(Too_few_fish) = c("SAMPLE_NUMBER", "Sample_size")

  if (remove_sparse) {

    Removed.sparse = length(Too_few_fish)
    atsea.data=atsea.data[!atsea.data$SAMPLE_NUMBER %in% Too_few_fish,]
    Too_few_fish = NULL

  } # End if

} # End if min_sample

##############################################################################
#
# Duplicate rows representing multiple fish for LENGTH comps only.
#
##############################################################################

Original.sum.freq = sum(atsea.data$FREQUENCY)

if (! BY_AGE) {

  atsea.data = atsea.data[rep(1:nrow(atsea.data), atsea.data$FREQUENCY),]

} # End if !BY_AGE

atsea.data$FREQUENCY = 1

N.Lengthed.fish = sum(atsea.data$FREQUENCY)

##############################################################################
#
# Data that come in without values can be assumed to approach the population
# median.  The next section of code gets these median values and adjusts
# for missing values.
#
##############################################################################

attach(atsea.data)

##############################################################################
#
# First get per-haul expansion factor by weight.  This gives us an
# estimate of the weight of the target species per haul.  This is
# the expansion factor used for age composition.  The length-composition
# expansion factor takes both this and the median individual weights
# into account.
#
# All WEIGHTS == 0 are Replaced with NA, to avoid skewing medians.
#
##############################################################################

atsea.data$WEIGHT[WEIGHT == 0] = NA

tmp.nlens = aggregate(FREQUENCY, list(Haul.ID), sum, na.rm=T)
names(tmp.nlens) = c("Haul.ID", "N_lens")

N_lens = find.matching.rows(atsea.data, tmp.nlens, "Haul.ID", "Haul.ID", "N_lens")

tmp_species_weight = aggregate(SPECIES_WEIGHT, list(Haul.ID), sum, na.rm=T)
tmp_total_sample_weight= aggregate(TOTAL_SAMPLE_WEIGHT, list(Haul.ID), sum, na.rm=T)

names(tmp_species_weight) = c("Haul.ID", "sum_sp_weight")
names(tmp_total_sample_weight) = c("Haul.ID", "sum_tot_weight")

# Replicate per-haul values per-row.

sum_byhaul_species_weight = find.matching.rows(atsea.data, tmp_species_weight,
                                               "Haul.ID", "Haul.ID", "sum_sp_weight")[[1]]
sum_byhaul_total_sample_weight = find.matching.rows(atsea.data, tmp_total_sample_weight,
                                              "Haul.ID", "Haul.ID", "sum_tot_weight")[[1]]

# Haul weight is either OBSVR_EST_CATCH or VESSEL_EST_CATCH*1000
# OBSVR_EST_CATCH is in mtonnes; VESSEL_EST_CATCH is in kg.

haul_weight = OBSVR_EST_CATCH

tmp.replace = replace.zeros(haul_weight, VESSEL_EST_CATCH * 1000)
haul_weight = tmp.replace[[1]]
Replaced.vessel.catch = tmp.replace[[2]]

# Create per-haul expansion factor

tmp.sum_total = aggregate(TOTAL_SAMPLE_WEIGHT, list(Haul.ID), sum, na.rm=T)
names(tmp.sum_total) = c("Haul.ID", "sum_total_sample_weight")

sum_total_sample_weight = find.matching.rows(atsea.data, tmp.sum_total, "Haul.ID", "Haul.ID",
                                             "sum_total_sample_weight")[[1]]

tmp.sum_spp = aggregate(SPECIES_WEIGHT, list(Haul.ID), sum, na.rm=T)
names(tmp.sum_total) = c("Haul.ID", "sum_species_weight")

sum_species_weight = find.matching.rows(atsea.data, tmp.sum_total, "Haul.ID", "Haul.ID",
                                        "sum_species_weight")[[1]]

haul_exp_factor = haul_weight/sum_total_sample_weight;

species_est_haul_weight = sum_species_weight*haul_exp_factor;

sp_wt_diff = sum(sum_species_weight, na.rm=T)  - sum(SPECIES_WEIGHT, na.rm=T);

# Determine the number of hauls with weight information 

tmp.wtsum = aggregate(WEIGHT, list(Haul.ID), sum, na.rm=T)
n.hauls.with.weight = length(tmp.wtsum)
names(tmp.wtsum) = c("Haul.ID", "wt.sum")
wtsum = find.matching.rows(atsea.data, tmp.wtsum, "Haul.ID", "Haul.ID", "wt.sum")[[1]]

# Staple it all together

detach(atsea.data)

atsea.data = cbind(atsea.data, sum_species_weight, sum_total_sample_weight,
                   haul_weight, haul_exp_factor, species_est_haul_weight, sp_wt_diff, wtsum, N_lens)


# Clean up

rm(sum_species_weight, sum_total_sample_weight, haul_weight, haul_exp_factor,
   species_est_haul_weight, sp_wt_diff, wtsum, N_lens)

if (BY_AGE) {

  # Determine the number of age samples taken

  n.age.samples = length(atsea.data$AGE[atsea.data$AGE != 0])

  # Determine the number of hauls with age information

  attach(atsea.data)

  tmp_agesum = aggregate(atsea.data$AGE, list(atsea.data$Haul.ID), sum, na.rm=T)
  names(tmp_agesum) = c("Haul.ID", "agesum")

  n.hauls.with.age = length(tmp_agesum$agesum[tmp_agesum$agesum!= 0])

  agesum = find.matching.rows(atsea.data, tmp_agesum,
                              "Haul.ID", "Haul.ID", "agesum")[[1]]

  # Get the weight of aged samples per haul

  tmp_aged_sample_weight = aggregate(WEIGHT, list(Haul.ID), sum, na.rm=T)
  names(tmp_aged_sample_weight) = c("Haul.ID", "aged_sample_weight")
  aged_sample_weight = find.matching.rows(atsea.data, tmp_aged_sample_weight,
                                 "Haul.ID", "Haul.ID", "aged_sample_weight")[[1]]

  ##############################################################################
  #
  # Quantiles for Age_exp_factor
  #
  ##############################################################################

  Age_exp_factor = species_est_haul_weight/aged_sample_weight
  pre_age_exp.factor.summary = summary(Age_exp_factor)

  # Reduces the impact of very-very large catches on the assessment

  pctl = quantile(Age_exp_factor, in.pctl, na.rm=T)

  Age_exp_factor[agesum < 40 & Age_exp_factor > pctl] = pctl

  # Deal with sparse samples

  if (!remove_sparse) {

    median_exp_factor = median(Age_exp_factor, na.rm=T)
    Age_exp_factor[SAMPLE_NUMBER %in% Too_few_fish] = median_exp_factor
    Replaced.age.exp.median = length(Too_few_fish)

  } # End if

  post_age_exp.factor.summary = summary(Age_exp_factor)

  detach(atsea.data)
  atsea.data = cbind(atsea.data, Age_exp_factor)
  rm(Age_exp_factor)

} # End if BY_AGE


# Clean up

rm(list=ls(pat="tmp*"))

##############################################################################
#
# Now get expansion for lengths
#
##############################################################################

attach(atsea.data)

# Get individual weight by haul, n_wts by haul

HAUL_WEIGHT = TOTAL_SAMPLE_WEIGHT/SAMPLE_NUMBER

tmp.mhw = aggregate(HAUL_WEIGHT, list(Haul.ID), median, na.rm=T)
tmp.n_wts = aggregate(FREQUENCY, list(Haul.ID), sum, na.rm=T)

names(tmp.mhw) = c("Haul.ID", "target")
names(tmp.n_wts) = c("Haul.ID", "target")

Median_Haul_Wt= find.matching.rows(atsea.data, tmp.mhw, "Haul.ID", "Haul.ID", "target")[[1]]
N_Wts = find.matching.rows(atsea.data, tmp.n_wts, "Haul.ID", "Haul.ID", "target")[[1]]

Median_Haul_Wt[is.na(Median_Haul_Wt)] = 0

Weighted_Medians = Median_Haul_Wt/N_Wts

# Get the by-trip, daily median.  Note that we weight the median by the number weighed.

Daily.Trip.ID = paste(Trip.ID, Date, sep=".")

detach(atsea.data)

atsea.data = cbind(atsea.data, Weighted_Medians, Median_Haul_Wt, N_Wts, Daily.Trip.ID)
rm(Weighted_Medians, Median_Haul_Wt, N_Wts, Daily.Trip.ID) 

attach(atsea.data)

By_Trip_Median_Indiv_Wt = aggregate(Weighted_Medians, list(Daily.Trip.ID), median, na.rm=T)
names(By_Trip_Median_Indiv_Wt) = c("Daily.Trip.ID", "target")

Daily_By_Trip_MIW = find.matching.rows(atsea.data, By_Trip_Median_Indiv_Wt,
                                       "Daily.Trip.ID", "Daily.Trip.ID", "target")[[1]]

# Accumulator for median by haul

Extrap_Median_Indiv_Wt = Median_Haul_Wt

# If there is no haul median, use the trip median for that {Trip.ID, Date}.

tmp.replace = replace.zeros(Extrap_Median_Indiv_Wt, Daily_By_Trip_MIW)
Extrap_Median_Indiv_Wt = tmp.replace[[1]]
Replaced.trip.median = tmp.replace[[2]]

atsea.data = cbind(atsea.data, Extrap_Median_Indiv_Wt, Daily_By_Trip_MIW)

detach(atsea.data)

#rm(list = ls(pat="tmp*"))
rm(Extrap_Median_Indiv_Wt, Daily_By_Trip_MIW)

attach(atsea.data)

# If the trip median doesn't exist, look for a median in the three days prior
# and following the trip.

to_fix = which(Extrap_Median_Indiv_Wt == 0)

for ( i in to_fix ) {

  # We may have fixed this one already, see below

  if ( Extrap_Median_Indiv_Wt[i] == 0 ) {

    offset = 1

    DONE = FALSE

    while (! DONE & offset < 4) {

      # Get median by trip and date +/- offset

      plus.median = Daily_By_Trip_MIW[Trip.ID == Trip.ID[i] & Date[i] == (Date[i]+offset)]
      
      minus.median = Daily_By_Trip_MIW[Trip.ID == Trip.ID[i] & Date[i] == (Date[i]-offset)]

      # Did we get a value?

      if ( length(plus.median) < 1 & length(minus.median) < 1 ) {

        # No.  Increment the offset and try the while loop again.

        offset = offset + 1
        next

      } # End if 

      # Yes!

      value = ifelse (plus.median > 0, plus.median, minus.median)

      # Find all the records for this trip, and set the Extrap value.

      which.ones = which(Trip.ID == Trip.ID[i])

      for ( j in which.ones) {

        # Don't overwrite a Haul median!

        if ( Extrap_Median_Indiv_Wt[j] == 0 ) {

          Extrap_Median_Indiv_Wt[j] = value

        } # End if

      } # End for j

      # Execution goes out to the "for i" loop

      DONE = TRUE

    } # End while

  } # End if

} # End for i

got_fixed = sum(Extrap_Median_Indiv_Wt == 0) - length(to_fix)

# Replace any remaining zeros with the daily median.

tmp.daily.miw = aggregate(Weighted_Medians, list(as.character(Date)), median, na.rm=T)
names(tmp.daily.miw) = c("Date", "target")

atsea.data$Date = as.character(atsea.data$Date)
Daily.miw = find.matching.rows(atsea.data, tmp.daily.miw, "Date", "Date", "target")[[1]]

tmp.replace = replace.zeros(Extrap_Median_Indiv_Wt, Daily.miw)
Extrap_Median_Indiv_Wt = tmp.replace[[1]]
Replaced.with.Daily.miw = tmp.replace[[2]]

# Replace any remaining zeros with the annual median.

Annual_MIW = median(Weighted_Medians)

tmp.replace = replace.zeros(Extrap_Median_Indiv_Wt, rep(Annual_MIW, length(Extrap_Median_Indiv_Wt)))
Extrap_Median_Indiv_Wt = tmp.replace[[1]]
Replaced.with.Annual.miw = tmp.replace[[2]]

Eff_sample_weight = Extrap_Median_Indiv_Wt * N_lens
Len_exp_factor = atsea.data$species_est_haul_weight/Eff_sample_weight

##############################################################################
#
# Quantiles for Len_exp_factor
#
##############################################################################

pre_Len_exp.factor.summary = summary(Len_exp_factor)

# Reduces the impact of very-very large catches on the assessment

pctl = quantile(Len_exp_factor, in.pctl, na.rm=T)

Len_exp_factor[N_lens < 40 & Len_exp_factor > pctl] = pctl

# Deal with sparse samples

if (!remove_sparse) {

  median_exp_factor = median(Len_exp_factor, na.rm=T)
  Len_exp_factor[Haul.ID %in% Too_few_fish] = median_exp_factor
  Replaced.Len.exp.median = length(Too_few_fish)

} # End if

post_Len_exp.factor.summary = summary(Len_exp_factor)

##############################################################################
#
# Now get means and sums of various quantities needed to create summary tables.
#
##############################################################################

# Get the median of the Eff_sample_weight and Extrap_Median_Indiv_Wt, and the sum 
# of FREQUENCY for each Haul.ID

toagg = cbind(Eff_sample_weight, Extrap_Median_Indiv_Wt, FREQUENCY)

tmp.agg = aggregate(toagg, list(Haul.ID), median, na.rm=T)
colnames(tmp.agg) = c("Haul.ID","MEAN_Eff_sample_weight", "MEAN_Extrap_Indiv_Wt")

tmp.agg2 = aggregate(toagg, list(Haul.ID), sum, na.rm=T)

tmp.agg = cbind(tmp.agg[,1:3],tmp.agg2[,4])
colnames(tmp.agg)[4] = "SUM_FREQUENCY"

targets = c("MEAN_Eff_sample_weight", "MEAN_Extrap_Indiv_Wt", "SUM_FREQUENCY")

Summary_data = find.matching.rows(atsea.data, tmp.agg, "Haul.ID", "Haul.ID", targets)

detach(atsea.data)

atsea.data = cbind(atsea.data, Summary_data, Len_exp_factor,
                   Eff_sample_weight)

rm(Summary_data, Len_exp_factor, Extrap_Median_Indiv_Wt,
   Eff_sample_weight)

##############################################################################
#
# Record and flag for removal those records for which we couldn't generate
# enough summary data
#
##############################################################################

attach(atsea.data)

atsea.data$MEAN_Eff_sample_weight[is.na(MEAN_Eff_sample_weight)] = 0
No.MEAN_Eff_sample_weight = sum(atsea.data$MEAN_Eff_sample_weight == 0, na.rm=T)
atsea.data$LENGTH_SIZE[atsea.data$MEAN_Eff_sample_weight == 0] = NA
                        
atsea.data$MEAN_Extrap_Indiv_Wt[is.na(MEAN_Extrap_Indiv_Wt)] = 0
No.MEAN_Extrap_Indiv_Wt = sum(atsea.data$MEAN_Extrap_Indiv_Wt == 0, na.rm=T)
atsea.data$LENGTH_SIZE[atsea.data$MEAN_Extrap_Indiv_Wt == 0] = NA
                         
atsea.data$SUM_FREQUENCY[is.na(SUM_FREQUENCY)] = 0
No.FREQUENCY = sum(atsea.data$SUM_FREQUENCY == 0, na.rm=T)
atsea.data$LENGTH_SIZE[atsea.data$SUM_FREQUENCY == 0] = NA
              
detach(atsea.data)

atsea.data = atsea.data[! is.na(atsea.data$LENGTH_SIZE), ] 

attach(atsea.data)

# How many records now?

N.lengthed.fish = nrow(atsea.data)

##############################################################################
#
# Create Use_Length:  the length-bins used in the assessment.
#
# If no length bins were read in, just use the rounded length
# in 2-cm bins
#
##############################################################################

if ( is.null(lbin.sizes) ) {
print("lbin.sizes is NULL")
  Use_Length = LENGTH_SIZE
  Use_Length[Use_Length <= 10] = 10
  Use_Length[Use_Length >= 92] = 92

  is.odd = Use_Length %% 2 > 0
  Use_Length[is.odd] = Use_Length[is.odd] - 1

  # Make sure all lengths show up in the table, even if not
  # in the data

  post_length.summary = summary(Use_Length)

  # Needs to be factored in order to have unused levels
  # show up in the tables.

  Use_Length = factor(Use_Length, levels=seq(20,70,2))

} else {
print("lbin.sizes is not NULL")

  Use_Length = rep(0, length(LENGTH_SIZE))

  # Iteratively replace zeros with next Lbin
  # Assumes Lbins are sorted least-greatest

  for (i in 2:length(Lbins)) {

    minsize = Lbins[i]  
    lesser = Lbins[i-1]

    Use_Length[Use_Length == 0 & LENGTH_SIZE < minsize] = lesser

  } # End for

  # Finish up largest bin
    
  Use_Length[Use_Length == 0] = Lbins[i]

  post_length.summary = summary(Use_Length * 10)

  # Make sure all lengths show up in the table, even if not
  # in the data

  Use_Length = factor(Use_Length, levels=Lbins)

} # End if

#print("USE_LENGTH")
#print(Use_Length)
#print("end Use_LENGTH")

detach(atsea.data)

atsea.data = cbind(atsea.data, Use_Length)
rm(Use_Length)


  # Weather or not we're doing age comps

  atsea.data$AGE[atsea.data$AGE < min(ages)] = min(ages)
  atsea.data$AGE[atsea.data$AGE > max(ages)] = max(ages)

  # Make sure all ages show up in the table, as zeros if not represented
  # in the data

  post_age.summary = summary(atsea.data$AGE)

  atsea.data$AGE = factor(atsea.data$AGE, levels=ages)

if (!BY_AGE) {

  atsea.data$AGE = addNA(atsea.data$AGE, ifany=T)

} # End if

##############################################################################
#
# Output data
#
##############################################################################

# Size comps by number and percentage

attach(atsea.data)

if ( BY_AGE ) {

  exp_factor = Age_exp_factor

} else {

  exp_factor = Len_exp_factor

} # End if-else

# The order in which factors are added matters!

by.factors = NULL

if (TEMPORAL) { by.factors =  "SEASON + " }

by.factors = paste(by.factors, "SEX_CODE + Use_Length + AGE", sep="")

form = as.formula(paste("exp_factor ~ ", by.factors, sep=""))
fish.form = as.formula(paste("FREQUENCY ~ ", by.factors, sep=""))
hauls.form = as.formula(paste("FREQUENCY ~ ", "Haul.ID + ", by.factors, sep=""))

atsea.comps.num = round(xtabs(form))
fish.num = round(xtabs(fish.form))
hauls.num = round(xtabs(hauls.form))

##############################################################################
#
#  Final SUMMARY STATISTICS
#
# Number of sample trips for length comps this year
# Total weight of fish sampled this year
# Number of individuals lengthed in sample hauls this year 
#
##############################################################################

# Helper function for counting trips
nfacts = function(x) { nlevels(factor(x)) }

# Collect Weight used for each Haul.ID
# NB:  Use of as.character and as.numeric because R insists on converting
# things to goddamn FACTORS!


SAMPLES = sort(unique(Haul.ID))

Weight = matrix(nrow = length(SAMPLES), ncol=2, 0.0)
colnames(Weight) = c("W","M_or_S")

for ( i in 1:length(SAMPLES) ) {

  Weight[i,1] = Eff_sample_weight[Haul.ID == SAMPLES[i]][1] 

  # Convert kg to metric tonnes for reporting.

  Weight[i,1] = as.numeric(Weight[i,1]) / 1000

  if (TEMPORAL) {

    Weight[i,2] = SEASON[Haul.ID == SAMPLES[i]][1]

  } # End if

} # End for
        
Weight = data.frame(Weight)

if (TEMPORAL) {

    Hauls = aggregate(Haul.ID, list(SEASON), FUN=nfacts)

    Rpt.Weight = aggregate(as.numeric(Weight$W),
                           list(Weight$M_or_S), sum, na.rm=T)

    N.Fish = aggregate(FREQUENCY, list(SEASON), sum, na.rm=T)

    sampled = cbind(Rpt.Weight, Hauls[,2], N.Fish[,2])
    colnames(sampled) = c("Season","Metric.Tonnes","Hauls","N.Fish")

} else {

  Hauls = nfacts(Haul.ID)

  Rpt.Weight = sum(Weight[,1])

  N.Fish = sum(FREQUENCY, na.rm=T)

  # Staple them together

  sampled = c(Rpt.Weight, Hauls, N.Fish)
  dim(sampled) = c(1,3)
  colnames(sampled) = c("Metric.Tonnes", "Hauls", "N.Fish")

} # End if-else

post_haul.summary = summary(atsea.data$Extrap_Median_Indiv_Wt)
post_total.wgt.summary = summary(atsea.data$Eff_sample_weight)

##############################################################################
#
#  Write out data
#
#  Each of these tables has to be converted to a matrix in order to be properly
#  formatted in output.  Data can be output as percentages, numbers or both.
#  Matrices may be 1-4 dimensional.
#
##############################################################################

options(warn=-1)

outfname = out.filename
if (BY_AGE) {

  cat(file=outfname, "\n\nAtsea age comps, Year:", YEAR[1], "\n\n", append=T)

  cat(file=outfname, "Aged Fish\n\n", append=T)

} else {

  cat(file=outfname, "\n\nAtsea size comps, Year:", YEAR[1], "\n\n", append=T)

  cat(file=outfname, "Lengthed Fish\n\n", append=T)

} # End if-else BY_AGE

# Write out the summary data

outdata = sampled

write.table(file=outfname, x=outdata, sep=",", col.names=T, row.names=F, append=T)

if (TEMPORAL) {

  Mynames = attributes(atsea.comps.num)$dimnames$SEASON

  if (BY_MONTH) {

    Mynames = month.abb[as.numeric(Mynames)]

  } # End if BY_MONTH

} # End if

while ( TRUE ) {

  target = atsea.comps.num

  # By percent?

  if (which == "num") {
   
    by_pct = FALSE
    cat(file=outfname, "\nIn numbers\n\n", append = T)

  } else {

    by_pct = TRUE
    cat(file=outfname, "\nAs percentages\n\n", append = T)

  } # End if-else which

  # Output summary control

  selection = "length.only"

  if (BY_AGE) {

    if (NO_LENGTH) {

      selection = "age.only"

    } else {

      selection = "age.and.length"

    } # End if-else NO_LENGTH

  } # End if BY_AGE

  
  if (TEMPORAL) {

    # Need to rebuild the tables for this case.

    outdata = NULL

    # Use m to index SEASON

    for ( m in 1:attributes(target)$dim[1] ) {

      outdata = get.age.or.length(target[m,,,], hauls.num[,m,,,], fish.num[m,,,],
                                  output=selection, by.sex=BY_GENDER, pct=by_pct)

     if (BY_AGE) {

          cat(file=outfname, "\n", Mynames[m], ",,,,,AGE\n", append=T)

      } else {

          cat(file=outfname, "\n", Mynames[m], ",,,,,LENGTH BIN\n", append=T)

      } # End if

      write.table(file=outfname, x=outdata, col.names=T, row.names=F, sep=",", append=T)

    } # End for m

  } else {

    # Not summarizing by Month/Season

    outdata = get.age.or.length(target, hauls.num, fish.num, output=selection,
                                by.sex=BY_GENDER, pct=by_pct)

    write.table(file=outfname, x=outdata, col.names=T, row.names=F, sep=",", append=T)

  } # End if-else TEMPORAL

  # Go through once more?

  if ( which == "both") {

    which = "num"

  } else {

    break

  } # End if-else which

} # End while


##############################################################################
#
#  Write out the error statistics
#
##############################################################################


outfile = file(report.filename, open="wt")

sink(outfile, split=T, append=T, type="output")

  cat("\n\nYear:", YEAR[1], "\n")
  cat("\nOriginal dataset:", Original.N.records, "  records\n")
  cat("\nSum of FREQUENCY:", Original.sum.freq, "  records\n\n")

  cat("Run conditions:\n\n")

  if (BY_AGE) {

    cat("\tAge comps generated\n")
 
  } else {

    cat("\tLength comps generated\n")

  } # End if-else


  if (BY_GENDER) { cat("\tConditioned on SEX\n") }
  if (BY_MONTH)  { cat("\tConditioned on SAMPLE_MONTH\n") }
  if (BY_SEASON) {

    cat("\tConditioned on SEASON\n")
    cat("\n\tSeasons were:\n")

    for ( i in 1:length(in.season) ) {

      cat("\t\t", "Season", i, " Months", in.season[[i]], "\n")

    } # End for

  } # End if

  if (min_sample > 0) {

    cat("\nMinimum sample size:", min_sample, "\n")
    if (nrow(Too_few_fish) > 0) {

      cat("Smaller samples:\n")
      print(Too_few_fish)

    } else {

      cat("No samples were smaller than the minimum\n")

    } # End if-else

    cat("\n\n")

  } # End if

  if (min_Haul > 0) {

    cat("\nMinimum haul weight:", min_Haul, "\n")
    cat("Replaced  ", Replaced.with.min.haul_weight,
        "  haul weights with minimum\n")

  } # End if

  if (min_T_weight > 0) {

    cat("\nMinimum total weight:", min_T_weight, "\n")
    cat("Replaced  ", Replaced.with.min_T_weight,
        "  total weights with minimum\n")

  } # End if


  cat("\nPre-processing lengths (in mm):\n")
  print(pre_length.summary)

  cat("\nPost-processing lengths (in mm):\n")
  print(post_length.summary)

  cat("\nPre-processing observer estimates (in kg):\n")
  print(pre_obs_est.summary)

  cat("\nPost-processing observer estimates (in kg):\n")
  print(summary(OBSVR_EST_CATCH))

  cat("\nPre-processing vessel estimates (in kg):\n")
  print(pre_ves_est.summary)

  cat("\nPost-processing vessel estimates (in mtonnes):\n")
  print(summary(VESSEL_EST_CATCH))

  cat("\nPre-processing total weights (in kg):\n")
  print(pre_total.wgt.summary)

  cat("\nPost-processing total Effective weights (in mtonnes):\n")
  print(post_total.wgt.summary)

  if (BY_AGE) {

    cat("\nPre-truncation expansion factor summary:\n")
    print(pre_age_exp.factor.summary)

    cat("\nPost-truncation expansion factor summary:\n")
    print(post_age_exp.factor.summary)

    cat("\nPre-processing ages:\n")
    print(pre_age.summary)

    cat("\nPost-processing ages:\n")
    print(post_age.summary)

    cat("\nNumber of hauls with aged fish:  ", n.hauls.with.age, "\n")

    cat("\nAged sample weight\n")
    print(summary(aged_sample_weight))

  } else {

    cat("\nPre-truncation expansion factor summary:\n")
    print(pre_Len_exp.factor.summary)

    cat("\nPost-truncation expansion factor summary:\n")
    print(post_Len_exp.factor.summary)

    cat("Number of lengthed fish:  ", N.Lengthed.fish, "\n")

  } # End if-else

  cat("\nNumber removed:\n")

  if (BY_AGE) {

    cat("\nRemoved  ", No.age, "  records without AGE.\n")

  } # End if

  cat("Removed  ", No.length, 
                       "  records without recorded LENGTH_SIZE.\n",
                        append=T)

  if (BY_GENDER) { 

    cat("Removed  ", No.gender, "  records without SEX_CODE.\n")

  } # End if

  if (TEMPORAL) {

    cat("Removed  ", No.month, "  records without SAMPLE_MONTH.\n")

  } # End if

  if (remove_sparse) {

    cat("Removed  ", Removed.sparse, "  samples without", min_sample, "fish\n")

  } # End if

  cat("\nWeights replaced in remaining data:\n")

  cat("Replaced  ", Replaced.vessel.catch, "  observer estimates with vessel", 
                       "estimates.\n")
  cat("Replaced  ", Replaced.trip.median, "  individual weights with median",
                       "by-trip individual weight\n")
  cat("Replaced  ", got_fixed, "  individual weights with median",
                       "by-trip individual weight from preceeding or following day(s).\n")
  cat("Replaced  ", Replaced.with.Daily.miw, "  individual weights with median",
                       "daily weights\n")
  cat("Replaced  ", Replaced.with.Annual.miw, "  individual weights with median",
                       "annual weights\n")

  cat("\n\nDifference between the sum of species weight and SPECIES_WEIGHT:", sp_wt_diff[1], "\n")

  if (!remove_sparse) {

    if (BY_AGE) {


      cat("\nReplaced  ", Replaced.age.exp.median, "  Age expansion factors with median",
                        "in sparse samples\n")


    } else {

      cat("\nReplaced  ", Replaced.Len.exp.median, "  Length expansion factors with median",
                        "in sparse samples\n")

    } # End if-else BY_AGE

  } # End if sparse

  cat("\nRecords removed because expansion factors couldn't be generated\n")
                       
  cat("Removed  ", No.MEAN_Eff_sample_weight, "  records with no MEAN_Eff_sample_weight\n")
  cat("Removed  ", No.MEAN_Extrap_Indiv_Wt, "  records with no MEAN_Extrap_Indiv_Wt\n")
  cat("Removed  ", No.FREQUENCY, "  records with no FREQUENCY\n")

  cat("\nRecords used in final analysis:", N.lengthed.fish, "\n")

sink()
close(outfile)

# Clean up

detach(atsea.data)
options(warn=0)



} # End function atsea.comps

#
# That's All, Folks!
#
##############################################################################
