############################################################################
#
# Atsea.comps.r
#
# Build the size- or age-compositions for the hake at-sea fishery.
# This version of the code uses the pre-2008 variables.
#
# Andi Stephens, August 2010
#
# The output data are (optionally) copied from the inputs, so that
# "original.data" can always be examined for comparison (and reality checks)
# with intermediate and final products.
#
# Input files:  
#
#    <year>.r.out.csv
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
#    The expansion factor (either age or length) is the normalized catch.
#    Tabulate these according to the summary-control arguments to the function.
#    Report summary statistics and write out results.
#
##############################################################################

 Old.Atsea.comps = function(x,BY_AGE=FALSE, BY_MONTH=FALSE, BY_GENDER=FALSE,
                        lbin.sizes=NULL, in.season=NULL, in.pctl=0.95,
                        which="pct", min_Haul=0, min_T_weight=0,
                        min_sample=15, remove_sparse=FALSE, NO_LENGTH=FALSE,
                        out.filename, report.filename) {

# Stop turning my character vectors into factors!
# Or the WRATH of ANDI shall rain down upon all the Ripleys of the earth!

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

  which = "pct"
  min_Haul      = 0
  min_T_weight  = 0
  min_sample    = 15
  remove_sparse = FALSE

  #in.filename     = "../Norpac_Annual/atsea.r.out.2005.csv"
  #in.filename     = "../Norpac_Annual/atsea.age.r.out.2005.csv"
  #out.filename    = "../results/NORPAC.Age.forSS.csv"
  #report.filename = "../results/NORPAC.2005.report"

  # Note input examples here

  in.pctl    = 0.95
  lbin.sizes = NULL
  in.season  = NULL
  #in.season = list(c(1,2,3),4,5,6,7,8,9,10,c(11,12))

  BY_AGE     = TRUE
  BY_MONTH   = FALSE
  BY_GENDER  = FALSE
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
# Note that the input file was created by merging two others,
# hence it is a "csv" file not an "sql" file as in other code.

# input.data = read.sql(in.filename)

input.data = x #read.csv(in.filename)

if ( ! is.null(lbin.sizes) ) { 

  Lbins = as.integer(lbin.sizes)

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

print("HERE")

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

# Can't use the "date" package on this format; have to do
# it manually.  Extract month in format "MAY", convert to
# numeric "5".

Month = gsub("\\d*-\\d*", "", input.data$HAUL_DATE)

Month = match(tolower(Month), tolower(month.abb))

YEAR[YEAR < 1945] = YEAR + 100

if (length(sort(unique(YEAR))) > 1 ) {

  stop("\n\n *** Fatal error: called Old.Atsea.comps with multi-year dataset.\n\n")

} # End if

print("HERE1")


if (DEBUG) {

  original.data = data.frame(HAUL_JOIN, HAUL_DATE, SAMPLE_WEIGHT, AGE,
                             SPECIES, YEAR, SEX, LENGTH, SAMPLE_NUMBER, 
                             FREQUENCY, EXTRAPOLATED_WEIGHT, Month,
                             stringsAsFactors=FALSE)
                             

} # End if DEBUG

atsea.data = data.frame(HAUL_JOIN, SAMPLE_WEIGHT, AGE,
                        SPECIES, YEAR, SEX, LENGTH, SAMPLE_NUMBER,
                        FREQUENCY, EXTRAPOLATED_WEIGHT,
                        Month, stringsAsFactors=FALSE)
                             

detach(input.data)

rm(input.data, Month, YEAR)

attach(atsea.data)

Original.N.records = nrow(atsea.data)

# Collect initial summary data

pre_length.summary = summary(atsea.data$LENGTH)
pre_EXTRAPOLATED_WEIGHT.summary = summary(atsea.data$EXTRAPOLATED_WEIGHT)
pre_SAMPLE_WEIGHT.summary = summary(atsea.data$SAMPLE_WEIGHT)
pre_age.summary = summary(atsea.data$AGE)


print("HERE2")

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

# Flag bad records by setting their LENGTH to 0 for deletion

# Data without length

atsea.data$LENGTH[is.na(LENGTH)] = 0
No.length = sum(atsea.data$LENGTH == 0, na.rm=T)

# Data without gender

if (BY_GENDER) {

  atsea.data$SEX[SEX == 0] = 0
  atsea.data$SEX[SEX == "U"] = 0
  atsea.data$SEX[is.na(SEX)] = 0

  No.gender = sum(atsea.data$SEX == 0, na.rm=T)
  atsea.data$LENGTH[atsea.data$SEX == 0] = 0

} else {

  #############################################################################%
  #
  # This is a terrible hack, but it helps when summing the sexes.
  # Note that if we are going by gender, then those not sexed have already
  # been deleted.
  #
  #############################################################################%

  atsea.data$SEX[atsea.data$SEX != "F"] = "M"

} # End if

print("HERE3")


# Data without Date

if (TEMPORAL) {

  atsea.data$Date[is.na(Month)] = 0
  No.month = sum(atsea.data$Month == 0, na.rm=T)
  atsea.data$LENGTH[atsea.data$Month == 0] = 0

} # End if TEMPORAL

# Unaged fish

if (BY_AGE) {

  atsea.data$AGE[is.na(AGE)] = 0
  atsea.data$AGE[AGE == -9] = 0
  No.age = sum(atsea.data$AGE == 0, na.rm=T)
  atsea.data$LENGTH[atsea.data$AGE == 0] = 0

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

atsea.data = atsea.data[atsea.data$LENGTH != 0 , ]

# Find too-sparse samples

Too_few_fish = NULL

if (min_sample > 0) {

  samples = sort(unique(atsea.data$HAUL_JOIN))
  sample_count = rep(0, length(samples))

  for ( i in 1:length(samples) ) {

    sample_count[i] = atsea.data$SAMPLE_NUMBER[atsea.data$HAUL_JOIN == samples[i]][1]

  } # End for

  Too_few_fish = cbind(samples[sample_count < min_sample], sample_count[sample_count < min_sample])
  colnames(Too_few_fish) = c("HAUL_JOIN", "SAMPLE_NUMBER")

  if (remove_sparse) {

    Removed.sparse = length(Too_few_fish)
    atsea.data=atsea.data[!atsea.data$HAUL_JOIN %in% Too_few_fish,]
    Too_few_fish = NULL

  } # End if

} # End if min_sample

##############################################################################
#
# Duplicate rows representing multiple lengthed (NOT AGED) fish
#
##############################################################################

atsea.data$FREQUENCY[is.na(atsea.data$FREQUENCY)] = 1
Original.sum.freq = sum(atsea.data$FREQUENCY, na.rm=T)

if (!BY_AGE) { 

  atsea.data = atsea.data[rep(1:nrow(atsea.data), atsea.data$FREQUENCY),]

}

atsea.data$FREQUENCY = 1

##############################################################################
#
# Data that come in without values can be assumed to approach the population
# median.  The next section of code gets these median values and adjusts
# for missing values.
#
# For these data, the expansion factor is the same for AGE and LENGTH comps.
#
##############################################################################

med_haul_weight = aggregate(atsea.data$EXTRAPOLATED_WEIGHT, list(atsea.data$HAUL_JOIN), median, na.rm=T)
med_sample_weight = aggregate(atsea.data$SAMPLE_WEIGHT, list(atsea.data$HAUL_JOIN), median, na.rm=T)
med_sample_number = aggregate(atsea.data$SAMPLE_NUMBER, list(atsea.data$HAUL_JOIN), median, na.rm=T)

sum_FREQ = aggregate(atsea.data$FREQUENCY, list(atsea.data$HAUL_JOIN), sum, na.rm=T)

tmp.stats = cbind(med_haul_weight, med_sample_weight[,2], med_sample_number[,2], sum_FREQ[,2])

names(tmp.stats) = c("HAUL_JOIN", "med_haul_weight", "med_sample_weight", "med_sample_number", "sum_FREQ")

targets = names(tmp.stats)[-1]

tmp.stats = find.matching.rows(atsea.data, tmp.stats, "HAUL_JOIN", "HAUL_JOIN", targets)

atsea.data = cbind(atsea.data, tmp.stats)

rm(med_haul_weight, med_sample_weight, med_sample_number, sum_FREQ, tmp.stats)

attach(atsea.data)

##############################################################################
#
# Now do replacements for missing values
#
##############################################################################

tmp.replace = replace.zeros(EXTRAPOLATED_WEIGHT, med_haul_weight)
atsea.data$EXTRAPOLATED_WEIGHT = tmp.replace[[1]]
Replaced.with.med.haul.weight = tmp.replace[[2]]

tmp.replace = replace.zeros(SAMPLE_WEIGHT, med_sample_weight)
atsea.data$SAMPLE_WEIGHT = tmp.replace[[1]]
Replaced.with.med.sample.weight = tmp.replace[[2]]

tmp.replace = replace.zeros(SAMPLE_NUMBER, med_sample_number)
atsea.data$SAMPLE_NUMBER = tmp.replace[[1]]
Replaced.with.med.sample.number = tmp.replace[[2]]

##############################################################################
#
# Calculate exp_factor
#
##############################################################################

indiv_weight = SAMPLE_WEIGHT / SAMPLE_NUMBER
eff_sample_weight = sum_FREQ * indiv_weight
exp_factor = EXTRAPOLATED_WEIGHT/eff_sample_weight

##############################################################################
#
# Quantiles for exp_factor
#
##############################################################################

pre_exp.factor.summary = summary(exp_factor)

# Reduces the impact of very-very large catches on the assessment

pctl = quantile(exp_factor, in.pctl, na.rm=T)

exp_factor[SAMPLE_NUMBER < 40 & exp_factor > pctl] = pctl

# Deal with sparse samples

if (!remove_sparse) {

  median_exp_factor = median(exp_factor, na.rm=T)
  exp_factor[HAUL_JOIN %in% Too_few_fish] = median_exp_factor
  Replaced.exp.median = length(Too_few_fish)

} # End if

post_exp.factor.summary = summary(exp_factor)

detach(atsea.data)
atsea.data = cbind(atsea.data, exp_factor)
rm(exp_factor)

if (BY_AGE) {

  # Determine the number of age samples taken

  n.age.samples = sum(atsea.data$FREQUENCY)

  # Determine the number of hauls with age information

  attach(atsea.data)

  tmp_agesum = aggregate(atsea.data$AGE, list(atsea.data$HAUL_JOIN), sum, na.rm=T)
  names(tmp_agesum) = c("HAUL_JOIN", "agesum")

  n.hauls.with.age = length(tmp_agesum$agesum[tmp_agesum$agesum != 0])

  agesum = find.matching.rows(atsea.data, tmp_agesum,
                              "HAUL_JOIN", "HAUL_JOIN", "agesum")[[1]]

  # Get the weight of aged samples per haul

  tmp_aged_sample_weight = aggregate(indiv_weight, list(HAUL_JOIN), sum, na.rm=T)
  names(tmp_aged_sample_weight) = c("HAUL_JOIN", "aged_sample_weight")
  aged_sample_weight = find.matching.rows(atsea.data, tmp_aged_sample_weight,
                                 "HAUL_JOIN", "HAUL_JOIN", "aged_sample_weight")[[1]]

  detach(atsea.data)

} # End if BY_AGE

# Clean up

rm(list=ls(pat="tmp*"))

attach(atsea.data)

##############################################################################
#
# Now get medians and sums of various quantities needed to create summary tables.
#
##############################################################################

# Get the median of the eff_sample_weight and indiv_weight, and the sum 
# of FREQUENCY for each HAUL_JOIN

toagg = cbind(indiv_weight, eff_sample_weight, FREQUENCY)

tmp.agg = aggregate(toagg, list(HAUL_JOIN), median, na.rm=T)
colnames(tmp.agg) = c("HAUL_JOIN","MEDIAN_indiv_weight", "MEDIAN_eff_sample_weight", "MEDIAN_FREQ") 

tmp.agg2 = aggregate(toagg, list(HAUL_JOIN), sum, na.rm=T)

tmp.agg = cbind(tmp.agg[,1:3],tmp.agg2[,4])
colnames(tmp.agg)[4] = "SUM_FREQUENCY"

targets = c("MEDIAN_eff_sample_weight", "MEDIAN_indiv_weight", "SUM_FREQUENCY")

Summary_data = find.matching.rows(atsea.data, tmp.agg, "HAUL_JOIN", "HAUL_JOIN", targets)

detach(atsea.data)

atsea.data = cbind(atsea.data, Summary_data, eff_sample_weight)

rm(Summary_data, indiv_weight, eff_sample_weight)

##############################################################################
#
# Record and flag for removal those records for which we couldn't generate
# enough summary data
#
##############################################################################

attach(atsea.data)

atsea.data$MEDIAN_eff_sample_weight[is.na(MEDIAN_eff_sample_weight)] = 0
No.MEDIAN_eff_sample_weight = sum(atsea.data$MEDIAN_eff_sample_weight == 0, na.rm=T)
atsea.data$LENGTH[atsea.data$MEDIAN_eff_sample_weight == 0] = NA
                        
atsea.data$MEDIAN_indiv_weight[is.na(MEDIAN_indiv_weight)] = 0
No.MEDIAN_indiv_weight = sum(atsea.data$MEDIAN_indiv_weight == 0, na.rm=T)
atsea.data$LENGTH[atsea.data$MEDIAN_indiv_weight == 0] = NA
                         
atsea.data$SUM_FREQUENCY[is.na(SUM_FREQUENCY)] = 0
No.FREQUENCY = sum(atsea.data$SUM_FREQUENCY == 0, na.rm=T)
atsea.data$LENGTH[atsea.data$SUM_FREQUENCY == 0] = NA
              
detach(atsea.data)

atsea.data = atsea.data[! is.na(atsea.data$LENGTH), ] 

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

  Use_Length = LENGTH
  Use_Length[Use_Length <= 20] = 20
  Use_Length[Use_Length >= 70] = 70

  is.odd = Use_Length %% 2 > 0
  Use_Length[is.odd] = Use_Length[is.odd] - 1

  # Make sure all lengths show up in the table, even if not
  # in the data

  post_length.summary = summary(Use_Length)

  # Needs to be factored in order to have unused levels
  # show up in the tables.

  Use_Length = factor(Use_Length, levels=seq(20,70,2))

} else {

  Use_Length = rep(0, length(LENGTH))

  # Iteratively replace zeros with next Lbin
  # Assumes Lbins are sorted least-greatest

  for (i in 2:length(Lbins)) {

    minsize = Lbins[i]  
    lesser = Lbins[i-1]

    Use_Length[Use_Length == 0 & LENGTH < minsize] = lesser

  } # End for

  # Finish up largest bin
    
  Use_Length[Use_Length == 0] = Lbins[i]

  post_length.summary = summary(Use_Length * 10)

  # Make sure all lengths show up in the table, even if not
  # in the data

  Use_Length = factor(Use_Length, levels=Lbins)

} # End if

detach(atsea.data)

atsea.data = cbind(atsea.data, Use_Length)
rm(Use_Length)

# Even if we're NOT doing age comps

atsea.data$AGE[atsea.data$AGE < 1] = 1
atsea.data$AGE[atsea.data$AGE > 15] = 15

post_age.summary = summary(atsea.data$AGE)

# Make sure all ages show up in the table, as zeros if not represented
# in the data

atsea.data$AGE = factor(atsea.data$AGE, levels=seq(1,15))

if (!BY_AGE) {

  atsea.data$AGE = addNA(atsea.data$AGE)

} # End if


##############################################################################
#
# Output data
#
##############################################################################

# Size comps by number and percentage

# The order in which factors are added matters!

attach(atsea.data)

by.factors = NULL

if (TEMPORAL) { by.factors =  "SEASON + " }

by.factors = paste(by.factors, "SEX + Use_Length + AGE", sep="")

form = as.formula(paste("exp_factor ~ ", by.factors, sep=""))
fish.form = as.formula(paste("FREQUENCY ~ ", by.factors, sep=""))
hauls.form = as.formula(paste("FREQUENCY ~ ", "HAUL_JOIN + ", by.factors, sep=""))

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

# Collect Weight used for each HAUL_JOIN
# NB:  Use of as.character and as.numeric because R insists on converting
# things to goddamn FACTORS!


SAMPLES = sort(unique(HAUL_JOIN))

Weight = matrix(nrow = length(SAMPLES), ncol=2, 0)
colnames(Weight) = c("W","M_or_S")

for ( i in 1:length(SAMPLES) ) {

  Weight[i,1] = eff_sample_weight[HAUL_JOIN == SAMPLES[i]][1]

  if (TEMPORAL) {

    Weight[i,2] = SEASON[HAUL_JOIN == SAMPLES[i]][1]

  } # End if

} # End for
        
Weight = data.frame(Weight)

if (TEMPORAL) {

    Hauls = aggregate(HAUL_JOIN, list(SEASON), FUN=nfacts)

    Rpt.Weight = aggregate(as.numeric(Weight$W),
                           list(Weight$M_or_S), sum, na.rm=T)

    N.Fish = aggregate(FREQUENCY, list(SEASON), sum, na.rm=T)

    sampled = cbind(Rpt.Weight, Hauls[,2], N.Fish[,2])
    colnames(sampled) = c("Season","Metric.Tonnes","Hauls","N.Fish")

} else {

  Hauls = nfacts(HAUL_JOIN)

  Rpt.Weight = sum(Weight[,1])

  N.Fish = sum(FREQUENCY, na.rm=T)

  # Staple them together

  sampled = c(Rpt.Weight, Hauls, N.Fish)
  dim(sampled) = c(1,3)
  colnames(sampled) = c("Metric.Tonnes", "Hauls", "N.Fish")

} # End if-else

post_haul.summary = summary(atsea.data$indiv_weight)
post_total.wgt.summary = summary(atsea.data$eff_sample_weight)

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

  cat("\nPre-processing observer estimates (in mtonnes):\n")
  print(pre_EXTRAPOLATED_WEIGHT.summary)

  cat("\nPost-processing observer estimates (in mtonnes):\n")
  print(summary(EXTRAPOLATED_WEIGHT))

  cat("\nPre-processing vessel estimates (in kg):\n")
  print(pre_SAMPLE_WEIGHT.summary)

  cat("\nPost-processing vessel estimates (in kg):\n")
  print(summary(SAMPLE_WEIGHT))

  if (BY_AGE) {

    cat("\nPre-processing ages:\n")
    print(pre_age.summary)

    cat("\nPost-processing ages:\n")
    print(post_age.summary)

    cat("\nNumber of hauls with aged fish:  ", n.hauls.with.age, "\n")

    cat("\nAged sample weight\n")
    print(summary(aged_sample_weight))

  } else {

    cat("\nNumber of lengthed fish:  ", sum(atsea.data$FREQUENCY), "\n")

  } # End if-else

  cat("\nPre-truncation expansion factor summary:\n")
  print(pre_exp.factor.summary)

  cat("\nPost-truncation expansion factor summary:\n")
  print(post_exp.factor.summary)

  cat("\nNumber removed:\n")

  if (BY_AGE) {

    cat("\nRemoved  ", No.age, "  records without AGE.\n")

  } # End if

  cat("Removed  ", No.length, 
                       "  records without recorded LENGTH.\n",
                        append=T)

  if (BY_GENDER) { 

    cat("Removed  ", No.gender, "  records without SEX.\n")

  } # End if

  if (TEMPORAL) {

    cat("Removed  ", No.month, "  records without SAMPLE_MONTH.\n")

  } # End if

  if (remove_sparse) {

    cat("Removed  ", Removed.sparse, "  samples without", min_sample, "fish\n")

  } # End if

  cat("\nWeights replaced in remaining data:\n")

  cat("Replaced  ", Replaced.with.med.haul.weight, "  EXTRAPOLATED_WEIGHTs with medians.\n") 
  cat("Replaced  ", Replaced.with.med.sample.weight, "  SAMPLE_WEIGHTs with medians.\n")
  cat("Replaced  ", Replaced.with.med.sample.number, "  SAMPLE_NUMBERs with medians.\n")

  if (!remove_sparse) {

      cat("\nReplaced  ", Replaced.exp.median, "  Length expansion factors with median",
                        "in sparse samples\n")

  } # End if sparse

  cat("\nRecords removed because expansion factors couldn't be generated\n")
                       
  cat("Removed  ", No.MEDIAN_eff_sample_weight, "  records with no MEDIAN_eff_sample_weight\n")
  cat("Removed  ", No.MEDIAN_indiv_weight, "  records with no MEDIAN_indiv_weight\n")
  cat("Removed  ", No.FREQUENCY, "  records with no FREQUENCY\n")

  cat("\nRecords used in final analysis:", N.lengthed.fish, "\n")

sink()
close(outfile)

# Clean up

detach(atsea.data)
options(warn=0)

} # End function old.atsea.comps

#
# That's All, Folks!
#
##############################################################################
