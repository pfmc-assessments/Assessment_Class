##################################################################################################
#
#	PacFIN Data Expansion for Dover sole 2021
# 		
#		Written by Chantel Wetzel
#
##################################################################################################
library(ggplot2)
# remotes::install_github('pfmc-assessments/PacFIN.Utilities')
# remotes::install_github('pfmc-assessments/nwfscSurvey')
library(PacFIN.Utilities)
library(nwfscSurvey)

dir = "C:/Assessments/2023/uw_class/data"
setwd(dir)

# Load in the PacFIN bds data pull and data should be 
# requested via the PacFIN.Utilities Issues using the 
# data request form.
bds_file = "PacFIN.DOVR.bds.12.Feb.2021.RData"
load(file.path(getwd(), "pacfin_bds", bds_file))
out = bds.pacfin 

# The cleanPacFIN function retains records that are
# randomly collected based on sampling protocols, removes
# any age reads that don't align with the keep_age_methods
# (e.g., sets the ages to NA), retains only records in the 
# U.S., and other various data factors
Pdata <- cleanPacFIN(
	Pdata = out,
	keep_age_method = c("B", "BB", 1),
  	CLEAN = TRUE, 
  	verbose = TRUE)

# N SAMPLE_TYPEs changed from M to S for special samples from OR: 0
# N not in keep_sample_type (SAMPLE_TYPE): 57790
# N with SAMPLE_TYPE of NA: 0
# N not in keep_sample_method (SAMPLE_METHOD): 24
# N with SAMPLE_NO of NA: 0
# N without length: 1265
# N without Age: 202318
# N without length and Age: 202350
# N sample weights not available for OR: 24126
# N records: 292378
# N remaining if CLEAN: 227326
# N removed if CLEAN: 65052

# Load in sex specific growth estimates from the survey
# These parameters are used to calculate missing weights
# for fish by sex within the first expansion.
load(file.path(dir, "biology", "growth_estimates_survey.Rdat"))
fa = as.numeric(est_growth$all_F[1])
fb = as.numeric(est_growth$all_F[2]) 
ma = as.numeric(est_growth$all_M[1])
mb = as.numeric(est_growth$all_M[2])
ua = (fa + ma) / 2
ub = (fb + mb) / 2

# Load in the catch file for expansion - the format of this
# catch file may not match your fleet structure in the model.
# Each state should be kept seperate for state based expansions
# for each gear fleet within that state.
catch = read.csv(file.path(dir, "pacfin_catch", "commercial_catch_by_state.csv"))
colnames(catch) = c("Year", "CA_ALL", "OR_ALL", "WA_ALL")

##########################################################################################
# Let do some data checking
##########################################################################################

# Check lengths by sex
# These are just starting values - estimates will be based on sex
k = 0.128; Linf = 47.8; L0 = 10.5; CV1 = 0.20; CV2 = 0.10
Pdata$Length_cm <- Pdata$lengthcm
Pdata$Sex <- Pdata$SEX

Pdata <- nwfscSurvey::est_growth(
	dat = Pdata,
  	Par = data.frame(K = k, Linf = Linf, L0 = L0, CV0 = CV1, CV1 = CV2),
  	sdFactor = 4)

# # See which records are outside 4*sd
remove = which(Pdata[,'lengthcm'] > Pdata[,'Lhat_high'] | Pdata[,'lengthcm'] < Pdata[,'Lhat_low'])
table(Pdata[remove, "lengthcm"], Pdata[remove, "Age"])

colors <- viridis::viridis(3)
f <- which(Pdata$SEX == "F")
plot(Pdata[f,"Age"], Pdata[f,"lengthcm"], type = 'pch', 
	col = colors[1], xlim = c(0, max(Pdata$Age, na.rm = TRUE)),
	ylim = c(0, max(Pdata$lengthcm, na.rm = TRUE)))
m <- which(Pdata$SEX == "M")
points(Pdata[m, "Age"], Pdata[m, "lengthcm"], pch = 1, col = colors[2], lwd = 3)
u <- which(Pdata$SEX == "U")
points(Pdata[u, "Age"], Pdata[u, "lengthcm"], pch = 1, col = colors[3], lwd = 3)
points(Pdata[remove, "Age"], Pdata[remove, "lengthcm"], 
	pch = 2, col = 'red', lwd = 3)

# Set these 144 records to NA for length and age
Pdata[remove, c("lengthcm", "Age")] = NA

quantile(Pdata$lengthcm, na.rm = TRUE) 
quantile(Pdata$age, na.rm = TRUE)
remove = which(Pdata$lengthcm > 70 | Pdata$age > 70)
# 3 total records being removed: 1 length 91 cm and 3 ages 71, 90, and 113
Pdata[remove, c("lengthcm", "Age")] = NA

# This shows only fish with length and age
ggplot(Pdata, aes(x = Age, y = lengthcm)) +
	geom_jitter() + 
	geom_point(aes(col = SEX), size = 2) +
	scale_colour_viridis_d()

####################################################################################################
# Specify fleets, the stratification for expansion, and calculate
# the expansions
####################################################################################################

# Set up the expected fleet structure
table(Pdata$geargroup)

Pdata$geargroup = "ALL"
Pdata$fleet[Pdata$state != "CA"] = paste0("WA_OR", "_", Pdata$geargroup)
Pdata$fleet[Pdata$state == "CA"] = paste0("CA", "_", Pdata$geargroup)

# First stage expansion: expand comps to the trip level
?getExpansion_1
Pdata_exp <- getExpansion_1(
	Pdata = Pdata,
	plot = file.path(dir, "pacfin_bds","plots"),
	fa = fa, fb = fb, ma = ma, mb = mb, ua = ua, ub = ub)

# Second stage expansion: expand comps up to the state and fleet
# The stratification.col input below needs to be the same
# as in the catch csv file
Pdata_exp <- getExpansion_2(
	Pdata = Pdata_exp, 
	Catch = catch, 
	Units = "MT",
  	stratification.cols = c("state", "geargroup"),
  	savedir = file.path(dir, "pacfin_bds", "plots"))

# Calculate the final expansion size
Pdata_exp$Final_Sample_Size <- capValues(Pdata_exp$Expansion_Factor_1_L * Pdata_exp$Expansion_Factor_2)

length_comps <- getComps(
	Pdata = Pdata_exp[!is.na(Pdata_exp$lengthcm), ], 
	Comps = "LEN")

table(Pdata$SOURCE_AGID, Pdata$SEX)

Pdata$count <- 1
ggplot(Pdata, aes(x = lengthcm, y = count, fill = SEX))  + 
	geom_histogram(aes(y = count), position="stack", stat="identity") +
	scale_fill_viridis_d()

####################################################################################################
# Create the length composition data
####################################################################################################

# Commenting out for now because I don't want to assign unsexed 
# due to the dimorphic growth
# There area a fair number of U in CA and in the early years of WA
# length_compSR <- doSexRatio(
#	CompData = length_comps, 
# 	ratioU = 0.5, 
# 	maxsizeU = 25, 
# 	savedir = file.path(dir, "commercial_comps"))

len_bins = c(seq(8, 60, 2))
out_name = sub(pattern = "(.*)\\..*$", replacement = "\\1", bds_file)

writeComps(
	inComps = length_comps, 
	fname = file.path(dir, "pacfin_bds", "forSS", paste0("Lengths_", out_name, ".csv")), 
	lbins = len_bins, 
	sum1 = TRUE, 
	partition = 2, 
	digits = 4,
	dummybins = FALSE)

##########################################################
# Calculate the expansion for age data
##########################################################

Pdata_exp$Final_Sample_Size = capValues(Pdata_exp$Expansion_Factor_1_A * Pdata_exp$Expansion_Factor_2)

age_comps <- getComps(
	Pdata_exp[!is.na(Pdata_exp$Age), ], 
	Comps = "Age")

##########################################################
# Create the age compositions
##########################################################

age_bins = 1:60

writeComps(
	inComps = age_comps, 
	fname = file.path(dir, "pacfin_bds", "forSS", paste0("Age_", out_name, ".csv")), 
	abins = age_bins, 
	sum1 = TRUE, 
	partition = 2, 
	digits = 4,
	dummybins = FALSE)

# Create condition-age-at-length compositions just in case
# you want to explore them
caal_comps <- getComps(
	Pdata = Pdata_exp[!is.na(Pdata_exp$age), ], 
	Comps = "AAL")

writeComps(
	inComps = caal_comps, 
	fname = file.path(dir, "pacfin_bds", "forSS", paste0("CAAL_", out_name, ".csv")), 
	lbins = len_bins, 
	abins = age_bins, 
	sum1 = TRUE, 
	partition = 2, 
	dummybins = FALSE)

###############################################################################################
# Let's format the csv files for direct use in SS3
#####################################################################################

out = read.csv(
	file.path(dir, "pacfin_bds", "forSS", paste0("Lengths_", out_name, ".csv")), 
	skip = 3, 
	header = TRUE)

start = 1 
end   = which(as.character(out[,1]) %in% c(" Females only ")) - 1 
cut_out = out[start:end,]

# format the length comps
cut_out$fleet[cut_out$fleet =="CA_ALL"] = 1
cut_out$fleet[cut_out$fleet =="WA_OR_ALL"] = 2

ind = which(colnames(cut_out) %in% paste0("F", min(len_bins))):
	  which(colnames(cut_out) %in% paste0("M", max(len_bins)))
format = cbind(cut_out$year, cut_out$month, cut_out$fleet, cut_out$sex, cut_out$partition, cut_out$InputN, cut_out[,ind])
colnames(format) = c("year", "month", "fleet", "sex", "part", "InputN", colnames(cut_out[ind]))
format = format[format$year != 2021, ]
write.csv(
	format, 
	file = file.path(dir, "pacfin_bds", "forSS", paste0("Lcomps_for_SS3_", out_name, ".csv")), 
	row.names = FALSE)


# Let's create the sample table
temp = Pdata[!is.na(Pdata$lengthcm) & Pdata$year < 2021,]
Nfish = table(temp$year, temp$state)
colnames(Nfish) = sort(unique(temp$state)) 

aa = sort(unique(temp$state)) 
yy = sort(unique(temp$year))
Ntows = matrix(0, length(yy), length(aa))
for(y in 1:length(yy)){
	for(a in 1:length(aa)){
		ind = which(temp$year == yy[y] & temp$state  == aa[a])
		if(length(ind) > 0) {Ntows[y, a] = length(unique(temp$SAMPLE_NO[ind])) }
	}
}
colnames(Ntows) = aa
rownames(Ntows) = yy

samples = cbind(rownames(Ntows), Ntows[,"CA"] , Nfish[,"CA"], 
			    Ntows[,"OR"] , Nfish[,"OR"],  Ntows[,"WA"] , Nfish[,"WA"])

colnames(samples) = c("Year", "CA Ntows", "CA Nfish",
	 				  "OR Ntows", "OR Nfish", "WA Ntows", "WA Nfish")
write.csv(
	samples, 
	file = file.path(dir, "pacfin_bds", "forSS", paste0("PacFIN_Length_Samples_by_State.csv")), 
	row.names = FALSE)

###############################################################################################
# Format Age Samples
###############################################################################################

out = read.csv(file.path(dir, "pacfin_bds", "forSS", paste0("Age_",out_name, ".csv")), 
				skip = 3, header = TRUE)
start = 1 
end   = which(as.character(out[,1]) %in% c(" Females only ")) - 1 
cut_out = out[start:end,]

# format the age comps
cut_out$fleet[cut_out$fleet =="CA_ALL"] = 1
cut_out$fleet[cut_out$fleet =="WA_OR_ALL"] = 2
cut_out$ageErr[cut_out$fleet == 1] = 2
cut_out$ageErr[cut_out$fleet == 2] = 1

ind = which(colnames(cut_out) %in% paste0("F",min(age_bins))):
	  which(colnames(cut_out) %in% paste0("M",max(age_bins)))
format = cbind(cut_out$year, cut_out$month, cut_out$fleet, cut_out$sex, cut_out$partition, 
			   cut_out$ageErr, cut_out$LbinLo, cut_out$LbinHi, cut_out$InputN, cut_out[,ind])
colnames(format) = c("year", "month", "fleet", "sex", "part", "AgeErr", "Low", "High", "InputN", colnames(cut_out[ind]))
format = format[format$year != 2021, ]
write.csv(
	format, 
	file = file.path(dir, "pacfin_bds", "forSS", paste0("Acomps_forSS3_", out_name, ".csv")), 
	row.names = FALSE)

# Let's create the sample table
temp = Pdata[!is.na(Pdata$age) & Pdata$year < 2021,]
Nfish = table(temp$year, temp$state)
colnames(Nfish) = unique(temp$state) 

aa = sort(unique(temp$state))
yy = sort(unique(temp$year))
Ntows = matrix(0, length(yy), length(aa))
for(y in 1:length(yy)){
	for(a in 1:length(aa)){
		ind = which(temp$year == yy[y] & temp$state  == aa[a])
		if(length(ind) > 0) {Ntows[y, a] = length(unique(temp$SAMPLE_NO[ind])) }
	}
}
colnames(Ntows) = aa
rownames(Ntows) = yy

samples = cbind(rownames(Ntows), Ntows[,"CA"] , Nfish[,"CA"], 
			    Ntows[,"OR"] , Nfish[,"OR"], Ntows[,"WA"] , Nfish[,"WA"])

colnames(samples) = c("Year", "CA Ntows", "CA Nfish",
	 				  "OR Ntows", "OR Nfish", "WA Ntows", "WA Nfish")
write.csv(
	samples, 
	file = file.path(dir, "pacfin_bds", "forSS", paste0("Com_Age_Samples_by_State.csv")), 
	row.names = FALSE)
######################################################################################################
######################################################################################################