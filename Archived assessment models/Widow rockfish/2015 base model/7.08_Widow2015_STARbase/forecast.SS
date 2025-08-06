# 2015 Widow Rockfish Assessment
# Allan C. Hicks and Chantel R. Wetzel
# NWFSC, NOAA, Seattle, WA

# for all year entries except rebuilder; enter either: actual year, -999 for styr, 0 for endyr, neg number for rel. endyr
1     # Benchmarks: 0=skip; 1=F(SPR); 2=F(MSY);3=F(Btarget); 4=F(endyr); 5=Ave recent F (not implemented); 6= read Fmult (not implemented)
2     # MSY: 1= set to F(SPR); 2=calc F(MSY); 3=set to F(Btgt); 4=set to F(endyr)
0.50  # SPR target (e.g. 0.40), 0.5 for west coast groundfish
0.40  # Biomass target (e.g. 0.40)
#_Bmark_years: beg_bio, end_bio, beg_selex, end_selex, beg_relF, end_relF (enter actual year, or values of 0 or -integer to be rel. endyr)
0 0 0 0 0 0

1   #Bmark_relF_Basis: 1 = use year range; 2 = set relF same as forecast below

1   # Forecast: 0=none; 1=F(SPR); 2=F(MSY) 3=F(Btgt); 4=Ave F (uses first-last relF yrs); 5=input annual F scalar
12  # N forecast years
1   # F scalar (only used for Do_Forecast==5)
#_Fcast_years:  beg_selex, end_selex, beg_relF, end_relF  (enter actual year, or values of 0 or -integer to be rel. endyr)
0 0 0 0

1       # Control rule method (1=catch=f(SSB) west coast; 2=F=f(SSB) )
0.40    # Control rule Biomass level for constant F (as frac of Bzero, e.g. 0.40)
0.10    # Control rule Biomass level for no F (as frac of Bzero, e.g. 0.10)
0.9548  # Control rule target as fraction of Flimit (e.g. 0.75)

3       #_N forecast loops (1-3) (fixed at 3 for now)
3       #_First forecast loop with stochastic recruitment
0       #_Forecast loop control #3 (reserved for future bells&whistles)
0       #_Forecast loop control #4 (reserved for future bells&whistles)
0       #_Forecast loop control #5 (reserved for future bells&whistles)

2027    #FirstYear for caps and allocations (should be after years with fixed inputs)
0       # stddev of log(realized catch/target catch) in forecast (set value>0.0 to cause active impl_error) (if=0, there will be N_forecase_years less parameters estimated)
0       # Do West Coast gfish rebuilder output (0/1)
2015    # Rebuilder:  first year catch could have been set to zero (Ydecl)(-1 to set to 1999)
2015    # Rebuilder:  year for current age structure (Yinit) (-1 to set to endyear+1)
1       # fleet relative F:  1=use first-last alloc year; 2=read seas(row) x fleet(col) below

# Note that fleet allocation is used directly as average F if Do_Forecast=4
2       # basis for fcast catch tuning and for fcast catch caps and allocation  (2=deadbio; 3=retainbio; 5=deadnum; 6=retainnum)
# Conditional input if relative F choice = 2
# Fleet relative F:  rows are seasons, columns are fleets
#_Fleet:  FISHERY
#  1
# max totalcatch by fleet (-1 to have no max) must enter value for each fleet
-1 -1 -1 -1 -1

# max totalcatch by area (-1 to have no max); must enter value for each fleet
-1

# fleet assignment to allocation group (enter group ID# for each fleet, 0 for not included in an alloc group)
0 0 0 0 0

#_Conditional on >1 allocation group
# allocation fraction for each of: 0 allocation groups
# no allocation groups
60      # Number of forecast catch levels to input (else calc catch from forecast F)
2       # basis for input Fcast catch:  2=dead catch; 3=retained catch; 99=input Hrate(F) (units are from fleetunits; note new codes in SSV3.20)
# Input fixed catch values
#Year Seas Fleet Catch(or_F)
2015	1	1	197.7855179
2015	1	2	850.7996
2015	1	3	949.4148821
2015	1	4	0
2015	1	5	2
2016	1	1	197.7855179
2016	1	2	850.7996
2016	1	3	949.4148821
2016	1	4	0
2016	1	5	2
2017	1	1	197.7855179
2017	1	2	850.7996
2017	1	3	949.4148821
2017	1	4	0
2017	1	5	2
2018	1	1	197.7855179
2018	1	2	850.7996
2018	1	3	949.4148821
2018	1	4	0
2018	1	5	2
2019	1	1	197.7855179
2019	1	2	850.7996
2019	1	3	949.4148821
2019	1	4	0
2019	1	5	2
2020	1	1	197.7855179
2020	1	2	850.7996
2020	1	3	949.4148821
2020	1	4	0
2020	1	5	2
2021	1	1	197.7855179
2021	1	2	850.7996
2021	1	3	949.4148821
2021	1	4	0
2021	1	5	2
2022	1	1	197.7855179
2022	1	2	850.7996
2022	1	3	949.4148821
2022	1	4	0
2022	1	5	2
2023	1	1	197.7855179
2023	1	2	850.7996
2023	1	3	949.4148821
2023	1	4	0
2023	1	5	2
2024	1	1	197.7855179
2024	1	2	850.7996
2024	1	3	949.4148821
2024	1	4	0
2024	1	5	2
2025	1	1	197.7855179
2025	1	2	850.7996
2025	1	3	949.4148821
2025	1	4	0
2025	1	5	2
2026	1	1	197.7855179
2026	1	2	850.7996
2026	1	3	949.4148821
2026	1	4	0
2026	1	5	2


999 # verify end of input
