# 2015 Widow Rockfish Assessment
# Allan C. Hicks and Chantel R. Wetzel
# NWFSC, NOAA, Seattle, WA

2015widow.dat
2015widow.ctl

1		# 0=use init values in control file; 1=use ss3.par
1		# run display detail (0,1,2)
1		# detailed age-structured reports in SS2.rep (0,1)
0		# write detailed checkup.sso file (0,1)
0		# write parm values to ParmTrace.sso (0=no,1=good,active; 2=good,all; 3=every_iter,all_parms)
0		# write to cumreport.sso (0=no,1=like&timeseries; 2=add survey fits)
1		# Include prior_like for non-estimated parameters (0,1)
1		# Use Soft Boundaries to aid convergence (0,1) (recommended)
0		# Number of bootstrap datafiles to produce
10		# Turn off estimation for parameters entering after this phase
1		# MCMC burn interval
1		# MCMC thin interval
0.0005		# jitter initial parm value by this fraction
-1		# begin annual SD report in start year
-2		# end annual SD report in end year (-2=end of annual SD report in last forecast year
0		# N individual STD years (0=none)
#vector of year values
0.0001	# final convergence criteria (e.g. 1.0e-04)
0		# retrospective year relative to end year (e.g. -4)
4		# min age for calc of summary biomass (age 4 is when all fleets significantly begin seeing fish)
1		# Depletion basis:  denom is: 0=skip; 1=rel X*B0; 2=rel X*Bmsy; 3=rel X*B_styr
1		# Fraction (X) for Depletion denominator (e.g. 0.4)
1		# (1-SPR)_reporting:  0=skip; 1=rel(1-SPR); 2=rel(1-SPR_MSY); 3=rel(1-SPR_Btarget); 4=no denominator (report actural 1-SPR values)
1		# F_report_units: 0=skip; 1=exploitation(Bio); 2=exploitation(Num); 3=sum(Frates); 4=true F for range of ages
0		# F_report_basis: 0=raw; 1=rel Fspr; 2=rel Fmsy ; 3=rel Fbtgt
999		# check value for end of file
