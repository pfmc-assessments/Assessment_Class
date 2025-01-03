 #V3.30.12.00-trans;_2018_08_01;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_11.6
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#_data_and_control_files: 2015widow.dat // 2015widow.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  2=main effects for GP, Settle timing, Area; 3=each Settle entity; 4=none, only when N_GP*Nsettle*pop==1
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
10 #_Nblock_Patterns
 3 2 1 1 1 1 3 1 1 1#_blocks_per_pattern 
# begin and end years of blocks
 1982 1989 1990 1997 1998 2010
 1982 1989 1990 2010
 1916 1982
 1916 2001
 1916 2002
 1995 2012
 1916 1982 1983 2001 2002 2010
 1915 1915
 1995 2004
 1991 1998
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#  autogen
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
# 
#
#_Available timevary codes
#_Block types: 0: Pblock=Pbase*exp(TVP); 1: Pblock=Pbase+TVP; 2: Pblock=TVP; 3: Pblock=Pblock(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=Pbase*exp(TVP*env(y));  2: P(y)=Pbase+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=env(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=not implemented
3 #_Age(post-settlement)_for_L1;linear growth below this
40 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
2 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
3 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
 0.01 0.3 0.10 -2.30 0.438 3 5 0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
 10 40 27.4948 27 99 0 3 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 35 60 50.0042 50 99 0 2 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.01 0.4 0.150077 0.15 99 0 2 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.01 0.4 0.0705642 0.07 99 0 3 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.01 0.4 0.041775 0.04 99 0 3 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
 -3 3 1.736e-05 0 99 0 -99 0 0 0 0 0 0 0 # Wtlen_1_Fem
 -3 10 2.962 2.962 99 0 -99 0 0 0 0 0 0 0 # Wtlen_2_Fem
 -3 50 5.47 7 99 0 -99 0 0 0 0 0 0 0 # Mat50%_Fem
 -3 3 -0.7747 -1 99 0 -99 0 0 0 0 0 0 0 # Mat_slope_Fem
 -1 1 1 1 99 0 -99 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem
 0 1 0 0 99 0 -99 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem
 0.01 0.3 0.10 -2.30 0.438 3 5 0 0 0 0 0 0 0 # NatM_p_1_Mal_GP_1
 10 40 26.0012 27 99 0 3 0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
 35 60 44.0029 45 99 0 2 0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
 0.01 0.4 0.210064 0.19 99 0 2 0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1
 0.01 0.4 0.0701206 0.07 99 0 3 0 0 0 0 0 0 0 # CV_young_Mal_GP_1
 0.01 0.4 0.0401227 0.04 99 0 3 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
 -3 3 1.484e-05 0 99 0 -99 0 0 0 0 0 0 0 # Wtlen_1_Mal
 -3 10 3.005 3.005 99 0 -99 0 0 0 0 0 0 0 # Wtlen_2_Mal
 0 2 1 1 99 0 -99 0 0 0 0 0 0 0 # RecrDist_GP_1
 0 2 1 1 99 0 -99 0 0 0 0 0 0 0 # RecrDist_Area_1
 0 2 1 1 99 0 -99 0 0 0 0 0 0 0 # RecrDist_timing_1
 0 2 1 1 99 0 -99 0 0 0 0 0 0 0 # CohortGrowDev
 0.000001 0.999999 0.5 0.5  0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             1            20       14.5006            10            99             0          2          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1         0.720         0.720         0.160             2         -5          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.6          0.65            99             0        -50          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             0        -99          0          0          0          0          0          0          0 # SR_regime
             0           0.5             0             0            99             0        -99          0          0          0          0          0          0          0 # SR_autocorr
1 #do_recdev:  0=none; 1=devvector; 2=simple deviations
1970 # first year of main recr_devs; early devs can preceed this era
2014 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase 
1 # (0/1) to read 13 advanced options
 1900 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 4 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1960 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1976 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2014 #_last_yr_fullbias_adj_in_MPD
 2017 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.85 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1900E 1901E 1902E 1903E 1904E 1905E 1906E 1907E 1908E 1909E 1910E 1911E 1912E 1913E 1914E 1915E 1916E 1917E 1918E 1919E 1920E 1921E 1922E 1923E 1924E 1925E 1926E 1927E 1928E 1929E 1930E 1931E 1932E 1933E 1934E 1935E 1936E 1937E 1938E 1939E 1940E 1941E 1942E 1943E 1944E 1945E 1946E 1947E 1948E 1949E 1950E 1951E 1952E 1953E 1954E 1955E 1956E 1957E 1958E 1959E 1960E 1961E 1962E 1963E 1964E 1965E 1966E 1967E 1968E 1969E 1970R 1971R 1972R 1973R 1974R 1975R 1976R 1977R 1978R 1979R 1980R 1981R 1982R 1983R 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011F 2012F 2013F 2014F 2015F 2016F 2017F 2018F 2019F 2020F 2021F 2022F 2023F 2024F 2025F 2026F
#  -0.00197699 0.00024231 -0.000610518 0.0009343 -0.000710894 0.00178719 0.000914743 0.000705787 0.000642114 -0.00190445 -0.0013381 0.00161387 0.0004615 -0.00171831 -0.00189711 0.00108002 0.00120945 0.000248233 -0.0019523 0.00197832 -0.000518339 -0.000714241 -0.000581909 -0.000712476 -0.00141057 0.000480737 0.00127825 -0.000387643 0.000714036 -0.000736057 -0.000613424 0.000875004 -0.00126539 -0.00022268 -0.00108961 -0.00190983 -0.000739032 0.00127379 4.57676e-05 -4.84407e-05 -0.000110513 0.00182372 -0.000181807 0.00174012 -2.80158e-05 -0.000848033 0.00094903 -0.00108725 -0.00152107 -0.00194648 -0.00166557 0.000493858 -0.000518935 0.0003545 0.000487365 -0.00062332 -0.00151522 0.000152153 -4.84485e-05 0.000156248 -5.67789e-05 -0.00157385 0.00114892 -0.000132076 0.00162723 -0.000888295 0.000817169 -0.000435262 0.000714717 0.000673137 0.00134954 6.75986e-05 -0.00162827 0.00133781 -0.00145496 -0.00137183 0.000665513 -0.000884725 -0.000407356 -0.00101511 0.000629278 -0.00159045 0.000781335 -0.00081777 0.00110797 -0.00161295 0.00192668 -0.000883153 0.000275409 -0.00111587 -0.00119724 0.000313546 -0.00107193 0.00076103 -0.000138742 -0.000308078 0.000684888 -0.00114981 1.97276e-05 -0.000342094 0.00143117 0.000541112 8.72455e-05 -0.00141298 0.000583378 0.000300142 0.000227363 0.00227038 0.00215869 0.00154362 -0.000660116 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.05 # F ballpark
-1982 # F ballpark year (neg value to disable)
1 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.9 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#2026 2035
# F rates by fleet
# Yr:  1916 1917 1918 1919 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024 2025 2026
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# BottomTrawl 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2.13771e-05 2.29696e-05 2.46785e-05 2.6512e-05 2.84786e-05 3.05862e-05 3.28416e-05 3.5251e-05 3.78206e-05 4.05568e-05 4.34658e-05 4.65533e-05
# MidwaterTrawl 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# Hake 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# Net 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# HnL 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         1         1         0         1         0         1  #  BottomTrawl
         3         1         0         1         1         0  #  Hake
         6         1         0         1         0         1  #  JuvSurvey
         7         1         0         1         1         0  #  Triennial
         8         1         0         1         0         1  #  NWFSC
         9         1         0         1         0         1  #  ForeignAtSea
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -25            25      -10.8418             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_BottomTrawl(1)
             0             2         0.002             0            99             0          2          0          0          0          0          0          0          0  #  Q_extraSD_BottomTrawl(1)
           -20             2      -10.0003             0            99             0          1          0          0          0          0          0         10         1  #  LnQ_base_Hake(3)
             0             2         0.002             0            99             0          2          0          0          0          0          0          0          0  #  Q_extraSD_Hake(3)
           -25            25      -8.38193             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_JuvSurvey(6)
             0             2         0.002             0            99             0          2          0          0          0          0          0          0          0  #  Q_extraSD_JuvSurvey(6)
            -4             4             0             0            99             0          2          0          0          0          0          0          9          1  #  LnQ_base_Triennial(7)
             0             2             0             0            99             0         -2          0          0          0          0          0          0          0  #  Q_extraSD_Triennial(7)
           -25            25      -6.64111             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_NWFSC(8)
             0             2             0             0            99             0         -2          0          0          0          0          0          0          0  #  Q_extraSD_NWFSC(8)
           -25            25      -15.9788             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_ForeignAtSea(9)
             0             2    0.00167016             0            99             0          2          0          0          0          0          0          0          0  #  Q_extraSD_ForeignAtSea(9)
# timevary Q parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type     PHASE  #  parm_name
        0.0001             2           0.5           0.5           0.5             6       3  # LnQ_base_Hake(3)_Block10_1991-1998
#        0.0001             2            99            99           0.5             6      -5  # LnQ_base_Hake(3)_dev_se
#         -0.99          0.99             0             0           0.5             6      -6  # LnQ_base_Hake(3)_dev_autocorr
        0.0001             2           0.5           0.5           0.5             6       3  # LnQ_base_Triennial(7)_Block9_1995-2004
#        0.0001             2            99            99           0.5             6      -5  # LnQ_base_Triennial(7)_dev_se
#         -0.99          0.99             0             0           0.5             6      -6  # LnQ_base_Triennial(7)_dev_autocorr
# info on dev vectors created for Q parms are reported with other devs after tag parameter section 
#
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline 
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 24 1 0 0 # 1 BottomTrawl
 24 1 0 0 # 2 MidwaterTrawl
 24 0 0 0 # 3 Hake
 24 0 0 0 # 4 Net
 24 1 0 0 # 5 HnL
 0 0 0 0 # 6 JuvSurvey
 27 0 0 3 # 7 Triennial
 27 0 0 3 # 8 NWFSC
 5 0 0 3 # 9 ForeignAtSea
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+nages+1; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 10 0 0 0 # 1 BottomTrawl
 10 0 0 0 # 2 MidwaterTrawl
 10 0 0 0 # 3 Hake
 10 0 0 0 # 4 Net
 10 0 0 0 # 5 HnL
 11 0 0 0 # 6 JuvSurvey
 10 0 0 0 # 7 Triennial
 11 0 0 0 # 8 NWFSC
 10 0 0 0 # 9 ForeignAtSea
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   BottomTrawl LenSelex
            10            59        37.995            45          0.05             0          1          0          0          0          0        0.5          4          2  #  SizeSel_P1_BottomTrawl(1)
            -5            10        2.4987             5          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P2_BottomTrawl(1)
            -4            12       3.99818             3          0.05             0          2          0          0          0          0        0.5          4          2  #  SizeSel_P3_BottomTrawl(1)
            -2            10             9            10          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P4_BottomTrawl(1)
            -9            10            -9           0.5          0.05             0         -3          0          0          0          0        0.5          0          0  #  SizeSel_P5_BottomTrawl(1)
            -9             9             8           0.5          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P6_BottomTrawl(1)
            -5            60      -2.93525             0            99             0          4          0          0          0          0          0          2          2  #  Retain_P1_BottomTrawl(1)
          0.01             8       1.20169             1            99             0          4          0          0          0          0          0          2          2  #  Retain_P2_BottomTrawl(1)
           -10            10       4.59512            10            99             0         -2          0          0          0          0          0          1          2  #  Retain_P3_BottomTrawl(1)
           -10            10             0             0            99             0        -99          0          0          0          0          0          0          0  #  Retain_P4_BottomTrawl(1)
# 2   MidwaterTrawl LenSelex
            10            59       37.9964            45          0.05             0          1          0          0          0          0        0.5          7          2  #  SizeSel_P1_MidwaterTrawl(2)
           -10            10    0.00388996             5          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P2_MidwaterTrawl(2)
            -4            12       2.99998             3          0.05             0          2          0          0          0          0        0.5          7          2  #  SizeSel_P3_MidwaterTrawl(2)
            -2            10       8.99092            10          0.05             0          4          0          0          0          0        0.5          7          2  #  SizeSel_P4_MidwaterTrawl(2)
            -9            10            -9           0.5          0.05             0         -3          0          0          0          0        0.5          0          0  #  SizeSel_P5_MidwaterTrawl(2)
            -9             9        7.9441           0.5          0.05             0          4          0          0          0          0        0.5          7          2  #  SizeSel_P6_MidwaterTrawl(2)
            -5            60            -5             0            99             0         -9          0          0          0          0          0          0          0  #  Retain_P1_MidwaterTrawl(2)
          0.01             8           1.2             1            99             0         -9          0          0          0          0          0          0          0  #  Retain_P2_MidwaterTrawl(2)
           -10            10       4.59512            10            99             0         -2          0          0          0          0          0          7          2  #  Retain_P3_MidwaterTrawl(2)
           -10            10             0             0            99             0        -99          0          0          0          0          0          0          0  #  Retain_P4_MidwaterTrawl(2)
# 3   Hake LenSelex
            10            59       39.9992            45          0.05             0          1          0          0          0          0        0.5          0          0  #  SizeSel_P1_Hake(3)
            -5            10       2.50126             5          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P2_Hake(3)
            -4            12       4.00279             3          0.05             0          2          0          0          0          0        0.5          0          0  #  SizeSel_P3_Hake(3)
            -2            10             9            10          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P4_Hake(3)
            -9            10            -9           0.5          0.05             0         -3          0          0          0          0        0.5          0          0  #  SizeSel_P5_Hake(3)
            -9             9             8           0.5          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P6_Hake(3)
# 4   Net LenSelex
            10            59       40.0001            45          0.05             0          1          0          0          0          0        0.5          0          0  #  SizeSel_P1_Net(4)
            -5            10       2.49866             5          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P2_Net(4)
            -4            12       4.00148             3          0.05             0          2          0          0          0          0        0.5          0          0  #  SizeSel_P3_Net(4)
            -2            10             9            10          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P4_Net(4)
            -9            10            -9           0.5          0.05             0         -3          0          0          0          0        0.5          0          0  #  SizeSel_P5_Net(4)
            -9             9             8           0.5          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P6_Net(4)
# 5   HnL LenSelex
            10            59       25            45          0.05             0          5          0          0          0          0        0.5          5          2  #  SizeSel_P1_HnL(5)
            -5            10        2.5005             5          0.05             0          3          0          0          0          0        0.5          0          0  #  SizeSel_P2_HnL(5)
            -5            12       4.00069             3          0.05             0          2          0          0          0          0        0.5          5          2  #  SizeSel_P3_HnL(5)
            -2            10             9            10          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P4_HnL(5)
            -9            10            -9           0.5          0.05             0         -3          0          0          0          0        0.5          0          0  #  SizeSel_P5_HnL(5)
            -9             9             8           0.5          0.05             0         -4          0          0          0          0        0.5          0          0  #  SizeSel_P6_HnL(5)
            -5            60       25.0099             0            99             0          2          0          0          0          0          0          3          2  #  Retain_P1_HnL(5)
          0.01             8      0.991006             1            99             0          3          0          0          0          0          0          3          2  #  Retain_P2_HnL(5)
           -10            10       2.19741            10            99             0          1          0          0          0          0          0          3          2  #  Retain_P3_HnL(5)
           -10            10             0             0            99             0        -99          0          0          0          0          0          0          0  #  Retain_P4_HnL(5)
# 6   JuvSurvey LenSelex
# 7   Triennial LenSelex
             0             2             0             0             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Code_Triennial(7)
        -0.001             1      0.148975             0             0             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_GradLo_Triennial(7)
            -1             1    -0.0300079             0             0             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_GradHi_Triennial(7)
             8            56            24           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_1_Triennial(7)
             8            56            34           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_2_Triennial(7)
             8            56            48           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_3_Triennial(7)
           -10            10      -3.00454           -10            99             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_Val_1_Triennial(7)
           -10            10            -1           -10            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Val_2_Triennial(7)
           -10            10    0.00205177           -10            99             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_Val_3_Triennial(7)
# 8   NWFSC LenSelex
             0             2             0             0             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Code_NWFSC(8)
        -0.001             1      0.150832             0             0             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_GradLo_NWFSC(8)
            -1             1    -0.0302647             0             0             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_GradHi_NWFSC(8)
             8            56            24           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_1_NWFSC(8)
             8            56            34           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_2_NWFSC(8)
             8            56            48           -10             0             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Knot_3_NWFSC(8)
           -10            10      -2.99769           -10            99             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_Val_1_NWFSC(8)
           -10            10            -1           -10            99             0        -99          0          0          0          0        0.5          0          0  #  SizeSpline_Val_2_NWFSC(8)
           -10            10    0.00335515           -10            99             0          2          0          0          0          0        0.5          0          0  #  SizeSpline_Val_3_NWFSC(8)
# 9   ForeignAtSea LenSelex
            -2            60             0             0           0.2             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P1_ForeignAtSea(9)
            -2            60             0             0           0.2             0        -99          0          0          0          0        0.5          0          0  #  SizeSel_P2_ForeignAtSea(9)
# 1   BottomTrawl AgeSelex
# 2   MidwaterTrawl AgeSelex
# 3   Hake AgeSelex
# 4   Net AgeSelex
# 5   HnL AgeSelex
# 6   JuvSurvey AgeSelex
             0             1             0             0            99             0        -99          0          0          0          0        0.5          0          0  #  AgeSel_P1_JuvSurvey(6)
             0             1             0             0            99             0        -99          0          0          0          0        0.5          0          0  #  AgeSel_P2_JuvSurvey(6)
# 7   Triennial AgeSelex
# 8   NWFSC AgeSelex
             0             1             0             0            99             0        -99          0          0          0          0        0.5          0          0  #  AgeSel_P1_NWFSC(8)
             0            50            40             0            99             0        -99          0          0          0          0        0.5          0          0  #  AgeSel_P2_NWFSC(8)
# 9   ForeignAtSea AgeSelex
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            10            59       34.0094            45          0.05             0      1  # SizeSel_P1_BottomTrawl(1)_BLK4repl_1916
            -4            12       5.70289             3          0.05             0      2  # SizeSel_P3_BottomTrawl(1)_BLK4repl_1916
            -5            50       34.9849            34            99             0      3  # Retain_P1_BottomTrawl(1)_BLK2repl_1982
            -5            50       35.0169            34            99             0      3  # Retain_P1_BottomTrawl(1)_BLK2repl_1990
          0.01             5       2.50005             1            99             0      3  # Retain_P2_BottomTrawl(1)_BLK2repl_1982
          0.01             5       2.49935             1            99             0      3  # Retain_P2_BottomTrawl(1)_BLK2repl_1990
           -10 		  10 	   4.59512 	      10            99             0      2  # Retain_P3_BottomTrawl(1)_BLK1repl_1982
           -10 		  10       4.59512 	      10            99             0      2  # Retain_P3_BottomTrawl(1)_BLK1repl_1990
           -10 		  10 	   4.59512 	      10            99             0      2  # Retain_P3_BottomTrawl(1)_BLK1repl_1998
            10            59       38.0042            45          0.05             0      1  # SizeSel_P1_MidwaterTrawl(2)_BLK7repl_1916
            10            59       37.9976            45          0.05             0      1  # SizeSel_P1_MidwaterTrawl(2)_BLK7repl_1983
            10            59       38.0034            45          0.05             0      1  # SizeSel_P1_MidwaterTrawl(2)_BLK7repl_2002
            -4            12       3.00242             3          0.05             0      2  # SizeSel_P3_MidwaterTrawl(2)_BLK7repl_1916
            -4            12        2.9987             3          0.05             0      2  # SizeSel_P3_MidwaterTrawl(2)_BLK7repl_1983
            -4            12       3.00111             3          0.05             0      2  # SizeSel_P3_MidwaterTrawl(2)_BLK7repl_2002
            -2            10       9.01133            10          0.05             0      4  # SizeSel_P4_MidwaterTrawl(2)_BLK7repl_1916
            -2            10       9.00781            10          0.05             0      4  # SizeSel_P4_MidwaterTrawl(2)_BLK7repl_1983
            -2            10       8.98781            10          0.05             0      4  # SizeSel_P4_MidwaterTrawl(2)_BLK7repl_2002
            -9             9       7.88418           0.5          0.05             0      4  # SizeSel_P6_MidwaterTrawl(2)_BLK7repl_1916
            -9             9       7.91362           0.5          0.05             0      4  # SizeSel_P6_MidwaterTrawl(2)_BLK7repl_1983
            -9             9       7.91958           0.5          0.05             0      4  # SizeSel_P6_MidwaterTrawl(2)_BLK7repl_2002
           -10 		  10 	   4.5912 	      10            99             0     -2  # Retain_P3_MidwaterTrawl(2)_BLK7repl_1916
           -10 		  10 	   4.59512 	      10            99             0      2  # Retain_P3_MidwaterTrawl(2)_BLK7repl_1983
           -10 		  10 	   4.59512 	      10            99             0      2  # Retain_P3_MidwaterTrawl(2)_BLK7repl_2002
            15            59       48.0017            45          0.05             0      1  # SizeSel_P1_HnL(5)_BLK5repl_1916
            -4            12       2.80146             3          0.05             0      2  # SizeSel_P3_HnL(5)_BLK5repl_1916
            -5            50            -5            34            99             0     -2  # Retain_P1_HnL(5)_BLK3repl_1916
           0.1             8           1.2             1            99             0     -3  # Retain_P2_HnL(5)_BLK3repl_1916
           -10 		  10 	   4.5912 	      10            99             0     -3  # Retain_P3_HnL(5)_BLK3repl_1916
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      3     3     1     0     0     0     0     1     1  1983  1998     5 0.000181693 0.00180972 0.001629 -0.00147081 -0.00193688 -0.00190652 -0.00192409 0.00191806 -0.000961358 -0.00177424 0.00153693 1.58565e-05 0.000139148 -0.00112142 0.000219548 -0.00152707
#      3     7     3     0     0     0     0     2     1  1980  2004     5 0.00104046 0.000577339 0.000989062 0.00156559 -0.000238446 -0.000608763 0.000150644 -0.000788755 0.00127865 -0.000634807 -0.00152794 -0.00122326 -0.00038408 0.000934731 0.000161483 -0.00150648 0.00177596 -0.000892146 0.000580528 0.00104247 0.000181536 -0.000538737 -0.000226183 0.000750919 -0.000152398
#      5     1     5     4     2     2     0     0     0     0     0     0
#      5     3     6     4     2     2     0     0     0     0     0     0
#      5     7     7     2     2     2     0     0     0     0     0     0
#      5     8     9     2     2     2     0     0     0     0     0     0
#      5     9    11     1     2     2     0     0     0     0     0     0
#      5    11    14     7     2     2     0     0     0     0     0     0
#      5    13    17     7     2     2     0     0     0     0     0     0
#      5    14    20     7     2     2     0     0     0     0     0     0
#      5    16    23     7     2     2     0     0     0     0     0     0
#      5    19    26     7     2     2     0     0     0     0     0     0
#      5    33    29     5     2     2     0     0     0     0     0     0
#      5    35    30     5     2     2     0     0     0     0     0     0
#      5    39    31     3     2     2     0     0     0     0     0     0
#      5    40    32     3     2     2     0     0     0     0     0     0
#      5    41    33     3     2     2     0     0     0     0     0     0
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 13 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 4 1 1 0.030 1 # BottomTrawl_Length_Comp
 4 2 1 0.095 1 # MidwaterTrawl_Length_Comp
 4 3 1 0.065 1 # Hake_Length_Comp
 4 4 1 0.237 1 # Net_Length_Comp
 4 5 1 0.138 1 # HnL_Length_Comp
 4 7 1 0.375 1 # Triennial_Length_Comp
 4 8 1 0.699 1 # NWFSC_Length_Comp
 5 1 1 0.081 1 # BottomTrawl_Marginal_Age_Comp
 5 2 1 0.130 1 # MidwaterTrawl_Marginal_Age_Comp
 5 3 1 0.110 1 # Hake_Marginal_Age_Comp
 5 4 1 0.240 1 # Net_Marginal_Age_Comp
 5 5 1 0.312 1 # HnL_Marginal_Age_Comp
 5 8 1 0.279 1 # NWFSC_CAAL
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  1 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  1 #_CPUE/survey:_3
#  0 #_CPUE/survey:_4
#  0 #_CPUE/survey:_5
#  1 #_CPUE/survey:_6
#  1 #_CPUE/survey:_7
#  1 #_CPUE/survey:_8
#  1 #_CPUE/survey:_9
#  1 #_discard:_1
#  1 #_discard:_2
#  0 #_discard:_3
#  0 #_discard:_4
#  1 #_discard:_5
#  0 #_discard:_6
#  0 #_discard:_7
#  0 #_discard:_8
#  0 #_discard:_9
#  0.035 #_lencomp:_1
#  0.13 #_lencomp:_2
#  0.06 #_lencomp:_3
#  0.23 #_lencomp:_4
#  0.2 #_lencomp:_5
#  0 #_lencomp:_6
#  0.38 #_lencomp:_7
#  0.73 #_lencomp:_8
#  0 #_lencomp:_9
#  0.08 #_agecomp:_1
#  0.16 #_agecomp:_2
#  0.11 #_agecomp:_3
#  0.23 #_agecomp:_4
#  0.31 #_agecomp:_5
#  0 #_agecomp:_6
#  0 #_agecomp:_7
#  0.33 #_agecomp:_8
#  0 #_agecomp:_9
#  1 #_init_equ_catch
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting 
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

