# 2015 Widow Rockfish Assessment
# Allan C. Hicks and Chantel R. Wetzel
# NWFSC, NOAA, Seattle, WA

1   #_N_Growth_Patterns
1   #_N_Morphs_Within_GrowthPattern

7         #_Nblock_Designs
3 2 1 1 1 1 3#_blocks_per_pattern
1982 1989 1990 1997 1998 2010 # Block Years for Bottom Trawl Retention    large trip limits 1982, smaller trip limits 1985, serious trip limits on bottom trawl, bottom trawl landings greatly declined (pre-1982 same as post-2010)
1982 1989 1990 2010           # Block Years for BT Retention    large trip limits 1982, smaller trip limits 1985, RCA in 2002 (does it matter for midwater?) & midwater trawl landings greatly declined (pre-1982 same as post-2010)
1916 1982                     # retention in HnL. Unknown pre 2004. 1983 is when trip limits went into affect for the entire year.
1916 2001                     # Block Years for before RCA's (Bottom Trawl)
1916 2002                     # HnL before RCA's
1995 2012 # Block Pattern for Triennial Selectivity when they started going deeper
1916 1982 1983 2001 2002 2010  #blocks for MWT selex and retention based on changes in catch

0.5 #_fracfemale
0   #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
#1  #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
#4  #_N_breakpoints
#3 4 25 26  # age(real) at M breakpoints
1   # GrowthModel: 1=vonBert with L1&L2; 2=vonBert with A0&Linf; 3=Richards; 4=readvector
3   #_Growth_Age_for_L1
40  #_Growth_Age_for_L2 (999 to use as Linf)
0   #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0   #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
2   #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity for each female GP; 4=read age-fecundity for each female GP
3   #_First_Mature_Age (from Barss & Echecverria 1987)
1   #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b;(4)eggs=a+b*L;(5)eggs=a+B*W
0   #_hermaphroditism   option: 0=none; 1=age-specific  fxn
1   #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
2   #_env/block/dev_adjust_method   (1=standard;    2=logistic  transform   keeps   in  base    parm    bounds; 3=standard  w/  no  bound   check)

#Biology parameters
#_LO    HI      INIT    PRIOR   PR_type SD      PHASE   env usdev   dminyr  dmaxyr  dev_std Block   Block_Fxn
0.01    0.3     0.081   -2.51    3      0.52    5       0   0       0       0       0       0       0   # NatM_p_1_Fem_GP_1
10      40      27.5    27      -1      99      3       0   0       0       0       0       0       0   # L_at_Amin_Fem_GP_1
35      60      50      50      -1      99      2       0   0       0       0       0       0       0   # L_at_Amax_Fem_GP_1
0.01    0.4     0.15    0.15    -1      99      2       0   0       0       0       0       0       0   # VonBert_K_Fem_GP_1
0.01    0.4     0.07    0.07    -1      99      3       0   0       0       0       0       0       0   # CV_young_Fem_GP_1
0.01    0.4     0.04    0.04    -1      99      3       0   0       0       0       0       0       0   # CV_old_Fem_GP_1

0.01    0.3     0.081   -2.51    3      0.52    5       0   0       0       0       0       0       0   # NatM_p_1_Male_GP_1
10      40      26      27      -1      99      3       0   0       0       0       0       0       0   # L_at_Amin_Male_GP_2
35      60      44      45      -1      99      2       0   0       0       0       0       0       0   # L_at_Amax_Male_GP_2
0.01    0.4     0.21    0.19    -1      99      2       0   0       0       0       0       0       0   # VonBert_K_Male_GP_2
0.01    0.4     0.07    0.07    -1      99      3       0   0       0       0       0       0       0   # CV_young_Male_GP_2
0.01    0.4     0.04    0.04    -1      99      3       0   0       0       0       0       0       0   # CV_old_Male_GP_2

-3      3   0.00001736  0.0     -1      99      -99     0   0       0       0       0       0       0   # Wtlen1_Fem
-3      10      2.962   2.962   -1      99      -99     0   0       0       0       0       0       0   # Wtlen2_Fem
#Maturity from CA and OR samples from Barss & Echeverria 1987
-3    50    5.47        7       -1      99      -99     0   0       0       0       0       0       0   # Mat50_Fem
-3    3   -0.7747      -1       -1      99      -99     0   0       0       0       0       0       0   # Mat_slope_Fem
#From Dick (2009) (Tables 10 & 11 median, converted to kg from g, intercept*1000, slope*1e6)
#-1   1000000  286200   0       -1      99      -99     0   0       0       0       0       0       0   # Eggs/kg_inter_Fem
#0      10000    2600   1       -1      99      -99     0   0       0       0       0       0       0   # Eggs/kg_slope_wt_Fem
#proportional using option 1 (Dick's results show non-significant slope)
-1       1       1       1      -1      99      -99     0   0       0       0       0       0       0   # Eggs/kg_inter_Fem
0        1       0       0      -1      99      -99     0   0       0       0       0       0       0   # Eggs/kg_slope_wt_Fem

-3      3   0.00001484  0.0     -1      99      -99     0   0       0       0       0       0       0   # Wtlen1_Mal
-3      10      3.005 3.005     -1      99      -99     0   0       0       0       0       0       0   # Wtlen2_Mal

# Unused recruitment interactions
0   2   1   1   -1  99  -99 0   0   0   0   0   0   0   #RecrDist_GP_1
0   2   1   1   -1  99  -99 0   0   0   0   0   0   0   #RecrDist_Area
0   2   1   1   -1  99  -99 0   0   0   0   0   0   0   #RecrDist_Seas
0   2   1   1   -1  99  -99 0   0   0   0   0   0   0   #CohortGrowDev

# 0  #custom_MG-env_setup (0/1)
# -2 2 0 0 -1 99 -2 #_placeholder for no MG-environ parameters

# 1  #custom_MG-block_setup (0/1)
# -2 2 0 0 -1 99 6 #_placeholder for no MG-block parameters
# -2 2 0 0 -1 99 6 #_placeholder for no MG-block parameters
# -2 2 0 0 -1 99 6 #_placeholder for no MG-block parameters

#_seasonal_effects_on_biology_parms
#_femwtlen1 femwtlen2   mat1    mat2    fec1    fec2    Malewtlen1  malewtlen2  L1  K
0           0           0       0       0       0       0           0           0   0

# -2 2 0 0 -1 99 -2 #_placeholder for no seasonal MG parameters
# -2 2 0 0 -1 99 -2 #_placeholder for no MG dev parameters
# if use Rick's recruit dist dev, active next line (phase for MGparm_dev)
#7  # placeholder for #_MGparm_Dev_Phase

#_Spawner-Recruitment
#_SR functions: 1=Beverton Holt with flat-top beyond Bzero; 2=Ricker; 3= Standard BH; 4=SCAA; 5=Hockey; 6=Shepard_3Parm
3   #_SR_function

#_LO    HI  INIT    PRIOR   PR_type     SD      PHASE
1       20  11.0     10      -1          99      1       # SR_R0
0.2     1   0.798    0.798   2       0.132      -5      # SR_steep: Thorson's new prior for 2015 assessments without widow
#0.2     1   0.773    0.773   2       0.147      -5      # SR_steep: Thorson's new prior for 2015 assessments
#0.2     1   0.76     0.76    2       0.17      -5       # SR_steep: Martin's prior (2010)
0       2   0.60    0.65    -1          99      -50     # SR_sigmaR
-5      5   0       0       -1          1       -99     # SR_envlink
-5      5   0       0       -1          1       -99     # SR_R1_offset
0       0.5 0       0       -1          99      -99     # SR_autocorr

0       #_SR_env_link
0       #_SR_env_target_0=none;1=devs;_2=R0;_3=steepness

1       # do_recdev:  0=none; 1=devvector; 2=simple deviations
1970    # first year of main recr_devs; early devs can preceed this era
2010   # last year of main recr_devs; forecast devs start in following year
2       #_recdev phase
1       # (0/1) to read 13 advanced options
    1900    #_recdev_early_start (0=none; neg value makes relative to recdev_start)
    4       #_recdev_early_phase
    0       #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
    1       #_lambda for forecast recr dev occurring before endyr+1
    1962    #_last_early_yr_nobias_adj_in_MPD
    1976    #_first_yr_fullbias_adj_in_MPD
    2010    #_last_yr_fullbias_adj_in_MPD
    2013    #_first_recent_yr_nobias_adj_in_MPD
    0.90    # Max bias adjustment
    0       #_period of cycles in recruitment (N parms read below)
    -5      #min rec_dev
    5       #max rec_dev
    0       #_read_recdevs
    #_end of advanced SR options


#Fishing Mortality info
0.05    # F ballpark for tuning early phases
-1982    # F ballpark year (neg value to disable)
1       # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.9     # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
# 4  # N iterations for tuning F in hybrid method (recommend 3 to 7)

#Fleet Year Seas F_value se phase (for detailed setup of F_Method=2)

#_initial_F_parms
#_LO    HI  INIT    PRIOR   PR_type SD  PHASE
0       0.5 0       0       -1      99  -99  # InitF_1 BottomTrawl
0       0.5 0       0       -1      99  -99  # InitF_2 MidwaterTrawl
0       0.5 0       0       -1      99  -99  # InitF_3 Hake
0       0.5 0       0       -1      99  -99  # InitF_4 Net
0       0.5 0       0       -1      99  -99  # InitF_5 HnL


#_Q_setup
# Q_type options:  <0=mirror, 0=median_float, 1=mean_float, 2=parameter, 3=parm_w_random_dev, 4=parm_w_randwalk, 5=mean_unbiased_float_assign_to_parm
#Den-dep  env-var  extra_se  Q_type
 0        0        1         0       # 1 BottomTrawl
 0        0        0         0       # 2 MidwaterTrawl
 0        0        1         4       # 3 Hake   #Q type 4 to break JV and Domestic
 0        0        0         0       # 4 Net
 0        0        0         0       # 5 HnL
 0        0        1         0       # 6 JuvSurvey
 0        0        1         4       # 7 Triennial
 0        0        1         0       # 8 NWFSCcombo
 0        0        1         0       # 9 ForeignAtSea

1 # Par setup: 0=read one parm for each fleet with random q; 1=read a parm for each year of index
#_Cond 0 #_If q has random component, then 0=read one parm for each fleet with random q; 1=read a parm for each year of index
#_Q_parms(if_any)

#Extra SD parameters for surveys
#Lo Hi  Init    Prior  Prior Prior Phase
0   2   0.01    0      -1    99    2   #BottomTrawl
0   2   0.01    0      -1    99    2   #Hake
0   2   0.01    0      -1    99    2   #JuvSurvey
0   2   0.00    0      -1    99    -2   #Triennial
0   2   0.00    0      -1    99    -2   #NWFSC_combo
0   2   0.01    0      -1    99    2   #ForeignAtSea

# Lo    Hi  Init Prior PrType PrSD Phase
# Early period
 -20    2  -10   0     -1      99   1    # Hake JV and Domestic (log) base parameter (1983)
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1985 deviation
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1986 deviation
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1987 deviation
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1988 deviation
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1989 deviation
 -4     4   0    0     -1      99   -50  # Hake JV and Domestic 1990 deviation
# Late period
  -4    4   0.4  0     -1      99   1    # Hake JV and Domestic 1991 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1992 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1993 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1994 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1995 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1996 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1997 deviation
  -4    4   0    0     -1      99   -50  # Hake JV and Domestic 1998 deviation

# Lo    Hi  Init Prior PrType PrSD Phase
# Early period
 -10    2  -2   0   -1  99  1   # Triennial (log) base parameter (1980)
 -4     4   0   0   -1  99  -50 # Triennial 1983 deviation
 -4     4   0   0   -1  99  -50 # Triennial 1986 deviation
 -4     4   0   0   -1  99  -50 # Triennial 1989 deviation
 -4     4   0   0   -1  99  -50 # Triennial 1992 deviation
# Late period
  -4    4   0   0   -1  99  -1   # Triennial 1995 deviation
  -4    4   0   0   -1  99  -50 # Triennial 1998 deviation
  -4    4   0   0   -1  99  -50 # Triennial 2001 deviation
  -4    4   0   0   -1  99  -50 # Triennial 2004 deviation


#_size_selex_Setup
#_SelPattern    Do_retain   Do_male Special
24 1 0 0    # 1 BottomTrawl
24 1 0 0    # 2 MidwaterTrawl
24 0 0 0    # 3 Hake
24 0 0 0    # 4 Net
24 1 0 0    # 5 HnL
 0 0 0 0    # 6 JuvSurvey
 27 0 0 3    # 7 Triennial 8-56 with
 27 0 0 3    # 8 NWFSCcombo
 5 0 0 3    # 9 ForeignAtSea mirrors Hake

#_age_selex_Setup
#_SelPattern    Do_retain   Do_male Special
 10 0 0 0   # 1 BottomTrawl
 10 0 0 0   # 2 MidwaterTrawl
 10 0 0 0   # 3 Hake
 10 0 0 0   # 4 Net
 10 0 0 0   # 5 HnL
 11 0 0 0   # 6 JuvSurvey (selects age 0)
 10 0 0 0   # 7 Triennial
 11 0 0 0   # 8 NWFSCcombo (type 11 to select age 0)
 10 0 0 0   # 9 ForeignAtSea

# double normal parameter comments
# P1=PEAL: beginning size for the plateau; P2=TOP: width of plateau, as logistic between PEAK and MAXLENG
# P3=ASC-WIDTH: ln(width); P4=DESC-WIDTH: ln(width); P5=INIT: logistive between 0 and 1; P6=FINAL: logistic between 0 and 1
# for initial P5 parameter: -999 or -1000: ignore the inital selectivity algorithm and simply decay small fish according to P3

#LO    HI       INIT    PRIOR  PRtype SD    PHASE    env-var use_dev dev_minyr dev_maxyr dev_stddev Block Block_Fxn
########### Bottom Trawl
10      59      38.0    45.0   -1      0.05    1       0   0   0   0   0.5   4   2  #PEAK
-5.0    10.0    2.5     5.0    -1      0.05    3       0   0   0   0   0.5   0   0  #TOP_WIDTH
-4.0    12.0    4.0     3.0    -1      0.05    2       0   0   0   0   0.5   4   2  #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   -4       0   0   0   0   0.5   0   0  #DESC_WIDTH
-9      10.0    -9      0.5    -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
-9      9.0      8      0.5    -1      0.05   -4       0   0   0   0   0.5   0   0  #FINAL
#-999    5.0    -999    -999   -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
#-999    10.0   -999     5.0   -1      0.05   -4       0   0   0   0   0.5   0   0  #FINAL
#RETENTION Bottom Trawl (discard comps may suggest Horizontal line/constant)
 -5      60      -3      0     -1      99      4       0   0   0   0   0     2   2    #inflection
0.01     8      1.2     1.0    -1      99      4       0   0   0   0   0     2   2    #slope
0.2      1      0.99    1      -1      99     -2       0   0   0   0   0     1   2    #asymptote
-10     10      0.0     0.0    -1      99    -99       0   0   0   0   0     0   0    #male offset to inflection (arithmetic)

############ Midwater Trawl
10      59      38.0    45.0   -1      0.05    1       0   0   0   0   0.5   7   2  #PEAK
-10.0    10.0    0.0     5.0    -1      0.05    3       0   0   0   0   0.5   0   0  #TOP_WIDTH
-4.0    12.0    3.0     3.0    -1      0.05    2       0   0   0   0   0.5   7   2  #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05    4       0   0   0   0   0.5   7   2  #DESC_WIDTH
-9      10.0    -9      0.5    -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
-9      9.0      8      0.5    -1      0.05    4       0   0   0   0   0.5   7   2  #FINAL
#RETENTION MIDWATER TRAWL (Horizontal line/constant because no comps)
 -5      60      -5      0     -1      99     -9       0   0   0   0   0     0   0    #inflection
0.01     8      1.2     1.0    -1      99     -9       0   0   0   0   0     0   0    #slope
0.2      1      0.99    1      -1      99     -2       0   0   0   0   0     7   2    #asymptote
-10     10      0.0     0.0    -1      99    -99       0   0   0   0   0     0   0    #male offset to inflection (arithmetic)

######### Hake
10      59      40.0    45.0   -1      0.05    1       0   0   0   0   0.5   0   0  #PEAK
-5.0    10.0    2.5     5.0    -1      0.05    3       0   0   0   0   0.5   0   0  #TOP_WIDTH
-4.0    12.0    4.0     3.0    -1      0.05    2       0   0   0   0   0.5   0   0  #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   -4       0   0   0   0   0.5   0   0  #DESC_WIDTH
-9      10.0    -9      0.5    -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
-9      9.0      8      0.5    -1      0.05   -4       0   0   0   0   0.5   0   0  #FINAL

########## Net
10      59      40.0    45.0   -1      0.05    1       0   0   0   0   0.5   0   0  #PEAK
-5.0    10.0    2.5     5.0    -1      0.05    3       0   0   0   0   0.5   0   0  #TOP_WIDTH
-4.0    12.0    4.0     3.0    -1      0.05    2       0   0   0   0   0.5   0   0  #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   -4       0   0   0   0   0.5   0   0  #DESC_WIDTH
-9      10.0    -9      0.5    -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
-9      9.0      8      0.5    -1      0.05   -4       0   0   0   0   0.5   0   0  #FINAL

#############  HnL
10      59      37.0    45.0   -1      0.05    1       0   0   0   0   0.5   5   2  #PEAK
-5.0    10.0    2.5     5.0    -1      0.05    3       0   0   0   0   0.5   0   0  #TOP_WIDTH
-5.0    12.0    4.0     3.0    -1      0.05    2       0   0   0   0   0.5   5   2  #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   -4       0   0   0   0   0.5   0   0  #DESC_WIDTH
-9      10.0    -9      0.5    -1      0.05   -3       0   0   0   0   0.5   0   0  #INIT
-9      9.0      8      0.5    -1      0.05   -4       0   0   0   0   0.5   0   0  #FINAL
#RETENTION HnL (discard comps suggest size based sorting)
 -5      60      25       0    -1      99      2       0   0   0   0   0     3   2    #inflection
0.01     8       1      1.0    -1      99      3       0   0   0   0   0     3   2    #slope
0.2      1     0.9      1      -1      99      1       0   0   0   0   0     3   2    #asymptote
-10     10      0.0     0.0    -1      99    -99       0   0   0   0   0     0   0    #male offset to inflection (arithmetic)

############### Juvenile Survey
#No parameters for length

############## Triennial Survey
0       2        0      0      -1      0     -99      0   0   0   0   0.5   0   0  #Spline code
-0.001  1        0.15   0      -1      0      2       0   0   0   0   0.5   0   0  #Spline GradLo
-1      1       -0.03   0      -1      0      2       0   0   0   0   0.5   0   0  #Spline GradHi
 8      56       24    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot1
 8      56       34    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot2
 8      56       48    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot3
-10     10       -3    -10     -1     99      2       0   0   0   0   0.5   0   0  #Spline Val1
-10     10       -1    -10     -1     99      -99     0   0   0   0   0.5   0   0  #Spline Val2
-10     10        0    -10     -1     99      2       0   0   0   0   0.5   0   0  #Spline Val3

############## NWFSCcombo Survey
0       2        0      0      -1      0     -99      0   0   0   0   0.5   0   0  #Spline code
-0.001  1        0.15   0      -1      0      2       0   0   0   0   0.5   0   0  #Spline GradLo
-1      1       -0.03   0      -1      0      2       0   0   0   0   0.5   0   0  #Spline GradHi
 8      56       24    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot1
 8      56       34    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot2
 8      56       48    -10     -1      0     -99      0   0   0   0   0.5   0   0  #Spline Knot3
-10     10       -3    -10     -1     99      2       0   0   0   0   0.5   0   0  #Spline Val1
-10     10       -1    -10     -1     99      -99     0   0   0   0   0.5   0   0  #Spline Val2
-10     10        0    -10     -1     99      2       0   0   0   0   0.5   0   0  #Spline Val3

############# Foreign At-sea fleet
-2      60       0      0      -1      0.2    -99      0   0   0   0   0.5   0   0 #MinBin (<=0 means first bin)
-2      60       0      0      -1      0.2    -99      0   0   0   0   0.5   0   0 #MaxBin (<=0 means last bin)

############### Juvenile Survey (Age paramters, select only age 0)
0      1        0       0      -1      99    -99        0   0   0   0   0.5  0   0    # Min
0      1        0       0      -1      99    -99        0   0   0   0   0.5  0   0    # Max

############### NWFSCcombo Survey (Age paramters, select age 0+)
0      1        0       0      -1      99    -99        0   0   0   0   0.5  0   0    # Min
0      50       40      0      -1      99    -99        0   0   0   0   0.5  0   0    # Max



1 #Custom Block Setup (0/1)
#LO     HI      INIT    PRIOR  PR_TYPE  SD    PHASE
#SELECTIVITY Bottom Trawl
10      59    34.0    45.0    -1      0.05   1  #PEAK (1930-2001)
-4.0    12.0   5.7     3.0    -1      0.05   2  #ASC_WIDTH (1930-2002)
#RETENTION Bottom trawl
 -5     50     35       34    -1      99     3  #inflection (1982-1989)
 -5     50     35       34    -1      99     3  #inflection (1990-2010)
0.01    5      2.5      1.0   -1      99     3  #slope (1982-1989)
0.01    5      2.5      1.0   -1      99     3  #slope (1990-2010)
0.2     1      0.8      1     -1      99     2  #asymptote (1982-1989)
0.2     1      0.5      1     -1      99     2  #asymptote (1990-1999)
0.2     1      0.5      1     -1      99     2  #asymptote (2000-2010)
#SELECTIVITY Midwtaer Trawl
10      59      38.0    45.0   -1      0.05    1   #PEAK
10      59      38.0    45.0   -1      0.05    1   #PEAK
10      59      38.0    45.0   -1      0.05    1   #PEAK
-4.0    12.0    3.0     3.0    -1      0.05    2   #ASC_WIDTH
-4.0    12.0    3.0     3.0    -1      0.05    2   #ASC_WIDTH
-4.0    12.0    3.0     3.0    -1      0.05    2   #ASC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   4   #DESC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   4   #DESC_WIDTH
-2.0    10.0    9.0    10.0    -1      0.05   4   #DESC_WIDTH
-9      9.0      8      0.5    -1      0.05   4       #FINAL
-9      9.0      8      0.5    -1      0.05   4       #FINAL
-9      9.0      8      0.5    -1      0.05   4       #FINAL
#RETENTION Midwater trawl
0.2     1     0.99     1      -1      99     -2  #asymptote (1916-1982)
0.2     1     0.8      1      -1      99     2  #asymptote (1983-2001)
0.2     1     0.8      1      -1      99     2  #asymptote (2002-2010)
#SELECTIVITY HnL
15      59    48.0    45.0    -1      0.05   1   #PEAK  (1916-2002)
-4.0    12.0  2.8     3.0     -1      0.05   2   #ASC_WIDTH (1916-2002)
#RETENTION HnL
 -5      50     -5     34    -1       99    -2   #inflection (1916-1983)
0.1      8      1.2    1.0   -1       99    -3   #slope  (1916-1983)
0.2      1      0.99   1     -1       99    -3   #asymptote  (1916-1983)
# #Selectivity Triennial
# -9      9.0      5     0.5   -1      0.05    3    #FINAL


#3 #selparm_dev_PH

1 #selparm_adjust_method: 1=standard; 2=logistic trans to keep in base parm bounds

# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
1 #_Variance_adjustments_to_input_values
#_fleet: 1    2    3     4     5     6    7     8     9
         0    0    0     0     0     0    0     0     0     #_add_to_survey_CV
         0    0    0     0     0     0    0     0     0     #_add_to_discard_stddev
         0    0    0     0     0     0    0     0     0     #_add_to_bodywt_CV
         1    1    1     1     1     1    1     1     1     #_mult_by_lencomp_N
         1    1    1     1     1     1    1     1     1    #_mult_by_agecomp_N
         1    1    1     1     1     1    1     1     1     #_mult_by_size-at-age_N

#
1 #_maxlambdaphase
1 #_sd_offset
#
13 # number of changes to make to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch;
# 9=init_equ_catch; 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin
#like_comp fleet/survey  phase  value  sizefreq_method
4           1            1      0.035    1    # lgth comps for bottom trawl
4           2            1      0.13    1    # lgth comps for midwater trawl
4           3            1      0.06    1    # lgth comps for hake
4           4            1      0.23    1    # lgth comps for net
4           5            1      0.20    1    # lgth comps for HnL
4           7            1      0.38    1    # lgth comps for Tri
4           8            1      0.73    1    # lgth comps for NWFSC
5           1            1      0.08    1    # age comps for bottom trawl
5           2            1      0.16    1    # age comps for midwater trawl
5           3            1      0.11    1    # age comps for hake
5           4            1      0.23    1    # age comps for net
5           5            1      0.31    1    # age comps for HnL
5           8            1      0.33    1    # age comps for NWFSC
# 4           1            1      0.07    1    # lgth comps for bottom trawl
# 4           2            1      0.26    1    # lgth comps for midwater trawl
# 4           3            1      0.12    1    # lgth comps for hake
# 4           4            1      0.47    1    # lgth comps for net
# 4           5            1      0.40    1    # lgth comps for HnL
# 4           7            1      0.38    1    # lgth comps for Tri
# 4           8            1      0.73    1    # lgth comps for NWFSC
# 5           1            1      0.16    1    # age comps for bottom trawl
# 5           2            1      0.32    1    # age comps for midwater trawl
# 5           3            1      0.21    1    # age comps for hake
# 5           4            1      0.47    1    # age comps for net
# 5           5            1      0.62    1    # age comps for HnL
# 5           8            1      0.33    1    # age comps for NWFSC

0 # (0/1) read specs for more stddev reporting

999
