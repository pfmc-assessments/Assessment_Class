#C file created using an r4ss function
#C file write time: 2025-04-17  11:43:43
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
2 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
11 #_Nblock_Patterns
3 2 1 1 1 1 3 1 1 1 1 #_blocks_per_pattern
#_begin and end years of blocks
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
1916 2018
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
0 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_no additional input for selected M option; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
3 #_Age(post-settlement)_for_L1;linear growth below this
40 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
2 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
3 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
 0.01	     0.3	 0.144401	 -2.3	0.438	3	  5	0	0	0	0	0	0	0	#_NatM_p_1_Fem_GP_1  
   10	      40	  20.8325	   27	   99	0	  3	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1 
   35	      60	  50.3914	   50	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1 
 0.01	     0.4	 0.171903	 0.15	   99	0	  2	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	     0.4	  0.10617	 0.07	   99	0	  3	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1  
 0.01	     0.4	0.0440323	 0.04	   99	0	  3	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	1.736e-05	    0	   99	0	-99	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	      10	    2.962	2.962	   99	0	-99	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1   
   -3	      50	     5.47	    7	   99	0	-99	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1    
   -3	       3	  -0.7747	   -1	   99	0	-99	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1 
   -1	       1	        1	    1	   99	0	-99	0	0	0	0	0	0	0	#_Eggs_alpha_Fem_GP_1
    0	       1	        0	    0	   99	0	-99	0	0	0	0	0	0	0	#_Eggs_beta_Fem_GP_1 
 0.01	     0.3	 0.154867	 -2.3	0.438	3	  5	0	0	0	0	0	0	0	#_NatM_p_1_Mal_GP_1  
   10	      40	  21.1828	   27	   99	0	  3	0	0	0	0	0	0	0	#_L_at_Amin_Mal_GP_1 
   35	      60	  44.1793	   45	   99	0	  2	0	0	0	0	0	0	0	#_L_at_Amax_Mal_GP_1 
 0.01	     0.4	 0.236074	 0.19	   99	0	  2	0	0	0	0	0	0	0	#_VonBert_K_Mal_GP_1 
 0.01	     0.4	0.0861627	 0.07	   99	0	  3	0	0	0	0	0	0	0	#_CV_young_Mal_GP_1  
 0.01	     0.4	0.0542116	 0.04	   99	0	  3	0	0	0	0	0	0	0	#_CV_old_Mal_GP_1    
   -3	       3	1.484e-05	    0	   99	0	-99	0	0	0	0	0	0	0	#_Wtlen_1_Mal_GP_1   
   -3	      10	    3.005	3.005	   99	0	-99	0	0	0	0	0	0	0	#_Wtlen_2_Mal_GP_1   
    0	       2	        1	    1	   99	0	-99	0	0	0	0	0	0	0	#_RecrDist_GP_1      
    0	       2	        1	    1	   99	0	-99	0	0	0	0	0	0	0	#_RecrDist_Area_1    
    0	       2	        1	    1	   99	0	-99	0	0	0	0	0	0	0	#_RecrDist_month_1   
    0	       2	        1	    1	   99	0	-99	0	0	0	0	0	0	0	#_CohortGrowDev      
1e-06	0.999999	      0.5	  0.5	  0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
  1	 20	10.8133	  10	  99	0	  2	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	  1	   0.72	0.72	0.16	2	 -5	0	0	0	0	0	0	0	#_SR_BH_steep
  0	  2	    0.6	0.65	  99	0	-50	0	0	0	0	0	0	0	#_SR_sigmaR  
 -5	  5	      0	   0	   1	0	-99	0	0	0	0	0	0	0	#_SR_regime  
  0	0.5	      0	   0	  99	0	-99	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1970 # first year of main recr_devs; early devs can preceed this era
2020 # last year of main recr_devs; forecast devs start in following year
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
#Fishing Mortality info
0.05 # F ballpark
-1982 # F ballpark year (neg value to disable)
1 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.9 # max F or harvest rate, depends on F_Method
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    1	1	0	1	0	1	#_BottomTrawl 
    3	1	0	1	1	0	#_Hake        
    6	1	0	1	0	1	#_JuvSurvey   
    7	1	0	1	1	0	#_Triennial   
    8	1	0	1	0	1	#_NWFSC       
    9	1	0	1	0	1	#_ForeignAtSea
-9999	0	0	0	0	0	#_terminator  
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-6.05923	0	 1	0	-1	0	0	0	0	0	 0	0	#_LnQ_base_BottomTrawl(1)  
  0	 2	0.157603	0	99	0	 2	0	0	0	0	0	 0	0	#_Q_extraSD_BottomTrawl(1) 
-20	 2	-11.1569	0	99	0	 1	0	0	0	0	0	10	1	#_LnQ_base_Hake(3)         
  0	 2	0.366788	0	99	0	 2	0	0	0	0	0	 0	0	#_Q_extraSD_Hake(3)        
-25	25	-1.02648	0	 1	0	-1	0	0	0	0	0	 0	0	#_LnQ_base_JuvSurvey(6)    
  0	 2	 0.83169	0	99	0	 2	0	0	0	0	0	 0	0	#_Q_extraSD_JuvSurvey(6)   
 -4	 4	-2.17803	0	99	0	 2	0	0	0	0	0	 9	1	#_LnQ_base_Triennial(7)    
  0	 2	       0	0	99	0	-2	0	0	0	0	0	 0	0	#_Q_extraSD_Triennial(7)   
-25	25	-3.16502	0	 1	0	-1	0	0	0	0	0	 0	0	#_LnQ_base_NWFSC(8)        
  0	 2	       0	0	99	0	-2	0	0	0	0	0	 0	0	#_Q_extraSD_NWFSC(8)       
-25	25	-11.4795	0	 1	0	-1	0	0	0	0	0	 0	0	#_LnQ_base_ForeignAtSea(9) 
  0	 2	0.579699	0	99	0	 2	0	0	0	0	0	 0	0	#_Q_extraSD_ForeignAtSea(9)
# timevary Q parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
1e-04	2	0.442829	0.5	0.5	6	3	#_LnQ_base_Hake(3)_BLK10add_1991    
1e-04	2	0.132873	0.5	0.5	6	3	#_LnQ_base_Triennial(7)_BLK9add_1995
# info on dev vectors created for Q parms are reported with other devs after tag parameter section
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	1	0	0	#_1 BottomTrawl  
24	1	0	0	#_2 MidwaterTrawl
24	0	0	0	#_3 Hake         
24	0	0	0	#_4 Net          
24	0	0	0	#_5 HnL          
 0	0	0	0	#_6 JuvSurvey    
27	0	0	3	#_7 Triennial    
27	0	0	3	#_8 NWFSC        
 5	0	0	3	#_9 ForeignAtSea 
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
10	0	0	0	#_1 BottomTrawl  
10	0	0	0	#_2 MidwaterTrawl
10	0	0	0	#_3 Hake         
10	0	0	0	#_4 Net          
10	0	0	0	#_5 HnL          
11	0	0	0	#_6 JuvSurvey    
10	0	0	0	#_7 Triennial    
11	0	0	0	#_8 NWFSC        
10	0	0	0	#_9 ForeignAtSea 
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
    10	59	   45.8061	 45	0.05	0	  1	0	0	0	0	0.5	 4	2	#_SizeSel_P_1_BottomTrawl(1)        
    -5	10	   2.49993	  5	0.05	0	  3	0	0	0	0	0.5	 0	0	#_SizeSel_P_2_BottomTrawl(1)        
    -4	12	   4.63698	  3	0.05	0	  2	0	0	0	0	0.5	 4	2	#_SizeSel_P_3_BottomTrawl(1)        
    -2	10	         9	 10	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_4_BottomTrawl(1)        
    -9	10	        -9	0.5	0.05	0	 -3	0	0	0	0	0.5	 0	0	#_SizeSel_P_5_BottomTrawl(1)        
    -9	 9	         8	0.5	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_6_BottomTrawl(1)        
    -5	60	   3.55191	  0	  99	0	  4	0	0	0	0	  0	 2	2	#_SizeSel_PRet_1_BottomTrawl(1)     
  0.01	 8	  0.955332	  1	  99	0	  4	0	0	0	0	  0	 2	2	#_SizeSel_PRet_2_BottomTrawl(1)     
   -10	10	   4.59512	 10	  99	0	 -2	0	0	0	0	  0	 1	2	#_SizeSel_PRet_3_BottomTrawl(1)     
   -10	10	         0	  0	  99	0	-99	0	0	0	0	  0	 0	0	#_SizeSel_PRet_4_BottomTrawl(1)     
    10	59	   37.4551	 45	0.05	0	  1	0	0	0	0	0.5	 7	2	#_SizeSel_P_1_MidwaterTrawl(2)      
   -10	10	  -9.30925	  5	0.05	0	  3	0	0	0	0	0.5	 0	0	#_SizeSel_P_2_MidwaterTrawl(2)      
    -4	12	   2.90638	  3	0.05	0	  2	0	0	0	0	0.5	 7	2	#_SizeSel_P_3_MidwaterTrawl(2)      
    -2	10	   4.51658	 10	0.05	0	  4	0	0	0	0	0.5	 7	2	#_SizeSel_P_4_MidwaterTrawl(2)      
    -9	10	        -9	0.5	0.05	0	 -3	0	0	0	0	0.5	 0	0	#_SizeSel_P_5_MidwaterTrawl(2)      
    -9	 9	 -0.614281	0.5	0.05	0	  4	0	0	0	0	0.5	 7	2	#_SizeSel_P_6_MidwaterTrawl(2)      
    -5	60	        -5	  0	  99	0	 -9	0	0	0	0	  0	 0	0	#_SizeSel_PRet_1_MidwaterTrawl(2)   
  0.01	 8	       1.2	  1	  99	0	 -9	0	0	0	0	  0	 0	0	#_SizeSel_PRet_2_MidwaterTrawl(2)   
   -10	10	   4.59512	 10	  99	0	 -2	0	0	0	0	  0	 7	2	#_SizeSel_PRet_3_MidwaterTrawl(2)   
   -10	10	         0	  0	  99	0	-99	0	0	0	0	  0	 0	0	#_SizeSel_PRet_4_MidwaterTrawl(2)   
    10	59	   43.2065	 45	0.05	0	  1	0	0	0	0	0.5	11	2	#_SizeSel_P_1_Hake(3)               
    -5	10	   2.50295	  5	0.05	0	  3	0	0	0	0	0.5	11	2	#_SizeSel_P_2_Hake(3)               
    -4	12	   3.72893	  3	0.05	0	  2	0	0	0	0	0.5	11	2	#_SizeSel_P_3_Hake(3)               
    -2	10	         9	 10	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_4_Hake(3)               
    -9	10	        -9	0.5	0.05	0	 -3	0	0	0	0	0.5	 0	0	#_SizeSel_P_5_Hake(3)               
    -9	 9	         8	0.5	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_6_Hake(3)               
    10	59	   42.7767	 45	0.05	0	  1	0	0	0	0	0.5	 0	0	#_SizeSel_P_1_Net(4)                
    -5	10	   2.50691	  5	0.05	0	  3	0	0	0	0	0.5	 0	0	#_SizeSel_P_2_Net(4)                
    -4	12	   3.58755	  3	0.05	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_P_3_Net(4)                
    -2	10	         9	 10	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_4_Net(4)                
    -9	10	        -9	0.5	0.05	0	 -3	0	0	0	0	0.5	 0	0	#_SizeSel_P_5_Net(4)                
    -9	 9	         8	0.5	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_6_Net(4)                
    10	59	   17.5003	 45	0.05	0	  5	0	0	0	0	0.5	 5	2	#_SizeSel_P_1_HnL(5)                
    -5	10	   -1.7446	  5	0.05	0	  3	0	0	0	0	0.5	 0	0	#_SizeSel_P_2_HnL(5)                
    -5	12	  -2.48616	  3	0.05	0	  2	0	0	0	0	0.5	 5	2	#_SizeSel_P_3_HnL(5)                
    -2	10	         9	 10	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_4_HnL(5)                
    -9	10	        -9	0.5	0.05	0	 -3	0	0	0	0	0.5	 0	0	#_SizeSel_P_5_HnL(5)                
    -9	 9	         8	0.5	0.05	0	 -4	0	0	0	0	0.5	 0	0	#_SizeSel_P_6_HnL(5)                
#    -5	60	   23.9094	  0	  99	0	  2	0	0	0	0	  0	 3	2	#_SizeSel_PRet_1_HnL(5)             
#  0.01	 8	  0.563213	  1	  99	0	  3	0	0	0	0	  0	 3	2	#_SizeSel_PRet_2_HnL(5)             
#   -10	10	    7.2008	 10	  99	0	  1	0	0	0	0	  0	 3	2	#_SizeSel_PRet_3_HnL(5)             
#   -10	10	         0	  0	  99	0	-99	0	0	0	0	  0	 0	0	#_SizeSel_PRet_4_HnL(5)             
     0	 2	         0	  0	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Code_Triennial(7)  
-0.001	 1	  0.124088	  0	   0	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_GradLo_Triennial(7)
    -1	 1	 0.0223155	  0	   0	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_GradHi_Triennial(7)
     8	56	        24	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_1_Triennial(7)
     8	56	        34	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_2_Triennial(7)
     8	56	        48	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_3_Triennial(7)
   -10	10	  -1.87994	-10	  99	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_1_Triennial(7)  
   -10	10	        -1	-10	  99	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_2_Triennial(7)  
   -10	10	  0.449632	-10	  99	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_3_Triennial(7)  
     0	 2	         0	  0	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Code_NWFSC(8)      
-0.001	 1	  0.493265	  0	   0	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_GradLo_NWFSC(8)    
    -1	 1	-0.0163988	  0	   0	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_GradHi_NWFSC(8)    
     8	56	        24	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_1_NWFSC(8)    
     8	56	        34	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_2_NWFSC(8)    
     8	56	        48	-10	   0	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spline_Knot_3_NWFSC(8)    
   -10	10	  -2.75945	-10	  99	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_1_NWFSC(8)      
   -10	10	        -1	-10	  99	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_2_NWFSC(8)      
   -10	10	-0.0378029	-10	  99	0	  2	0	0	0	0	0.5	 0	0	#_SizeSel_Spine_Val_3_NWFSC(8)      
    -2	60	         0	  0	 0.2	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_P_1_ForeignAtSea(9)       
    -2	60	         0	  0	 0.2	0	-99	0	0	0	0	0.5	 0	0	#_SizeSel_P_2_ForeignAtSea(9)       
#_AgeSelex
0	 1	 0	0	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_1_JuvSurvey(6)
0	 1	 0	0	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_2_JuvSurvey(6)
0	 1	 0	0	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_1_NWFSC(8)    
0	50	40	0	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_2_NWFSC(8)    
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
  10	59	  39.3021	 45	0.05	0	 1	#_SizeSel_P_1_BottomTrawl(1)_BLK4repl_1916     
  -4	12	  3.46297	  3	0.05	0	 2	#_SizeSel_P_3_BottomTrawl(1)_BLK4repl_1916     
  -5	50	  26.7803	 34	  99	0	 3	#_SizeSel_PRet_1_BottomTrawl(1)_BLK2repl_1982  
  -5	50	  25.7623	 34	  99	0	 3	#_SizeSel_PRet_1_BottomTrawl(1)_BLK2repl_1990  
0.01	 5	 0.998554	  1	  99	0	 3	#_SizeSel_PRet_2_BottomTrawl(1)_BLK2repl_1982  
0.01	 5	  1.79122	  1	  99	0	 3	#_SizeSel_PRet_2_BottomTrawl(1)_BLK2repl_1990  
 -10	10	  1.71493	 10	  99	0	 2	#_SizeSel_PRet_3_BottomTrawl(1)_BLK1repl_1982  
 -10	10	 0.673116	 10	  99	0	 2	#_SizeSel_PRet_3_BottomTrawl(1)_BLK1repl_1990  
 -10	10	 0.085133	 10	  99	0	 2	#_SizeSel_PRet_3_BottomTrawl(1)_BLK1repl_1998  
  10	59	  38.7217	 45	0.05	0	 1	#_SizeSel_P_1_MidwaterTrawl(2)_BLK7repl_1916   
  10	59	  38.0921	 45	0.05	0	 1	#_SizeSel_P_1_MidwaterTrawl(2)_BLK7repl_1983   
  10	59	  37.5949	 45	0.05	0	 1	#_SizeSel_P_1_MidwaterTrawl(2)_BLK7repl_2002   
  -4	12	  3.35403	  3	0.05	0	 2	#_SizeSel_P_3_MidwaterTrawl(2)_BLK7repl_1916   
  -4	12	  3.07871	  3	0.05	0	 2	#_SizeSel_P_3_MidwaterTrawl(2)_BLK7repl_1983   
  -4	12	  2.83281	  3	0.05	0	 2	#_SizeSel_P_3_MidwaterTrawl(2)_BLK7repl_2002   
  -2	10	  4.30674	 10	0.05	0	 4	#_SizeSel_P_4_MidwaterTrawl(2)_BLK7repl_1916   
  -2	10	  3.15164	 10	0.05	0	 4	#_SizeSel_P_4_MidwaterTrawl(2)_BLK7repl_1983   
  -2	10	 -1.63262	 10	0.05	0	 4	#_SizeSel_P_4_MidwaterTrawl(2)_BLK7repl_2002   
  -9	 9	 -2.32627	0.5	0.05	0	 4	#_SizeSel_P_6_MidwaterTrawl(2)_BLK7repl_1916   
  -9	 9	-0.454742	0.5	0.05	0	 4	#_SizeSel_P_6_MidwaterTrawl(2)_BLK7repl_1983   
  -9	 9	  1.74863	0.5	0.05	0	 4	#_SizeSel_P_6_MidwaterTrawl(2)_BLK7repl_2002   
 -10	10	   4.5912	 10	  99	0	-2	#_SizeSel_PRet_3_MidwaterTrawl(2)_BLK7repl_1916
 -10	10	  1.65362	 10	  99	0	 2	#_SizeSel_PRet_3_MidwaterTrawl(2)_BLK7repl_1983
 -10	10	   1.8541	 10	  99	0	 2	#_SizeSel_PRet_3_MidwaterTrawl(2)_BLK7repl_2002
  10	59	  43.2065	 45	0.05	0	 1	#_SizeSel_P_1_Hake(3)_BLK11repl_1916           
  -5	10	  2.50295	  5	0.05	0	 3	#_SizeSel_P_2_Hake(3)_BLK11repl_1916           
  -4	12	  3.72893	  3	0.05	0	 2	#_SizeSel_P_3_Hake(3)_BLK11repl_1916           
  15	59	  37.9869	 45	0.05	0	 1	#_SizeSel_P_1_HnL(5)_BLK5repl_1916             
  -4	12	  3.86062	  3	0.05	0	 2	#_SizeSel_P_3_HnL(5)_BLK5repl_1916             
#  -5	50	       -5	 34	  99	0	-2	#_SizeSel_PRet_1_HnL(5)_BLK3repl_1916          
# 0.1	 8	      1.2	  1	  99	0	-3	#_SizeSel_PRet_2_HnL(5)_BLK3repl_1916          
# -10	10	   4.5912	 10	  99	0	-3	#_SizeSel_PRet_3_HnL(5)_BLK3repl_1916          
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_factor	fleet	value
    4	1	0.067399	#_Variance_adjustment_list1 
    4	2	0.213207	#_Variance_adjustment_list2 
    4	3	0.124653	#_Variance_adjustment_list3 
    4	4	0.481985	#_Variance_adjustment_list4 
    4	5	0.239654	#_Variance_adjustment_list5 
    4	7	0.373267	#_Variance_adjustment_list6 
    4	8	 0.70329	#_Variance_adjustment_list7 
    5	1	0.162182	#_Variance_adjustment_list8 
    5	2	0.266212	#_Variance_adjustment_list9 
    5	3	0.221241	#_Variance_adjustment_list10
    5	4	0.489099	#_Variance_adjustment_list11
    5	5	0.625465	#_Variance_adjustment_list12
    5	8	0.241542	#_Variance_adjustment_list13
-9999	0	       0	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 13 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    4	1	1	0.5	1	#_length_BottomTrawl_sizefreq_method_1_Phz1  
    4	2	1	0.5	1	#_length_MidwaterTrawl_sizefreq_method_1_Phz1
    4	3	1	0.5	1	#_length_Hake_sizefreq_method_1_Phz1         
    4	4	1	0.5	1	#_length_Net_sizefreq_method_1_Phz1          
    4	5	1	0.5	1	#_length_HnL_sizefreq_method_1_Phz1          
    4	7	1	  1	1	#_length_Triennial_sizefreq_method_1_Phz1    
    4	8	1	  1	1	#_length_NWFSC_sizefreq_method_1_Phz1        
    5	1	1	0.5	1	#_age_BottomTrawl_Phz1                       
    5	2	1	0.5	1	#_age_MidwaterTrawl_Phz1                     
    5	3	1	0.5	1	#_age_Hake_Phz1                              
    5	4	1	0.5	1	#_age_Net_Phz1                               
    5	5	1	0.5	1	#_age_HnL_Phz1                               
    5	8	1	  1	1	#_age_NWFSC_Phz1                             
-9999	0	0	  0	0	#_terminator                                 
#
0 # 0/1 read specs for more stddev reporting
#
999
