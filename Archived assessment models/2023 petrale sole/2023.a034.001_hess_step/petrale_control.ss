#V3.30
#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2023-06-21 16:38:18
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
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
4 #_Nblock_Patterns
6 4 2 1 #_blocks_per_pattern
#_begin and end years of blocks
1973 1982 1983 1992 1993 2002 2003 2010 2011 2017 2018 2022
2002 2002 2003 2008 2009 2010 2011 2022
2010 2010 2011 2022
1995 2004
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
1 #_Age(post-settlement)_for_L1;linear growth below this
17 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
3 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
0.005	0.5	 0.142675	  -1.7793	0.31	3	 2	0	0	0	0	0	0	0	#_NatM_p_1_Fem_GP_1  
    5	 45	  8.87713	    17.18	  10	0	 3	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1 
   35	 80	  47.4413	     54.2	  10	0	 3	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1 
 0.04	0.5	 0.197191	    0.157	  99	0	 3	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1 
  0.5	 15	   1.2682	        3	  99	0	 3	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1  
  0.5	 15	   4.8766	        3	  99	0	 4	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1    
   -3	  3	2.035e-06	2.035e-06	  99	0	-3	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1   
    1	  5	    3.478	    3.478	  99	0	-3	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1   
   10	 50	    35.45	    35.45	  99	0	-3	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1    
   -3	  3	 -0.48921	 -0.48921	  99	0	-3	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1 
   -3	  1	  3.2e-11	        1	   1	0	-3	0	0	0	0	0	0	0	#_Eggs_alpha_Fem_GP_1
   -3	  5	     4.55	        0	   1	0	-3	0	0	0	0	0	0	0	#_Eggs_beta_Fem_GP_1 
0.005	0.6	 0.156902	  -1.6809	0.31	3	 2	0	0	0	0	0	0	0	#_NatM_p_1_Mal_GP_1  
    0	 45	        0	    17.18	  10	0	-3	0	0	0	0	0	0	0	#_L_at_Amin_Mal_GP_1 
   35	 80	  39.9318	     41.1	  10	0	 3	0	0	0	0	0	0	0	#_L_at_Amax_Mal_GP_1 
 0.04	0.5	 0.245868	    0.247	  99	0	 3	0	0	0	0	0	0	0	#_VonBert_K_Mal_GP_1 
  0.5	 15	  1.34386	        3	  99	0	 3	0	0	0	0	0	0	0	#_CV_young_Mal_GP_1  
  0.5	 15	  3.41291	        3	  99	0	 4	0	0	0	0	0	0	0	#_CV_old_Mal_GP_1    
   -3	  3	3.043e-06	3.043e-06	  99	0	-3	0	0	0	0	0	0	0	#_Wtlen_1_Mal_GP_1   
   -3	  5	    3.359	    3.359	  99	0	-3	0	0	0	0	0	0	0	#_Wtlen_2_Mal_GP_1   
    0	  1	        1	        1	   0	0	-4	0	0	0	0	0	0	0	#_CohortGrowDev      
  0.3	0.7	      0.5	      0.4	  99	0	-5	0	0	0	0	0	0	0	#_FracFemale_GP_1    
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
  5	20	9.64716	  9	  10	0	  1	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	 1	    0.8	0.8	0.09	0	 -5	0	0	0	0	0	0	0	#_SR_BH_steep
  0	 2	    0.5	0.9	   5	0	-99	0	0	0	0	0	0	0	#_SR_sigmaR  
 -5	 5	      0	  0	 0.2	0	 -2	0	0	0	0	0	0	0	#_SR_regime  
  0	 0	      0	  0	   0	0	-99	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1959 # first year of main recr_devs; early devs can preceed this era
2020 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
1845 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
3 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1935.3   #_last_early_yr_nobias_adj_in_MPD 
2002.0   #_first_yr_fullbias_adj_in_MPD 
2015.7   #_last_yr_fullbias_adj_in_MPD 
2021.8   #_first_recent_yr_nobias_adj_in_MPD 
0.8405  #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models)
0 #_period of cycles in recruitment (N parms read below)
-4 #min rec_dev
4 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.3 # F ballpark
-2001 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 0
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    3	1	0	1	0	1	#_Triennial 
    4	1	0	0	0	1	#_WCGBTS    
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  -15	15	-0.687865	   0	 1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_Triennial(3) 
0.001	 2	 0.259618	0.22	-1	0	 5	0	0	0	0	0	0	0	#_Q_extraSD_Triennial(3)
  -15	15	  1.40805	   0	 1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_WCGBTS(4)    
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	1	3	0	#_1 North    
24	1	3	0	#_2 South    
24	0	3	0	#_3 Triennial
24	0	3	0	#_4 WCGBTS   
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
10	0	0	0	#_1 North    
10	0	0	0	#_2 South    
10	0	0	0	#_3 Triennial
10	0	0	0	#_4 WCGBTS   
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  15	75	   61.478	43.1	99	0	 2	0	0	0	0	0	1	2	#_SizeSel_P_1_North(1)          
 -15	 4	      -15	 -15	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_North(1)          
  -4	12	  5.29517	3.42	99	0	 3	0	0	0	0	0	1	2	#_SizeSel_P_3_North(1)          
  -2	20	       20	0.21	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_4_North(1)          
-999	 9	     -999	-8.9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_North(1)          
-999	 9	     -999	0.15	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_6_North(1)          
  10	40	  28.5514	  15	99	0	 2	0	0	0	0	0	2	2	#_SizeSel_PRet_1_North(1)       
 0.1	10	    1.406	   3	99	0	 4	0	0	0	0	0	2	2	#_SizeSel_PRet_2_North(1)       
 -10	10	  9.71481	  10	99	0	 4	0	0	0	0	0	2	2	#_SizeSel_PRet_3_North(1)       
 -10	10	        0	   0	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_PRet_4_North(1)       
 -25	15	 -17.9766	   0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_1_North(1)    
 -15	15	 -1.72602	   0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_2_North(1)    
 -15	15	        0	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_3_North(1)    
 -15	15	        0	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_4_North(1)    
 -15	15	        1	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_5_North(1)    
  15	75	  54.5977	43.1	99	0	 2	0	0	0	0	0	1	2	#_SizeSel_P_1_South(2)          
 -15	 4	      -15	 -15	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_2_South(2)          
  -4	12	  6.02843	3.42	99	0	 3	0	0	0	0	0	1	2	#_SizeSel_P_3_South(2)          
  -2	20	       20	0.21	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_P_4_South(2)          
-999	 9	     -999	-8.9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_South(2)          
-999	 9	     -999	0.15	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_6_South(2)          
  10	40	  28.2307	  15	99	0	 2	0	0	0	0	0	3	2	#_SizeSel_PRet_1_South(2)       
 0.1	10	  1.16452	   3	99	0	 3	0	0	0	0	0	3	2	#_SizeSel_PRet_2_South(2)       
 -10	10	  6.62751	  10	99	0	 4	0	0	0	0	0	3	2	#_SizeSel_PRet_3_South(2)       
 -10	10	        0	   0	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_PRet_4_South(2)       
 -25	15	 -16.3837	   0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_1_South(2)    
 -15	15	 -2.07396	   0	99	0	 4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_2_South(2)    
 -15	15	        0	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_3_South(2)    
 -15	15	        0	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_4_South(2)    
 -15	15	        1	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_5_South(2)    
  15	61	  37.7261	43.1	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_Triennial(3)      
 -15	 4	      -15	  -1	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_2_Triennial(3)      
  -4	12	   4.7366	3.42	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_3_Triennial(3)      
  -2	20	       20	0.21	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_4_Triennial(3)      
-999	 9	     -999	-8.9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_Triennial(3)      
-999	 9	     -999	0.15	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_6_Triennial(3)      
 -15	15	 -4.46666	   0	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_1_Triennial(3)
 -15	15	-0.368962	   0	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_2_Triennial(3)
 -15	15	        0	   0	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_3_Triennial(3)
 -15	15	        0	   0	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_4_Triennial(3)
 -15	15	        1	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_5_Triennial(3)
  15	61	  51.3245	43.1	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_WCGBTS(4)         
 -15	 4	      -15	  -1	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_2_WCGBTS(4)         
  -4	12	  5.69654	3.42	99	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_3_WCGBTS(4)         
  -2	20	       20	0.21	99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_4_WCGBTS(4)         
-999	 9	     -999	-8.9	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_5_WCGBTS(4)         
-999	 9	     -999	0.15	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_6_WCGBTS(4)         
 -15	15	 -11.5553	   0	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_1_WCGBTS(4)   
 -15	15	-0.901813	   0	99	0	 3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_2_WCGBTS(4)   
 -15	15	        0	   0	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_3_WCGBTS(4)   
 -15	15	        0	   0	99	0	-3	0	0	0	0	0	0	0	#_SizeSel_PMalOff_4_WCGBTS(4)   
 -15	15	        1	   0	99	0	-4	0	0	0	0	0	0	0	#_SizeSel_PMalOff_5_WCGBTS(4)   
#_AgeSelex
#_No age_selex_parm
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
 15	75	 61.8126	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_1973   
 15	75	  58.425	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_1983   
 15	75	 57.2809	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_1993   
 15	75	 58.3231	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_2003   
 15	75	  58.573	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_2011   
 15	75	 58.7676	43.1	99	0	5	#_SizeSel_P_1_North(1)_BLK1repl_2018   
 -4	12	 5.52079	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_1973   
 -4	12	 5.42032	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_1983   
 -4	12	 5.60774	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_1993   
 -4	12	 5.49584	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_2003   
 -4	12	 5.42821	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_2011   
 -4	12	 5.30087	3.42	99	0	6	#_SizeSel_P_3_North(1)_BLK1repl_2018   
 10	40	 31.2756	  15	99	0	5	#_SizeSel_PRet_1_North(1)_BLK2repl_2002
 10	40	 30.0364	  15	99	0	5	#_SizeSel_PRet_1_North(1)_BLK2repl_2003
 10	40	 31.6662	  15	99	0	5	#_SizeSel_PRet_1_North(1)_BLK2repl_2009
 10	40	 27.1721	  15	99	0	5	#_SizeSel_PRet_1_North(1)_BLK2repl_2011
0.1	10	0.887097	   3	99	0	5	#_SizeSel_PRet_2_North(1)_BLK2repl_2002
0.1	10	 1.33968	   3	99	0	5	#_SizeSel_PRet_2_North(1)_BLK2repl_2003
0.1	10	 1.72013	   3	99	0	5	#_SizeSel_PRet_2_North(1)_BLK2repl_2009
0.1	10	 1.66993	   3	99	0	5	#_SizeSel_PRet_2_North(1)_BLK2repl_2011
-10	10	 9.76406	   0	99	0	5	#_SizeSel_PRet_3_North(1)_BLK2repl_2002
-10	10	 6.18735	   0	99	0	5	#_SizeSel_PRet_3_North(1)_BLK2repl_2003
-10	10	 3.64149	   0	99	0	5	#_SizeSel_PRet_3_North(1)_BLK2repl_2009
-10	10	 7.44097	   0	99	0	5	#_SizeSel_PRet_3_North(1)_BLK2repl_2011
 15	75	 51.1682	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_1973   
 15	75	 49.2949	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_1983   
 15	75	 50.3237	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_1993   
 15	75	 51.2959	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_2003   
 15	75	 50.9801	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_2011   
 15	75	 51.7918	43.1	99	0	5	#_SizeSel_P_1_South(2)_BLK1repl_2018   
 -4	12	 6.27617	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_1973   
 -4	12	 5.29108	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_1983   
 -4	12	 4.98497	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_1993   
 -4	12	 5.01588	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_2003   
 -4	12	 5.00364	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_2011   
 -4	12	 4.95519	3.42	99	0	6	#_SizeSel_P_3_South(2)_BLK1repl_2018   
 10	40	 30.8509	  15	99	0	5	#_SizeSel_PRet_1_South(2)_BLK3repl_2010
 10	40	 25.1312	  15	99	0	5	#_SizeSel_PRet_1_South(2)_BLK3repl_2011
0.1	10	 1.84447	   3	99	0	5	#_SizeSel_PRet_2_South(2)_BLK3repl_2010
0.1	10	 1.60882	   3	99	0	5	#_SizeSel_PRet_2_South(2)_BLK3repl_2011
-10	10	  9.3848	   0	99	0	5	#_SizeSel_PRet_3_South(2)_BLK3repl_2010
-10	10	 8.14009	   0	99	0	5	#_SizeSel_PRet_3_South(2)_BLK3repl_2011
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Data_type	Fleet	Value
    4	1	0.277766	#_Variance_adjustment_list1
    4	2	0.135516	#_Variance_adjustment_list2
    4	3	0.289957	#_Variance_adjustment_list3
    4	4	0.101749	#_Variance_adjustment_list4
    5	1	0.279775	#_Variance_adjustment_list5
    5	2	0.081769	#_Variance_adjustment_list6
    5	4	0.040054	#_Variance_adjustment_list7
-9999	0	       0	#_terminator               
#
15 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
-9999 0 0 0 0 # terminator
#
0 # 0/1 read specs for more stddev reporting
#
999
