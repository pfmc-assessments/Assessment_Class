#C Yelloweye 2017 control file
#_data_and_control_files: yelloweye_data.SS // yelloweye_control.SS
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
2 # recr_dist_method for parameters:  1=like 3.24; 2=main effects for GP, Settle timing, Area; 3=each Settle entity; 4=none when N_GP*Nsettle*pop==1
1 # Recruitment: 1=global; 2=by area (future option)
2 #  number of recruitment settlement assignments 
0 # year_x_area_x_settlement_event interaction requested (only for recr_dist_method=1)
#GPat month  area age (for each settlement assignment)
 1 1 1 0
 1 1 2 0
# 1 1 3 0
#
0 #_N_movement_definitions
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) if do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, GP=1, source=1 dest=2, age1=4, age2=10
#
4 #_Nblock_Patterns
1	1	1	1	 #_blocks_per_pattern 
# begin and end years of blocks
1992 	2004 	# allow change in Triennial catchability
2005	2016	# to allow chance in OR_Rec index catchability
2002 	2016 	#	to allow change in selex with IFQ
2002	2016	#	No yelloweye retention allowed in rec fisheries
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
1 1 1 1 1 # autogen
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if min=-12345
# 1st element for biology, 2nd for SR, 3rd for Q, 5th for selex, 4th reserved
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
1 		#_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
1 		#_N_breakpoints
4 		# age(real) at M breakpoints
1 		# GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K; 4=not implemented
1 		#_Growth_Age_for_L1
70 		#_Growth_Age_for_L2 (999 to use as Linf)
-999 	#_exponential decay for growth above maxage (fixed at 0.2 in 3.24; value should approx initial Z; -999 replicates 3.24)
0 		#_placeholder for future growth feature
0 		#_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 		#_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 		#_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2 		#_First_Mature_Age
2 		#_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 		#_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 		#_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO 			HI 				INIT 					PRIOR 					PR_SD 			PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
 0.01 			0.15 			0.04390336 	 			-3.12576436 		0.4384383 				3 			-1 		0 	0 	0 	0 	0 	0 	0 	#	Hamel prior, maxage = 123
 #0.01 			0.15 			0.06 	 			    -3.12576436 		0.4384383 				3 			-2 		0 	0 	0 	0 	0 	0 	0 	#	Hamel prior, maxage = 123
 1  			35 				23 						30 						99 					0 			2 		0 	0 	0 	0 	0 	0 	0 	# L_at_Amin_Fem_GP_1
 40 			120 			61 						66 						99 					0 			2 		0 	0 	0 	0 	0 	0 	0 	# L_at_Amax_Fem_GP_1
 0.01 			0.2 			0.05 					0.05 					99 					0 			1 		0 	0 	0 	0 	0 	0 	0 	# VonBert_K_Fem_GP_1
 0.01 			0.5				0.1 					0.1 					99 					0 			3 		0 	0 	0 	0 	0 	0 	0 	# CV_young_Fem_GP_1
 0.01			0.5				0.1 					0.1 					99 					0 		    7 		0 	0 	0 	0 	0 	0 	0 	# CV_old_Fem_GP_1
 -3 			3 				7.312807e-06			7.312807e-06 			99 					0 			-50  	0 	0 	0 	0 	0 	0 	0 	# Wtlen_1_Fem
 -3 			4 				3.242482 				3.242482 				99 					0 		 	-50  	0 	0 	0 	0 	0 	0 	0 	# Wtlen_2_Fem

#Maturity updated with Melissa Head data 
 38 			45		 		41.765		 			41.765 			   		99 					0 			-50  	0 	0 	0 	0 	0 	0 	0 	# Mat50%_Fem
 -3 			3 		        -0.36886 	 	    	-0.36886 	 			99 					0  	 	 	-50  	0 	0 	0 	0 	0 	0 	0 	# Mat_slope_Fem
# Fecundity updated with E.J.Dick 2017 results 
 -3 				300000 		7.21847E-08		7.21847E-08 					1 					6 		 	-6 		0 	0 	0 	0 	0 	0 	0 	# Million_Eggs/cm_a_Fem
 -3 				39000 		4.043 				4.043 						1 					6 		 	-6 		0 	0 	0 	0 	0 	0 	0 	# Million_Eggs/cm_b_Fem
# Male parameters
#		0.01 			0.15 			0.04390336 	 -3.12576436 		0.4384383 	3 			2 		0 	0 	0 	0 	0 	0 	0 	#	Hamel prior, maxage = 123
# 	-1 				1 				0 						0 						99 					0 		 -50  	0 	0 	0 	0 	0 	0 	0 	# L_at_Amin_Mal_GP_1
# 	40 				120 			63 						66 						99 					0 			2 		0 	0 	0 	0 	0 	0 	0 	# L_at_Amax_Mal_GP_1
# 	0.01 			0.2 			0.05 					0.05 					99 					0 			2 		0 	0 	0 	0 	0 	0 	0 	# VonBert_K_Mal_GP_1
# 	-3 				3 				0 						0 						99 					0 		 -3 		0 	0 	0 	0 	0 	0 	0 	# CV_young_Mal_GP_1
# 	0.5 			15 				2.5 					2.5 					99 					0 			5 		0 	0 	0 	0 	0 	0 	0 	# CV_old_Mal_GP_1
# 	-3 				3 				7.314623e-06 	7.314623e-06 	99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# Wtlen_1_Mal
# 	-3 				4 				3.237174			3.237174 			99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# Wtlen_2_Mal
#
  0 				2 				1 						1 						99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# RecrDist_GP_1
 -4 				4 				0 						0 						99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# RecrDist_Area_1
 -4 				4 		 	 -0.1 						0 						99 					0 		  3 		0 	0 	0 	0 	0 	0 	0 	# RecrDist_Area_2
  0 				2 				1 						1 						99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# RecrDist_Bseas_1
  0 				2 				1 						1 						99 					0 		-50 		0 	0 	0 	0 	0 	0 	0 	# CohortGrowDev
  0.000001 			0.999999  		0.5        				0.5  					0.5 				0 		-99 		0 	0 	0 	0 	0 	0 	0 	# FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 	0 	0 	0 	0 	0 	0 	0 	0 	0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
#_Spawner-Recruitment
3 #_SR_function: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepard_3Parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_LO          HI         INIT      PRIOR        PR_SD      PR_type   PHASE   env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn #  parm_name
  3           15         6.1          5           99            0      3       0      	0      	0       	0         0       0      0 			# SR_LN(R0)
  0.2          1         0.718      0.718      	0.158         2     -3       0      	0      	0       	0         0       0      0 			# SR_BH_steep
  0            5         0.5     	0.5       	99            0     -2       0      	0      	0       	0         0       0      0 			# SR_sigmaR
 -5            5         0          0           99            0     -50      0      	0      	0       	0         0       0      0 			# SR_regime
 -1            2         0          1           99            0     -50      0      	0      	0       	0         0       0      0 			# SR_autocorr

 1 		#do_recdev:  0=none; 1=devvector; 2=simple deviations
 1980 	# first year of main recr_devs; early devs can preceed this era
 2015 	# last year of main recr_devs; forecast devs start in following year
 7 		#_recdev phase 
 1 		# (0/1) to read 13 advanced options
 1889 	#_recdev_early_start (0=none; neg value makes relative to recdev_start)
 7			#_recdev_early_phase
 0			#_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 		#_lambda for Fcast_recr_like occurring before endyr+1
 1917 	#_last_early_yr_nobias_adj_in_MPD
 1975 	#_first_yr_fullbias_adj_in_MPD
 2013 	#_last_yr_fullbias_adj_in_MPD
 2014 	#_first_recent_yr_nobias_adj_in_MPD
 0.36 	#_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
 0 		#_period of cycles in recruitment (N parms read below)
 -5 		#min rec_dev
 5 		#max rec_dev
 0 		#_read_recdevs
# #
#1 #do_recdev:  0=none; 1=devvector; 2=simple deviations
#1916 # first year of main recr_devs; early devs can preceed this era
#1916 # last year of main recr_devs; forecast devs start in following year
#-8 #_recdev phase 
#1 # (0/1) to read 13 advanced options
# 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
# -8 #_recdev_early_phase
# -8 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
# 1 #_lambda for Fcast_recr_like occurring before endyr+1
# -1965 #_last_early_yr_nobias_adj_in_MPD
# -1970 #_first_yr_fullbias_adj_in_MPD
# -1990 #_last_yr_fullbias_adj_in_MPD
# -1995 #_first_recent_yr_nobias_adj_in_MPD
# 1 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
# 0 #_period of cycles in recruitment (N parms read below)
# -4 #min rec_dev
# 4 #max rec_dev
# 0 #_read_recdevs

#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1916R 1917F 1918F 1919F 1920F 1921F 1922F 1923F 1924F 1925F 1926F 1927F 1928F 1929F 1930F 1931F 1932F 1933F 1934F 1935F 1936F 1937F 1938F 1939F 1940F 1941F 1942F 1943F 1944F 1945F 1946F 1947F 1948F 1949F 1950F 1951F 1952F 1953F 1954F 1955F 1956F 1957F 1958F 1959F 1960F 1961F 1962F 1963F 1964F 1965F 1966F 1967F 1968F 1969F 1970F 1971F 1972F 1973F 1974F 1975F 1976F 1977F 1978F 1979F 1980F 1981F 1982F 1983F 1984F 1985F 1986F 1987F 1988F 1989F 1990F 1991F 1992F 1993F 1994F 1995F 1996F 1997F 1998F 1999F 2000F 2001F 2002F 2003F 2004F 2005F 2006F 2007F 2008F 2009F 2010F 2011F 2012F 2013F 2014F 2015F 2016F 2017F 2018F 2019F 2020F 2021F 2022F
#  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# implementation error by year in forecast:  0 0 0 0 0 0 0 0 0 0 0 0
#
#Fishing Mortality info 
0.09 # F ballpark
1999 # F ballpark year (neg value to disable)
1 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
0.9 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
#
#_initial_F_parms; count = 0
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
#
# F rates by fleet
# Yr:  1916 1917 1918 1919 1920 1921 1922 1923 1924 1925 1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# 1_CA_REC 0 0 0 0 0 0 0 0 0 0 0 0 0 4.81194e-005 7.78088e-005 0.000116104 0.000155086 0.000194134 0.000233179 0.000272284 0.000310817 0.000371248 0.000364186 0.000318687 0.000454094 0.000414544 0.000449942 0.000484663 0.000520173 0.000556401 0.00059335 0.000336733 0.000679471 0.000880362 0.00107377 0.00120696 0.0010767 0.000944318 0.00126808 0.00163023 0.00184205 0.00168416 0.00244479 0.00207537 0.00150757 0.00121018 0.00151226 0.0015835 0.00142901 0.0021646 0.00243174 0.00252311 0.00295523 0.0031126 0.00360407 0.00321969 0.00417701 0.00529263 0.00571588 0.00582153 0.00638788 0.00580109 0.00545396 0.00634743 0.00601608 0.00328856 0.00608016 0.00297015 0.00553786 0.00952133 0.0060757 0.00456315 0.0047589 0.00433399 0.00319963 0.00219636 0.00132728 0.000579289 0.00111868 0.000818205 0.00103327 0.00118398 0.000393324 0.000761123 0.000460896 0.00051168 0.000199075 0.000297472 4.74855e-005 7.08751e-005 0.000321336 0.000624155 0.000131336 0.000297076 0.000207981 0.000370495
# 2_CA_TWL 0.000160912 0.000264811 0.000310968 0.000158088 0.000174212 0.00016838 0.00015083 0.000161831 0.000206524 0.000282734 0.000356793 0.000433842 0.000404668 0.00041506 0.000495905 0.000412464 0.000596923 0.000326926 0.000424803 0.000587526 0.000594547 0.000447701 0.000468612 0.000474067 0.000337128 0.000394895 0.00024889 0.000435207 0.00183956 0.00433807 0.004295 0.00121591 0.00174199 0.000740735 0.000602064 0.00127525 0.00106399 0.000886303 0.000888051 0.000526937 0.000786226 0.000997439 0.00101745 0.000779665 0.000676751 0.000400601 0.000414612 0.000830007 0.000575514 0.000718535 0.000688331 0.000603558 0.000590082 0.00198432 0.00214652 0.0036185 0.00497891 0.00389979 0.00447262 0.0048168 0.00467257 0.00467081 0.012634 0.00825547 0.00353516 0.0143275 0.0132201 0.00544801 0.00468274 0.00107261 0.00295996 0.00481856 0.00583177 0.00459158 0.00731211 0.0133079 0.0101267 0.00484757 0.00512853 0.00470592 0.00700279 0.00629274 0.0020074 0.00214616 0.000366104 0.00039416 8.02477e-005 6.27736e-005 0.000232823 0.000304358 0.000164206 0.000422576 0.000150468 5.31057e-005 0.000159444 0.000219971
# 3_OR_REC 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000492478 0.000516763 0.00053325 0.000558051 0.000574624 0.000599316 0.000616276 0.000641261 0.000658463 0.000682749 0.00069869 0.000726298 0.000742374 0.000766938 0.00078325 0.000810484 0.00106179 0.00132265 0.00060113 0.00104024 0.000504269 0.00157774 0.00162045 0.00199792 0.00317554 0.0025951 0.00298561 0.00476271 0.00541933 0.00376489 0.00232682 0.00240355 0.00262487 0.000832595 0.0015071 0.00168805 0.00155762 0.00268736 0.00256687 0.00202453 0.00158767 0.000779167 0.00147257 0.00180971 0.0017395 0.000918972 0.000464255 0.000300463 0.000287527 0.000349533 0.000405346 0.000233594 0.000266082 0.000302008 0.000189638 0.000257795 0.000266715
# 4_OR_TWL 0.000131238 0.000138203 0.00014517 0.000152139 0.000158339 0.000165313 0.000172289 0.000179269 0.000185478 0.000192463 0.00019945 0.00020644 0.000344122 0.00057705 0.000514647 0.000403376 0.000138632 0.000212221 0.00024091 0.000220028 0.000514498 0.00058679 0.000566117 0.000296368 0.000865199 0.00124084 0.00185383 0.005142 0.00367333 0.00490481 0.00333884 0.00201678 0.00163994 0.00126153 0.00140975 0.00119749 0.00118067 0.000898635 0.00112987 0.00119466 0.00143581 0.00206057 0.00145891 0.00164725 0.00202329 0.00197128 0.00226686 0.000664887 0.000155827 0.00563449 0.000311501 0.000915863 0.000640191 0.00447018 0.000511415 0.00145651 0.00120256 0.00121688 0.00120765 0.000802763 0.00120852 0.00139489 0.00341767 0.00551464 0.00627932 0.00881337 0.0131046 0.0121281 0.0071128 0.0115429 0.00499044 0.00649348 0.0098078 0.0151789 0.00551732 0.0124736 0.0151626 0.0169563 0.00959141 0.0140105 0.00882357 0.011051 0.00399619 0.00591239 0.000351371 0.000598821 0.000149275 8.75909e-005 0.000252914 0.000159311 0.000273933 0.000306227 0.000360552 0.000148935 0.000207157 0.000286878
# 5_WA_REC 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000113114 0.000113146 0.000113178 0.000226421 0.000226507 0.000226593 0.000226679 0.000226764 0.000340275 0.000340497 0.000340716 0.000340933 0.000341148 0.000341361 0.000341571 0.000455704 0.000456067 0.00045648 0.000456925 0.000457444 0.000458002 0.000492796 0.00100974 0.000517545 0.000403517 0.000277454 0.000394313 0.000394644 0.000778537 0.00142021 0.00102646 0.0010528 0.00123043 0.000975408 0.00172032 0.00117185 0.00213654 0.00192858 0.00215156 0.00123584 0.0011902 0.00130071 0.00137535 0.00174066 0.00128253 0.001226 0.00151836 0.000450467 0.000316348 0.000449841 0.000631785 0.000206443 0.000302127 0.000291017 0.000197505 0.000314801 0.000281921
# 6_WA_TWL 0 0 0 0 0 0 0 0 0 0.000173573 0.000173591 0.000173609 0.000173627 0.000173645 0.000173662 0.000173679 0.000173695 0.000173711 0.000173727 0.000173743 0.000173758 0.000173773 0.000173787 0.000173802 0.000173816 0.000173829 0.000173843 0.000173856 0.000173869 0.000173882 0.000173895 0.000173907 0.00017392 0.000173932 0.000173944 0.000173957 0.000173969 0.000173982 0.000173995 0.000348016 0.000348114 0.000348213 0.000348311 0.000348444 0.000348577 0.00034871 0.000348842 0.000697951 0.000698435 0.000698916 0.000699394 0.000699868 0.00070034 0.000700807 0.00089412 0.00112477 0.00128398 0.00161947 0.00181521 0.00125178 0.0018178 0.0031599 0.0042349 0.00506422 0.00624966 0.00173554 0.00225663 0.00304677 0.00241156 0.00475631 0.0026994 0.00454224 0.00464222 0.00719552 0.00481022 0.00373973 0.00623787 0.00550922 0.00364005 0.00336681 0.00315317 0.003494 0.00104408 0.00617659 0.00148014 0.00411561 0.000292832 0.000185007 0.000124482 0.000139456 0.000143134 0.00030292 0.000146641 0.000216014 0.000170781 0.000236279
#
#_Q_setup
#_    fleet       link  link_info   extra_se    biasadj      float  #  fleetname
          3          1          0          1          0          1  #  3_CA_REC
          6          1          0          1          0          1  #  6_OR_REC
          7          1          0          1          0          1  #  7_WA_REC
          8          1          0          1          0          1  #  8_CACPFV
          9          1          0          1          0          1  #  9_OR_RECOB
         10          1          0          1          0          1  #  10_TRI_ORWA
         11          1          0          1          0          1  #  11_NWFSC_ORWA
         12          1          0          1          0          1  #  12_IPHC_ORWA
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_         LO           HI         INIT        PRIOR        PR_SD      PR_type     PHASE   env-var   use_dev  dev_mnyr  dev_mxyr    dev_PH     Block   Blk_Fxn  #  parm_name
           -15          15           0            0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_CA_REC(3)
          	0            5    		 0            0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_CA_REC(3)
           -15          15           0            0           99            0        -1         0         0         0         0         0         2         1 	#  LnQ_base_OR_REC(6)
          	0            5    		 0 		      0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_OR_REC(6)
           -20          15           0            0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_WA_REC(7)
          	0            5    		 0            0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_WA_REC(7)
 #          -15           15     -2.21024          0           99            0        -1         0         0         0         0         0         1         1 	#  LnQ_base_CACPFV(8)
           -15          15      -2.21024          0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_CACPFV(10)
          	0            5    		 0            0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_CACPFV(8)
           -15          15      -16.1653          0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_OR_RECOB(9)
          	0            5    		 0            0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_OR_RECOB(9)
#           -15           15     -2.21024          0           99            0        -1         0         0         0         0         0         1         1 	#  LnQ_base_TRI_ORWA(10)
           -15          15      -2.21024          0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_TRI_ORWA(10)
           	0            5    		 0            0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_TRI_ORWA(10)
           -15          15       -3.1832          0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_NWFSC_ORWA(11)
            0            5    		 0       	  0.01        99            0        -5         0         0         0         0         0         0         0 	#  Q_extraSD_NWFSC_ORWA(11)
           -15          15      -12.0696          0           99            0        -1         0         0         0         0         0         0         0 	#  LnQ_base_IPHC_ORWA(12)
            0            5    		 0         	  0.01        99            0         5         0         0         0         0         0         0         0 	#  Q_extraSD_IPHC_ORWA(12)
#_timevary Q parameters
#	HI	LO		init			PRIOR		PR_SD		PRIOR_TYPE PHASE
#  -4    4       -0.6    0       99      -1      1   # OR_rec deviation, to accomodate two rec dockside surveys
#  -4    4       -0.6    0       99      -1      1   # Triennial 1995 deviation
  -4    4       -0.6    0       99      -1      1   # Triennial 1995 deviation

#
#_size_selex_types
#discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead
#_Pattern Discard Male Special
24 0 0 0 #  1_CA_TWL  		
24 0 0 0 #  2_CA_NONTWL   
24 0 0 0 #  3_CA_REC  		
24 0 0 1 #  4_ORWA_TWL  	
24 0 0 2 #  5_ORWA_NONTWL 
24 0 0 0 #  6_OR_REC  		
24 0 0 0 #  7_WA_REC  		
15 0 0 3 #  8_CACPFV  		
24 0 0 0 #  9_OR_RECOB 	  
24 0 0 0 #  10_TRI_ORWA   
24 0 0 0 #  11_NWFSC_ORWA 
24 0 0 0 #  12_IPHC_ORWA  	  
#
#_age_selex_types
#_Pattern Discard Male Special
 10 0 0 0 #	1_CA_TWL  		  
 10 0 0 0 #	2_CA_NONTWL    
 10 0 0 0 #	3_CA_REC  		  
 10 0 0 0 #	4_ORWA_TWL  	  
 10 0 0 0 #	5_ORWA_NONTWL  
 10 0 0 0 #	6_OR_REC  		  
 10 0 0 0 #	7_WA_REC  		  
 10 0 0 0 #	8_CACPFV  		  
 10 0 0 0 #	9_OR_RECOB 	  
 10 0 0 0 #	10_TRI_ORWA    
 10 0 0 0 #	11_NWFSC_ORWA  
 10 0 0 0 #	12_IPHC_ORWA  	  

#_size_selex_settings
#_         LO           HI         INIT        PRIOR        PR_SD      PR_type     PHASE   env-var   use_dev  dev_mnyr  dev_mxyr    dev_PH     Block   Blk_Fxn  #  parm_name
#CA_TWL
           20           60           45           40          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           20           9           99            0         5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999         9           -999    	 -999  		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999         9           -999       	  9        	  99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#CA_NONTWL       
           20           60           45           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           20           9           99            0         5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999   		 -5     	  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999         9          	-999       	  9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#CA_REC        
           20           60           44           40          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30         	 20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999   		 -5   		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999         9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#ORWA_TWL       
           20           60           40           40          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30         	 20           9           99            0         5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999         9           -999     	 -999   	  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9          	-999       	  9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#ORWA_NONTWL       
           20           60           50           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999       	9           -999     	 -5    		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9          	-999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
##OR_REC        
           20           60           40           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           12           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999         -5    		  99            0      	 -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#WA_REC        
           20           60           50           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30         	 20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999     	 -5    		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#OR_observer
           20           60           40           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999     	 -5 		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#TRI_ORWA        
           20           80           50           30          99            0         4         0         0         0         0         0         0         0  #  PEAK     
           -15          4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1           9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1           30           12           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999     	 -5    		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#NWFSC_ORWA                                                                                      
           20           60           50           40          99            0         4         0         0         0         0         0         0         0  #  PEAK     
          -15           4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1        	9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1        	30         	 20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999    	 -5    	      99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999          9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
#IPHC_ORWA                                                             
           20           60           50           40          99            0         4         0         0         0         0         0         0         0  #  PEAK     
          -15           4           -15          -15          99            0        -5         0         0         0         0         0         0         0  #  TOP      
           -1        	9            6            6           99            0         4         0         0         0         0         0         0         0  #  ASC-WIDTH
           -1        	30         	 20           9           99            0        -5         0         0         0         0         0         0         0  #  DSC-WIDTH
           -999        	9           -999   		 -5    		  99            0        -4         0         0         0         0         0         0         0  #  INIT     
           -999        	9           -999       	  9           99            0        -5         0         0         0         0         0         0         0  #  FINAL    
   
#_timevary selex parameters
#CA_TWL
#           20           60           40           40          99            0         4         #  PEAK     
#           -15          4           -15          -15          99            0        -5         #  TOP      
#           -1           9            6            6           99            0         4         #  ASC-WIDTH
#           -1           30         	 20            9           99            0         -5         #  DSC-WIDTH
#           -5         	9           -5     		  -5    	   99            0        	4         #  INIT     
#           -5       	9            9         	  9        	  99            0         5         #  FINAL    
#CA_NONTWL       
#           20           60           40           30          99            0         4         #  PEAK     
#           -15          4           -15          -15          99            0        -5         #  TOP      
#           -1           9            6            6           99            0         4         #  ASC-WIDTH
#           -1           30         	 20            9           99            0         -5         #  DSC-WIDTH
#           -5         	9          	 9        	  9        	  99            0         5         #  FINAL    
##CA_REC                                                                                         
#           20           60           40           40          99            0         4         #  PEAK     
##           -15         4           -15          -15          99            0        -5         #  TOP      
#           -1           9            6            6           99            0         4         #  ASC-WIDTH
#           -1           9            9            9           99            0         5         #  DSC-WIDTH
#           -5     		 9            9            9      	   99            0         5         #  FINAL    
#ORWA_TWL
#            20           60           40           40          99            0         4         #  PEAK     
# #           -15          4           -15          -15         99            0        -5         #  TOP      
#            -1           9            6            6           99            0         4         #  ASC-WIDTH
#            -1           30         	 20            9           99            0         5         #  DSC-WIDTH
# #          -5       	9           -5     		 -5    		  99            0         4         #  INIT     
# #           -5       	9            9         	  9        	  99            0         5         #  FINAL    
#ORWA_NONTWL       
#            20           60           40           30          99            0         4         #  PEAK     
# #           -15          4           -15          -15          99            0        -5         #  TOP      
#            -1           9            6            6           99            0         4         #  ASC-WIDTH
#            -1           30         	 20            9           99            0         5         #  DSC-WIDTH
# #           -5       	9            9         	  9        	  99            0         5         #  FINAL    
###OR_REC        
#           20           60           30           30          99            0         4         #  PEAK     
##           -15         4           -15          -15          99            0        -5         #  TOP      
#           -1           9            6            6           99            0         4         #  ASC-WIDTH
#           -1           9            9            9           99            0         5         #  DSC-WIDTH
#           -5         	 9           -5     	  -5    	   99            0         4         #  INIT     
#           -5       	 9            9            9           99            0         5         #  FINAL    
##WA_REC        
#           20           60           30           30          99            0         4         #  PEAK     
##           -15         4           -15          -15          99            0        -5         #  TOP      
#           -1           9            6            6           99            0         4         #  ASC-WIDTH
#           -1           9            9            9           99            0         5         #  DSC-WIDTH
#           -5       	 9            9            9           99            0         5         #  FINAL    
#
#
0  #_ 0/1 to request experimental 2D_AR selectivity smoother options
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      3     1     1     0     0     0     0     0     0     0     0     0
#      3     3     1     0     0     0     0     0     0     0     0     0
#      3     5     1     0     0     0     0     0     0     0     0     0
#      3     8     1     0     0     0     0     0     0     0     0     0
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
4	1	0.8594
4	2	0.3079
4	3	0.5901
4	4	0.1938
4	5	0.5066
4	6	0.3287
4	7	0.8445
4	8	0.75
4	9	0.7836
4	10	0.4054
4	11	0.6474
4	12	1
5	2	0.55541
5	3	0.57354
5	4	0.48938
5	5	0.45
5	6	0.37193
5	7	1
5	11	1
5	12	0.025


 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset
# read 0 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark
#like_comp fleet  phase  value  sizefreq_method
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  1 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  1 #_CPUE/survey:_3
#  0 #_CPUE/survey:_4
#  1 #_CPUE/survey:_5
#  0 #_CPUE/survey:_6
#  1 #_CPUE/survey:_7
#  1 #_CPUE/survey:_8
#  1 #_CPUE/survey:_9
#  1 #_CPUE/survey:_10
#  1 #_CPUE/survey:_11
#  1 #_CPUE/survey:_12
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_lencomp:_5
#  1 #_lencomp:_6
#  1 #_lencomp:_7
#  1 #_lencomp:_8
#  1 #_lencomp:_9
#  1 #_lencomp:_10
#  1 #_lencomp:_11
#  1 #_lencomp:_12
#  1 #_agecomp:_1
#  1 #_agecomp:_2
#  1 #_agecomp:_3
#  1 #_agecomp:_4
#  1 #_agecomp:_5
#  1 #_agecomp:_6
#  0 #_agecomp:_7
#  0 #_agecomp:_8
#  1 #_agecomp:_9
#  0 #_agecomp:_10
#  1 #_agecomp:_11
#  0 #_agecomp:_12
#  1 #_init_equ_catch
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  1 # F_ballpark_lambda
0 # (0/1) read specs for more stddev reporting 
 # 0 1 -1 5 1 5 1 -1 5 # placeholder for selex type, len/age, year, N selex bins, Growth pattern, N growth ages, NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

