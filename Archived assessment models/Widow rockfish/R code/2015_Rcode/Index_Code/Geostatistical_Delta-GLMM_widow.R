
library(devtools)
install_github("nwfsc-assess/geostatistical_delta-GLMM", ref="3.2.0")  # This is the developement version.  Please check GitHub for the latest release number.
devtools::install_github("nwfsc-assess/nwfscDeltaGLM")
# Load libraries
library(TMB)
library(INLA)
library(SpatialDeltaGLMM)
library(nwfscDeltaGLM)
#############################
#############################
# File structure
my.wd <- "C:/Assessments/Widow2015/Data/Index standardization/4_strata_geo_strat_nwfsc"
setwd(my.wd)

# This is where all runs will be located
DateFile = paste(getwd(),'/',Sys.Date(),'-withpass/',sep='')
dir.create(DateFile)


# Settings
  Data_Set = "NWFSC"
  Version = "geo_index_v3b"
  n_x = c(250, 500, 1000, 2000)[3] # Number of stations
  FieldConfig = c("Omega1"=1, "Epsilon1"=1, "Omega2"=1, "Epsilon2"=1) # 1=Presence-absence; 2=Density given presence
  CovConfig = c("Depth_km"=0, "Depth_km2"=0, "N_km"=0) # DON'T USE DURING REAL-WORLD DATA FOR ALL SPECIES (IT IS UNSTABLE FOR SOME)
  Q_Config = c("Pass"=1)
  VesselConfig = c("Vessel"=0, "VesselYear"=1)
  ObsModel_Set = c(1,2,11,12)  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
  Aniso = 1 # 0=No; 1=Yes
  Kmeans_Config = list( "Locs"=c("Samples","Domain")[1], "nstart"=1000, "iter.max"=1e3)     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid
  CovConception = FALSE # Include Point Conception as sampling intensity stratum
  Options = c( 0, 0 )  # First slot: Encounter; Second slot: Positive density; 0: Don't ADREPORT() linear predictors; 1: Do ADREPORT() them
  ConvergeTol = c(1,2) # 1:Default; 2:Increased; 3:High
  BiasCorr = FALSE    # TRUE crashes R

# Save settings
  Settings = c("Data_Set"=Data_Set, "Version"=Version, "n_x"=n_x, "FieldConfig"=FieldConfig, "CovConfig"=CovConfig, "Q_Config"=Q_Config, "VesselConfig"=VesselConfig, "ObsModel_Set"=ObsModel_Set, "Aniso"=Aniso, "Kmeans_Config"=Kmeans_Config, "CovConception"=CovConception, "Options"=Options, "ConvergeTol"=ConvergeTol, "BiasCorr"=BiasCorr)
  save( Settings, file=paste0(DateFile,"Settings.RData"))
  capture.output(Settings, file=paste0(DateFile,"Settings.txt"))

# Compile TMB software
  setwd( system.file("executables", package="SpatialDeltaGLMM") )
  #dyn.unload( dynlib(Version) )
  #file.remove( paste0(Version,c(".o", ".dll")) )
  compile( paste(Version,".cpp",sep="") )

  if(Data_Set=="NWFSC"){
    #data( Example_NWFSC_shelf_slope_trawl )
    #NWFSC_Trawl <- Example_NWFSC_shelf_slope_trawl
    NWFSC_Trawl <- read.csv("C:/Assessments/Widow2015/Data/TrawlSurvey/NWFSCsurvey/Widow2015_HaulCatchWt&Effort - Copy.csv", skip = 8, header = T)
    NWFSC_Trawl[,ncol(NWFSC_Trawl)+1]<-as.numeric(substr(as.character(NWFSC_Trawl$PROJECT_CYCLE),7,10))
    colnames(NWFSC_Trawl)[ncol(NWFSC_Trawl)]<-"YEAR"
    NWFSC_Trawl= subset(NWFSC_Trawl,YEAR>2002)
    NWFSC_Trawl[,'VESSEL'] = as.factor( as.character(NWFSC_Trawl[,'VESSEL']) )
    Data_Geostat = cbind( "Catch_KG"=NWFSC_Trawl[,'HAUL_WT_KG'], "Year"=as.numeric(sapply(NWFSC_Trawl[,'PROJECT_CYCLE'],FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]})), "Vessel"=NWFSC_Trawl[,"VESSEL"], "AreaSwept_km2"=NWFSC_Trawl[,"AREA_SWEPT_HA"]/1e2, "Lat"=NWFSC_Trawl[,'BEST_LAT_DD'], "Lon"=NWFSC_Trawl[,'BEST_LON_DD'],"Pass"=NWFSC_Trawl[,'SURVEY_PASS']-1.5)
    # 
  }

# Load data and strata
#strata.limits <- nwfscDeltaGLM::readIn(ncol=5,nlines=3)
#  STRATA  NLat SLat MinDepth MaxDepth
#  shallow 49.0 34.5  55       183
#  deep    49.0 34.5  183      400

strata.limits <- readIn(ncol=5,nlines=6)
  STRATA    NLat SLat MinDepth MaxDepth
  Coastwide 49.0 34.5 55       183
  Sshallow 40.5 34.5  55       183
  Sdeep    40.5 34.5  183      400
  Nshallow 49.0 40.5  55       183
  Ndeep    49.0 40.5  183      400

################
# Prepare data
################

# Read extrapolation data
  data( extrapolation_data )
  Data_Extrap <- extrapolation_data

# Augment with strata for each extrapolation cell
  Tmp = cbind("BEST_DEPTH_M"=(-1000)*Data_Extrap[,'Depth_km'], "BEST_LAT_DD"=Data_Extrap[,'Lat'], "propInWCGBTS"=Data_Extrap[,'propInWCGBTS'])
  a_el = as.data.frame(matrix(NA, nrow=nrow(Data_Extrap), ncol=nrow(strata.limits)))
  for(l in 1:ncol(a_el)){
    a_el[,l] = apply(Tmp , MARGIN=1, FUN=nwfscDeltaGLM::strata.fn, Strata.df=strata.limits[l,])
    a_el[,l] = ifelse( is.na(a_el[,l]), 0, 4*Data_Extrap[,'propInWCGBTS'])
  }

# Convert extrapolation-data to an Eastings-Northings coordinate system
  Tmp = cbind('PID'=1,'POS'=1:nrow(Data_Extrap),'X'=Data_Extrap[,'Lon'],'Y'=Data_Extrap[,'Lat'])
  attr(Tmp,"projection") = "LL"
  attr(Tmp,"zone") = "10"
  tmpUTM = convUL(Tmp)                                                         #$
  Data_Extrap = cbind( Data_Extrap, 'E_tmp'=tmpUTM[,'X'], 'N_tmp'=tmpUTM[,'Y'], 'Include'=(Data_Extrap[,'Cowcod']==0 & Data_Extrap[,'Ngdc_m']<(-35)))

# Read or simulate trawl data
  #NWFSC_Trawl = read.csv( paste0(DataFile,"FisheryIndices2015_Canary_V2--HaulCatchWt&Effort.csv"), header=TRUE, skip=8 )
  NWFSC_Trawl = cbind(NWFSC_Trawl, "AREA_SWEPT_MSQ"=NWFSC_Trawl[,'AREA_SWEPT_HA']*1e4)

# Add PASS
  Date = as.character(NWFSC_Trawl[,'CAPTURE_DATE'])
  NWFSC_Trawl = cbind(NWFSC_Trawl, "Date"=as.Date(sapply( Date, FUN=function(Char){ paste( c('2000',strsplit(Char,'/')[[1]][1:2]),collapse="-") })) )
  NWFSC_Trawl = cbind(NWFSC_Trawl, "PASS"=NWFSC_Trawl[,'SURVEY_PASS'])
  range( NWFSC_Trawl[which(NWFSC_Trawl[,'PASS']==1),'Date'] )

# Make simplified data frame
  Data_Geostat = cbind( "Catch_KG"=NWFSC_Trawl[,'HAUL_WT_KG'], "Year"=as.numeric(sapply(NWFSC_Trawl[,'PROJECT_CYCLE'],FUN=function(Char){strsplit(as.character(Char)," ")[[1]][2]})), "Vessel"=NWFSC_Trawl[,"VESSEL"], "AreaSwept_km2"=NWFSC_Trawl[,"AREA_SWEPT_HA"]/1e2, "Lat"=NWFSC_Trawl[,'BEST_LAT_DD'], "Lon"=NWFSC_Trawl[,'BEST_LON_DD'], "Pass"=NWFSC_Trawl[,'PASS']-1.5)
  Year_Set = sort(unique(Data_Geostat[,'Year']))

# Convert to an Eastings-Northings coordinate system 
  Tmp = cbind('PID'=1,'POS'=1:nrow(Data_Geostat),'X'=Data_Geostat[,'Lon'],'Y'=Data_Geostat[,'Lat'])
  attr(Tmp,"projection") = "LL"
  attr(Tmp,"zone") = "10"
  tmpUTM = convUL(Tmp)                                                         #$ 
  Data_Geostat = cbind( Data_Geostat, 'E_tmp'=tmpUTM[,'X'], 'N_tmp'=tmpUTM[,'Y'])
  # Find nearest knot in extrapolation data, and define covariates via nearest knot
  NN_Tmp = nn2( data=Data_Extrap[,c('E_tmp','N_tmp')], query=Data_Geostat[,c('E_tmp','N_tmp')], k=1 )
  Data_Geostat = cbind( Data_Geostat, Data_Extrap[NN_Tmp$nn.idx,c('Depth_km','Depth_km2','E_km','N_km')] )

# Calculate k-means centroids (but only once for all species)
  setwd( DateFile )
  Kmeans = Calc_Kmeans(n_x=n_x, Kmeans_Config=Kmeans_Config, Data_Geostat=Data_Geostat, Data_Extrap=Data_Extrap)
  save( Kmeans, file=paste0("Kmeans-",n_x,".RData") )
  loc_x = Kmeans$centers

# Calculate areas and average characteristics
  Voronoi = calcVoronoi( xydata=cbind("X"=loc_x[,'E_km'],"Y"=loc_x[,'N_km']), xlim=range(Data_Extrap[,'E_km']), ylim=range(Data_Extrap[,'N_km']))
  NN = nn2( data=loc_x[,c('E_km','N_km')], query=Data_Geostat[,c('E_km','N_km')], k=1 )
  NN_Extrap = nn2( data=loc_x[,c('E_km','N_km')], query=Data_Extrap[,c('E_km','N_km')], k=1 )
  a_xl = matrix(NA, ncol=ncol(a_el), nrow=n_x)
  for(l in 1:ncol(a_xl)) a_xl[,l] = tapply(a_el[,l], INDEX=NN_Extrap$nn.idx, FUN=sum)
  
# Make design matrix (X_xj)
  X_xj = NULL
  for( j in 1:length(CovConfig) ){
    Which = names(CovConfig)[j]
    if(Which=="Depth_km2") X_xj = cbind(X_xj, tapply( Data_Extrap[,"Depth_km"], INDEX=NN_Extrap$nn.idx, FUN=mean, MARGIN=2)^2)
    if(Which!="Depth_km2") X_xj = cbind(X_xj, tapply( Data_Extrap[,Which], INDEX=NN_Extrap$nn.idx, FUN=mean, MARGIN=2))
    colnames(X_xj)[j] = Which
  }
  X_xj = apply(X_xj, MARGIN=2, FUN=function(Vec){(Vec-mean(Vec))/(sd(Vec))})
  # Drop elements from X_xj
  X_xj = X_xj[,which(CovConfig==1)]
  # Add Point Conception if necessary
  if(CovConception==TRUE) X_xj = cbind( X_xj, 'S_of_Concep'=ifelse( tapply( Data_Extrap[,'Lat'], INDEX=NN_Extrap$nn.idx, FUN=mean, MARGIN=2)<34.5, 1, 0) )
  # Make filler matrix if necessary
  if(ncol(X_xj)==0) X_xj = cbind( "Dummy"=rep(0,n_x) )

# Make catchability matrix (Q_i)
  if( sum(Q_Config)==0 ){
    Q_ik = matrix(0, ncol=1, nrow=nrow(Data_Extrap))
  }else{
    Q_ik = as.matrix(Data_Geostat[,names(Q_Config)[which(Q_Config==1)],drop=FALSE])
  }
    
# Make mesh and info for anisotropy
  MeshList = Calc_Anisotropic_Mesh(loc_x=loc_x)

################
# Make and Run TMB model
################

for(ConfigI in 1:length(ObsModel_Set)){
  # Settings
  ObsModel = ObsModel_Set[ConfigI]
  
  # Config file
  ConfigFile = paste0(DateFile,"ObsModel=",ObsModel,"/")
  dir.create(ConfigFile)
  
  # Covariates
  TmbData = Data_Fn("Aniso"=Aniso, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=NN$nn.idx[,1]-1, "t_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "MeshList"=MeshList)

  # Parameters
  Parameters = Param_Fn( TmbData )

  # Which are random
  Random = c("Epsiloninput1_st", "Omegainput1_s", "Epsiloninput2_st", "Omegainput2_s", "nu1_v", "nu2_v", "nu1_vt", "nu2_vt")
  # Random = c( Random, "logSigmaM")

  # Which parameters are turned off
  Map = Make_Map( VesselConfig=VesselConfig, TmbData=TmbData, FieldConfig=FieldConfig, CovConfig=CovConfig, CovConception=CovConception, ObsModel=ObsModel, Aniso=Aniso)
  #Map[["logSigmaM"]] = factor( c(1,NA,2,NA) )

  # Build object                                                              
  dyn.load( paste0(system.file("executables", package="SpatialDeltaGLMM"),"/",dynlib(Version)) )
  #dyn.load( paste0("C:/Users/James.Thorson/Desktop/Project_git/geostatistical_delta-GLMM/inst/executables/",Version) )
  if(any(FieldConfig!=0)|any(VesselConfig!=0)){
    Obj <- MakeADFun(data=TmbData, parameters=Parameters, inner.method=c("newton","BFGS")[1], random=Random, hessian=FALSE, map=Map)
  }else{
    Obj <- MakeADFun(data=TmbData, parameters=Parameters, hessian=FALSE, map=Map)
  }
  Obj$control <- list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100)
  #obj$env$inner.control <- c( obj$env$inner.control, list("power"=0.01, "u0"=0.1, "ustep"=1, "silent"=FALSE))  
  
  # Reconfigure Obj$gr
  if(TRUE){
    Obj$gr_mod <- function(x,...){
      save_gr = list( "x"=x )
      save( save_gr, file=paste0(DateFile,"save_gr.RData"))
      Obj$gr(x=x,...)
    }
  }else{
    Obj$gr_mod <- Obj$gr
  }
  
  # Run first time
  Obj$fn(Obj$par)
  Obj$gr_mod(Obj$par)

  # Declare upper and lower bounds for parameter search
  Lower = rep(-50, length(Obj$par))
  Lower[grep("logsigmaV",names(Obj$par))] = log(0.01)
  Upper = rep( 50, length(Obj$par))
  Upper[grep("logtau",names(Obj$par))] = 10   # Version < v2i
  Upper[grep("logeta",names(Obj$par))] = log(1/(1e-2*sqrt(4*pi))) # Version >= v2i: Lower bound on margSD = 1e-4
  Upper[grep("SigmaM",names(Obj$par))] = 10 # ZINB can crash if it gets > 20
  if( "gamma1" %in% names(Obj$par) ){
    Lower[grep("gamma1",names(Obj$par))] = -20
    Upper[grep("gamma1",names(Obj$par))] = 20
  }
  if( "gamma2" %in% names(Obj$par) ){
    Lower[grep("gamma2",names(Obj$par))] = -20
    Upper[grep("gamma2",names(Obj$par))] = 20
  }
  if( "lambda1" %in% names(Obj$par) ){
    Lower[grep("lambda1",names(Obj$par))] = -20
    Upper[grep("lambda1",names(Obj$par))] = 20
  }
  if( "lambda2" %in% names(Obj$par) ){
    Lower[grep("lambda2",names(Obj$par))] = -20
    Upper[grep("lambda2",names(Obj$par))] = 20
  }

  # Change convergence tolerance
  Obj$env$inner.control$step.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol[1]] # Default : 1e-8  # Change in parameters limit inner optimization
  Obj$env$inner.control$tol10 <- c(1e-6,1e-8,1e-12)[ConvergeTol[1]]  # Default : 1e-3     # Change in pen.like limit inner optimization
  Obj$env$inner.control$grad.tol <- c(1e-8,1e-12,1e-15)[ConvergeTol[1]] # # Default : 1e-8  # Maximum gradient limit inner optimization

  # Run model
  Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr_mod, lower=Lower, upper=Upper, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=c(1e-8,1e-10,1e-14)[ConvergeTol[2]]))  # , rel.tol=1e-20
  Opt[["final_diagnostics"]] = data.frame( "Name"=names(Opt$par), "Lwr"=Lower, "Est"=Opt$par, "Upr"=Upper, "Gradient"=Obj$gr(Opt$par) )
  Opt[["num_par"]] = length(Opt$par)
  Opt[["AIC"]] = 2*Opt$objective + 2*Opt[["num_par"]]
  #numGr = numDeriv::grad( func=Obj$fn, x=Opt$par )
  # load(file=paste0(DateFile,"save_gr.RData")); par<-save_gr$x
  
  # Debugging
  if(FALSE){
    Report = Obj$report( Obj$env$last.par )
  }
    
  # Reports
  Report = Obj$report()                                      
  Sdreport = sdreport(Obj, bias.correct=BiasCorr)
  #Hess = optimHess( par=Opt$par, fn=Obj$fn, gr=Obj$gr )
  
  # Save stuff
  Save = list("Opt"=Opt, "Report"=Report, "Sdreport"=Sdreport)
  save(Save, file=paste0(ConfigFile,"Save.RData"))
  capture.output( Opt, file=paste0(ConfigFile,"Opt.txt"))
  capture.output( summary(Sdreport), file=paste0(ConfigFile,"summary-Sdreport.txt"))
  file.copy( from=paste0(system.file("executables", package="SpatialDeltaGLMM"),"/",dynlib(Version)), to=paste0(ConfigFile,"Version.cpp"), overwrite=TRUE)
  
################
# Make diagnostic plots
################

  # Plot Anisotropy  
  if(Aniso==1){
    PlotAniso_Fn( FileName=paste0(ConfigFile,"Aniso.png"), Report=Report )
  }

  # Plot surface
  PlotMap_Fn(MappingDetails=list("state", c("Oregon","Washington","California")), Report=Report, MapSizeRatio=c("Height(in)"=4,"Width(in)"=1.55), Xlim=c(-126,-117), Ylim=c(32,49), FileName=paste0(ConfigFile,"Field_"), Year_Set=Year_Set, Rotate=20, mfrow=c(3,4), mar=c(0,0,2,0), oma=c(3.5,3.5,0,0))

  # Covariate effect
  PlotCov_Fn(Report=Report, NN_Extrap=NN_Extrap, X_xj=X_xj, FileName=paste0(ConfigFile,"Cov_"))
  
  # Time series measures
  Timeseries_Fn(Report=Report, Year_Set=Year_Set, FileName=paste0(ConfigFile,"Summaries.png"))
  
  # Positive catch rate Q-Q plot
  Q = QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(ConfigFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(ConfigFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(ConfigFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(ConfigFile,"Q-Q_hist.jpg"))

  # Vessel effects
  Return = Vessel_Fn(TmbData=TmbData, Sdreport=Sdreport, FileName_VYplot=paste0(ConfigFile,"VY-effect.jpg"))

  # Plot index
  png( file=paste0(ConfigFile,"Index.png"), width=4, height=4, res=200, units="in")
    par( mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02, yaxs="i")
    log_Index = array( summary(Sdreport)[which(rownames(summary(Sdreport))=="ln_Index_tl"),], dim=c(unlist(TmbData[c('n_t','n_l')]),2), dimnames=list(NULL,NULL,c('Estimate','Std. Error')) )
    Index = array( summary(Sdreport)[which(rownames(summary(Sdreport))=="Index_tl"),], dim=c(unlist(TmbData[c('n_t','n_l')]),2), dimnames=list(NULL,NULL,c('Estimate','Std. Error')) )
    plot(1, type="n", xlim=range(Year_Set), ylim=1.05*c(0,max(exp(log_Index[,,'Estimate']+1*log_Index[,,'Std. Error']))), xlab="Year", ylab="Abundance" )
    for(l in 1:dim(Index)[2]){
      lines( y=Index[,l,'Estimate'], x=Year_Set+seq(-0.1,0.1,length=dim(Index)[2])[l], type="b", col=rainbow(TmbData[['n_l']])[l] )
      for(t in 1:dim(Index)[1]){
        lines( x=rep(Year_Set[t],2)+seq(-0.1,0.1,length=dim(Index)[2])[l], y=exp(log_Index[t,l,'Estimate']+c(-1,1)*log_Index[t,l,'Std. Error']), col=rainbow(TmbData[['n_l']])[l])
      }
    }
    # Write to file
    Names = strata.limits[,'STRATA']
    Table = data.frame( "Year"=Year_Set, "Unit"=1, "Fleet"=rep(Names,each=dim(Index)[1]), "Estimate"=as.vector(Index[,,'Estimate']), "SD"=as.vector(log_Index[,,'Std. Error']) )
    write.csv( Table, file=paste0(ConfigFile,"Table_for_SS3.csv"), row.names=FALSE)
  dev.off()
}
