#####
## JJCode: 22/02/2018.
## LRE: Mar2018
##  MERGED code to run the functions in: 
##   "/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_Calculate_ImplausibilityMeasure_SingleObs.r"
##   "/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_Calculate_NMAEF_NMBF_SingleObs.r"
## Over a set of observations for a given variable/month
##  -> Initially assume that the obs table is generated and the emulation has been done prior to running this code
##  
## JJ: 19/11/18: Code tidied for sharing...
##  The "OnJASMIN=" only applying to PPE output file... *** Other filepaths would need adapting for new user areas/JASMIN... ***
###
## source("/nfs/a173/earjsj/constraint/R_code_for_constraint/RunCode_Calculate_ImplausMeasure_NMAEF_NMBF_2627_LSUnifLSPrior_code_embedded.r",echo=TRUE)
## source("/nfs/a173/earjsj/constraint/R_code_for_constraint/RunCode_Calculate_ImplausMeasure_NMAEF_NMBF_2627_LSUnifLSPrior_code_embedded.r")
## system.time(source("/nfs/a173/earjsj/constraint/R_code_for_constraint/RunCode_Calculate_ImplausMeasure_NMAEF_NMBF_2627_LSUnifLSPrior_code_embedded.r"))
## system.time(source("/nfs/a173/earjsj/constraint/R_code_for_constraint/RunCode_Calculate_ImplausMeasure_NMAEF_NMBF_2627_LSUnifLSPrior_code_SLembedded.r",echo=TRUE))
###
## DSEm method used consistently throughout
###############################################################################################
###############################################################################################

#################################
## Source in the R librarys and JJcodes to use:
#################################

  library(ncdf4)
  library(DiceKriging)
  source("/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_EmPred_LargeSample.r")
  source("/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_GetN48GBvalue_FromN96OutputMat.r")
  source("/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_MakeAve_GBplusSurrNeighbours.r")


######################################################
## Remove old netCDF files?
## Avoids network clashes with removed files seeming to still exist
######################################################

  RemNetCDF = TRUE
#  RemNetCDF = FALSE


################################
## PPE Source definition. Many dependencies
################################

  PPEnameIn = "26aer"  ## Alternatives "26aer" / "27aeratm"
  PPE_EmsYear = 2008 ## Emissions year of output. e.g. 2008, 19782008, 1850, 18502008, etc.
  MonthVec = c("jan") ## ONE MONTH AT A TIME IS PLENTY FOR A SINGLE JOB 
  GridResolutionIn = "N48"  ## could be "N96"?

  if(GridResolutionIn=="N48"){
    Resolution_text="Lores"
  }else{
    Resolution_text="Hires"
  }


##################################
## Source code for emulation
##################################

## LRE May2018. Replaced JJs original ExtractEmulate file.
## Can now accomodate 26aert/27aeratm structures and paths
  source("/nfs/a134/mm11lr/constraint/R_code_for_constraint/JJCode_N48_ExtractEmulate_ObservationGBs_unif_switch.r")

## DS method of fitting linear models prior to emulation. Used consistenly. Always.
  source("/nfs/see-fs-01_users/earjsj/CSSP/UKCA_LargePPE_Jan15/EmulationSA_Codes/JJCode_RealObsConstraint_CSSP/JJCode_FitEmulators_DSEmMethod_Validate_SA.r")

  if(PPEnameIn=="27aeratm"){
    Tr_Val_text= "UKCA27aeratmos"
    PPE_DataFileNameMidIn = "_xmcta-xmdbi_pm2006"
    PPE_DataFileNameEndIn = paste("_191_",GridResolutionIn,".nc",sep="")
  }else{
    Tr_Val_text= "UKCA26aerPPE"
    PPE_DataFileNameMidIn = "_tebaa-tebiz_teafw_pm2008"
    PPE_DataFileNameEndIn = paste("_",GridResolutionIn,".nc",sep="")
  }


###################################################
## Identify the case to run the code for and set up the specifications for it.
## Remember to make directories for VariableIn
###################################################

#  VariableIn = "N50"
#  VariableIn = "N100"
#  VariableIn = "CCN0p2"
#  VariableIn = "Ntot" ## OR "aer_corsol_conc" ## OR "N700"  Ntot no good b/c no .nc files made; aero_corsol_conc no good b/c entire column average in kg/m3
#  VariableIn = "N3"
#  VariableIn = "nconc_corsol"  # Folders marked _ACESPACE with symbolic links
  VariableIn = "PM2p5"    
#  VariableIn="AODs_TOTAL"
#  VariableIn= "H2SO4_TOTAL_conc"
#  VariableIn= "OC_TOTAL_conc"

  if(VariableIn=="N50" | VariableIn=="Ntot" | VariableIn=="N100" | VariableIn=="N3" | VariableIn=="CCN0p2" | VariableIn=="nconc_corsol" | VariableIn=="H2SO4_TOTAL_conc"){
    VariableUnits = "cm^-3"
  }else if(VariableIn=="AODs_TOTAL"){
    VariableUnits = ""
  }else if(VariableIn=="PM2p5" |  VariableIn=="OC_TOTAL_conc" | VariableIn=="BC_TOTAL_conc"){
    VariableUnits = "ug m^-3"
  }else{
    VariableUnits = "!needs_units!"
  }


################################
## Observation source text
##  -> JJ: Always use "" here...
##    -> This only applies to the name of the variable folder where the "obsID table" is stored? This is not feeding through to everything.
##    -> 'ObsTab_MultiSetIDIn' below allows different observation sets of the same variable...
################################

  Obs_text= ""   


############################
## HPC, res/olution and vertical level
##  -> JJ: 19/11/18: The "OnJASMIN" only applying to PPE output file... *** Other filepaths would need adapting for new user areas/JASMIN... ***
############################

  OnJASMIN = FALSE   ## Is the code run on JASMIN
  SelectVLIn = TRUE  ## Does a Vertical Level need to be selected; alternatively the netCDF file will be a single VL
  ReScaleOutputIn = FALSE  ## e.g. for converting m-3 to cm-3
  ReScaleDivisorIn = 1
 #      ReScaleOutputIn = TRUE
 #      ReScaleDivisorIn = 1000000  ## take from per m^3 to per cm^3 ## with KP data?? 
  UseSurrNeighbourAveIn = FALSE ## This is for N96 resolution, so FALSE for N48

  if(VariableIn=="AODs_TOTAL"){
    VLnoIn = 2
    VLIDIn = "VL2"
  }else{
    VLnoIn = 1
    VLIDIn = "VL1"
  }

  if(SelectVLIn==FALSE){
     VLnoIn = NA
     VLIDIn = "NA"
  }


#######################################################
## Make the percentage uncertainty specifications for the implausibility calculations:
#######################################################

  ObsMeasUnc_PlusMinusPercIn = 10
  SpatialCoLocUnc_PlusMinusPercIn = 20
  TemporalCoLocUnc_PlusMinusPercIn = 10  ##10 ##50 for Ntot_ACESPACE  ##27 for CCN  ##10 LRE: the 27% comes from the ratio of s.d. to mean. see 'CCN_convert_to_netcdf_and_dat.py'


################################
## Structural uncertainty
## Switch for including extra uncertainty term
## within embedded code.
## Both PlusMinus and ValueIn need to be set.
## Both used in embedded code
################################

  StructVarTerm_AsPercOnModelAveIn = FALSE     ## IS STRUCTURAL UNCERTAINTY QUANTIFIED AND INCLUDED
#  StructVarTerm_AsPercOnModelAveIn = TRUE
  StructVarTerm_PlusMinusPercIn = 0
#    StructVarTerm_PlusMinusPercIn = 10
  StructVarTermValueIn = 0


####################
## Switch for Large Sample (LS) type
## 1= uses prior pdfs
## 2= uniform
####################

  LS_switch=2 ## See comment above

  if(LS_switch==1){
    EmPredLS_text= "EmPredLSPrior"
    LS_text= "LSPrior"
  }else{
    EmPredLS_text="EmPredLSUnif"
    LS_text= "LSUnif"
  }


#############################
## Paths for storing output
## and locating observations
## NOTE: PPE extension for storage
## when PPE_specific output req'd
#############################

#  user_text= "a134/mm11lr"  ## This would equate to /nfs/a134/mm11lr/constraint/ etc.
  user_text= "a173/earjsj"  ## This would equate to /nfs/a173/earjsj/constraint/ etc.

  FnOutputStoragePathStartInPPE = paste("/nfs/",user_text,"/constraint/",PPEnameIn,"/",sep="")

  FnOutputStoragePathStartIn = paste("/nfs/",user_text,"/constraint/",sep="")
##  FnOutputStoragePathStartIn = "/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/27aeratm/real_obs_constraint/"  ## JJ version for CSSP

## Observations may be in a different folder to other output. If so, use symbolic links

#  Local_Obs_Path_Start = "/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/27aeratm/real_obs_constraint/"
  Local_Obs_Path_Start = paste("/nfs/",user_text,"/constraint/",sep="")


###################################################################
## Define the storage file path and the start of the ObsIDTab filename that will work from
## the ObsIDTab filename is of form: "N48_ObsLocTab_PM2p5Test1_jul_VL1.dat". 
## Start is the part before 'ObsLocTab'; the rest comes from the set-up inputs above
## ObsTab can have rows removed if want to run over fewer than total observations
## BUT, read in 'EmulatorList' or 'MeanPredTab' is assumed to have all observations
## Fine if making within this file
##
## NOTE: Obs data needs to be copied (or symbolically linked) into a directory with write access
## b/c this directory is used elsewhere as a starting point for writing to file
###################################################################

  ObsTab_MultiSetIDIn = "_GASSP" ## For N50, PM2p5 data from Jo.
#  ObsTab_MultiSetIDIn = "_440nm_ave0515" ## AOD 440nm wavelength, averaged over the period 2005-2015.
##  ObsTab_MultiSetIDIn = "_GASSP_CSSP"  ##For changing filenames in input. Measures in place to account for empty.
##  ObsTab_MultiSetIDIn = "_Test_mini" ## May run different obs files for same variable ('surface'? or 'TOA' etc? )
#  ObsTab_MultiSetIDIn = "_GASSPT4" ## May run different obs files for same variable ('surface'? or 'TOA' etc? )


  ObsIDLocTab_FileNameStartIn = GridResolutionIn
##  ObsIDLocTab_FileNameStartIn = paste(GridResolutionIn,"_MergedObs",sep="")


#####################################################################
##  FLAG: Does the observation Table have "all observations" (and match EmList or EmPredTab for size)?
##   If a subset, which are removed?
##   If starting from scratch and emulating here, then this is TRUE
##
##   If TRUE, and reading in EmList or EmPred, assume 'EmulatorList' positions or 'MeanPredTab' 
##       positions correspond in exact order to rows of tab
##  -> If FALSE, then assumes that 'EmulatorList' or 'MeanPredTab' has all obs and selects Emulator or 
##       PredTab Column according to obsID in list
#####################################################################

  AllObsInTab = TRUE
#  AllObsInTab = FALSE


#########################################################
## FLAG: Emulator prediction required in observation gridbox locations?
## i.e. sample from emulator?
## Could be previously created and read in, so that percentage uncertainties can change
##  ->  only takes around 10 mins to run if already predicted
#########################################################

  EmPredRequired = TRUE
#  EmPredRequired = FALSE


#########################################################
## FLAG: Retain emulator prediction (sample) for observational gridboxes? 
## Only used if EmPredRequired==TRUE
#########################################################

  if(EmPredRequired==TRUE){
    SaveEmPredToFile = TRUE  ## Makes netcdf filename and file for sumary statistics of prediction (and predictions themselves - big!)
    RetainEmSDPred = TRUE # Save the StDev on each emulator prediction as well at the mean prediction...
    MakeEmPredMeanHistPlots = TRUE ## make histogram plots looking at the mean emulator predictions over the Large sample uncertainty for each observation.
  }else{
    SaveEmPredToFile = FALSE ## Already saved to file if being read in!
    RetainEmSDPred = FALSE
    MakeEmPredMeanHistPlots = FALSE
  }

#########################################################################
## Emulators/Emulation:  
##   Either a) Set all of the required inputs/parameters for emulation
##   b) set up to read in the fitted emulators from file (e.g. "Emulators_N48_PM2p5_DSEm_jul_VL1.RData")
##   set by 'FitEmulators'.
##  'EmMultiSetIDIn' is a fname extension for netcdf output files
##  'UseDSEmIn' sets emulation method to David Sexton's linear regression
##  alternative is the dicekriging fit
## NOTE: JJ and LR agreed on 6thJun2018 to always use DSEm
## Therefore 'UseDSEmIn' is always TRUE. Options retained in code
## but EmMultiSetIDIn now always empty "" because it's currently only applied to implausibility filenames
## and if used should be applied to all other output files. Removing the best option.
##
##  'ReEmulateAllRunsIn' TRUE if validation runs to be used for new emulator
##  'SaveEmulatorsToFileIn' FLAG: Save emulators? Good idea generally.
#########################################################################

  FitEmulators = TRUE    ## Do new emulators need to be fitted
#  FitEmulators = FALSE
  UseDSEmIn = TRUE       ## Use David Sexton's linear regression proceedure pre-emulation
  EmMultiSetIDIn = ""    ## Remember: DSEm method applied consistently, but no extension used
  ReEmulateAllRunsIn = TRUE     ## Combine validation and training simulations to make an emulator with presumably space-filling properties
  SaveEmulatorsToFileIn = TRUE  ## Are you wanting to save the emulators for later use


################################
## Want to calculate NMAEF and NMBF also? 
################################

  Calc_NMAEF_NMBF = TRUE


##################################################################################
##################################################################################
################################### DO NOT CHANGE FROM HERE ######################
##################################################################################
##  NOTHING BELOW HERE NEEDS DIRECT ALTERATION FOR ROUTINE RUNNING OF CODE
##      HOWEVER, MORE DEPENDENT VARIABLES ARE DEFINED HERE 
##################################################################################
##################################################################################

############################
## Make all directories
############################

  cat("\nStarting main routine. Creating directories as required.\n",sep="")

  dir.create(paste(FnOutputStoragePathStartInPPE,"mean_sd_constrained_and_unconstrained/Variables/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"mean_sd_constrained_and_unconstrained/Variables/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"predict_emulator_varsapplyconstraintto/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"implausibility_measure_outputdata/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"implausibility_measure_outputdata/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"apply_implausibility_tholdtol_RRinddata/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"apply_implausibility_tholdtol_RRinddata/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"samples_from_constraint/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"samples_from_constraint/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_implausibility_values/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_implausibility_values/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"NMAEF_NMBF_outputdata/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"NMAEF_NMBF_outputdata/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_appliedthreshtolcon/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_appliedthreshtolcon/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_NMAEF_NMBF/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_NMAEF_NMBF/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",VariableIn,"/varsapplyconstraintto/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",VariableIn,"/varsapplyconstraintto/validation_plots/",sep='')) ## Used?
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",VariableIn,"/varsapplyconstraintto/Validation_Tables/",sep='')) ## Used??
  dir.create(paste(FnOutputStoragePathStartInPPE,"fitted_emulator_models/",LS_text,"/",VariableIn,"/validation_plots/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"emulator_validation/",LS_text,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"emulator_validation/",LS_text,"/",VariableIn,"/",sep=''))
  dir.create(paste(FnOutputStoragePathStartInPPE,"emulator_validation/",LS_text,"/",VariableIn,"/validation_plots/",sep=''))

  for(kp in 1:length(MonthVec)){
    dir.create(paste(FnOutputStoragePathStartInPPE,"implausibility_measure_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"apply_implausibility_tholdtol_RRinddata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"samples_from_constraint/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"NMAEF_NMBF_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_implausibility_values/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_appliedthreshtolcon/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"plotting_NMAEF_NMBF/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
    dir.create(paste(FnOutputStoragePathStartInPPE,"emulator_validation/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",sep=''))
  }


##################################
## Remove .nc files if pre-existing
## and would be replaced by this file only
##################################

  if(RemNetCDF==TRUE){
    cat("\nRemoving all existing netCDF files related to this file.\n",sep="")
    for(kp in 1:length(MonthVec)){
      if(Calc_NMAEF_NMBF==TRUE){
        file.remove(paste(FnOutputStoragePathStartInPPE,"NMAEF_NMBF_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/NMAEF_NMBF_",GridResolutionIn,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,".nc",sep=''))
      }
      if(SaveEmPredToFile==TRUE){
        file.remove(paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",VariableIn,"/",EmPredLS_text,"_ConstraintObs_",GridResolutionIn,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,".nc",sep=''))
      }
      file.remove(paste(FnOutputStoragePathStartInPPE,"implausibility_measure_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/ImplausMeas_",GridResolutionIn,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,"_1Mmv_Ob",ObsMeasUnc_PlusMinusPercIn,"_Sp",SpatialCoLocUnc_PlusMinusPercIn,"_T",TemporalCoLocUnc_PlusMinusPercIn,"_St",StructVarTerm_PlusMinusPercIn,EmMultiSetIDIn,".nc",sep=''))
    }
  }


############################
## Emulation switches and paths
############################

  if(EmPredRequired==TRUE){
    if(FitEmulators==TRUE){
      PPE_TrainingInputs = read.csv(paste("/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/",PPEnameIn,"/design/Unit_",Tr_Val_text,"_TrainingDesign.csv",sep=''),header=TRUE,stringsAsFactors=FALSE)
      PPE_ValidationInputs = read.csv(paste("/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/",PPEnameIn,"/design/Unit_",Tr_Val_text,"_ValidationDesign.csv",sep=''),header=TRUE,stringsAsFactors=FALSE)
      InputNames_PPE = names(PPE_TrainingInputs)
      if(OnJASMIN==TRUE){
        PPE_DataFilePathStartIn = paste("/group_workspaces/jasmin2/gassp/myoshioka/um/UKCA_",PPEnameIn,"_PPE/",PPE_EmsYear,"/Variables/",VariableIn,"/",Resolution_text,"_",GridResolutionIn,"/",sep="")
      }else{
        PPE_DataFilePathStartIn = paste(Local_Obs_Path_Start,"/",PPEnameIn,"/observations_PPEmembers_datafiles/PPEoutput_",PPEnameIn,"_",GridResolutionIn,"/Variables/",VariableIn,"/",sep="")
      }
    }else{
      EmulatorList_FileNameStartIn = paste("Emulators_",GridResolutionIn,sep="")
    }
## If retaining the emulator predictons to file...
    if(SaveEmPredToFile==TRUE){
      OutputPred_FileNameStart = paste(EmPredLS_text,"_ConstraintObs_",GridResolutionIn,sep="")
      HistsFileNameStart = paste("HistsEmMeanPred_ConstraintObs_",GridResolutionIn,sep="")
    }
  }


#########################################################################
## Is emulation prediction (sample) required within all gridboxes or will it be read in?
## THIS IS THE COMPUTATIONALLY EXPENSIVE PART, SO TAKE CARE
## Large sample of parameter combinations to be read in.
## This is JJs 1 million member sample made using prior pdfs - consistency is priority.
## If LS_UseNetCDF == FALSE then data read in from .dat table.
##
## LS_switch defines the type of LargeSample to use 1=Prior / 2=Unif
## Additional switch for PPE type (26/27) embedded in LS_FilePath
#########################################################################

  if(EmPredRequired==TRUE){  ## define the file with the input combinations for the lrge sample want to use.
    LS_UseNetCDF = TRUE
    if(LS_switch==1){
      LS_FilePath = paste("/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/",PPEnameIn,"/regional_synthetic_constraint/largeinputsamples/SingleRLHS_size1M/RLHSamp_1Mcombs_MargDistsApp_Sample1.nc",sep='')
    }else{
      LS_FilePath = paste("/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/",PPEnameIn,"/regional_synthetic_constraint/largeinputsamples/SingleRLHS_size1M/RLHSamp_1Mcombs_MargUnif_Sample1.nc",sep='')
    }
  }else{  ## if already predicted, make the start of the filename where the mean predictions and their StDev's etc are stored.
    EmPredLS_FileNameStartIn = paste(EmPredLS_text,"_ConstraintObs_",GridResolutionIn,sep="")
  }


#############################################################################################
## Implausibility output filename for specified observation/variable/month
## ' EmMultiSetIDIn' for testing smaller datasets
## Current filename format to state large sample size (1Million) and then an optional test # (Test1) with the chosen Percentage uncertainty 
## specifications used for the 'ObsMeasUnc', 'SpatialCoLocUnc' and 'TemporalCoLocUnc' in that order
## and StructuralUnc on end (zero until formally specified) 
##
## 'SaveImplausTabToR' saves implausibility output to a table in the background R that you are running in.
##   -> Table is over-written if do multiple months in same run of code...
#############################################################################################

  SaveImplausTabToR = TRUE
#  SaveImplausTabToR = FALSE

  if(StructVarTerm_AsPercOnModelAveIn==TRUE){  
    Implaus_MultiSetIDStore = paste("_1Mmv_Ob",ObsMeasUnc_PlusMinusPercIn,"_Sp",SpatialCoLocUnc_PlusMinusPercIn,"_T",TemporalCoLocUnc_PlusMinusPercIn,"_St",StructVarTerm_PlusMinusPercIn,EmMultiSetIDIn,sep="")
  }else{
    Implaus_MultiSetIDStore = paste("_1Mmv_Ob",ObsMeasUnc_PlusMinusPercIn,"_Sp",SpatialCoLocUnc_PlusMinusPercIn,"_T",TemporalCoLocUnc_PlusMinusPercIn,"_St",StructVarTermValueIn,EmMultiSetIDIn,sep="")
  }

  ImplausOutput_FileNameStart = paste("ImplausMeas_",GridResolutionIn,sep="")
  ncdf_OutPrecision = "float" ## "float" saves to 7sf? use "double" if want ncdf to be same as when retain in R space. With "float", ncdf and R same to 6dp, but float used less memory.
  ncdf_MissingValueSpec = NA


#################################################################
## NMAEF and NMBF metric filenames for the specified observation/variable/month
## 'EmMultiSetIDIn' defines filename extension for specific test cases
## 'SaveMetricsToR' keeps information in background of R for ongoing use.
## If saving to file, not needed.
## SaveMetricsToR will continually overwrite variables to only have the most recent month
#################################################################

  if(Calc_NMAEF_NMBF==TRUE){
    MetricOut_MultiSetIDStore = EmMultiSetIDIn ## see comment above
    OutputMetrics_FileNameStart = paste("NMAEF_NMBF_",GridResolutionIn,sep="")
    SaveMetricsToR = TRUE
  }


################################################################################################
################################################################################################
##   MAIN CODE BEGINS HERE   
################################################################################################
################################################################################################

###########################################
## 1. If required, read in the large sample of input combinations
###########################################

  if(EmPredRequired==TRUE){
    cat("\nReading in large sample from ",LS_FilePath," \n",sep="")
    if(LS_UseNetCDF==TRUE){
      ncLIS = nc_open(LS_FilePath)
      LargeInputsSample = ncvar_get(ncLIS,"InputsSample") ##[1:2,]  ## LRE TEMPORARY!!! [1:2,] SHOULD BE REMOVED ASAP
    }else{
      LargeInputsSample = read.table(LS_FilePath,header = TRUE, stringsAsFactors=FALSE)
    }
    dim(LargeInputsSample)
    NoVariants = dim(LargeInputsSample)[1]
    NoPar = dim(LargeInputsSample)[2]
  }


##########################################################################
## Loop over months, for this variable:
## NOTE: All filenames have month included, but SaveMetricsToR == TRUE will store information in R for ongoing 
## use, but will rewrite variables using current month.
##########################################################################

  cat("\nVariable ",VariableIn,"\n",sep="")
  for(kp in 1:length(MonthVec)){
    cat("\nMonth ",MonthVec[kp]," \n",sep="")


################################################
## 2. Read in the Observations Info Table for the given Variable/month 
################################################

    ObsIDLocInfoTab_File = paste(Local_Obs_Path_Start,"obsID_gblocation_infotables/",VariableIn,Obs_text,"/",ObsIDLocTab_FileNameStartIn,"_ObsLocTab_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,".dat",sep="")

    ObsIDLocInfoTab = read.table(ObsIDLocInfoTab_File,header=TRUE,stringsAsFactors=FALSE)

    TotalObs = dim(ObsIDLocInfoTab)[1]
    NoObs = dim(ObsIDLocInfoTab)[1]
    rm(ObsIDLocInfoTab_File)


#####################################################################
## 3. Emulators, if need to do the predictions:
##    Either: Fit the emulators over the observation locations (takes a while), OR, 
##            Read in the emulator list for the observations of 'Variable' in 'Month'
##  -> Emulators can be pre-fitted by running "/nfs/a173/earjsj/UM-UKCA_vn8.4_PPEs/27aeratm/real_obs_constraint/RunCode_N48_EmulateValidateEmulate.r"
##            Then, set up file if want to save them:
## ELSE: read in the emulator predictions output that has already been generated and saved:
#####################################################################

    if(EmPredRequired==TRUE){ 
      EmMeanPredTab_OverObs = NULL
      EmSDPredTab_OverObs = NULL
      if(FitEmulators==TRUE){
        cat("\nFitting emulators.\n",sep="")
        EmulatorList = JJCode_N48_EmulateValidateEmulate_DSEm_SingleVariableMonth(Variable=VariableIn, Month=MonthVec[kp], ObsTabMultiSetID=ObsTab_MultiSetIDIn, EmMultiSetID=EmMultiSetIDIn, VLID=VLIDIn, SelectVL=SelectVLIn, VLno=VLnoIn, PPEname=PPEnameIn, PPE_TrainingInputs_Table=PPE_TrainingInputs, PPE_ValidationInputs_Table=PPE_ValidationInputs, PPE_DataFilePathStart=PPE_DataFilePathStartIn, PPE_DataFileNameMid=PPE_DataFileNameMidIn, PPE_DataFileNameEnd=PPE_DataFileNameEndIn, ObsPath= Local_Obs_Path_Start, StoragePath=FnOutputStoragePathStartInPPE, ObsIDLocTab_FileNameStart=ObsIDLocTab_FileNameStartIn, SourceLibs=FALSE, UseDSEm=UseDSEmIn, ReEmulateAllRuns=ReEmulateAllRunsIn, ReScaleOutput=ReScaleOutputIn, ReScaleDivisor=ReScaleDivisorIn, SaveEmulatorsToFile=SaveEmulatorsToFileIn, ReturnEmListToR=TRUE,LS_textIn=LS_text)
      }else{
        cat("\nReading in emulators.\n",sep="")
        EmModelsFile = paste(FnOutputStoragePathStartIn,"/",LS_text,"/fitted_emulator_models/",VariableIn,"/",EmulatorList_FileNameStartIn,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,EmMultiSetIDIn,".RData",sep="")
        load(EmModelsFile) ## loads the emulators as a list object called "EmulatorList"
        rm(EmModelsFile)
      }
      if(SaveEmPredToFile==TRUE){ ## setting up the NetCDF file to store the emulator predictions:
        ncStorageFileEmPred = paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",VariableIn,"/",OutputPred_FileNameStart,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,EmMultiSetIDIn,".nc",sep="")
        Predxdim = ncdim_def(name="variants",units="InputCombRowID",vals=1:NoVariants,unlim=TRUE)
        Predydim = ncdim_def(name="obs",units="ObsID",vals=ObsIDLocInfoTab[,1],unlim=TRUE)
        PredInfodim = ncdim_def(name="Infodim",units="",vals=1:20,unlim=FALSE)
        ObsEmMeanPred = ncvar_def(name="ObsEmMeanPred",units=VariableUnits,dim=list(Predxdim,Predydim),missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision,chunksizes=c(1024,NoObs))
        if(RetainEmSDPred==TRUE){
          ObsEmSDPred = ncvar_def(name="ObsEmSDPred",units=VariableUnits,dim=list(Predxdim,Predydim),missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision,chunksizes=c(1024,NoObs))
        }
        Min_EmPred = ncvar_def(name="Min_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        Max_EmPred = ncvar_def(name="Max_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        L95_EmPred = ncvar_def(name="L95_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        U95_EmPred = ncvar_def(name="U95_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        Median_EmPred = ncvar_def(name="Median_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        Mean_EmPred = ncvar_def(name="Mean_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        SD_EmPred = ncvar_def(name="SD_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        LQ_EmPred = ncvar_def(name="LQ_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        UQ_EmPred = ncvar_def(name="UQ_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        L90_EmPred = ncvar_def(name="L90_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        U90_EmPred = ncvar_def(name="U90_EmPred",units=VariableUnits,dim=Predydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
        VariableUsed = ncvar_def(name="VariableUsed",units="",dim=PredInfodim,prec="char")
        MonthUsed = ncvar_def(name="MonthUsed",units="",dim=PredInfodim,prec="char")
        VLIDUsed = ncvar_def(name="VLIDUsed",units="",dim=PredInfodim,prec="char")
        if(RetainEmSDPred==TRUE){ ## Might choose not to keep $sd, as is big?? Default is to store it...
          Pred_ncout = nc_create(filename=ncStorageFileEmPred,vars=list(ObsEmMeanPred,ObsEmSDPred,Min_EmPred,Max_EmPred,L95_EmPred,U95_EmPred,Median_EmPred,Mean_EmPred,SD_EmPred,LQ_EmPred,UQ_EmPred,L90_EmPred,U90_EmPred,VariableUsed,MonthUsed,VLIDUsed))
        }else{
          Pred_ncout = nc_create(filename=ncStorageFileEmPred,vars=list(ObsEmMeanPred,Min_EmPred,Max_EmPred,L95_EmPred,U95_EmPred,Median_EmPred,Mean_EmPred,SD_EmPred,LQ_EmPred,UQ_EmPred,L90_EmPred,U90_EmPred,VariableUsed,MonthUsed,VLIDUsed))
        }
        if(MakeEmPredMeanHistPlots==TRUE){ ## Initiate file for histograms...
          HistPlotFile = paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",VariableIn,"/",HistsFileNameStart,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,EmMultiSetIDIn,".pdf",sep="")
          pdf(HistPlotFile,width=14,height=14,paper="a4",onefile=TRUE)
          par(mfrow=c(3,2))
        }
      }
    }else{   ## If not doing emulator predicting in this run:
      cat("\nNo emulation prediction for this run.\n",sep="")
      EmulatorList = NULL
      LargeInputsSample = NULL
      EmMeanPredLSTab_OverObs_File = paste(FnOutputStoragePathStartInPPE,"predict_emulator_obsconstrain/",LS_text,"/",VariableIn,"/",EmPredLS_FileNameStartIn,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,EmMultiSetIDIn,".nc",sep="")
      ncEmPred = nc_open(EmMeanPredLSTab_OverObs_File) ## Read in the predictions for the large sample...
      EmMeanPredTab_OverObs = ncvar_get(ncEmPred,"ObsEmMeanPred")
      EmSDPredTab_OverObs = ncvar_get(ncEmPred,"ObsEmSDPred")
      NoVariants = dim(EmMeanPredTab_OverObs)[1]
    }


##############################################################
## 4. Read in the HiRes longitude/latitude GridBox ID conversion tables 
##    Read in the corresponding HiRes grid areas matrix
##    Masaru's interannual variability trends are all on HiRes - this converts.
##############################################################

    LonFilePath = paste(Local_Obs_Path_Start,"obsID_gblocation_infotables/LonLatConverterFiles/N96_LongitudeGB_Converter.dat",sep="")
    LatFilePath = paste(Local_Obs_Path_Start,"obsID_gblocation_infotables/LonLatConverterFiles/N96_LatitudeGB_Converter.dat",sep="")
    N96_LonIDTable = read.table(LonFilePath,header=TRUE,stringsAsFactors=FALSE)
    N96_LatIDTable = read.table(LatFilePath,header=TRUE,stringsAsFactors=FALSE)
    GridAreaFilepath = paste(FnOutputStoragePathStartIn,"obsID_gblocation_infotables/MO_grid_areas_N96.nc",sep="")
    GAncid = nc_open(GridAreaFilepath)
###    print(GAncid)
    N96_GridAreaMat = ncvar_get(GAncid,"land_area_fraction")
###    dim(GridAreaMat) ## MO area mat uses LonLat at midpoints... only has 144 cols for lat... MY files have 145... duplicate final lat col to match
    N96_GridAreaMat = cbind(N96_GridAreaMat,N96_GridAreaMat[,144])
###   dim(GridAreaMat)
    rm(LonFilePath,LatFilePath,GridAreaFilepath,GAncid)


##################################################################
## 5. Read in the "Inter-annual uncertainty Relative StDev Output Grid" from MY analysis
## on HiRes grid
##################################################################

    if(VariableIn=="AODs_TOTAL"){
      InterAnnUncFile = paste("/nfs/a68/earmy/PEGASOS_Hindcast/trend_stats/trend_stats_AOD440_TOTAL_xjum_1980-2009_",MonthVec[kp],".nc",sep="")
    }else{
      InterAnnUncFile = paste("/nfs/a68/earmy/PEGASOS_Hindcast/trend_stats/trend_stats_",VariableIn,"_sfc_xjum_1980-2009_",MonthVec[kp],".nc",sep="")
    }
##    InterAnnUncFile = paste("/nfs/a134/mm11lr/constraint/trend_stats/temporary_link_for_",VariableIn,"_",MonthVec[kp],".nc",sep="")  # TEMPORARY FOR nconc_corsol ONLY
    cat("\n.Interannual variability being read in from ",InterAnnUncFile,".\n",sep="")
    IAncid<-nc_open(InterAnnUncFile)
    InterAnn_RelSDData = ncvar_get(IAncid,"rel_stdev_res") 
##    dim(InterAnn_RelSDData)
    rm(InterAnnUncFile,IAncid)


###############################################################################################
## Set up the NetCDF file to write the NMAEF and NMBF output to.
###############################################################################################

    if(Calc_NMAEF_NMBF==TRUE){
      Metricydim = ncdim_def(name="obs",units="obsID",vals=ObsIDLocInfoTab[,1],unlim=TRUE)
      chardim = ncdim_def(name="nchar",units="",vals=1:20,create_dimvar=FALSE)
      NMAEF = ncvar_def(name="NMAEF",units="No units",dim=Metricydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
      NMBF = ncvar_def(name="NMBF",units="No units",dim=Metricydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
      S_ErrorDirection = ncvar_def(name="S_ErrorDirection",units="No units",dim=Metricydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
      M_bar = ncvar_def(name="M_bar",units="VariableUnits",dim=Metricydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
      ObsVal = ncvar_def(name="ObsVal",units="VariableUnits",dim=Metricydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
      VariableUse = ncvar_def(name="Variable",units="",dim=chardim,prec="char")
      MonthUse = ncvar_def(name="Month",units="",dim=chardim,prec="char")
      VLUse = ncvar_def(name="VLID",units="",dim=chardim,prec="char")
      ncOutputFile = paste(FnOutputStoragePathStartInPPE,"NMAEF_NMBF_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",OutputMetrics_FileNameStart,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,MetricOut_MultiSetIDStore,".nc",sep="")
      Metric_ncout = nc_create(filename=ncOutputFile,vars=list(NMAEF,NMBF,S_ErrorDirection,M_bar,ObsVal,VariableUse,MonthUse,VLUse))
      if(SaveMetricsToR==TRUE){
        NMAEF_Vec = vector("numeric",NoObs)
        NMBF_Vec = vector("numeric",NoObs)
        S_ErrorDirectionVec = vector("numeric",NoObs)
      }
    }


###############################################################################################
## Set up the NetCDF file to write the implausibility output to.
##  -> Also create a table to save the output to, if want to retain it in the R-Space as well.
##  -> **** Need to think how would save output if only use a subset of the observations?? ****
###############################################################################################

    xdim = ncdim_def(name="variants",units="InputCombRowID",vals=1:NoVariants,unlim=TRUE)
    ydim = ncdim_def(name="obs",units="obsID",vals=ObsIDLocInfoTab[,1],unlim=TRUE)
    Percdim = ncdim_def(name="percSet",units="",vals=c(1),unlim=FALSE)
    ImplausMeasure = ncvar_def(name="ImplausMeasure",units="No units",dim=list(xdim,ydim),missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision,chunksizes=c(1024,NoObs)) 
    MeanEmScaler_IAU_Str = ncvar_def(name="MeanEmScaler_IAU_Str",units="variable units",dim=ydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    RelSDToScaleIAU = ncvar_def(name="RelSDToScaleIAU",units="No units",dim=ydim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    ObsMeasUnc_PMPercOnObs = ncvar_def(name="ObsMeasUnc_PMPercOnObs",units="percentage",dim=Percdim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    SpatialCoLocUnc_PMPercOnObs = ncvar_def(name="SpatialCoLocUnc_PMPercOnObs",units="percentage",dim=Percdim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    TemporalCoLocUnc_PMPercOnObs = ncvar_def(name="TemporalCoLocUnc_PMPercOnObs",units="percentage",dim=Percdim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    if(StructVarTerm_AsPercOnModelAveIn==TRUE){
      StructUnc_PMPercOnModAve = ncvar_def(name="StructUnc_PMPercOnModAve",units="percentage",dim=Percdim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    }else{
      StructUnc_VarTermValueApplied = ncvar_def(name="StructUnc_VarTermValueApplied",units="variable units",dim=Percdim,missval=ncdf_MissingValueSpec,prec=ncdf_OutPrecision)
    }
    ncOutputFile = paste(FnOutputStoragePathStartInPPE,"implausibility_measure_outputdata/",LS_text,"/",VariableIn,"/",MonthVec[kp],"/",ImplausOutput_FileNameStart,"_",VariableIn,ObsTab_MultiSetIDIn,"_",MonthVec[kp],"_",VLIDIn,Implaus_MultiSetIDStore,".nc",sep="")
    if(StructVarTerm_AsPercOnModelAveIn==TRUE){
      Implaus_ncout = nc_create(filename=ncOutputFile,vars=list(ImplausMeasure,MeanEmScaler_IAU_Str,RelSDToScaleIAU,ObsMeasUnc_PMPercOnObs,SpatialCoLocUnc_PMPercOnObs,TemporalCoLocUnc_PMPercOnObs,StructUnc_PMPercOnModAve))
    }else{
      Implaus_ncout = nc_create(filename=ncOutputFile,vars=list(ImplausMeasure,MeanEmScaler_IAU_Str,RelSDToScaleIAU,ObsMeasUnc_PMPercOnObs,SpatialCoLocUnc_PMPercOnObs,TemporalCoLocUnc_PMPercOnObs,StructUnc_VarTermValueApplied))
    }
    if(SaveImplausTabToR==TRUE){
      ImplausMeasure_OverObsTable = NULL
    }


#########################################################################
## Loop over the observations and apply the implausibility measure calculation
## saving to the NetCDF as go along and retaining output to Rspace Table if want to
## If want to retain the EmMeanPred and EmSDPred, output it back from the Implausibility 
## calculation run then put that to the EmPred NetCDF
##
## Predict the output variable with the emulator at all input combinations. Retain the mean
##   predictions, the corresponding StDev on these mean predictions, and then take the mean 
##    value to get an "average emulator prediction" for use with the InterAnnUnc term
##
#########################################################################
#########################################################################
## LRE May2018
## Equivalent to calls of the functions 'JJCode_Calculate_ImplausibilityMeasure_SingleObs.r' and
##  'JJCode_Calculate_NMAEF_NMBF_SingleObs_OverInSamp'... These codes now embedded in here...
##
## THIS IS INTENDED TO MAKE IT SIMPLE TO ADAPT THE CODE HERE FOR OBSERVATION-SPECIFIC UNCERTAINTIES
## (rather than using single percentage values, as is currently the case)
## Code variable names are consistent b/w files, meaning there are duplicates here.
## Seperated code for setting emulated quantities, calculating implausibility and errors
## to avoid repeat prediction
#########################################################################

    for(jt in 1:NoObs){
      cat("\nOn Observation ",jt," out of ",NoObs,".\n",sep="")

      ObsValue = ObsIDLocInfoTab[jt,2]
      ObsLon = ObsIDLocInfoTab[jt,5]
      ObsLat = ObsIDLocInfoTab[jt,6]
      ObsIDtoUse = ObsIDLocInfoTab[jt,1] 


######################################################################
## CODES PREVIOUSLY IN FUNCTIONS IN JJ FILES GO HERE:
## -> from 'JJCode_CalculateImplausibility_SingleObs_OverInSamp':
######################################################################
######################################################################
## If needs to make predictions from the Emulator for the large sample 
##  of input combinations, do that here... 
##  -> Otherwise, set EmMeanPred and EmSDPred from file read in...
######################################################################

      if(EmPredRequired==TRUE){

        cat("\n.Making sample predictions from emulator.\n",sep="")
        nInSamp = dim(LargeInputsSample)[1]


############################################
## Select the emulator from the emulator list:
##  -> position select is dependent on if "AllObsInTab"...
############################################

        if(FitEmulators==TRUE | AllObsInTab==TRUE){
          EmulatorModel = EmulatorList[[jt]]
        }else{
          EmulatorModel = EmulatorList[[ObsIDtoUse]]
        }


############################################
## Predict the output at the large sample input combinations using the emulator:
############################################      

        if(nInSamp>10000){
          PredOut = JJCode_PredictFromEm_UsingLargeSample(EmModIn=EmulatorModel,LargeSampInputCombs=LargeInputsSample,nPredBreakVal=10000,PredMean=TRUE,PredSD=TRUE,Pred95CIs=FALSE)
        }else{
          PredOut = predict.km(EmulatorModel,LargeInputsSample,"UK",se.compute=TRUE,light.return=TRUE,checkNames=FALSE)
        }
        EmMeanPred = PredOut$mean
        EmSDPred = PredOut$sd
        rm(PredOut)
      }else{

############################################
## Otherwise, if not doing prediction, use the read in emulator predictions for the large sample:
##  -> Column select from read in tables again dependent on if "AllObsInTab"
############################################

        if(AllObsInTab==TRUE){
          EmMeanPred = EmMeanPredTab_OverObs[,jt]
          EmSDPred = EmSDPredTab_OverObs[,jt]
        }else{
          EmMeanPred = EmMeanPredTab_OverObs[,ObsIDtoUse]
          EmSDPred = EmSDPredTab_OverObs[,ObsIDtoUse]
        }
      }
      MeanEmOut_IAU_Str_scaling = mean(EmMeanPred)


##################################################################
## Make the denominator variance error term that comes from the 'Emulator Uncertainty'
##  -> Different for each input combinations in the InputsSample
##  -> Taken directly from DiceKriging prediction$sd
##################################################################

      EmVarVec = (EmSDPred)^2


####################################################################################
## Make the denominator variance error term that comes from the 'Observational Measurement Uncertainty'
##  -> Take as a percentage error on the observation value
##  -> Obs measured to within +/- %age of their value
##  -> Assume Gaussian distribution and equate (+/- %age/100)*obs to 95% Gaussian Confidence Interval, of ~ +/-2*StDev 
##  -> Is currently the same for all input combinations (model variants)
####################################################################################

      ObsMeasSD = (ObsValue*(ObsMeasUnc_PlusMinusPercIn/100))/2
      ObsMeasVar = (ObsMeasSD)^2

##########################################################################
## Make the denominator variance error term that comes from the 'Spatial Co-location Uncertainty'
##  -> Treat in same way as the 'Observational Measurement Uncertainty' term above
##  -> Take as a percentage error on the observation value
##  -> Again, is currently the same for all input combinations (model variants)
##########################################################################

      SpatialCoLocSD = (ObsValue*(SpatialCoLocUnc_PlusMinusPercIn/100))/2
      SpatialCoLocVar = (SpatialCoLocSD)^2


###########################################################################
## Make the denominator variance error term that comes from the 'Temporal Co-location Uncertainty'
##  -> Treat in same way as the 'Observational Measurement Uncertainty' term above
##  Same idea as above. Currently single value for all model variants
###########################################################################

      TempCoLocSD = (ObsValue*(TemporalCoLocUnc_PlusMinusPercIn/100))/2
      TempCoLocVar = (TempCoLocSD)^2


#################################################################################
## Make the denominator variance error term that comes from the 'Inter-annual variability Uncertainty'
##  -> Use the Relative StDev (StDev/mean value) from Masaru's analysis of Steve Turnock's 40 year runs
##  and scale back to variable scale using the mean of the full sample of emulator predictions that cover the full space
##  -> Method to extract the Relative StDev from MY's array is dependent on resolution used
##     -> The array is N96 (HiRes): code allows for making an average with surrounding neigbour GBs (for N96) 
##         or for extracting an N48 value from the N96 matrix
#################################################################################

      if(GridResolutionIn=="N48"){
        InterAnn_RelSDVal = JJCode_ExtractN48GB_FromN96Mat(N48Lon=ObsLon,N48Lat=ObsLat,N96LonIDTab=N96_LonIDTable,N96LatIDTab=N96_LatIDTable,N96_OutputMat=InterAnn_RelSDData,N96_AreaMat=N96_GridAreaMat)
      }else{
        if(UseSurrNeighbourAveIn==TRUE){
          N96_LonID = N96_LonIDTable[N96_LonIDTable[,2]==ObsLon,1]
          N96_LatID = N96_LatIDTable[N96_LatIDTable[,2]==ObsLat,1]
          InterAnn_RelSDVal = JJCode_MakeAve_GBplusSurrNeighbours(GBLonID=N96_LonID,GBLatID=N96_LatID,ModelOutputGBMat=InterAnn_RelSDData,GBAreaMat=N96_GridAreaMat)      
        }else{
          N96_LonID = N96_LonIDTable[N96_LonIDTable[,2]==ObsLon,1]
          N96_LatID = N96_LatIDTable[N96_LatIDTable[,2]==ObsLat,1]
          InterAnn_RelSDVal = InterAnn_RelSDData[N96_LonID,N96_LatID]
        }
      }
      InterAnnSD = MeanEmOut_IAU_Str_scaling*(InterAnn_RelSDVal)
      InterAnnVar = (InterAnnSD)^2


##########################################################################
## Make the denominator variance error term that comes from the 'Structural Uncertainty/Model Discrepancy'
##  -> Allow this to be a value that is just entered, or calculated as a percentage error on the 
##      mean from the emulator predictions over the full sample (using 'MeanEmOut_IAU_Str_scaling')
##########################################################################

      if(StructVarTerm_AsPercOnModelAveIn==TRUE){
        StructSD = (MeanEmOut_IAU_Str_scaling*(StructVarTerm_PlusMinusPercIn/100))/2   
        StructVar = (StructSD)^2 
      }else{
        StructVar = StructVarTermValueIn
      }


############################################
## Calculate the Implausibility measure over the InputsSample:
############################################

      ImplausTopVec = abs(ObsValue-EmMeanPred)
      SumVarTermsVec = EmVarVec+ObsMeasVar+SpatialCoLocVar+TempCoLocVar+InterAnnVar+StructVar
      ImplausBottomVec = sqrt(SumVarTermsVec)
      ImplausibilityVec = ImplausTopVec/ImplausBottomVec


######################################################################
## CODES PREVIOUSLY IN FUNCTIONS IN JJ FILES GO HERE:
## -> from 'JJCode_Calculate_NMAEF_NMBF_SingleObs.r':
######################################################################
###############################################################
## Calculate the 'normalised mean absolute error factor' (NMAEF) and the 'normalised mean bias 
##  factor' (NMBF) over the InputsSample:
##  -> From the formula in the paper by "S. Yu et al, 2006", the observed values O_i are the 
##      same value in all comparisons (same observation each time!) and this is compared to all 
##       the different model output variants over the uncertainty, M_i. (1 milion of them)
###############################################################

###############################################################
## 1. Calculate N, O_bar and M_bar; the mean of the observations and model predictions:
## -> O_bar is the mean of the observations, which is just equal to the observed value as 
##       the same value is used for the observation in all comparisons over the uncertainty
###############################################################

      N = length(EmMeanPred)
      O_Vec = rep(ObsValue,times=N)
      O_bar = ObsValue
      M_bar_val = mean(EmMeanPred)


#################################################
## 2. Calculate the "S" indicator; S = (M_bar - O_bar)/|M_bar - O_bar|
#################################################

      Sval = (M_bar_val-O_bar)/(abs(M_bar_val-O_bar)) 


#################################################
## 3. Calculate the 'normalised mean absolute error factor': E_NMAEF
##   -> (sum(|M_i - O_i|))/(((sum(O_i))^((1+S)/2))*((sum(M_i))^((1-S)/2)))
#################################################

      TopEq = sum(abs(EmMeanPred - O_Vec))
      BottomEq = ((sum(O_Vec))^((1+Sval)/2))*((sum(EmMeanPred))^((1-Sval)/2))
      E_NMAEF = TopEq/BottomEq


################################################
## 4. Calculate the 'normalised mean bias factor': B_NMBF
##   -> S*[exp(|ln((sum(M_i))/(sum(O_i)))|)-1]
################################################

      LogPart = log(((sum(EmMeanPred))/(sum(O_Vec)))) 
      B_NMBF = Sval*((exp(abs(LogPart)))-1)


###############################################################
###############################################################
## Store the Implausibility output to NetCDF:
###############################################################
###############################################################

      ncvar_put(nc=Implaus_ncout,varid=ImplausMeasure,vals=ImplausibilityVec,start=c(1,jt),count=c(-1,1))
      ncvar_put(nc=Implaus_ncout,varid=MeanEmScaler_IAU_Str,vals=MeanEmOut_IAU_Str_scaling,start=jt,count=1)
      ncvar_put(nc=Implaus_ncout,varid=RelSDToScaleIAU,vals=InterAnn_RelSDVal,start=jt,count=1)
      if(SaveImplausTabToR==TRUE){
        ImplausMeasure_OverObsTable = cbind(ImplausMeasure_OverObsTable,ImplausibilityVec)
      }

###############################################################
###############################################################
## If retaining the emulator predictions, write these to their NetCDF file:
###############################################################
###############################################################

      if(SaveEmPredToFile==TRUE){
        ncvar_put(nc=Pred_ncout,varid=ObsEmMeanPred,vals=EmMeanPred,start=c(1,jt),count=c(-1,1))
        if(RetainEmSDPred==TRUE){
          ncvar_put(nc=Pred_ncout,varid=ObsEmSDPred,vals=EmSDPred,start=c(1,jt),count=c(-1,1))
        }
        ncvar_put(nc=Pred_ncout,varid=Min_EmPred,vals=min(EmMeanPred),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=Max_EmPred,vals=max(EmMeanPred),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=L95_EmPred,vals=quantile(EmMeanPred,probs=c(0.025)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=U95_EmPred,vals=quantile(EmMeanPred,probs=c(0.975)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=Median_EmPred,vals=quantile(EmMeanPred,probs=c(0.5)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=Mean_EmPred,vals=mean(EmMeanPred),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=SD_EmPred,vals=sd(EmMeanPred),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=LQ_EmPred,vals=quantile(EmMeanPred,probs=c(0.25)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=UQ_EmPred,vals=quantile(EmMeanPred,probs=c(0.75)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=L90_EmPred,vals=quantile(EmMeanPred,probs=c(0.05)),start=jt,count=1)
        ncvar_put(nc=Pred_ncout,varid=U90_EmPred,vals=quantile(EmMeanPred,probs=c(0.95)),start=jt,count=1)
        if(MakeEmPredMeanHistPlots==TRUE){
          PlotTitle = paste(VariableIn,": Observation ",ObsIDtoUse,"\n Lon = ",ObsLon,"; Lat = ",ObsLat,sep="")
          hist(EmMeanPred,freq=FALSE,main=PlotTitle,breaks=20,xlab=paste(VariableIn," (",VariableUnits,")",sep=""),cex.axis=1.1,cex.lab=1.2,cex.main=1.4)
          abline(v=ObsValue,col=2,lwd=2)
          abline(v=mean(EmMeanPred),col=4,lwd=2)
        }
      }


###############################################################
###############################################################
## If generating the NMAEF and NMBF metrics, write these to their NetCDF file here:
###############################################################
###############################################################

      if(Calc_NMAEF_NMBF==TRUE){
        ncvar_put(nc=Metric_ncout,varid=NMAEF,vals=E_NMAEF,start=jt,count=1)
        ncvar_put(nc=Metric_ncout,varid=NMBF,vals=B_NMBF,start=jt,count=1)
        ncvar_put(nc=Metric_ncout,varid=S_ErrorDirection,vals=Sval,start=jt,count=1)
        ncvar_put(nc=Metric_ncout,varid=M_bar,vals=M_bar_val,start=jt,count=1)
        ncvar_put(nc=Metric_ncout,varid=ObsVal,vals=ObsValue,start=jt,count=1)
        if(SaveMetricsToR==TRUE){
          NMAEF_Vec[jt] = E_NMAEF
          NMBF_Vec[jt] = B_NMBF
          S_ErrorDirectionVec[jt] = Sval
        }
      }


###############################################
## End the jt loop over the observations:
###############################################

    } 


###############################################
## Add in extra "meta-data" outputs to the Implausibility and Metrics NetCDF files:
###############################################

    ncvar_put(nc=Implaus_ncout,varid=ObsMeasUnc_PMPercOnObs,vals=ObsMeasUnc_PlusMinusPercIn)
    ncvar_put(nc=Implaus_ncout,varid=SpatialCoLocUnc_PMPercOnObs,vals=SpatialCoLocUnc_PlusMinusPercIn)
    ncvar_put(nc=Implaus_ncout,varid=TemporalCoLocUnc_PMPercOnObs,vals=TemporalCoLocUnc_PlusMinusPercIn)
    if(StructVarTerm_AsPercOnModelAveIn==TRUE){
      ncvar_put(nc=Implaus_ncout,varid=StructUnc_PMPercOnModAve,vals=StructVarTerm_PlusMinusPercIn)
    }else{
      ncvar_put(nc=Implaus_ncout,varid=StructUnc_VarTermValueApplied,vals=StructVarTermValueIn)
    }
    if(Calc_NMAEF_NMBF==TRUE){
      ncvar_put(nc=Metric_ncout,varid=VariableUse,vals=VariableIn)
      ncvar_put(nc=Metric_ncout,varid=MonthUse,vals=MonthVec[kp])
      ncvar_put(nc=Metric_ncout,varid=VLUse,vals=VLIDIn)
    }


##############################################################
## Close the NetCDF output files we are writing to: Implausibility and prediction (if using)
## Close the plotting file for the histograms of the predicted values (if plotting them)
## Tidy up the Output table if retained to R Space:
##############################################################

    nc_close(Implaus_ncout)
    if(SaveEmPredToFile==TRUE){
      nc_close(Pred_ncout)
      if(MakeEmPredMeanHistPlots==TRUE){
        dev.off()
      }
    }
    if(Calc_NMAEF_NMBF==TRUE){
      nc_close(Metric_ncout)
    }
    if(SaveImplausTabToR==TRUE){
      ColNameVec = paste("Obs",ObsIDLocInfoTab[,1],sep="")
      colnames(ImplausMeasure_OverObsTable) = ColNameVec
      rownames(ImplausMeasure_OverObsTable) = NULL
      rm(ColNameVec)
    }


###############################################
##  End the kp loop over the months:
###############################################

  }  


###############################################################################################
###############################################################################################
## END OF CODE
###############################################################################################
###############################################################################################
