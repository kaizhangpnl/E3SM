#!/bin/csh
date

set echo verbose


set fetch_code    = 0   # 0 = No, >0 = Yes
set compile_model = 1   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG E3SM_20190418
setenv CCSMROOT /compyfs/${USER}/model/${CCSMTAG}

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv COMPSET FC5AV1C-04Z1
setenv RESOLUTION ne30_ne30
setenv MRES ne30
setenv MACH      compy
setenv PTMP      /compyfs/${USER}/bld

setenv ntasks 1600
setenv nthrds 1

setenv MYSRC     ${CCSMROOT}/mods_v1_p3_cmdv
setenv MYCLM     ${CCSMROOT}/mods_clm

setenv CASE     ${MACH}_${COMPSET}_${MRES}_${CCSMTAG}_v54a_clean_1600p
setenv COMCASE  ${MACH}_${COMPSET}_${MRES}_${CCSMTAG}_v54a_clean_1600p

setenv CASEROOT  ${CCSMROOT}/cases/$CASE
setenv RUNDIR    /compyfs/${USER}/csmruns/$CASE

####################################################################
# Compile model
####################################################################
if ($compile_model > 0) then

   rm -rf $CASEROOT
   cd  $CCSMROOT/cime/scripts

   ./create_newcase --case $CASEROOT --project e3sm --mach $MACH \
                    --res $RESOLUTION --compset $COMPSET

#====================================================================
# set up case
#====================================================================

   cd $CASEROOT

   ./xmlchange -file env_run.xml   -id RUNDIR  -val $RUNDIR

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ROF -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ROF -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_WAV -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_WAV -val '0'

   ./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $ntasks
   ./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $nthrds
   ./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val '0'

   ./case.setup

#====================================================================
# my mods of source code
#====================================================================
cd $CASEROOT

   ln -s ${MYSRC}/* SourceMods/src.cam    # put your mods in here
   ln -s ${MYCLM}/* SourceMods/src.clm    # put your mods in here

   ./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append -val ' -cosp'

   # Build the model

   cd $CASEROOT

   ./case.build

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then

#------------------
## set environment
#------------------

cd $CASEROOT

./xmlchange  -file env_run.xml  -id  RUN_STARTDATE   -val '0000-01-01'
./xmlchange  -file env_run.xml  -id  RESUBMIT        -val '0'
./xmlchange  -file env_run.xml  -id  STOP_N          -val '5'
./xmlchange  -file env_run.xml  -id  STOP_OPTION     -val 'ndays'
./xmlchange  -file env_run.xml  -id  REST_N          -val '5'
./xmlchange  -file env_run.xml  -id  REST_OPTION     -val 'ndays'
./xmlchange  -file env_run.xml  -id  DOUT_S          -val 'FALSE'

cat << EOF
cld_macmic_num_steps = 6
nucleate_ice_subgrid = 1.2
micro_p3_l_mg2_qidep = .true.
micro_p3_l_satadj = .false.
micro_p3_l_crconevp = .false.
micro_p3_l_massclip  = .true.
micro_p3_l_limit_qidep_qinuc = .true.
micro_p3_l_limit_qisub_qrevp = .true.
micro_p3_num_steps = 6
micro_p3_scale_berg = 1.0
micro_p3_scale_qidep = 1.0
micro_p3_l_cshd = .true.
micro_p3_l_imlt = .true.
micro_p3_l_ccol = .true.
micro_p3_opt_cheti = 1
micro_p3_opt_inuc  = 1
use_hetfrz_classnuc = .true.
EOF

cp ${CCSMROOT}/work/p3_lookup_table_1.dat ${RUNDIR}
cp ${CCSMROOT}/work/p3_lookup_table_2.dat ${RUNDIR}

# goto the case directory, make changes, and submit the job 
./$CASE.submit

endif

