subroutine micro_p3_addfld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! register output field for P3 microphysics
!!
!! Author: Kai Zhang 
!!
!! Last updated: 2017-11-28
!!
!! Current version: 0.52
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use shr_kind_mod,   only: r8=>shr_kind_r8
   use phys_control,   only: phys_getopts 
   use constituents,   only: cnst_add,      &
                             cnst_get_ind,  &
                             cnst_name,     &
                             cnst_longname, &
                             sflxnam,       &
                             apcnst,        &
                             bpcnst,        &
                             pcnst
   use micro_p3_acme,  only: ncnst,      &
                             cnst_names, & 
                             ixcldliq, &
                             ixcldice, & 
                             ixrain,   &
                             ixcldrim, &
                             ixbvrim,  &
                             ixnumliq, &
                             ixnumice, &
                             ixnumrain
   use cam_history,    only: addfld,      &
                             horiz_only,  &
                             add_default, &
                             outfld
   use cam_logfile,    only: iulog
   use cam_abortutils, only: endrun
   use error_messages, only: handle_errmsg
                          
   implicit none

   logical :: &
      history_amwg,  &
      history_budget

   character(32) ::  &
      tnvar,         & 
      tunit,         &
      tflag 
      
   character(128) :: & 
      tname  
   
   integer ::            &
      budget_histfile,   & ! output history file number for budget fields
      m,                 & 
      mm,                & 
      ierr

   !!.....................................................................................
   !! surface flux
   !!.....................................................................................
   
   do m = 1, ncnst
      call cnst_get_ind(cnst_names(m), mm)
      if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixcldrim /)) ) then
         ! mass mixing ratios
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', &
           cnst_longname(mm) )
         call addfld(sflxnam(mm), horiz_only, 'A', 'kg/m2/s', &
           trim(cnst_name(mm))//' surface flux')
      else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain /)) ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', &
           cnst_longname(mm) )
         call addfld(sflxnam(mm), horiz_only, 'A', '1/m2/s', &
           trim(cnst_name(mm))//' surface flux')
      else if ( mm == ixbvrim ) then
         ! number concentrations
         call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'm3/kg', &
           cnst_longname(mm) )
         call addfld(sflxnam(mm), horiz_only, 'A', 'm3/m2/s', &
           trim(cnst_name(mm))//' surface flux')
      else
         call endrun( "micro_p3_acme_init: &
              &Could not call addfld for constituent with unknown units.")
      endif
   end do
   
   !!.....................................................................................
   !!
   !! mass before and after physics 
   !!
   !! apcnst(m) = trim(cnst_name(m))//'AP'
   !! bpcnst(m) = trim(cnst_name(m))//'BP'
   !!.....................................................................................

   call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixcldliq))//' after physics'  )
   call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixcldice))//' after physics'  )
   call addfld(apcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixrain))//' after physics'  )
   call addfld(apcnst(ixcldrim), (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixcldrim))//' after physics'  )

   call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixcldliq))//' before physics' )
   call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixcldice))//' before physics' )
   call addfld(bpcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', &
           trim(cnst_name(ixrain))//' before physics' )
   call addfld(bpcnst(ixcldrim),  (/ 'lev' /), 'A', 'm3/kg', &
           trim(cnst_name(ixcldrim))//' before physics' )



   !!.....................................................................................
   !! time scales 
   !!.....................................................................................
   
   tnvar = 'P3_PRAIN'
   tname = 'Total precipitation (rain+snow) production'
   tunit = 'kg/kg'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NEVAPR'
   tname = 'Evaporation of total precipitation (rain + snow)'
   tunit = 'kg/kg'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   !!.....................................................................................
   !! time scales 
   !!.....................................................................................
      
   tnvar = 'P3_EPSC'
   tname = 'supersaturation relaxation time scales for droplets '
   tunit = 's'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
      
   tnvar = 'P3_EPSI'
   tname = 'supersaturation relaxation time scales for ice '
   tunit = 's'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
      
   tnvar = 'P3_EPSR'
   tname = 'supersaturation relaxation time scales for rain '
   tunit = 's'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
      
   tnvar = 'P3_REPS'
   tname = '1/(EPSC+EPSI+EPSR)'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
 
   tnvar = 'P3_AC1'
   tname = 'Ac term in P3'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
 
   tnvar = 'P3_AC2'
   tname = 'Ac term in P3'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
 
   tnvar = 'P3_AC3'
   tname = 'Ac term in P3'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCCON1'
   tname = 'QCCON1'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCCON2'
   tname = 'QCCON2'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QICON1'
   tname = 'QICON1'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QICON2'
   tname = 'QICON2'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QICON3'
   tname = 'QICON3'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   !!.....................................................................................
   !! liquid microphysics process rates 
   !!.....................................................................................
   
   tnvar = 'P3_QCCON'
   tname = 'qctend due to condensation on droplets'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
  
   tnvar = 'P3_QRCON'
   tname = 'qrtend due to condensation on rain'
   tunit = 'kg/kg/s'
   tflag = 'A'

   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
 
   tnvar = 'P3_QCAUT'
   tname = 'qctend due to autoconversion to rain'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCAUTC'
   tname = 'nctend due to autoconversion to rain (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCAUTR'
   tname = 'nrtend due to autoconversion of cloud water'
   tunit = '#/kg/s'
   tflag = 'A'

   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCACC'
   tname = 'qctend due to accretion by rain'
   tunit = 'kg/kg/s'
   tflag = 'A'

   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCACC'
   tname = 'nctend due to accretion by rain (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCSLF'
   tname = 'nctend due to cloud droplet self-collection'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NRSLF'
   tname = 'nrtend due to rain self-collection (negative) '
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCNUC'
   tname = 'qctend due to activation of CCN'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCNUC'
   tname = 'nctend due to activation of CCN (positive or negative, due to diffusion)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QREVP'
   tname = 'qrtend due to rain evaporation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NREVP'
   tname = 'nrtend due to rain evaporation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCEVP'
   tname = 'qctend due to cloud droplet evaporation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QBERG'
   tname = 'qctend/qitend due to Bergeron process'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 


   !!.....................................................................................
   !! ice microphysics process rates
   !!.....................................................................................
   
   tnvar = 'P3_QCCOL'
   tname = 'qc/qitend due to collection cloud water'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_QIDEP'
   tname = 'qv/qitend due to vapor deposition'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_QRCOL'
   tname = 'qi/qrtend due to collection rain mass by ice'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_QINUMC'
   tname = 'qitend due to deposition/condensation freezing nuc'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_NCCOL'
   tname = 'change in cloud droplet number from collection by ice (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_NRCOL'
   tname = 'change in rain number from collection by ice (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )

   tnvar = 'P3_QINUC'
   tname = 'qitend due to cirrus ice nucleation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NINUC'
   tname = 'change in ice number from deposition/cond-freezing nucleation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_QISUB'
   tname = 'sublimation of ice'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_QIMLT'
   tname = 'qitend due to melting of ice'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_NIMLT'
   tname = 'nitend due to melting of ice (negative)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_NISUB'
   tname = 'nitend due to sublimation (negative)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )
   
   tnvar = 'P3_NISLF'
   tname = 'nitend due to collection within a category (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) )

   tnvar = 'P3_QCHETI'
   tname = 'qctend/qitend due to immersion freezing of droplets (positive)'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCHETI'
   tname = 'nctend/nitend due to immersion freezing of droplets (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCHETC'
   tname = 'qctend/qitend due to contact freezing of droplets (if not zero, positive)'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NCHETC'
   tname = 'nctend/nitend due to contact freezing of droplets'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QRHETI'
   tname = 'qctend/qitend due to immersion freezing of rain (positive)'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NRHETI'
   tname = 'nrtend/nitend due to immersion freezing of rain (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QRHETC'
   tname = 'qrtend/qitend due to contact freezing of rain'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NRHETC'
   tname = 'nrtend/nitend due to contact freezing of rain'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QCSHD'
   tname = 'qctend due to cloud water/ice collision above freezing and shedding or wet growth and shedding'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NCSHDC'
   tname = 'nrtend due to cloud water/ice collision above freezing and shedding (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_NRSHDR'
   tname = 'nrtend due to collision of rain/ice above freezing and shedding (positive)'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_QRMUL'
   tname = 'qrtend/qitend due to ice multiplication from rime-splitnering of rain'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NIMUL'
   tname = 'nitend due to ice multiplication from rime-splintering'
   tunit = '#/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 


   !!.....................................................................................
   !! sedimentation 
   !!.....................................................................................

   tnvar = 'P3_QCSED'
   tname = 'qctend due to sedimentation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QISED'
   tname = 'qitend due to sedimentation' 
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_QRSED'
   tname = 'qrtend due to sedimentation'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 


   !!.....................................................................................
   !! diag  
   !!.....................................................................................

 
   tnvar = 'P3_CME'
   tname = 'Rate of cond-evap within the cloud'
   tunit = 'kg/kg/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_CMEIOUT'
   tname = 'qitend due to deposition/sublimation of cloud ice'
   tunit = 'kg/kg/s'
   tflag = 'A'
   !! diag_cmei  (i,k,iice) = (qidep(iice) - qisub(iice) + qinuc(iice)) * icldm(i,k)   
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_ICWMRST'
   tname = 'Prognostic in-stratus qc'
   tunit = 'kg/kg'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_ICIMRST'
   tname = 'Prognostic in-stratus qi'
   tunit = 'kg/kg'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_UMC'
   tname = 'Mass-weighted qc fallspeed'
   tunit = 'm/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_UMR'
   tname = 'Mass-weighted qr fallspeed'
   tunit = 'm/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_UMI'
   tname = 'Mass-weighted qi fallspeed'
   tunit = 'm/s'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_ICWNC'
   tname = 'Prognostic in-cloud nc'
   tunit = 'm-3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_ICINC'
   tname = 'Prognostic in-cloud ni'  
   tunit = 'm-3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_AWNC'
   tname = 'Average cloud ice number conc'  
   tunit = 'm-3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_AWNI'
   tname = 'Prognostic in-cloud ni'  
   tunit = 'm-3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_AREL'
   tname = 'Average droplet effective radius'  
   tunit = 'um'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_AREI'
   tname = 'Average ice effective radius'  
   tunit = 'um'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_FREQL'
   tname = 'Fractional occurrence of liquid'  
   tunit = 'm-3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_FREQI'
   tname = 'Fractional occurrence of ice'  
   tunit = '1'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_FREQR'
   tname = 'Fractional occurrence of rain'  
   tunit = '1'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_REL'
   tname = 'stratiform cloud effective radius liquid'  
   tunit = 'um'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_REI'
   tname = 'stratiform cloud effective radius ice'  
   tunit = 'um'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NCAL'
   tname = 'Number Concentation Activated for Liquid'  
   tunit = '1/m3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   
   tnvar = 'P3_NCAI'
   tname = 'Number Concentation Activated for Ice'  
   tunit = '1/m3'
   tflag = 'A'
   
   call addfld (trim(tnvar), (/'lev'/), trim(tflag), trim(tunit), trim(tname) ) 
   

   !!
   !! vertically-integrated fields   
   !!

   tnvar = 'P3_MPICLWPI'
   tname = 'in-cloud LWP before P3'
   tunit = 'kg/m2'
   tflag = 'A'
   
   call addfld (trim(tnvar), horiz_only, trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_MPICIWPI'
   tname = 'in-cloud IWP before P3'
   tunit = 'kg/m2'
   tflag = 'A'
   
   call addfld (trim(tnvar), horiz_only, trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_CDNUMC'
   tname = 'Vertically-integrated droplet concentration'
   tunit = '1/m2'
   tflag = 'A'
   
   call addfld (trim(tnvar), horiz_only, trim(tflag), trim(tunit), trim(tname) ) 

   tnvar = 'P3_ADRAIN'
   tname = 'Average rain effective Diameter'
   tunit = 'um'
   tflag = 'A'
   
   call addfld (trim(tnvar), horiz_only, trim(tflag), trim(tunit), trim(tname) ) 

        
!!!   call addfld ('P3_FICE', (/ 'lev' /), 'A', 'fraction', 'Fractional ice content within cloud'                     )
!!!
!!!   ! History variables for CAM5 microphysics
!!!   call addfld ('P3_MPDT', (/ 'lev' /), 'A', 'W/kg', 'Heating tendency - Morrison microphysics'                )
!!!   call addfld ('P3_MPDQ', (/ 'lev' /), 'A', 'kg/kg/s', 'Q tendency - Morrison microphysics'                      )
!!!   call addfld ('P3_MPDLIQ', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDLIQ tendency - Morrison microphysics'                 )
!!!   call addfld ('P3_MPDICE', (/ 'lev' /), 'A', 'kg/kg/s', 'CLDICE tendency - Morrison microphysics'                 )
!!!   call addfld ('P3_MPDW2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Vapor tendency - Morrison microphysics'       )
!!!   call addfld ('P3_MPDW2I', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Ice tendency - Morrison microphysics'         )
!!!   call addfld ('P3_MPDW2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Water <--> Precip tendency - Morrison microphysics'      )
!!!   call addfld ('P3_MPDI2V', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Vapor tendency - Morrison microphysics'         )
!!!   call addfld ('P3_MPDI2W', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Water tendency - Morrison microphysics'         )
!!!   call addfld ('P3_MPDI2P', (/ 'lev' /), 'A', 'kg/kg/s', 'Ice <--> Precip tendency - Morrison microphysics'
!!!   call addfld ('P3_EFFLIQ_IND', (/ 'lev' /), 'A','Micron', 'Prognostic droplet effective radius (indirect effect)'   )

!!!   ! This is provided as an example on how to write out subcolumn output
!!!   ! NOTE -- only 'I' should be used for sub-column fields as subc-columns could shift from time-step to time-step
!!!   if (use_subcol_microp) then
!!!      call addfld('P3_FICE_SCOL', (/'psubcols','lev     '/), 'I', 'fraction', &
!!!           'Sub-column fractional ice content within cloud', flag_xyfill=.true., fill_value=1.e30_r8)
!!!   end if
!!!
   
   
   
   
!!!   ! Average cloud top particle size and number (liq, ice) and frequency
!!!   call addfld ('P3_ACTREL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet effective radius'              )
!!!   call addfld ('P3_ACTREI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice effective radius'                  )
!!!   call addfld ('P3_ACTNL', horiz_only,    'A', 'Micron', 'Average Cloud Top droplet number'                        )
!!!   call addfld ('P3_ACTNI', horiz_only,    'A', 'Micron', 'Average Cloud Top ice number'                            )
!!!
!!!   call addfld ('P3_FCTL', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top liquid'                )
!!!   call addfld ('P3_FCTI', horiz_only,    'A', 'fraction', 'Fractional occurrence of cloud top ice'                   )
!!!
!!!   call addfld ('P3_LS_FLXPRC', (/ 'ilev' /), 'A', 'kg/m2/s', 'ls stratiform gbm interface rain+snow flux')
!!!
!!!   call addfld ('P3_LS_REFFRAIN', (/ 'lev' /), 'A', 'micron', 'ls stratiform rain effective radius')
!!!   call addfld ('P3_CV_REFFLIQ', (/ 'lev' /), 'A', 'micron', 'convective cloud liq effective radius')
!!!   call addfld ('P3_CV_REFFICE', (/ 'lev' /), 'A', 'micron', 'convective cloud ice effective radius')
!!!
!!!   ! diagnostic precip
!!!
!!!   ! diagnostic radar reflectivity, cloud-averaged
!!!   call addfld ('P3_REFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity'       )
!!!   call addfld ('P3_AREFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity'       )
!!!   call addfld ('P3_FREFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity'       )
!!!
!!!   call addfld ('P3_CSRFL',(/ 'lev' /), 'A','DBz','94 GHz radar reflectivity (CloudSat thresholds)'       )
!!!   call addfld ('P3_ACSRFL',(/ 'lev' /), 'A','DBz','Average 94 GHz radar reflectivity (CloudSat thresholds)'       )
!!!   call addfld ('P3_FCSRFL',(/ 'lev' /), 'A','fraction','Fractional occurrence of radar reflectivity (CloudSat thresholds)' &
!!!        )
!!!
!!!   call addfld ('P3_AREFLZ',(/ 'lev' /), 'A','mm^6/m^3','Average 94 GHz radar reflectivity'       )



!!!   ! precipitation efficiency & other diagnostic fields
!!!   call addfld('P3_PE'    ,     horiz_only, 'A', '1', 'Stratiform Precipitation Efficiency  (precip/cmeliq)' )
!!!   call addfld('P3_APRL'  ,     horiz_only, 'A', 'm/s', 'Average Stratiform Precip Rate over efficiency calculation' )
!!!   call addfld('P3_PEFRAC',     horiz_only, 'A', '1', 'Fraction of timesteps precip efficiency reported' )
!!!   call addfld('P3_QCRAT', (/ 'lev' /), 'A', 'fraction', 'Qc Limiter: Fraction of qc tendency applied')


   !!.....................................................................................
   !! determine the add_default fields
   !!.....................................................................................
   
   call phys_getopts(history_amwg_out           = history_amwg         , &
                     history_budget_out         = history_budget       , &
                     history_budget_histfile_num_out = budget_histfile)


   !!.....................................................................................
   !! amwg diag 
   !!.....................................................................................
   
   if (history_amwg) then
!!!      call add_default ('P3_FICE    ', 1, ' ')
      call add_default ('P3_ADRAIN   ', 1, ' ')
      call add_default ('P3_AREI     ', 1, ' ')
      call add_default ('P3_AREL     ', 1, ' ')
      call add_default ('P3_AWNC     ', 1, ' ')
      call add_default ('P3_AWNI     ', 1, ' ')
      call add_default ('P3_CDNUMC   ', 1, ' ')
      call add_default ('P3_FREQR    ', 1, ' ')
      call add_default ('P3_FREQL    ', 1, ' ')
      call add_default ('P3_FREQI    ', 1, ' ')
      
      do m = 1, ncnst
         call cnst_get_ind(cnst_names(m), mm)
         call add_default(cnst_name(mm), 1, ' ')
         ! call add_default(sflxnam(mm),   1, ' ')
      end do
   end if

   !!.....................................................................................
   !! cloud water/ice budget
   !!.....................................................................................

      if ( history_budget ) then
    
      !!
      !! QC budget 
      !!
!! 
!!    qctend = + qcnuc  &nucleation 
!!             - qcaut  &autoconversion 
!!             - qcacc  &accretion 
!!             + qccon  &condensation 
!!             - qcevp  &evaporation 
!!             - qchetc &contact 
!!             - qcheti &immersion
!!             - qccol  &collection 
!!             - qcshd  &shedding 
!!             - qcsed  &sedimentation 

      call add_default ('P3_QCCON ', budget_histfile, ' ')
      call add_default ('P3_QCNUC ', budget_histfile, ' ')
      call add_default ('P3_QCAUT ', budget_histfile, ' ') !! qrtend + 
      call add_default ('P3_QCACC ', budget_histfile, ' ') !! qrtend + 
      call add_default ('P3_QCEVP ', budget_histfile, ' ')
      call add_default ('P3_QBERG ', budget_histfile, ' ')
      call add_default ('P3_QCHETC', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QCHETI', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QCCOL ', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QCSHD ', budget_histfile, ' ') !! qrtend + 
      call add_default ('P3_QCSED ', budget_histfile, ' ')
          
      !!
      !! QR budget 
      !!
!!    qrtend = + qcaut  &accretion
!!             + qcacc  &autoconversion
!!             + qcshd  &shedding 
!!             + qrcon  &condensation 
!!             - qrevp  &evaporation 
!!             - qrcol  &collection 
!!             + qimlt  &melting of ice/snow 
!!             - qrhetc &contact 
!!             - qrheti &immersion
!!             - qrmul  &rime-splitnering of rain 
!!             - qrsed  &sedimentation 

      call add_default ('P3_QRCON ', budget_histfile, ' ')
!     call add_default ('P3_QCAUT ', budget_histfile, ' ') !! qrtend + 
!     call add_default ('P3_QCACC ', budget_histfile, ' ') !! qrtend +  
!     call add_default ('P3_QCSHD ', budget_histfile, ' ') !! qrtend + 
      call add_default ('P3_QREVP ', budget_histfile, ' ')
      call add_default ('P3_QRCOL ', budget_histfile, ' ') !! qitend + 
!     call add_default ('P3_QIMLT ', budget_histfile, ' ') !! qitend - 
      call add_default ('P3_QRHETC', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QRHETI', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QRMUL ', budget_histfile, ' ') !! qitend + 
      call add_default ('P3_QRSED ', budget_histfile, ' ')
      
      !!
      !! NC budget 
      !!

      call add_default ('P3_NCNUC ', budget_histfile, ' ')
      call add_default ('P3_NCACC ', budget_histfile, ' ')
      call add_default ('P3_NCAUTC', budget_histfile, ' ')
      call add_default ('P3_NCAUTR', budget_histfile, ' ')
      call add_default ('P3_NCHETI', budget_histfile, ' ')
      call add_default ('P3_NCHETC', budget_histfile, ' ')
      call add_default ('P3_NCSLF ', budget_histfile, ' ')
      call add_default ('P3_NCSHDC', budget_histfile, ' ')
      call add_default ('P3_NCCOL ', budget_histfile, ' ')
      
      !!
      !! NR budget 
      !!
      call add_default ('P3_NRSLF ', budget_histfile, ' ')
      call add_default ('P3_NRHETI', budget_histfile, ' ')
      call add_default ('P3_NRHETC', budget_histfile, ' ')
      call add_default ('P3_NRSHDR', budget_histfile, ' ')
      call add_default ('P3_NRCOL ', budget_histfile, ' ')
      
      !!
      !! QI budget 
      !!
!!    qitend = - qisub  &sublimation 
!!             - qimlt  &melting
!!             + qidep  &deposition 
!!             + qinuc  &nucleation 
!!             + qrcol  &collection qr
!!             + qccol  &collection qc
!!             + qrhetc &contact qr
!!             + qrheti &immersion qr
!!             + qchetc &contact qc
!!             + qcheti &immersion qc
!!             + qrmul  &rime-splitnering of rain 

      call add_default ('P3_QISUB ', budget_histfile, ' ')
      call add_default ('P3_QIMLT ', budget_histfile, ' ') ! qrtend - 
      call add_default ('P3_QIDEP ', budget_histfile, ' ')
      call add_default ('P3_QINUC ', budget_histfile, ' ')
      call add_default ('P3_QINUMC', budget_histfile, ' ')
      call add_default ('P3_QISED ', budget_histfile, ' ')

      call add_default ('P3_NINUC ', budget_histfile, ' ')
      call add_default ('P3_NIMLT ', budget_histfile, ' ')
      call add_default ('P3_NISLF ', budget_histfile, ' ')
      call add_default ('P3_NIMUL ', budget_histfile, ' ')
      
      
      !!
      !! others 
      !!
      
      call add_default ('P3_CMEIOUT',  budget_histfile, ' ')
      call add_default ('P3_ICWMRST',  budget_histfile, ' ')
      call add_default ('P3_ICIMRST',  budget_histfile, ' ')
      call add_default ('P3_MPICLWPI', budget_histfile, ' ')
      call add_default ('P3_MPICIWPI', budget_histfile, ' ')
      
!!!      call add_default ('P3_MPDW2V   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDW2P   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDW2I   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDT     ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDQ     ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDLIQ   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDICE   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDI2W   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDI2V   ', budget_histfile, ' ')
!!!      call add_default ('P3_MPDI2P   ', budget_histfile, ' ')

      call add_default(cnst_name(ixcldliq), budget_histfile, ' ')
      call add_default(cnst_name(ixcldice), budget_histfile, ' ')
      call add_default(cnst_name(ixrain),   budget_histfile, ' ')
      call add_default(cnst_name(ixcldrim), budget_histfile, ' ')

      call add_default(apcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(apcnst   (ixcldice), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldliq), budget_histfile, ' ')
      call add_default(bpcnst   (ixcldice), budget_histfile, ' ')
      call add_default(apcnst   (ixrain),   budget_histfile, ' ')
      call add_default(apcnst   (ixcldrim), budget_histfile, ' ')
      call add_default(bpcnst   (ixrain),   budget_histfile, ' ')
      call add_default(bpcnst   (ixcldrim), budget_histfile, ' ')

   end if

   return 
   
end subroutine micro_p3_addfld
