module micro_p3_acme

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! E3SM interface for P3 microphysics
!!
!! Author: Kai Zhang 
!!
!! Features: 
!!
!!   1. support two kinds of subgrid treatment 
!!   2. can be run along with MG2 without impact on met fields 
!!   3. all output variables with 'P3_' in name 
!! 
!! Current version: 0.54ac
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use shr_kind_mod,   only: r8=>shr_kind_r8
use spmd_utils,     only: masterproc
use ppgrid,         only: pcols, &
                          pver,  &
                          pverp, &
                          psubcols
use physconst,      only: gravit, &
                          rair,   &
                          tmelt,  &
                          cpair,  &
                          rh2o,   &
                          rhoh2o, &
                          latvap, &
                          latice, &
                          mwh2o,  &
                          mwdry
use phys_control,   only: phys_getopts, &
                          use_hetfrz_classnuc
use physics_types,  only: physics_state, &
                          physics_ptend, &
                          physics_ptend_init, &
                          physics_state_copy, &
                          physics_state_dealloc, &
                          physics_ptend_sum, &
                          physics_ptend_scale
use physics_update_mod, only: physics_update
use physics_buffer, only: physics_buffer_desc, &
                          pbuf_add_field, &
                          dyn_time_lvls, &
                          pbuf_old_tim_idx, &
                          pbuf_get_index, &
                          dtype_r8, &
                          dtype_i4, &
                          pbuf_get_field, &
                          pbuf_set_field, &
                          col_type_subcol, &
                          pbuf_register_subcol
use constituents,   only: cnst_add, &
                          cnst_get_ind, &
                          qmin, &
                          cnst_name, &
                          cnst_longname, &
                          sflxnam, &
                          apcnst, &
                          bpcnst, &
                          pcnst
use cldfrc2m,       only: rhmini=>rhmini_const
use cam_history,    only: addfld, &
                          horiz_only, &
                          add_default, &
                          outfld
use cam_logfile,    only: iulog
use cam_abortutils, only: endrun
use error_messages, only: handle_errmsg
use ref_pres,       only: top_lev=>trop_cloud_top_lev
use subcol_utils,   only: subcol_get_scheme
use perf_mod,       only: t_startf, &
                          t_stopf
use cam_debug,      only: l_debug, l_summary_debug

implicit none
private
save


public :: &
   micro_p3_acme_readnl,    & 
   micro_p3_acme_register,  & 
   micro_p3_acme_init_cnst, & 
   micro_p3_acme_implements_cnst, & 
   micro_p3_acme_init,      & 
   micro_p3_acme_tend

!!........................................................................................
!! parameters 
!!........................................................................................

   
logical, parameter :: &
   microp_uniform = .False. 

character(len=16) :: &
   micro_p3_precip_frac_method = 'max_overlap' ! type of precipitation fraction method

real(r8) :: &
   ice_sed_ai = 700.0_r8      ! Fall speed parameter for cloud ice

logical, public :: &
   do_cldliq = .True.,   &
   do_cldice = .True.,   &
   do_nccons = .False.,  &
   do_nicons = .False.

logical :: & 
   l_mg2_qidep = .False. , & 
   l_massclip  = .True. , & 
   l_satadj    = .false. , & 
   l_crconevp  = .True. , & 
   l_limit_qidep_qinuc = .True. , &
   l_limit_qisub_qrevp = .True. , &
   l_cshd      = .True. , &
   l_imlt      = .True. , &
   l_ccol      = .True. 

integer, public :: & 
   opt_inuc   = 1,    &
   opt_cheti  = 1   

integer :: num_steps = 1 ! Number of P3 substeps

real(r8) :: scale_berg  = 1.0_r8 
real(r8) :: scale_qidep = 1.0_r8 

logical, parameter :: &
   l_pnog_nc = .True. 

integer, parameter, public :: &
   ncnst = 8 ! Number of constituents

!! CLDRIM : ice mass from rime growth
!! BVRIM  : bulk rime volume

character(len=8), public :: &      ! Constituent names
!  cnst_names(8) = (/'P3CLDLIQ','P3CLDICE','P3NUMLIQ','P3NUMICE', &
!                    'P3RAINQM','P3CLDRIM','P3NUMRAI','P3BVRIM'/)
   cnst_names(8) = (/'CLDLIQ','CLDICE','NUMLIQ','NUMICE', &
                     'RAINQM','CLDRIM','NUMRAI','BVRIM'/)

!! tracer index 
integer, public :: &
   ixcldliq = -1,      &
   ixcldice = -1,      &
   ixnumliq = -1,      &
   ixnumice = -1,      &
   ixrain   = -1,      &
   ixnumrain= -1,      &
   ixcldrim = -1,      &
   ixbvrim  = -1 

!! pbuf 
integer :: &
   cldo_idx,           & !! needed by microp_aero 
   qme_idx,            & !! needed by wetdep 
   prain_idx,          & !! needed by wetdep 
   nevapr_idx,         & !! needed by wetdep 
   rate1_cw2pr_st_idx, & !! needed by aero_model 
   dei_idx,            & !! needed by cloud_rad_props 
   mu_idx,             & !! needed by cloud_rad_props  
   lambdac_idx,        & !! needed by cloud_rad_props 
   rei_idx,            & !! needed by cosp  
   rel_idx,            & !! needed by cosp  
   ls_flxprc_idx,      & !! needed by cosp 
   ls_flxsnw_idx,      & !! needed by cosp 
   ls_reffrain_idx,    & !! needed by cosp 
   ls_reffsnow_idx,    & !! needed by cosp 
   cv_reffliq_idx,     & !! needed by cosp 
   cv_reffice_idx,     & !! needed by cosp 
   prer_evap_idx,      & !! needed by clubb 
   cmeliq_idx,         & !! needed by clubb 
   relvar_idx,         & !! input from clubb 
   accre_enhan_idx,    & !! input from clubb 
   iciwpst_idx,        & !! internal  
   iclwpst_idx,        & !! internal 
   cc_t_idx,           & !! needed by macrop_driver 
   cc_qv_idx,          & !! needed by macrop_driver 
   cc_ql_idx,          & !! needed by macrop_driver 
   cc_qi_idx,          & !! needed by macrop_driver 
   cc_nl_idx,          & !! needed by macrop_driver 
   cc_ni_idx,          & !! needed by macrop_driver 
   cc_qlst_idx


! Index fields for precipitation efficiency.
integer :: &
     acpr_idx = -1, &
     acgcme_idx = -1, &
     acnum_idx = -1

! Physics buffer indices for fields registered by other modules
integer :: &
   ast_idx = -1,            &
   cld_idx = -1,            &
   concld_idx = -1

! Pbuf fields needed for subcol_SILHS
integer :: &
     qrain_idx=-1, &
     nrain_idx=-1

integer :: &
   naai_idx = -1,           &
   naai_hom_idx = -1,       &
   npccn_idx = -1,          &
   rndst_idx = -1,          &
   nacon_idx = -1,          &
   prec_str_idx = -1,       &
   prec_pcw_idx = -1,       &
   prec_sed_idx = -1,       &
   snow_str_idx = -1,       &
   snow_pcw_idx = -1,       &
   snow_sed_idx = -1

! pbuf fields for heterogeneous freezing
integer :: &
   frzimm_idx = -1, &
   frzcnt_idx = -1, &
   frzdep_idx = -1

logical :: &
   allow_sed_supersat  ! allow supersaturated conditions after sedimentation loop



real(r8) :: &
   micro_mg_accre_enhan_fac = huge(1.0_r8), & !Accretion enhancement factor from namelist
   prc_coef1_in             = huge(1.0_r8), &
   prc_exp_in               = huge(1.0_r8), &
   prc_exp1_in              = huge(1.0_r8), &
   cld_sed_in               = huge(1.0_r8), & !scale fac for cloud sedimentation velocity
   nccons                   = huge(1.0_r8), &
   nicons                   = huge(1.0_r8)

interface p
   module procedure p1
   module procedure p2
end interface p



contains



subroutine micro_p3_acme_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Namelist variables
  integer :: micro_p3_num_steps = 1      ! Number of substepping iterations 
  logical :: micro_p3_l_mg2_qidep
  logical :: micro_p3_l_satadj
  logical :: micro_p3_l_massclip
  logical :: micro_p3_l_crconevp
  logical :: micro_p3_l_limit_qidep_qinuc
  logical :: micro_p3_l_limit_qisub_qrevp
  real(r8) :: micro_p3_scale_berg
  real(r8) :: micro_p3_scale_qidep
  logical :: micro_p3_l_cshd
  logical :: micro_p3_l_imlt
  logical :: micro_p3_l_ccol
  integer :: micro_p3_opt_cheti
  integer :: micro_p3_opt_inuc

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_p3_acme_readnl'

  namelist /micro_p3_nl/ micro_p3_num_steps, &
                         micro_p3_scale_berg, &
                         micro_p3_scale_qidep, &
                         micro_p3_l_mg2_qidep, &
                         micro_p3_l_satadj, &
                         micro_p3_l_massclip, &
                         micro_p3_l_crconevp, &
                         micro_p3_l_limit_qidep_qinuc, &
                         micro_p3_l_limit_qisub_qrevp, & 
                         micro_p3_l_cshd, &
                         micro_p3_l_imlt, &
                         micro_p3_l_ccol, &
                         micro_p3_opt_cheti, &
                         micro_p3_opt_inuc

  !-----------------------------------------------------------------------------

  if (masterproc) then
  
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_p3_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_p3_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     ! set local variables
     num_steps   = micro_p3_num_steps 
     scale_berg  = micro_p3_scale_berg
     scale_qidep  = micro_p3_scale_qidep
     l_mg2_qidep  = micro_p3_l_mg2_qidep 
     l_satadj     = micro_p3_l_satadj
     l_massclip  = micro_p3_l_massclip
     l_crconevp  = micro_p3_l_crconevp
     l_limit_qidep_qinuc = micro_p3_l_limit_qidep_qinuc
     l_limit_qisub_qrevp = micro_p3_l_limit_qisub_qrevp
     l_cshd      = micro_p3_l_cshd
     l_imlt      = micro_p3_l_imlt 
     l_ccol      = micro_p3_l_ccol
     opt_cheti     = micro_p3_opt_cheti
     opt_inuc      = micro_p3_opt_inuc
     
  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(num_steps,               1, mpiint, 0, mpicom)
  call mpibcast(scale_berg,              1, mpir8,  0, mpicom)
  call mpibcast(scale_qidep,             1, mpir8,  0, mpicom)
  call mpibcast(l_mg2_qidep,             1, mpilog, 0, mpicom)
  call mpibcast(l_satadj,                1, mpilog, 0, mpicom)
  call mpibcast(l_massclip,              1, mpilog, 0, mpicom)
  call mpibcast(l_crconevp,              1, mpilog, 0, mpicom)
  call mpibcast(l_limit_qidep_qinuc,     1, mpilog, 0, mpicom)
  call mpibcast(l_limit_qisub_qrevp,     1, mpilog, 0, mpicom)
  call mpibcast(l_cshd,                  1, mpilog, 0, mpicom)
  call mpibcast(l_imlt,                  1, mpilog, 0, mpicom)
  call mpibcast(l_ccol,                  1, mpilog, 0, mpicom)
  call mpibcast(opt_cheti,            1, mpiint, 0, mpicom)
  call mpibcast(opt_inuc,             1, mpiint, 0, mpicom)

#endif

end subroutine micro_p3_acme_readnl

!================================================================================================

subroutine micro_p3_acme_register

  ! Register microphysics constituents and fields in the physics buffer.
  !-----------------------------------------------------------------------

  logical :: prog_modal_aero
  logical :: use_subcol_microp  ! If true, then are using subcolumns in microphysics
  logical :: save_subcol_microp ! If true, then need to store sub-columnized fields in pbuf

  if(l_summary_debug) write(6,*) 'micro_p3_acme_register - 001 -' 
  
  call phys_getopts(use_subcol_microp_out = use_subcol_microp, &
                    prog_modal_aero_out   = prog_modal_aero, &
                    micro_mg_accre_enhan_fac_out = micro_mg_accre_enhan_fac)

  ! Register microphysics constituents and save indices.

  call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
                longname='Grid box averaged cloud liquid amount', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
                longname='Grid box averaged cloud ice amount', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(3), mwh2o, cpair, 0._r8, ixnumliq, &
                longname='Grid box averaged cloud liquid number', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(4), mwh2o, cpair, 0._r8, ixnumice, &
                longname='Grid box averaged cloud ice number', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(5), mwh2o, cpair, 0._r8, ixrain, &
                longname='Grid box averaged rain amount', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(6), mwh2o, cpair, 0._r8, ixcldrim, &
                longname='Grid box averaged riming amount', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(7), mwh2o, cpair, 0._r8, ixnumrain, &
                longname='Grid box averaged rain number', &
                is_convtran1=.true.)
  call cnst_add(cnst_names(8), mwh2o, cpair, 0._r8, ixbvrim, &
                longname='Grid box averaged riming volume', &
                is_convtran1=.true.)

  if(l_summary_debug) write(6,*) 'micro_p3_acme_register - 002 -' 
  
  !!  
  !! OUTPUT   
  !! 
  !!    RATE1_CW2PR_ST: needed by aero_model
  !!    CLDO     : needed by microp_aero 
  !!    PRER_EVAP: needed by clubb 
  !!    DEI      : needed by radiation 
  !!    MU       : needed by radiation 
  !!    LAMBDAC  : needed by radiation 
  !!    REL      : needed by cosp 
  !!    REI      : needed by cosp 
  !!    LS_FLXPRC   : needed by cosp 
  !!    LS_FLXSNW   : needed by cosp 
  !!    LS_REFFRAIN : needed by cosp 
  !!    LS_REFFSNOW : needed by cosp 
  !!    CV_REFFLIQ  : needed by cosp 
  !!    CV_REFFICE  : needed by cosp 
  !!    CC_xx       : needed by macro 
  !! 
  !! LOCAL
  !!
  !!    ICLWPST
  !!    ICIWPST
  !! 

  !! module microp_aero
  call pbuf_add_field('CLDO','global', dtype_r8,(/pcols,pver,dyn_time_lvls/), cldo_idx)
  
  !! module wetdep 
  call pbuf_add_field('QME',  'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
  call pbuf_add_field('PRAIN','physpkg',dtype_r8,(/pcols,pver/), prain_idx)
  call pbuf_add_field('NEVAPR','physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)

  !! module aero_model
  if (prog_modal_aero) then
     call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/), rate1_cw2pr_st_idx)
  endif
  
  !! module clubb_intr
  call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), prer_evap_idx)
  
  !! module radiation_data & module cloud_rad_props
  call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
  call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
  call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)
  
  !! module cospsimulator_intr
  call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)
  call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
  call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
  call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)
  call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
  call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
  call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
  call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)
  
  !! module macrop_driver
  call pbuf_add_field('CC_T',     'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_t_idx)
  call pbuf_add_field('CC_qv',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qv_idx)
  call pbuf_add_field('CC_ql',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ql_idx)
  call pbuf_add_field('CC_qi',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qi_idx)
  call pbuf_add_field('CC_nl',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_nl_idx)
  call pbuf_add_field('CC_ni',    'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_ni_idx)
  call pbuf_add_field('CC_qlst',  'global',  dtype_r8, (/pcols,pver,dyn_time_lvls/), cc_qlst_idx)

  !! (internal) Stratiform only in cloud liquid/ice water path for radiation
  call pbuf_add_field('ICLWPST',    'physpkg',dtype_r8,(/pcols,pver/), iclwpst_idx)
  call pbuf_add_field('ICIWPST',    'physpkg',dtype_r8,(/pcols,pver/), iciwpst_idx)


  if(l_summary_debug) write(6,*) 'micro_p3_acme_register - 003 -' 
  
  ! Register subcolumn pbuf fields
  if (use_subcol_microp) then
  
    call pbuf_register_subcol('CLDO',        'micro_p3_acme_register', cldo_idx)
    call pbuf_register_subcol('QME',         'micro_mg_cam_register', qme_idx)
    call pbuf_register_subcol('PRAIN',       'micro_mg_cam_register', prain_idx)
    call pbuf_register_subcol('NEVAPR',      'micro_mg_cam_register', nevapr_idx)
    call pbuf_register_subcol('PRER_EVAP',   'micro_p3_acme_register', prer_evap_idx)
    call pbuf_register_subcol('DEI',         'micro_p3_acme_register', dei_idx)
    call pbuf_register_subcol('MU',          'micro_p3_acme_register', mu_idx)
    call pbuf_register_subcol('LAMBDAC',     'micro_p3_acme_register', lambdac_idx)
    call pbuf_register_subcol('REL',         'micro_p3_acme_register', rel_idx)
    call pbuf_register_subcol('REI',         'micro_p3_acme_register', rei_idx)
    call pbuf_register_subcol('LS_FLXPRC',   'micro_p3_acme_register', ls_flxprc_idx)
    call pbuf_register_subcol('LS_FLXSNW',   'micro_p3_acme_register', ls_flxsnw_idx)
    call pbuf_register_subcol('LS_REFFRAIN', 'micro_p3_acme_register', ls_reffrain_idx)
    call pbuf_register_subcol('LS_REFFSNOW', 'micro_p3_acme_register', ls_reffsnow_idx)
    call pbuf_register_subcol('CV_REFFLIQ',  'micro_p3_acme_register', cv_reffliq_idx)
    call pbuf_register_subcol('CV_REFFICE',  'micro_p3_acme_register', cv_reffice_idx)
    call pbuf_register_subcol('CC_T',        'micro_p3_acme_register', cc_t_idx)
    call pbuf_register_subcol('CC_qv',       'micro_p3_acme_register', cc_qv_idx)
    call pbuf_register_subcol('CC_ql',       'micro_p3_acme_register', cc_ql_idx)
    call pbuf_register_subcol('CC_qi',       'micro_p3_acme_register', cc_qi_idx)
    call pbuf_register_subcol('CC_nl',       'micro_p3_acme_register', cc_nl_idx)
    call pbuf_register_subcol('CC_ni',       'micro_p3_acme_register', cc_ni_idx)
    call pbuf_register_subcol('CC_qlst',     'micro_p3_acme_register', cc_qlst_idx)
    call pbuf_register_subcol('ICIWPST',     'micro_p3_acme_register', iciwpst_idx)
    call pbuf_register_subcol('ICLWPST',     'micro_p3_acme_register', iclwpst_idx)
    
    if (prog_modal_aero) then
      call pbuf_register_subcol('RATE1_CW2PR_ST', 'micro_p3_acme_register', rate1_cw2pr_st_idx)
    end if
    
  end if

  !! (internal) Precipitation efficiency fields across timesteps.
  call pbuf_add_field('ACPRECL',    'global',dtype_r8,(/pcols/), acpr_idx)   ! accumulated precip
  call pbuf_add_field('ACGCME',     'global',dtype_r8,(/pcols/), acgcme_idx) ! accumulated condensation
  call pbuf_add_field('ACNUM',      'global',dtype_i4,(/pcols/), acnum_idx)  ! counter for accumulated # timesteps

  !! module clubb_intr
  call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/), relvar_idx)
  call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

  ! Diagnostic fields needed for subcol_SILHS, need to be grid-only
  if (subcol_get_scheme() == 'SILHS') then
     call pbuf_add_field('QRAIN',   'global',dtype_r8,(/pcols,pver/), qrain_idx)
     call pbuf_add_field('NRAIN',   'global',dtype_r8,(/pcols,pver/), nrain_idx)
  end if

  if(l_summary_debug) write(6,*) 'micro_p3_acme_register - 004 -' 

end subroutine micro_p3_acme_register


!!
!!........................................................................................
!!
!!


function micro_p3_acme_implements_cnst(name)

   ! Return true if specified constituent is implemented by the
   ! microphysics package

   character(len=*), intent(in) :: name        ! constituent name
   logical :: micro_p3_acme_implements_cnst    ! return value

   micro_p3_acme_implements_cnst = any(name == cnst_names)

end function micro_p3_acme_implements_cnst


!!
!!........................................................................................
!!


subroutine micro_p3_acme_init_cnst(name, q, gcid)

   ! Initialize the microphysics constituents, if they are
   ! not read from the initial file.

   character(len=*), intent(in) :: name     ! constituent name
   real(r8), intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
   integer,  intent(in)  :: gcid(:)  ! global column id

   if (micro_p3_acme_implements_cnst(name)) q = 0.0_r8

end subroutine micro_p3_acme_init_cnst


!!
!!........................................................................................
!!


subroutine micro_p3_acme_init(pbuf2d)

   use time_manager,   only: is_first_step, is_first_restart_step 
   use micro_p3_utils, only: micro_p3_utils_init
   use micro_p3,       only: micro_p3_init, &
                             micro_p3_lookuptable_init

   !!
   !! Initialization for P3 microphysics
   !!

   type(physics_buffer_desc), pointer :: &
      pbuf2d(:,:)

   integer :: &
      m,      &
      mm
      
   logical :: &
      history_amwg,      & ! output the variables used by the AMWG diag package
      history_budget,    & ! for cloud water budgets.
      use_subcol_microp, &
      do_clubb_sgs 
   
   integer ::            &
      budget_histfile,   & ! output history file number for budget fields
      ierr
      
   character(128) ::     &
      errstring,         & ! return status (non-blank for error return)
      lookup_file_1,     & 
      lookup_file_2
   

   character(32) ::  &
      tnvar,         & 
      tunit,         &
      tflag 
      
   character(128) :: & 
      tname 

!!=== 

   if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 001 -' 


   call phys_getopts(use_subcol_microp_out= use_subcol_microp, &
                     do_clubb_sgs_out     = do_clubb_sgs,      &
                     prc_coef1_out        = prc_coef1_in,      &
                     prc_exp_out          = prc_exp_in,        &
                     prc_exp1_out         = prc_exp1_in,       &
                     cld_sed_out          = cld_sed_in)

   if (do_clubb_sgs) then
     allow_sed_supersat = .false.
   else
     allow_sed_supersat = .true.
   endif

   call micro_p3_utils_init(r8, rh2o, cpair, tmelt, latvap, latice, &
                            ice_sed_ai, errstring) 

   call handle_errmsg(errstring, subname="micro_p3_utils_init")

if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 002 -' 


   call micro_p3_init( &
        r8, gravit, rair, rh2o, cpair, &
        tmelt, latvap, latice, rhmini, &
        microp_uniform, do_cldice, use_hetfrz_classnuc, &
        do_nccons, do_nicons, nccons, nicons, &
        allow_sed_supersat, ice_sed_ai, prc_coef1_in,prc_exp_in, &
        prc_exp1_in, cld_sed_in, errstring)

   call handle_errmsg(errstring, subname="micro_p3_init")


if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 003 -' 


   !! register output 
   
   call micro_p3_addfld()


if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 004 -' 

  !!  
  !! INPUT   
  !! 
  !!    AST      : from clubb 
  !!    CLD      : from clubb 
  !!    CONCLD   : from clubb 
  !!    CMELIQ   : from clubb 
  !! 
  !! LOCAL
  !!
  !! 
  
   !! 
   !! cloud macro 
   !!
   
   ast_idx      = pbuf_get_index('AST')    !! from CLUBB 
   cld_idx      = pbuf_get_index('CLD')    !! from CLUBB 
   concld_idx   = pbuf_get_index('CONCLD') !! from CLUBB 
   cmeliq_idx   = pbuf_get_index('CMELIQ') !! from CLUBB Rate of cond-evap of liq within the cloud

   !!
   !! for ice nucleation 
   !!
   
   naai_idx     = pbuf_get_index('NAAI')    !! from microp 
   naai_hom_idx = pbuf_get_index('NAAI_HOM')!! from microp 
   npccn_idx    = pbuf_get_index('NPCCN')   !! from microp 
   rndst_idx    = pbuf_get_index('RNDST')   !! from microp 
   nacon_idx    = pbuf_get_index('NACON')   !! from microp 

   frzimm_idx   = pbuf_get_index('FRZIMM',ierr)  !! from microp 
   frzcnt_idx   = pbuf_get_index('FRZCNT',ierr)  !! from microp 
   frzdep_idx   = pbuf_get_index('FRZDEP',ierr)  !! from microp 

   prec_str_idx = pbuf_get_index('PREC_STR') !! from physpkg 
   snow_str_idx = pbuf_get_index('SNOW_STR') !! from physpkg 
   prec_sed_idx = pbuf_get_index('PREC_SED') !! from physpkg 
   snow_sed_idx = pbuf_get_index('SNOW_SED') !! from physpkg 
   prec_pcw_idx = pbuf_get_index('PREC_PCW') !! from physpkg 
   snow_pcw_idx = pbuf_get_index('SNOW_PCW') !! from physpkg 

   !! not used 
   qrain_idx    = pbuf_get_index('QRAIN', ierr) !! local 
   nrain_idx    = pbuf_get_index('NRAIN', ierr) !! local 

   if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 005 -' 

   if (is_first_step() .or. is_first_restart_step()) then

      lookup_file_1 = 'p3_lookup_table_1.dat'
      lookup_file_2 = 'p3_lookup_table_2.dat'

      call micro_p3_lookuptable_init(lookup_file_1,lookup_file_2)

   end if

  ! Initialize physics buffer grid fields for accumulating precip and condensation
   if (is_first_step()) then
   
!!!      lookup_file_1 = 'p3_lookup_table_1.dat'
!!!      lookup_file_2 = 'p3_lookup_table_2.dat'
!!!      
!!!      call micro_p3_lookuptable_init(lookup_file_1,lookup_file_2)
   
      call pbuf_set_field(pbuf2d, cldo_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8)
      call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8)
      call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8)
      call pbuf_set_field(pbuf2d, acpr_idx,   0._r8)
      call pbuf_set_field(pbuf2d, acgcme_idx, 0._r8)
      call pbuf_set_field(pbuf2d, acnum_idx,  0)
      call pbuf_set_field(pbuf2d, relvar_idx, 2._r8)
      call pbuf_set_field(pbuf2d, accre_enhan_idx, micro_mg_accre_enhan_fac)
      call pbuf_set_field(pbuf2d, prer_evap_idx,  0._r8)

      !! not used 
      if (qrain_idx > 0)   call pbuf_set_field(pbuf2d, qrain_idx, 0._r8)
      if (nrain_idx > 0)   call pbuf_set_field(pbuf2d, nrain_idx, 0._r8) 

      ! If sub-columns turned on, need to set the sub-column fields as well
      if (use_subcol_microp) then
         call pbuf_set_field(pbuf2d, cldo_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_t_idx,   0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qv_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ql_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qi_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_nl_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_ni_idx,  0._r8, col_type=col_type_subcol)
         call pbuf_set_field(pbuf2d, cc_qlst_idx,0._r8, col_type=col_type_subcol)
      end if

   end if

   if(l_summary_debug) write(6,*) 'micro_p3_acme_init - 006 -' 
   
end subroutine micro_p3_acme_init




!!
!!........................................................................................
!!


subroutine micro_p3_acme_tend(state, ptend, dtime, pbuf)

   use time_manager,   only: is_first_step, is_first_restart_step
   use micro_p3_utils, only: size_dist_param_basic, &
                             size_dist_param_liq, &
                             mg_liq_props, &
                             mg_ice_props, &
                             avg_diameter, &
                             rhoi, &
                             rhosn, &
                             rhow, &
                             rhows, &
                             qsmall, &
                             mincld
   use micro_mg_data,  only: MGPacker, &
                             MGPostProc, &
                             accum_null, &
                             accum_mean
   use micro_p3,       only: micro_p3_tend, &
                             micro_p3_get_cols
   use physics_buffer, only: pbuf_col_type_index
   use subcol,         only: subcol_field_avg

   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtime
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! Local variables
   integer :: lchnk, ncol, psetcols, ngrdcol

   integer :: i, k, itim_old, it

   logical :: lq(pcnst)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! subcol values  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   real(r8), pointer :: naai(:,:)      ! ice nucleation number
   real(r8), pointer :: naai_hom(:,:)  ! ice nucleation number (homogeneous)
   real(r8), pointer :: npccn(:,:)     ! liquid activation number tendency
   real(r8), pointer :: rndst(:,:,:)
   real(r8), pointer :: nacon(:,:,:)

   !! derived 
   
   real(r8), pointer :: prec_str(:)          ! [Total] Sfc flux of precip from stratiform [ m/s ]
   real(r8), pointer :: prec_sed(:)          ! Surface flux of total cloud water from sedimentation
   real(r8), pointer :: prec_pcw(:)          ! Sfc flux of precip from microphysics [ m/s ]
   real(r8), pointer :: snow_str(:)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
   real(r8), pointer :: snow_pcw(:)          ! Sfc flux of snow from microphysics [ m/s ]
   real(r8), pointer :: snow_sed(:)          ! Surface flux of cloud ice from sedimentation

  !!  
  !! OUTPUT   
  !! 
  !!    RATE1_CW2PR_ST: needed by aero_model
  !!    CLDO     : needed by microp_aero 
  !!    PRER_EVAP: needed by clubb 
  !!    DEI      : needed by radiation 
  !!    MU       : needed by radiation 
  !!    LAMBDAC  : needed by radiation 
  !!    REL      : needed by cosp 
  !!    REI      : needed by cosp 
  !!    LS_FLXPRC   : needed by cosp 
  !!    LS_FLXSNW   : needed by cosp 
  !!    LS_REFFRAIN : needed by cosp 
  !!    LS_REFFSNOW : needed by cosp 
  !!    CV_REFFLIQ  : needed by cosp 
  !!    CV_REFFICE  : needed by cosp 
  !!    CC_xx       : needed by macro 
  !! 
  !! LOCAL
  !!
  !!    ICLWPST
  !!    ICIWPST
  !! 
  
   real(r8), pointer :: cldo(:,:)         ! Old cloud fraction
   real(r8), pointer :: prer_evap(:,:)    ! precipitation evaporation rate
   
   !! radiation 
   real(r8), pointer :: dei(:,:)          ! Ice effective diameter (um)
   real(r8), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
   real(r8), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
   
   !! wetdep 
   real(r8), pointer :: qme(:,:)          ! Total condensation rate 
   real(r8), pointer :: prain(:,:)        ! Total precipitation (rain + snow)
   real(r8), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
   
   !! COSP simulator
   real(r8), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
   real(r8), pointer :: rei(:,:)          ! Ice effective drop size (microns)
   real(r8), pointer :: mgflxprc(:,:)     ! P3 grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgflxsnw(:,:)     ! P3 grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
   real(r8), pointer :: mgreffrain_grid(:,:)   ! P3 diagnostic rain effective radius (um)
   real(r8), pointer :: mgreffsnow_grid(:,:)   ! P3 diagnostic snow effective radius (um)
   real(r8), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
   real(r8), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)
   
   !! macro 
   real(r8), pointer :: CC_T(:,:)         ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qv(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ql(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qi(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_nl(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_ni(:,:)        ! Grid-mean microphysical tendency
   real(r8), pointer :: CC_qlst(:,:)      ! In-liquid stratus microphysical tendency
   
   real(r8), pointer :: ast(:,:)          ! Relative humidity cloud fraction
   real(r8), pointer :: alst_mic(:,:)
   real(r8), pointer :: aist_mic(:,:)
   real(r8), pointer :: relvar(:,:)       ! relative variance of cloud water
   real(r8), pointer :: accre_enhan(:,:)  ! optional accretion enhancement for experimentation

   real(r8) :: rho(state%psetcols,pver)
   real(r8) :: cldmax(state%psetcols,pver)

   real(r8) :: th  (state%psetcols,pver)       ! potential temperature            K
   real(r8) :: uzpl(state%psetcols,pver)       ! vertical air velocity            m s-1
   real(r8) :: pres(state%psetcols,pver)       ! pressure                         Pa
   real(r8) :: dzq (state%psetcols,pver)       ! vertical grid spacing            m

   real(r8) :: ztop (state%psetcols,pver)       ! vertical grid spacing            m
   real(r8) :: zbot (state%psetcols,pver)       ! vertical grid spacing            m
   
   real(r8) :: pcprt_liq(pver)                 ! precipitation rate, liquid       m s-1
   real(r8) :: pcprt_sol(pver)                 ! precipitation rate, solid        m s-1
   real(r8) :: pcprt_tot(pver)                 ! precipitation rate, solid        m s-1
   
   real(r8), target :: tlat(state%psetcols,pver)
   real(r8), target :: qvlat(state%psetcols,pver)
   real(r8), target :: qcten(state%psetcols,pver)
   real(r8), target :: qiten(state%psetcols,pver)
   real(r8), target :: ncten(state%psetcols,pver)
   real(r8), target :: niten(state%psetcols,pver)
   real(r8), target :: qrten(state%psetcols,pver)
   real(r8), target :: qsten(state%psetcols,pver)
   real(r8), target :: nrten(state%psetcols,pver)
   real(r8), target :: nsten(state%psetcols,pver)
   real(r8), target :: qirimten(state%psetcols,pver)
   real(r8), target :: bvrimten(state%psetcols,pver)

   real(r8), target :: epsc(state%psetcols,pver)
   real(r8), target :: epsi(state%psetcols,pver)
   real(r8), target :: epsr(state%psetcols,pver)
   real(r8), target :: reps(state%psetcols,pver)
   real(r8), target :: ac1(state%psetcols,pver)
   real(r8), target :: ac2(state%psetcols,pver)
   real(r8), target :: ac3(state%psetcols,pver)
   real(r8), target :: qccon1(state%psetcols,pver)
   real(r8), target :: qccon2(state%psetcols,pver)
   real(r8), target :: qicon1(state%psetcols,pver)
   real(r8), target :: qicon2(state%psetcols,pver)
   real(r8), target :: qicon3(state%psetcols,pver)
   
   real(r8), target :: qcaut(state%psetcols,pver)
   real(r8), target :: ncautc(state%psetcols,pver)
   real(r8), target :: qccon(state%psetcols,pver)
   real(r8), target :: qrcon(state%psetcols,pver)
   real(r8), target :: ncautr(state%psetcols,pver)
   real(r8), target :: ncacc(state%psetcols,pver)
   real(r8), target :: qcacc(state%psetcols,pver)
   real(r8), target :: ncslf(state%psetcols,pver)
   real(r8), target :: nrslf(state%psetcols,pver)
   real(r8), target :: ncnuc(state%psetcols,pver)
   real(r8), target :: qcnuc(state%psetcols,pver)
   real(r8), target :: qcevp(state%psetcols,pver)
   real(r8), target :: qberg(state%psetcols,pver)
   real(r8), target :: qrevp(state%psetcols,pver)
   real(r8), target :: nrevp(state%psetcols,pver)
   real(r8), target :: qccol(state%psetcols,pver)
   real(r8), target :: qidep(state%psetcols,pver)
   real(r8), target :: qrcol(state%psetcols,pver)
   real(r8), target :: qinuc(state%psetcols,pver)
   real(r8), target :: nccol(state%psetcols,pver)
   real(r8), target :: nrcol(state%psetcols,pver)
   real(r8), target :: ninuc(state%psetcols,pver)
   
   real(r8), target :: qisub(state%psetcols,pver)
   real(r8), target :: qimlt(state%psetcols,pver)
   real(r8), target :: nimlt(state%psetcols,pver)
   real(r8), target :: nisub(state%psetcols,pver)
   real(r8), target :: nislf(state%psetcols,pver)
   real(r8), target :: qchetc(state%psetcols,pver)
   real(r8), target :: qcheti(state%psetcols,pver)
   real(r8), target :: qrhetc(state%psetcols,pver)
   real(r8), target :: qrheti(state%psetcols,pver)
   real(r8), target :: nchetc(state%psetcols,pver)
   real(r8), target :: ncheti(state%psetcols,pver)
   real(r8), target :: nrhetc(state%psetcols,pver)
   real(r8), target :: nrheti(state%psetcols,pver)
   real(r8), target :: nrshdr(state%psetcols,pver)
   real(r8), target :: qcshd(state%psetcols,pver)
   real(r8), target :: qrmul(state%psetcols,pver)
   real(r8), target :: nimul(state%psetcols,pver)
   real(r8), target :: ncshdc(state%psetcols,pver)
   
   real(r8), target :: rate1cld(state%psetcols,pver) ! array to hold rate1ord_cw2pr_st from microphysics

   real(r8), target  :: ze   (state%psetcols,pver) ! equivalent reflectivity     dBZ
   real(r8), target  :: umc  (state%psetcols,pver) ! mass-weighted qc fallspeed  m/s 
   real(r8), target  :: umi  (state%psetcols,pver) ! mass-weighted qi fallspeed  m/s 
   real(r8), target  :: umr  (state%psetcols,pver) ! mass-weighted qr fallspeed  m/s 
   real(r8), target  :: di   (state%psetcols,pver) ! mean size of ice            m
   real(r8), target  :: rhopo(state%psetcols,pver) ! bulk ice density            kg m-3

   real(r8), target :: prect(state%psetcols)
   real(r8), target :: preci(state%psetcols)
   
   real(r8), target :: cmeice(state%psetcols,pver)     ! Rate of cond-evap of ice within the cloud
   real(r8), target :: rflx(state%psetcols,pverp)      ! grid-box average rain flux (kg m^-2 s^-1)
   real(r8), target :: sflx(state%psetcols,pverp)      ! grid-box average snow flux (kg m^-2 s^-1) 
   real(r8), target :: cmeiout(state%psetcols,pver)    ! Deposition/sublimation rate of cloud ice
   real(r8), target :: qcsedten(state%psetcols,pver)   ! Cloud water mixing ratio tendency from sedimentation
   real(r8), target :: qisedten(state%psetcols,pver)   ! Cloud ice mixing ratio tendency from sedimentation
   real(r8), target :: qrsedten(state%psetcols,pver)   ! Rain mixing ratio tendency from sedimentation
   real(r8), target :: ncsedten(state%psetcols,pver)   ! Cloud droplet num tendency from sedimentation
   real(r8), target :: nisedten(state%psetcols,pver)   ! Cloud ice num tendency from sedimentation
   real(r8), target :: nrsedten(state%psetcols,pver)   ! Rain num tendency from sedimentation

   real(r8), target :: prao(state%psetcols,pver)
   real(r8), target :: prco(state%psetcols,pver)
   real(r8), target :: mnuccco(state%psetcols,pver)
   real(r8), target :: mnuccto(state%psetcols,pver)
   real(r8), target :: msacwio(state%psetcols,pver)
   real(r8), target :: psacwso(state%psetcols,pver)
   real(r8), target :: bergso(state%psetcols,pver)
   real(r8), target :: bergo(state%psetcols,pver)
   real(r8), target :: melto(state%psetcols,pver)
   real(r8), target :: homoo(state%psetcols,pver)
   real(r8), target :: qcreso(state%psetcols,pver)
   real(r8), target :: praio(state%psetcols,pver)
   real(r8), target :: qireso(state%psetcols,pver)
   real(r8), target :: mnuccro(state%psetcols,pver)
   real(r8), target :: frzrdt (state%psetcols,pver)
   real(r8), target :: mnuccdo(state%psetcols,pver)
   real(r8), target :: refl(state%psetcols,pver)    ! analytic radar reflectivity
   real(r8), target :: arefl(state%psetcols,pver)   ! average reflectivity will zero points outside valid range
   real(r8), target :: areflz(state%psetcols,pver)  ! average reflectivity in z.
   real(r8), target :: frefl(state%psetcols,pver)
   real(r8), target :: csrfl(state%psetcols,pver)   ! cloudsat reflectivity
   real(r8), target :: acsrfl(state%psetcols,pver)  ! cloudsat average
   real(r8), target :: fcsrfl(state%psetcols,pver)
   real(r8), target :: ncai(state%psetcols,pver)    ! output number conc of ice nuclei available (1/m3)
   real(r8), target :: ncal(state%psetcols,pver)    ! output number conc of CCN (1/m3)
   real(r8), target :: freqr(state%psetcols,pver)
   real(r8), target :: nfice(state%psetcols,pver)
   real(r8), target :: qcrat(state%psetcols,pver)   ! qc limiter ratio (1=no limit)




   real(r8), pointer :: rate1ord_cw2pr_st(:,:) ! 1st order rate for direct conversion of
                                               ! strat. cloud water to precip (1/s)    ! rce 2010/05/01



  ! variables for heterogeneous freezing
   real(r8), pointer :: frzimm(:,:)
   real(r8), pointer :: frzcnt(:,:)
   real(r8), pointer :: frzdep(:,:)

   
   ! A local copy of state is used for diagnostic calculations
   type(physics_state) :: state_loc
   type(physics_ptend) :: ptend_loc

   real(r8) :: icecldf(state%psetcols,pver) ! Ice cloud fraction
   real(r8) :: liqcldf(state%psetcols,pver) ! Liquid cloud fraction (combined into cloud)


   real(r8), pointer :: cmeliq(:,:)

   real(r8), pointer :: cld(:,:)          ! Total cloud fraction
   real(r8), pointer :: concld(:,:)       ! Convective cloud fraction
   real(r8), pointer :: iciwpst(:,:)      ! Stratiform in-cloud ice water path for radiation
   real(r8), pointer :: iclwpst(:,:)      ! Stratiform in-cloud liquid water path for radiation

   real(r8) :: icimrst(state%psetcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst(state%psetcols,pver) ! In stratus water mixing ratio
   real(r8) :: icinc(state%psetcols,pver)   ! In cloud ice number conc
   real(r8) :: icwnc(state%psetcols,pver)   ! In cloud water number conc

   real(r8) :: iclwpi(state%psetcols)       ! Vertically-integrated in-cloud Liquid WP before microphysics
   real(r8) :: iciwpi(state%psetcols)       ! Vertically-integrated in-cloud Ice WP before microphysics
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! packer 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Object that packs columns with clouds/precip.
   type(MGPacker) :: packer

   ! Output field post-processing.
   type(MGPostProc) :: post_proc
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! packed input fields 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   real(r8), allocatable :: packed_t(:,:)
   real(r8), allocatable :: packed_q(:,:)
   real(r8), allocatable :: packed_qc(:,:)
   real(r8), allocatable :: packed_nc(:,:)
   real(r8), allocatable :: packed_qi(:,:)
   real(r8), allocatable :: packed_ni(:,:)
   real(r8), allocatable :: packed_qr(:,:)
   real(r8), allocatable :: packed_nr(:,:)
   real(r8), allocatable :: packed_qirim(:,:)
   real(r8), allocatable :: packed_bvrim(:,:)
   
   real(r8), allocatable :: packed_th(:,:)
   real(r8), allocatable :: packed_uzpl(:,:)
   real(r8), allocatable :: packed_pres(:,:)
   real(r8), allocatable :: packed_dzq(:,:)
   real(r8), allocatable :: packed_relvar(:,:)
   real(r8), allocatable :: packed_accre_enhan(:,:)
   real(r8), allocatable :: packed_p(:,:)
   real(r8), allocatable :: packed_pdel(:,:)
   real(r8), allocatable :: packed_pint(:,:)
   real(r8), allocatable :: packed_cldn(:,:)
   real(r8), allocatable :: packed_liqcldf(:,:)
   real(r8), allocatable :: packed_icecldf(:,:)
   real(r8), allocatable :: packed_naai(:,:)
   real(r8), allocatable :: packed_npccn(:,:)
   real(r8), allocatable :: packed_rndst(:,:,:)
   real(r8), allocatable :: packed_nacon(:,:,:)

   real(r8), pointer :: packed_frzimm(:,:)
   real(r8), pointer :: packed_frzcnt(:,:)
   real(r8), pointer :: packed_frzdep(:,:)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! packed output fields 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   real(r8), allocatable, target :: packed_tlat(:,:)
   real(r8), allocatable, target :: packed_qvlat(:,:)
   real(r8), allocatable, target :: packed_qctend(:,:)
   real(r8), allocatable, target :: packed_nctend(:,:)
   real(r8), allocatable, target :: packed_qitend(:,:)
   real(r8), allocatable, target :: packed_nitend(:,:)
   real(r8), allocatable, target :: packed_qrtend(:,:)
   real(r8), allocatable, target :: packed_nrtend(:,:)
   real(r8), allocatable, target :: packed_qirimtend(:,:)
   real(r8), allocatable, target :: packed_bvrimtend(:,:)

   real(r8), allocatable, target :: packed_epsc(:,:)
   real(r8), allocatable, target :: packed_epsi(:,:)
   real(r8), allocatable, target :: packed_epsr(:,:)
   real(r8), allocatable, target :: packed_reps(:,:)
   real(r8), allocatable, target :: packed_ac1(:,:)
   real(r8), allocatable, target :: packed_ac2(:,:)
   real(r8), allocatable, target :: packed_ac3(:,:)
   real(r8), allocatable, target :: packed_qccon1(:,:)
   real(r8), allocatable, target :: packed_qccon2(:,:)
   real(r8), allocatable, target :: packed_qicon1(:,:)
   real(r8), allocatable, target :: packed_qicon2(:,:)
   real(r8), allocatable, target :: packed_qicon3(:,:)
   
   real(r8), allocatable, target :: packed_qccon(:,:)
   real(r8), allocatable, target :: packed_qrcon(:,:)
   real(r8), allocatable, target :: packed_qcaut(:,:)
   real(r8), allocatable, target :: packed_ncautc(:,:)
   real(r8), allocatable, target :: packed_ncautr(:,:)
   real(r8), allocatable, target :: packed_qcacc(:,:)
   real(r8), allocatable, target :: packed_ncacc(:,:)
   real(r8), allocatable, target :: packed_ncslf(:,:)
   real(r8), allocatable, target :: packed_nrslf(:,:)
   real(r8), allocatable, target :: packed_ncnuc(:,:)
   real(r8), allocatable, target :: packed_qcnuc(:,:)
   real(r8), allocatable, target :: packed_prer_evap(:,:)
   real(r8), allocatable, target :: packed_qcevp(:,:)
   real(r8), allocatable, target :: packed_qberg(:,:)
   real(r8), allocatable, target :: packed_qrevp(:,:)
   real(r8), allocatable, target :: packed_nrevp(:,:)
   
   real(r8), allocatable, target :: packed_qccol(:,:)
   real(r8), allocatable, target :: packed_qidep(:,:)
   real(r8), allocatable, target :: packed_qrcol(:,:)
   real(r8), allocatable, target :: packed_qinuc(:,:)
   real(r8), allocatable, target :: packed_nccol(:,:)
   real(r8), allocatable, target :: packed_nrcol(:,:)
   real(r8), allocatable, target :: packed_ninuc(:,:)
   real(r8), allocatable, target :: packed_qisub(:,:)
   real(r8), allocatable, target :: packed_qimlt(:,:)
   real(r8), allocatable, target :: packed_nimlt(:,:)
   real(r8), allocatable, target :: packed_nisub(:,:)
   real(r8), allocatable, target :: packed_nislf(:,:)
   real(r8), allocatable, target :: packed_qchetc(:,:)
   real(r8), allocatable, target :: packed_qcheti(:,:)
   real(r8), allocatable, target :: packed_qrhetc(:,:)
   real(r8), allocatable, target :: packed_qrheti(:,:)
   real(r8), allocatable, target :: packed_nchetc(:,:)
   real(r8), allocatable, target :: packed_ncheti(:,:)
   real(r8), allocatable, target :: packed_nrhetc(:,:)
   real(r8), allocatable, target :: packed_nrheti(:,:)
   real(r8), allocatable, target :: packed_nrshdr(:,:)
   real(r8), allocatable, target :: packed_qcshd(:,:)
   real(r8), allocatable, target :: packed_qrmul(:,:)
   real(r8), allocatable, target :: packed_nimul(:,:)
   real(r8), allocatable, target :: packed_ncshdc(:,:)
   
   
   real(r8), allocatable, target :: packed_rate1ord_cw2pr_st(:,:)
   
   real(r8), allocatable, target :: packed_rel(:,:)
   real(r8), allocatable, target :: packed_rei(:,:)
   real(r8), allocatable, target :: packed_lambdac(:,:)
   real(r8), allocatable, target :: packed_mu(:,:)
   real(r8), allocatable, target :: packed_dei(:,:)
   real(r8), allocatable, target :: packed_ze(:,:) 
   real(r8), allocatable, target :: packed_di(:,:)
   real(r8), allocatable, target :: packed_rhopo(:,:)
   real(r8), allocatable, target :: packed_umc(:,:)
   real(r8), allocatable, target :: packed_umi(:,:)
   real(r8), allocatable, target :: packed_umr(:,:)


   real(r8), allocatable, target :: packed_prect(:)
   real(r8), allocatable, target :: packed_preci(:)
   real(r8), allocatable, target :: packed_nevapr(:,:)
   real(r8), allocatable, target :: packed_prain(:,:)
   real(r8), allocatable, target :: packed_rflx(:,:)
   real(r8), allocatable, target :: packed_sflx(:,:)
   real(r8), allocatable, target :: packed_cmei(:,:)
   real(r8), allocatable, target :: packed_qcsedten(:,:)
   real(r8), allocatable, target :: packed_qisedten(:,:)
   real(r8), allocatable, target :: packed_qrsedten(:,:) 
   real(r8), allocatable, target :: packed_ncsedten(:,:)
   real(r8), allocatable, target :: packed_nisedten(:,:)
   real(r8), allocatable, target :: packed_nrsedten(:,:) 
   
   real(r8), allocatable, target :: packed_refl(:,:)
   real(r8), allocatable, target :: packed_arefl(:,:)
   real(r8), allocatable, target :: packed_areflz(:,:)
   real(r8), allocatable, target :: packed_frefl(:,:)
   real(r8), allocatable, target :: packed_csrfl(:,:)
   real(r8), allocatable, target :: packed_acsrfl(:,:)
   real(r8), allocatable, target :: packed_fcsrfl(:,:)
   real(r8), allocatable, target :: packed_ncai(:,:)
   real(r8), allocatable, target :: packed_ncal(:,:)
   real(r8), allocatable, target :: packed_freqr(:,:)
   real(r8), allocatable, target :: packed_nfice(:,:)
   real(r8), allocatable, target :: packed_qcrat(:,:)

         
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! grid mean values 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   real(r8) :: efiout_grid(pcols,pver)
   real(r8) :: efcout_grid(pcols,pver)
   real(r8) :: ncout_grid(pcols,pver)
   real(r8) :: niout_grid(pcols,pver)
   real(r8) :: freqi_grid(pcols,pver)
   real(r8) :: freql_grid(pcols,pver)

   real(r8) :: cdnumc_grid(pcols)           ! Vertically-integrated droplet concentration
   real(r8) :: icimrst_grid_out(pcols,pver) ! In stratus ice mixing ratio
   real(r8) :: icwmrst_grid_out(pcols,pver) ! In stratus water mixing ratio

   ! Cloud fraction used for precipitation.
   real(r8) :: cldmax_grid(pcols,pver)

   ! Average cloud top radius & number
   real(r8) :: ctrel_grid(pcols)
   real(r8) :: ctrei_grid(pcols)
   real(r8) :: ctnl_grid(pcols)
   real(r8) :: ctni_grid(pcols)
   real(r8) :: fcti_grid(pcols)
   real(r8) :: fctl_grid(pcols)

   real(r8) :: ftem_grid(pcols,pver)

   ! Variables for precip efficiency calculation
   real(r8) :: minlwp        ! LWP threshold

   real(r8), pointer, dimension(:) :: acprecl_grid ! accumulated precip across timesteps
   real(r8), pointer, dimension(:) :: acgcme_grid  ! accumulated condensation across timesteps
   integer,  pointer, dimension(:) :: acnum_grid   ! counter for # timesteps accumulated



 
   !!..................................................................................... 
   !! derived variables
   !!..................................................................................... 
   
   real(r8) :: tgliqwp_grid(pcols)   ! column liquid
   real(r8) :: tgcmeliq_grid(pcols)  ! column condensation rate (units)
   real(r8) :: pe_grid(pcols)        ! precip efficiency for output
   real(r8) :: pefrac_grid(pcols)    ! fraction of time precip efficiency is written out
   real(r8) :: tpr_grid(pcols)       ! average accumulated precipitation rate in pe calculation


   real(r8) :: icimrst_grid(pcols,pver) ! stratus ice mixing ratio - on grid
   real(r8) :: icwmrst_grid(pcols,pver) ! stratus water mixing ratio - on grid

   real(r8), pointer :: lambdac_grid(:,:)
   real(r8), pointer :: mu_grid(:,:)
   real(r8), pointer :: rel_grid(:,:)
   real(r8), pointer :: rei_grid(:,:)
   real(r8), pointer :: dei_grid(:,:)
   real(r8), pointer :: iclwpst_grid(:,:)

   real(r8) :: rho_grid(pcols,pver)
   real(r8) :: liqcldf_grid(pcols,pver)
   real(r8) :: ncic_grid(pcols,pver)
   real(r8) :: niic_grid(pcols,pver)
   real(r8) :: rel_fn_grid(pcols,pver)    ! Ice effective drop size at fixed number (indirect effect) (microns) - on grid

   real(r8) :: drout2_grid(pcols,pver)
   
   !! time scales 
   real(r8) :: epsc_grid(pcols,pver)
   real(r8) :: epsi_grid(pcols,pver)
   real(r8) :: epsr_grid(pcols,pver)
   real(r8) :: reps_grid(pcols,pver)
   real(r8) :: ac1_grid(pcols,pver)
   real(r8) :: ac2_grid(pcols,pver)
   real(r8) :: ac3_grid(pcols,pver)
   real(r8) :: qccon1_grid(pcols,pver)
   real(r8) :: qccon2_grid(pcols,pver)
   real(r8) :: qicon1_grid(pcols,pver)
   real(r8) :: qicon2_grid(pcols,pver)
   real(r8) :: qicon3_grid(pcols,pver)
   
   !! liquid process rates 
   real(r8) :: qcaut_grid(pcols,pver)
   real(r8) :: ncautc_grid(pcols,pver)
   real(r8) :: qccon_grid(pcols,pver)
   real(r8) :: qrcon_grid(pcols,pver)
   real(r8) :: ncautr_grid(pcols,pver)
   real(r8) :: ncacc_grid(pcols,pver)
   real(r8) :: qcacc_grid(pcols,pver)
   real(r8) :: ncslf_grid(pcols,pver)
   real(r8) :: nrslf_grid(pcols,pver)
   real(r8) :: ncnuc_grid(pcols,pver)
   real(r8) :: qcnuc_grid(pcols,pver)
   real(r8) :: qcevp_grid(pcols,pver)
   real(r8) :: qberg_grid(pcols,pver)
   real(r8) :: qrevp_grid(pcols,pver)
   real(r8) :: nrevp_grid(pcols,pver)
   
   !! ice process rates 
   real(r8) :: qccol_grid(pcols,pver)
   real(r8) :: qidep_grid(pcols,pver)
   real(r8) :: qrcol_grid(pcols,pver)
   real(r8) :: qinuc_grid(pcols,pver)
   real(r8) :: nccol_grid(pcols,pver)
   real(r8) :: nrcol_grid(pcols,pver)
   real(r8) :: ninuc_grid(pcols,pver)
   real(r8) :: qisub_grid(pcols,pver)
   real(r8) :: qimlt_grid(pcols,pver)
   real(r8) :: nimlt_grid(pcols,pver)
   real(r8) :: nisub_grid(pcols,pver)
   real(r8) :: nislf_grid(pcols,pver)
   real(r8) :: qchetc_grid(pcols,pver)
   real(r8) :: qcheti_grid(pcols,pver)
   real(r8) :: qrhetc_grid(pcols,pver)
   real(r8) :: qrheti_grid(pcols,pver)
   real(r8) :: nchetc_grid(pcols,pver)
   real(r8) :: ncheti_grid(pcols,pver)
   real(r8) :: nrhetc_grid(pcols,pver)
   real(r8) :: nrheti_grid(pcols,pver)
   real(r8) :: nrshdr_grid(pcols,pver)
   real(r8) :: qcshd_grid(pcols,pver)
   real(r8) :: qrmul_grid(pcols,pver)
   real(r8) :: nimul_grid(pcols,pver)
   real(r8) :: ncshdc_grid(pcols,pver)

   real(r8) :: prer_evap_grid(pcols,pver)
   
   real(r8) :: reff_rain_grid(pcols,pver)
   real(r8) :: cld_grid(pcols,pver)
   real(r8) :: pdel_grid(pcols,pver)
   real(r8) :: prco_grid(pcols,pver)
   real(r8) :: prao_grid(pcols,pver)
   real(r8) :: icecldf_grid(pcols,pver)
   real(r8) :: icwnc_grid(pcols,pver)
   real(r8) :: icinc_grid(pcols,pver)
   real(r8) :: qcreso_grid(pcols,pver)
   real(r8) :: melto_grid(pcols,pver)
   real(r8) :: mnuccco_grid(pcols,pver)
   real(r8) :: mnuccto_grid(pcols,pver)
   real(r8) :: bergo_grid(pcols,pver)
   real(r8) :: homoo_grid(pcols,pver)
   real(r8) :: msacwio_grid(pcols,pver)
   real(r8) :: psacwso_grid(pcols,pver)
   real(r8) :: bergso_grid(pcols,pver)
   real(r8) :: cmeiout_grid(pcols,pver)
   real(r8) :: qireso_grid(pcols,pver)
   real(r8) :: praio_grid(pcols,pver)

   real(r8) :: nc_grid(pcols,pver)
   real(r8) :: ni_grid(pcols,pver)
   real(r8) :: qr_grid(pcols,pver)
   real(r8) :: nr_grid(pcols,pver)
   real(r8) :: qirim_grid(pcols,pver)
   real(r8) :: bvrim_grid(pcols,pver)

   real(r8), pointer :: cmeliq_grid(:,:)

   real(r8), pointer :: prec_str_grid(:)
   real(r8), pointer :: snow_str_grid(:)
   real(r8), pointer :: prec_pcw_grid(:)
   real(r8), pointer :: snow_pcw_grid(:)
   real(r8), pointer :: prec_sed_grid(:)
   real(r8), pointer :: snow_sed_grid(:)
   real(r8), pointer :: cldo_grid(:,:)
   real(r8), pointer :: nevapr_grid(:,:)
   real(r8), pointer :: prain_grid(:,:)
   real(r8), pointer :: mgflxprc_grid(:,:) !! for cosp 
   real(r8), pointer :: mgflxsnw_grid(:,:) !! for cosp  
   real(r8), pointer :: cvreffliq_grid(:,:)
   real(r8), pointer :: cvreffice_grid(:,:)
   real(r8), pointer :: rate1ord_cw2pr_st_grid(:,:)
   real(r8), pointer :: CC_t_grid(:,:)
   real(r8), pointer :: CC_qv_grid(:,:)
   real(r8), pointer :: CC_ql_grid(:,:)
   real(r8), pointer :: CC_qi_grid(:,:)
   real(r8), pointer :: CC_nl_grid(:,:)
   real(r8), pointer :: CC_ni_grid(:,:)
   real(r8), pointer :: CC_qlst_grid(:,:)
   real(r8), pointer :: iciwpst_grid(:,:)
   real(r8), pointer :: ast_grid(:,:)
   real(r8), pointer :: qme_grid(:,:)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! others 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   integer :: nlev   ! number of levels where cloud physics is done
   integer :: mgncol ! size of mgcols
   integer, allocatable :: mgcols(:) ! Columns with microphysics performed

   logical :: use_subcol_microp
   integer :: col_type ! Flag to store whether accessing grid or sub-columns in pbuf_get_field

   character(128) :: errstring   ! return status (non-blank for error return)

   ! For rrtmg optics. specified distribution.
   real(r8), parameter :: dcon   = 25.e-6_r8         ! Convective size distribution effective radius (um)
   real(r8), parameter :: mucon  = 5.3_r8            ! Convective size distribution shape parameter
   real(r8), parameter :: deicon = 50._r8            ! Convective ice effective diameter (um)

   real(r8), pointer :: pckdptr(:,:)

   logical :: first_step 
   


   call t_startf('micro_p3_acme_tend')
   
   call t_startf('micro_p3_acme_tend_init')
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! micro_p3_acme_tend starts !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 001 -' 

   ! Find the number of levels used in the microphysics.
   nlev  = pver - top_lev + 1

   lchnk = state%lchnk
   ncol  = state%ncol
   psetcols = state%psetcols
   ngrdcol  = state%ngrdcol

   itim_old = pbuf_old_tim_idx()

   call phys_getopts(use_subcol_microp_out=use_subcol_microp)

   ! Set the col_type flag to grid or subcolumn dependent on the value of use_subcol_microp
   call pbuf_col_type_index(use_subcol_microp, col_type=col_type)

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 002 -' 
   
   !!.....................................................................................
   !! input 
   !!.....................................................................................

   call pbuf_get_field(pbuf, naai_idx, &
                             naai,     &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, naai_hom_idx, &
                             naai_hom,     &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, npccn_idx,   &
                             npccn,       &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, rndst_idx,   &
                             rndst,       &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, nacon_idx,   &
                             nacon,       &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, relvar_idx,      &
                             relvar,      &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, accre_enhan_idx, &
                             accre_enhan, &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, cmeliq_idx,  &
                             cmeliq,      &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)

   call pbuf_get_field(pbuf, cld_idx, &
                             cld,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, concld_idx, &
                             concld,  &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)
   call pbuf_get_field(pbuf, ast_idx, &
                             ast,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type, &
                             copy_if_needed=use_subcol_microp)

   if (use_hetfrz_classnuc) then
      call pbuf_get_field(pbuf, frzimm_idx, &
                                frzimm, &
                                col_type=col_type, &
                                copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzcnt_idx, &
                                frzcnt, &
                                col_type=col_type, &
                                copy_if_needed=use_subcol_microp)
      call pbuf_get_field(pbuf, frzdep_idx, &
                                frzdep, &
                                col_type=col_type, &
                                copy_if_needed=use_subcol_microp)
   end if

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 003 -' 
   
   !!.....................................................................................
   !! output 
   !!.....................................................................................

   call pbuf_get_field(pbuf, prec_str_idx,    prec_str,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_str_idx,    snow_str,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw,    col_type=col_type)
   call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed,    col_type=col_type)
   call pbuf_get_field(pbuf, nevapr_idx,      nevapr,      col_type=col_type)
   call pbuf_get_field(pbuf, prer_evap_idx,   prer_evap,   col_type=col_type)
   call pbuf_get_field(pbuf, prain_idx,       prain,       col_type=col_type)
   call pbuf_get_field(pbuf, dei_idx,         dei,         col_type=col_type)
   call pbuf_get_field(pbuf, mu_idx,          mu,          col_type=col_type)
   call pbuf_get_field(pbuf, lambdac_idx,     lambdac,     col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc,    col_type=col_type)
   call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw,    col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq,   col_type=col_type)
   call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice,   col_type=col_type)
   call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst,     col_type=col_type)
   call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst,     col_type=col_type) 
   call pbuf_get_field(pbuf, rel_idx,         rel,         col_type=col_type)
   call pbuf_get_field(pbuf, rei_idx,         rei,         col_type=col_type)
   call pbuf_get_field(pbuf, qme_idx,         qme,         col_type=col_type)

   call pbuf_get_field(pbuf, cldo_idx, &
                             cldo,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_t_idx, &
                             CC_t,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_qv_idx, &
                             CC_qv,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_ql_idx, &
                             CC_ql,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_qi_idx, &
                             CC_qi,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_nl_idx, &
                             CC_nl,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_ni_idx, &
                             CC_ni,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)
   call pbuf_get_field(pbuf, cc_qlst_idx, &
                             CC_qlst,     &
                             start=(/1,1,itim_old/), &
                             kount=(/psetcols,pver,1/), &
                             col_type=col_type)

   if (rate1_cw2pr_st_idx > 0) then
      call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, &
                                rate1ord_cw2pr_st, &
                                col_type=col_type)
   end if

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 004 -' 
   
   !!.....................................................................................
   !! grid values 
   !!.....................................................................................

   if (use_subcol_microp) then
   
      call pbuf_get_field(pbuf, prec_str_idx,    prec_str_grid)
      call pbuf_get_field(pbuf, snow_str_idx,    snow_str_grid)
      call pbuf_get_field(pbuf, prec_pcw_idx,    prec_pcw_grid)
      call pbuf_get_field(pbuf, snow_pcw_idx,    snow_pcw_grid)
      call pbuf_get_field(pbuf, prec_sed_idx,    prec_sed_grid)
      call pbuf_get_field(pbuf, snow_sed_idx,    snow_sed_grid)
      call pbuf_get_field(pbuf, nevapr_idx,      nevapr_grid)
      call pbuf_get_field(pbuf, prain_idx,       prain_grid)
      call pbuf_get_field(pbuf, dei_idx,         dei_grid)
      call pbuf_get_field(pbuf, mu_idx,          mu_grid)
      call pbuf_get_field(pbuf, lambdac_idx,     lambdac_grid)
      call pbuf_get_field(pbuf, ls_flxprc_idx,   mgflxprc_grid)
      call pbuf_get_field(pbuf, ls_flxsnw_idx,   mgflxsnw_grid)
      call pbuf_get_field(pbuf, cv_reffliq_idx,  cvreffliq_grid)
      call pbuf_get_field(pbuf, cv_reffice_idx,  cvreffice_grid)
      call pbuf_get_field(pbuf, iciwpst_idx,     iciwpst_grid)
      call pbuf_get_field(pbuf, iclwpst_idx,     iclwpst_grid)
      call pbuf_get_field(pbuf, rel_idx,         rel_grid)
      call pbuf_get_field(pbuf, rei_idx,         rei_grid)
      call pbuf_get_field(pbuf, qme_idx,         qme_grid)

      call pbuf_get_field(pbuf, cldo_idx,     &
                                cldo_grid,     &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_t_idx,     &
                                CC_t_grid,     &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qv_idx,    &
                                CC_qv_grid,    &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ql_idx,    &
                                CC_ql_grid,    &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qi_idx,    &
                                CC_qi_grid,    &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_nl_idx,    &
                                CC_nl_grid,    &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_ni_idx,    &
                                CC_ni_grid,    &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))
      call pbuf_get_field(pbuf, cc_qlst_idx,  &
                                CC_qlst_grid,  &
                                start=(/1,1,itim_old/), &
                                kount=(/pcols,pver,1/))

      if (rate1_cw2pr_st_idx > 0) then
         call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, &
                                   rate1ord_cw2pr_st_grid)
      end if

   end if
   
   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 005 -' 
   
   !!.....................................................................................
   !! grid only values 
   !!.....................................................................................
   
   call pbuf_get_field(pbuf, ls_reffrain_idx, mgreffrain_grid)
   call pbuf_get_field(pbuf, ls_reffsnow_idx, mgreffsnow_grid)
   call pbuf_get_field(pbuf, acpr_idx,        acprecl_grid)
   call pbuf_get_field(pbuf, acgcme_idx,      acgcme_grid)
   call pbuf_get_field(pbuf, acnum_idx,       acnum_grid)
   call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq_grid)
   call pbuf_get_field(pbuf, ast_idx,  &
                             ast_grid, &
                             start=(/1,1,itim_old/), &
                             kount=(/pcols,pver,1/))

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! some init 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !! liquid stratus frac 
   !!    = ice stratus frac
   !!    = max( liquid stratus frac, ice stratus frac

   alst_mic => ast
   aist_mic => ast

   ! Output initial in-cloud LWP (before microphysics)

   iclwpi = 0._r8
   iciwpi = 0._r8

   do i = 1, ncol
      do k = top_lev, pver
         iclwpi(i) = iclwpi(i) + &
                     min(state%q(i,k,ixcldliq) / max(mincld,ast(i,k)),0.005_r8) &
                   * state%pdel(i,k) / gravit
         iciwpi(i) = iciwpi(i) + &
                     min(state%q(i,k,ixcldice) / max(mincld,ast(i,k)),0.005_r8) &
                   * state%pdel(i,k) / gravit
      end do
   end do

   cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! local state 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   call physics_state_copy(state, state_loc)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! init ptend  
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   lq           = .false.
   lq(1)        = .true.
   lq(ixcldliq) = .true.
   lq(ixcldice) = .true.
   lq(ixnumliq) = .true.
   lq(ixnumice) = .true.
   lq(ixrain)   = .true.
   lq(ixcldrim) = .true.
   lq(ixnumrain)= .true.
   lq(ixbvrim)  = .true.

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 006 -' 
   
   ! the name 'cldwat' triggers special tests on cldliq
   ! and cldice in physics_update
   
   call physics_ptend_init(ptend, psetcols, "cldwat", ls=.true., lq=lq)


   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 007 -' 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! Determines which columns microphysics should operate over by
   !! checking for non-zero cloud water/ice.
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call micro_p3_get_cols( &
           ncol, &
           nlev, &
           top_lev, &
           state%q(:,:,ixcldliq), &
           state%q(:,:,ixcldice), &
           state%q(:,:,ixrain), &
           mgncol, & !! number of columns scheme will use 
           mgcols)   !! column indices


   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 008 -' 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! prepare packer 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   packer = MGPacker(psetcols, pver, mgcols, top_lev)
   
   post_proc = MGPostProc(packer)


   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 009 -' 
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! pack data 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !!
   !! tendency of T, Q, and tracers 
   !!
   
   allocate(packed_tlat(mgncol,nlev))
      call post_proc%add_field(p(tlat), p(packed_tlat))
      
   allocate(packed_qvlat(mgncol,nlev))
      call post_proc%add_field(p(qvlat), p(packed_qvlat))
      
   allocate(packed_qctend(mgncol,nlev))
      call post_proc%add_field(p(qcten), p(packed_qctend))
   
   allocate(packed_nctend(mgncol,nlev))
      call post_proc%add_field(p(ncten), p(packed_nctend))
   
   allocate(packed_qitend(mgncol,nlev))
      call post_proc%add_field(p(qiten), p(packed_qitend))
   
   allocate(packed_nitend(mgncol,nlev))
      call post_proc%add_field(p(niten), p(packed_nitend))

   allocate(packed_qrtend(mgncol,nlev))
      call post_proc%add_field(p(qrten), p(packed_qrtend))
   
   allocate(packed_nrtend(mgncol,nlev))
      call post_proc%add_field(p(nrten), p(packed_nrtend))
   
   allocate(packed_qirimtend(mgncol,nlev))
      call post_proc%add_field(p(qirimten), p(packed_qirimtend))

   allocate(packed_bvrimtend(mgncol,nlev))
      call post_proc%add_field(p(bvrimten), p(packed_bvrimtend))

   !! time scales 
   allocate(packed_epsc(mgncol,nlev))
      call post_proc%add_field(p(epsc), p(packed_epsc))
   allocate(packed_epsi(mgncol,nlev))
      call post_proc%add_field(p(epsi), p(packed_epsi))
   allocate(packed_epsr(mgncol,nlev))
      call post_proc%add_field(p(epsr), p(packed_epsr))
   allocate(packed_reps(mgncol,nlev))
      call post_proc%add_field(p(reps), p(packed_reps))
   allocate(packed_ac1(mgncol,nlev))
      call post_proc%add_field(p(ac1), p(packed_ac1))
   allocate(packed_ac2(mgncol,nlev))
      call post_proc%add_field(p(ac2), p(packed_ac2))
   allocate(packed_ac3(mgncol,nlev))
      call post_proc%add_field(p(ac3), p(packed_ac3))
   allocate(packed_qccon1(mgncol,nlev))
      call post_proc%add_field(p(qccon1), p(packed_qccon1))
   allocate(packed_qccon2(mgncol,nlev))
      call post_proc%add_field(p(qccon2), p(packed_qccon2))
   allocate(packed_qicon1(mgncol,nlev))
      call post_proc%add_field(p(qicon1), p(packed_qicon1))
   allocate(packed_qicon2(mgncol,nlev))
      call post_proc%add_field(p(qicon2), p(packed_qicon2))
   allocate(packed_qicon3(mgncol,nlev))
      call post_proc%add_field(p(qicon3), p(packed_qicon3))
   
   !!
   !! process rates liquid 
   !!
   
   allocate(packed_qccon(mgncol,nlev))
      call post_proc%add_field(p(qccon), p(packed_qccon))
      
   allocate(packed_qrcon(mgncol,nlev))
      call post_proc%add_field(p(qrcon), p(packed_qrcon))
   
   allocate(packed_qcaut(mgncol,nlev))
      call post_proc%add_field(p(qcaut), p(packed_qcaut))

   allocate(packed_ncautc(mgncol,nlev))
      call post_proc%add_field(p(ncautc), p(packed_ncautc))

   allocate(packed_ncautr(mgncol,nlev))
      call post_proc%add_field(p(ncautr), p(packed_ncautr))

   allocate(packed_ncacc(mgncol,nlev))
      call post_proc%add_field(p(ncacc), p(packed_ncacc))

   allocate(packed_qcacc(mgncol,nlev))
      call post_proc%add_field(p(qcacc), p(packed_qcacc))

   allocate(packed_ncslf(mgncol,nlev))
      call post_proc%add_field(p(ncslf), p(packed_ncslf))

   allocate(packed_nrslf(mgncol,nlev))
      call post_proc%add_field(p(nrslf), p(packed_nrslf))

   allocate(packed_ncnuc(mgncol,nlev))
      call post_proc%add_field(p(ncnuc), p(packed_ncnuc))

   allocate(packed_qcnuc(mgncol,nlev))
      call post_proc%add_field(p(qcnuc), p(packed_qcnuc))

   allocate(packed_qcevp(mgncol,nlev))
      call post_proc%add_field(p(qcevp), p(packed_qcevp))

   allocate(packed_qberg(mgncol,nlev))
      call post_proc%add_field(p(qberg), p(packed_qberg))

   allocate(packed_qrevp(mgncol,nlev))
      call post_proc%add_field(p(qrevp), p(packed_qrevp))

   allocate(packed_nrevp(mgncol,nlev))
      call post_proc%add_field(p(nrevp), p(packed_nrevp))

   allocate(packed_qccol(mgncol,nlev))
      call post_proc%add_field(p(qccol), p(packed_qccol))

   !!
   !! process rates ice 
   !!
   
   allocate(packed_qidep(mgncol,nlev))
      call post_proc%add_field(p(qidep), p(packed_qidep))

   allocate(packed_qrcol(mgncol,nlev))
      call post_proc%add_field(p(qrcol), p(packed_qrcol))

   allocate(packed_qinuc(mgncol,nlev))
      call post_proc%add_field(p(qinuc), p(packed_qinuc))

   allocate(packed_nccol(mgncol,nlev))
      call post_proc%add_field(p(nccol), p(packed_nccol))

   allocate(packed_nrcol(mgncol,nlev))
      call post_proc%add_field(p(nrcol), p(packed_nrcol))

   allocate(packed_ninuc(mgncol,nlev))
      call post_proc%add_field(p(ninuc), p(packed_ninuc))

   allocate(packed_qisub(mgncol,nlev))
      call post_proc%add_field(p(qisub), p(packed_qisub))

   allocate(packed_qimlt(mgncol,nlev))
      call post_proc%add_field(p(qimlt), p(packed_qimlt))

   allocate(packed_nimlt(mgncol,nlev))
      call post_proc%add_field(p(nimlt), p(packed_nimlt))

   allocate(packed_nisub(mgncol,nlev))
      call post_proc%add_field(p(nisub), p(packed_nisub))

   allocate(packed_nislf(mgncol,nlev))
      call post_proc%add_field(p(nislf), p(packed_nislf))

   allocate(packed_qchetc(mgncol,nlev))
      call post_proc%add_field(p(qchetc), p(packed_qchetc))

   allocate(packed_qcheti(mgncol,nlev))
      call post_proc%add_field(p(qcheti), p(packed_qcheti))

   allocate(packed_qrhetc(mgncol,nlev))
      call post_proc%add_field(p(qrhetc), p(packed_qrhetc))

   allocate(packed_qrheti(mgncol,nlev))
      call post_proc%add_field(p(qrheti), p(packed_qrheti))

   allocate(packed_nchetc(mgncol,nlev))
      call post_proc%add_field(p(nchetc), p(packed_nchetc))

   allocate(packed_ncheti(mgncol,nlev))
      call post_proc%add_field(p(ncheti), p(packed_ncheti))

   allocate(packed_nrhetc(mgncol,nlev))
      call post_proc%add_field(p(nrhetc), p(packed_nrhetc))

   allocate(packed_nrheti(mgncol,nlev))
      call post_proc%add_field(p(nrheti), p(packed_nrheti))

   allocate(packed_nrshdr(mgncol,nlev))
      call post_proc%add_field(p(nrshdr), p(packed_nrshdr))

   allocate(packed_qcshd(mgncol,nlev))
      call post_proc%add_field(p(qcshd), p(packed_qcshd))

   allocate(packed_qrmul(mgncol,nlev))
      call post_proc%add_field(p(qrmul), p(packed_qrmul))

   allocate(packed_nimul(mgncol,nlev))
      call post_proc%add_field(p(nimul), p(packed_nimul))

   allocate(packed_ncshdc(mgncol,nlev))
      call post_proc%add_field(p(ncshdc), p(packed_ncshdc))

   allocate(packed_qcsedten(mgncol,nlev))
      call post_proc%add_field(p(qcsedten), p(packed_qcsedten))
   
   allocate(packed_qisedten(mgncol,nlev))
      call post_proc%add_field(p(qisedten), p(packed_qisedten))

   allocate(packed_qrsedten(mgncol,nlev))
      call post_proc%add_field(p(qrsedten), p(packed_qrsedten))

   allocate(packed_ncsedten(mgncol,nlev))
      call post_proc%add_field(p(ncsedten), p(packed_ncsedten))
   
   allocate(packed_nisedten(mgncol,nlev))
      call post_proc%add_field(p(nisedten), p(packed_nisedten))

   allocate(packed_nrsedten(mgncol,nlev))
      call post_proc%add_field(p(nrsedten), p(packed_nrsedten))
      
   !!
   !! diag fields 
   !!

   allocate(packed_rate1ord_cw2pr_st(mgncol,nlev))
      pckdptr => packed_rate1ord_cw2pr_st 
      call post_proc%add_field(p(rate1cld), pckdptr)
      
   allocate(packed_prect(mgncol))
      call post_proc%add_field(p(prect), p(packed_prect))
      
   allocate(packed_preci(mgncol))
      call post_proc%add_field(p(preci), p(packed_preci))

   allocate(packed_nevapr(mgncol,nlev))
      call post_proc%add_field(p(nevapr), p(packed_nevapr))
    
   allocate(packed_prain(mgncol,nlev))
      call post_proc%add_field(p(prain), p(packed_prain))
   
   allocate(packed_umc(mgncol,nlev))
      call post_proc%add_field(p(umc), p(packed_umc))
      
   allocate(packed_umi(mgncol,nlev))
      call post_proc%add_field(p(umi), p(packed_umi))
   
   allocate(packed_umr(mgncol,nlev))
      call post_proc%add_field(p(umr), p(packed_umr))

   allocate(packed_rflx(mgncol,nlev+1))
      call post_proc%add_field(p(rflx), p(packed_rflx))
      
   allocate(packed_sflx(mgncol,nlev+1))
      call post_proc%add_field(p(sflx), p(packed_sflx))
      
   allocate(packed_cmei(mgncol,nlev))
      call post_proc%add_field(p(cmeiout), p(packed_cmei))

   

!!!
!!!   allocate(packed_refl(mgncol,nlev))
!!!   call post_proc%add_field(p(refl), p(packed_refl), fillvalue=-9999._r8)
!!!   allocate(packed_arefl(mgncol,nlev))
!!!   call post_proc%add_field(p(arefl), p(packed_arefl))
!!!   allocate(packed_areflz(mgncol,nlev))
!!!   call post_proc%add_field(p(areflz), p(packed_areflz))
!!!   allocate(packed_frefl(mgncol,nlev))
!!!   call post_proc%add_field(p(frefl), p(packed_frefl))
!!!   allocate(packed_csrfl(mgncol,nlev))
!!!   call post_proc%add_field(p(csrfl), p(packed_csrfl), fillvalue=-9999._r8)
!!!   allocate(packed_acsrfl(mgncol,nlev))
!!!   call post_proc%add_field(p(acsrfl), p(packed_acsrfl))
!!!   allocate(packed_fcsrfl(mgncol,nlev))
!!!   call post_proc%add_field(p(fcsrfl), p(packed_fcsrfl))
!!!
   allocate(packed_ncai(mgncol,nlev))
      call post_proc%add_field(p(ncai), p(packed_ncai))
      
   allocate(packed_ncal(mgncol,nlev))
      call post_proc%add_field(p(ncal), p(packed_ncal))

!!!   allocate(packed_freqr(mgncol,nlev))
!!!     call post_proc%add_field(p(freqr), p(packed_freqr))
!!!   allocate(packed_nfice(mgncol,nlev))
!!!     call post_proc%add_field(p(nfice), p(packed_nfice))
!!!   allocate(packed_qcrat(mgncol,nlev))
!!!     call post_proc%add_field(p(qcrat), p(packed_qcrat), fillvalue=1._r8)
!!!
!!!
!!!   ! The following are all variables related to sizes, where it does not
!!!   ! necessarily make sense to average over time steps. Instead, we keep
!!!   ! the value from the last substep, which is what "accum_null" does.
!!!   

   allocate(packed_rel(mgncol,nlev))
      call post_proc%add_field(p(rel), p(packed_rel), &
           fillvalue=10._r8, accum_method=accum_null)

   allocate(packed_rei(mgncol,nlev))
      call post_proc%add_field(p(rei), p(packed_rei), &
           fillvalue=25._r8, accum_method=accum_null)

   allocate(packed_mu(mgncol,nlev))
      call post_proc%add_field(p(mu), p(packed_mu), &
           accum_method=accum_null)
        
   allocate(packed_lambdac(mgncol,nlev))
      call post_proc%add_field(p(lambdac), p(packed_lambdac), &
           accum_method=accum_null)
   
   allocate(packed_dei(mgncol,nlev))
      call post_proc%add_field(p(dei), p(packed_dei), &
           accum_method=accum_null)
      
   allocate(packed_prer_evap(mgncol,nlev))
      call post_proc%add_field(p(prer_evap), p(packed_prer_evap), &
           accum_method=accum_null)
        
   allocate(packed_ze(mgncol,nlev))
      call post_proc%add_field(p(ze), p(packed_ze), &
           accum_method=accum_null)
        
   allocate(packed_di(mgncol,nlev))
      call post_proc%add_field(p(di), p(packed_di), &
           accum_method=accum_null)
        
   allocate(packed_rhopo(mgncol,nlev))
      call post_proc%add_field(p(rhopo), p(packed_rhopo), &
           accum_method=accum_null)

   !!
   !! Pack input variables that are not updated during substeps.
   !!
   
   allocate(packed_relvar(mgncol,nlev))
      packed_relvar = packer%pack(relvar)
      
   allocate(packed_accre_enhan(mgncol,nlev))
      packed_accre_enhan = packer%pack(accre_enhan)

   allocate(packed_p(mgncol,nlev))
      packed_p = packer%pack(state_loc%pmid)
      
   allocate(packed_pdel(mgncol,nlev))
      packed_pdel = packer%pack(state_loc%pdel)

   allocate(packed_pint(mgncol,nlev+1))
      packed_pint = packer%pack_interface(state_loc%pint)

   allocate(packed_cldn(mgncol,nlev))
      packed_cldn = packer%pack(ast)
      
   allocate(packed_liqcldf(mgncol,nlev))
      packed_liqcldf = packer%pack(alst_mic)
      
   allocate(packed_icecldf(mgncol,nlev))
      packed_icecldf = packer%pack(aist_mic)

   allocate(packed_naai(mgncol,nlev))
      packed_naai = packer%pack(naai)
   
   allocate(packed_npccn(mgncol,nlev))
      packed_npccn = packer%pack(npccn)

   allocate(packed_rndst(mgncol,nlev,size(rndst, 3)))
      packed_rndst = packer%pack(rndst)
      
   allocate(packed_nacon(mgncol,nlev,size(nacon, 3)))
      packed_nacon = packer%pack(nacon)
 
   if (use_hetfrz_classnuc) then
      allocate(packed_frzimm(mgncol,nlev))
         packed_frzimm = packer%pack(frzimm)
      allocate(packed_frzcnt(mgncol,nlev))
         packed_frzcnt = packer%pack(frzcnt)
      allocate(packed_frzdep(mgncol,nlev))
         packed_frzdep = packer%pack(frzdep)
   else
      nullify(packed_frzimm)
      nullify(packed_frzcnt)
      nullify(packed_frzdep)
   end if


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! prepare needed fields and pack them 
   !!
   !!    density 
   !!    vertical velocity in m/s 
   !!    potential temperature 
   !!    layer height 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   rho(:ncol,top_lev:) = state_loc%pmid(:ncol,top_lev:) &
                         / (rair*state_loc%t(:ncol,top_lev:))
                         
   pres(:ncol,top_lev:) = state_loc%pmid(:ncol,top_lev:) 

   uzpl(:ncol,top_lev:) = -1._r8*state_loc%omega(:ncol,top_lev:) &
                          / (rho(:ncol,top_lev:)*gravit)

   th(:ncol,top_lev:) = state_loc%t(:ncol,top_lev:) &
                        * state_loc%exner(:ncol,top_lev:)
 
   do k=top_lev,pver
        ztop(1:ncol,k) = state_loc%zi(1:ncol,k)
        zbot(1:ncol,k) = state_loc%zi(1:ncol,k+1)
         dzq(1:ncol,k) = ztop(1:ncol,k) - zbot(1:ncol,k)
   end do
   
   allocate(packed_th(mgncol,nlev))
      packed_th = packer%pack(th)
   allocate(packed_uzpl(mgncol,nlev))
      packed_uzpl = packer%pack(uzpl)
   allocate(packed_pres(mgncol,nlev))
      packed_pres = packer%pack(pres)
   allocate(packed_dzq(mgncol,nlev))
      packed_dzq = packer%pack(dzq)
      
   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 010 -' 
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! Allocate input variables that are updated during substeps
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   allocate(packed_t(mgncol,nlev))
   allocate(packed_q(mgncol,nlev))
   allocate(packed_qc(mgncol,nlev))
   allocate(packed_nc(mgncol,nlev))
   allocate(packed_qi(mgncol,nlev))
   allocate(packed_ni(mgncol,nlev))
   allocate(packed_qr(mgncol,nlev))
   allocate(packed_nr(mgncol,nlev))
   allocate(packed_qirim(mgncol,nlev))
   allocate(packed_bvrim(mgncol,nlev))

   call t_stopf('micro_p3_acme_tend_init')

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 011 -' 



   call t_startf('micro_p3_acme_tend_loop')
   
   first_step = is_first_step()
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! P3 microphysics subcycle 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   do it = 1, num_steps

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 012 -' 

      ! Pack input variables that are updated during substeps.
      
      packed_t  = packer%pack(state_loc%t)
      packed_q  = packer%pack(state_loc%q(:,:,1))
      packed_qc = packer%pack(state_loc%q(:,:,ixcldliq))
      packed_nc = packer%pack(state_loc%q(:,:,ixnumliq))
      packed_qi = packer%pack(state_loc%q(:,:,ixcldice))
      packed_ni = packer%pack(state_loc%q(:,:,ixnumice))
      packed_qr = packer%pack(state_loc%q(:,:,ixrain))
      packed_nr = packer%pack(state_loc%q(:,:,ixnumrain))
      packed_qirim = packer%pack(state_loc%q(:,:,ixcldrim))
      packed_bvrim = packer%pack(state_loc%q(:,:,ixbvrim))

      call t_startf('micro_p3')

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 013 -'
   if(l_summary_debug) write(6,*) 'mgncol, top_lev : ', mgncol, top_lev 

      call micro_p3_tend( &
!!
!! in 
!!
         mgncol,          &! number of columns 
         top_lev,         &! 
         pver,            &! 
         first_step,      &!
         dtime/num_steps, &! model time step (s)
         scale_qidep,     & 
         scale_berg,      & 
         l_mg2_qidep,     &
         l_satadj, & 
         l_crconevp,      &
         l_massclip ,     &
         l_limit_qidep_qinuc, &
         l_limit_qisub_qrevp, &
         l_cshd,          &
         l_imlt,          &
         l_ccol,          &
         opt_cheti,    &
         opt_inuc,     & 
         packed_pres,     &! pressure (Pa) 
         packed_pdel,     &! pressure thickness 
         packed_uzpl,     &! vertical air velocity (m s-1)
         packed_dzq,      &! vertical grid spacing (m) 
         packed_cldn,     &! cloud fraction 
         packed_liqcldf,  &! liquid cloud fraction 
         packed_icecldf,  &! ice cloud fraction 
         packed_relvar,   & 
         packed_accre_enhan, &
         packed_naai,     &
         packed_npccn,    &
         packed_rndst,    &
         packed_nacon,    &
         packed_frzimm,   & 
         packed_frzcnt,   & 
         packed_frzdep,   & 
!!
!! in/out 
!!
         packed_t,        &! temperature (K) 
         packed_q,        &! humidity (kg/kg) 
         packed_qc,       &! 
         packed_nc,       &!
         packed_qr,       &!
         packed_nr,       &!
         packed_qi,       &!
         packed_ni,       &!
         packed_qirim,    &! ice, rime mass mixing ratio      kg kg-1
         packed_bvrim,    &! ice, rime volume mixing ratio    m3 kg-1
!!
!! out tendency 
!!
         packed_tlat,     &
         packed_qvlat,    &
         packed_qctend,   &
         packed_nctend,   &
         packed_qrtend,   &
         packed_nrtend,   &
         packed_qitend,   &
         packed_nitend,   &
         packed_qirimtend,&
         packed_bvrimtend,&
!!
!! out diag 
!!
         packed_prect,    &! surface precip rate (m/s)
         packed_preci,    &! cloud ice/snow precip rate (m/s)
         packed_prain,    &! Total precipitation (rain + snow)
         packed_prer_evap, &! qrvap qrtend due to rain evaporation
         packed_nevapr,   &! Evaporation of total precipitation (rain + snow)
         packed_ze,       &!
         packed_mu,       &! diag_mu_c 
         packed_lambdac,  &! diag_lamc
         packed_rel,      &! effc
         packed_rei,      &! diag_effi
         packed_dei,      &! diag_deffi 
         packed_rflx,     &! rflx 
         packed_sflx,     &! sflx 
         packed_umc,      &!
         packed_umr,      &!
         packed_umi,      &!
         packed_ncal,     &!
         packed_ncai,     &!
         packed_di,       &!
         packed_rhopo,    &!
!!....................................................................
!! timescale 
!!....................................................................
         packed_epsc,      & ! 
         packed_epsi,      & ! 
         packed_epsr,      & ! 
         packed_reps,      & ! 
         packed_ac1,       & ! 
         packed_ac2,       & ! 
         packed_ac3,       & ! 
         packed_qccon1,    & ! 
         packed_qccon2,    & ! 
         packed_qicon1,    & ! 
         packed_qicon2,    & ! 
         packed_qicon3,    & ! 
!!
!! liquid process 
!!
         packed_rate1ord_cw2pr_st,    & ! 1stqc2qr
         packed_qccon,    &! qctend due to condensation of cloud droplets
         packed_qrcon,    &! qrtend due to condensation of rain 
         packed_qcaut,    &! qctend due to autoconversion to rain
         packed_ncautc,   &! nctend due to autoconversion to rain
         packed_ncautr,   &! nrtend due to autoconversion of cloud water
         packed_qcacc,    &! qctend due to accretion by rain
         packed_ncacc,    &! nctend due to accretion by rain
         packed_ncslf,    &! nctend due to cloud droplet self-collection
         packed_nrslf,    &! nrtend due to rain self-collection
         packed_ncnuc,    &! nctend due to activation of CCN
         packed_qcnuc,    &! qctend due to activation of CCN
         packed_qcevp,    &! qctend due to cloud droplet evaporation
         packed_qberg,    &! qctend/qitend due to Bergeron process 
         packed_qrevp,    &! qrtend due to rain droplet evaporation
         packed_nrevp,    &! nrtend due to rain droplet evaporation
!!
!! ice process
!!
         packed_qccol,    &
         packed_qidep,    &
         packed_qrcol,    &
         packed_qinuc,    &
         packed_nccol,    &
         packed_nrcol,    &
         packed_ninuc,    &
         packed_qisub,    &
         packed_qimlt,    &
         packed_nimlt,    &
         packed_nisub,    &
         packed_nislf,    &
         packed_qchetc,   &
         packed_qcheti,   &
         packed_qrhetc,   &
         packed_qrheti,   &
         packed_nchetc,   &
         packed_ncheti,   &
         packed_nrhetc,   &
         packed_nrheti,   &
         packed_nrshdr,   &
         packed_qcshd,    &
         packed_qrmul,    &
         packed_nimul,    &
         packed_ncshdc,   &
!!
!! sedimenation  
!!
         packed_qcsedten, &
         packed_ncsedten, &
         packed_qisedten, &
         packed_nisedten, &
         packed_qrsedten, &
         packed_nrsedten, &
!!
!! mg  
!!
         packed_cmei,     & 
!!
!! last 
!!
         errstring) 

!!         packed_cflx,            &
!!         packed_iflx,            &
!!         packed_refl,            &
!!         packed_arefl,           &
!!         packed_areflz,          &
!!         packed_frefl,           &
!!         packed_csrfl,           &
!!         packed_acsrfl,          &
!!         packed_fcsrfl,          &
!!         packed_freqr,           &
!!         packed_nfice,           &

      call t_stopf('micro_p3')

      call handle_errmsg(errstring, subname="micro_p3_tend")

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 014 -' 
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! ptend_init 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      call physics_ptend_init(ptend_loc,  &
                              psetcols,   &
                              "micro_p3", &
                              ls=.true.,  &
                              lq=lq)

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 015 -' 
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! Set local tendencies
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ptend_loc%s = &
        packer%unpack(packed_tlat, 0._r8)
        
      ptend_loc%q(:,:,1) = &
        packer%unpack(packed_qvlat, 0._r8)
        
      ptend_loc%q(:,:,ixcldliq) = &
        packer%unpack(packed_qctend, 0._r8)
        
      ptend_loc%q(:,:,ixcldice) = &
        packer%unpack(packed_qitend, 0._r8)
        
      ptend_loc%q(:,:,ixnumliq) = &
        packer%unpack(packed_nctend, -state_loc%q(:,:,ixnumliq)/(dtime/num_steps))
        
      ptend_loc%q(:,:,ixnumice) = &
        packer%unpack(packed_nitend, -state_loc%q(:,:,ixnumice)/(dtime/num_steps))
      
      ptend_loc%q(:,:,ixcldrim)  = &
        packer%unpack(packed_qirimtend, 0._r8)
      
      ptend_loc%q(:,:,ixbvrim)   = &
        packer%unpack(packed_bvrimtend, 0._r8)
      
      ptend_loc%q(:,:,ixrain)    = &
        packer%unpack(packed_qrtend, 0._r8)
      
      ptend_loc%q(:,:,ixnumrain) = &
        packer%unpack(packed_nrtend, -state_loc%q(:,:,ixnumrain)/(dtime/num_steps))

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 016 -' 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! ptend_sum
   !! physics_update
   !! sum all outputs for averaging
   !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call physics_ptend_sum(ptend_loc, ptend, ncol)

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 017 -' 

      call physics_update(state_loc, ptend_loc, dtime/num_steps)

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 018 -' 

      call post_proc%accumulate()

if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 019 -' 

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! end of subcycle 
   !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   end do
   
   
   
   call t_stopf('micro_p3_acme_tend_loop')

   call t_startf('micro_p3_acme_tend_fini')
   
   !!
   !! Divide ptend by substeps.
   !!
   
   call physics_ptend_scale(ptend, 1._r8/num_steps, ncol)

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 020 -' 
   
   !!
   !! Use summed outputs to produce averages
   !!
   
   call post_proc%process_and_unpack()

   call post_proc%finalize()

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 021 -' 
   
   if (associated(packed_frzimm)) deallocate(packed_frzimm)
   if (associated(packed_frzcnt)) deallocate(packed_frzcnt)
   if (associated(packed_frzdep)) deallocate(packed_frzdep)


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! derived fields 
   !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   ! array must be zeroed beyond trop_cloud_top_pre otherwise undefined values will be used in cosp.
   mgflxprc(:ncol,1:top_lev) = 0.0_r8
   mgflxsnw(:ncol,1:top_lev) = 0.0_r8

   mgflxprc(:ncol,top_lev:pverp) = rflx(:ncol,top_lev:pverp) + sflx(:ncol,top_lev:pverp)
   mgflxsnw(:ncol,top_lev:pverp) = sflx(:ncol,top_lev:pverp)

   !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
   !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
   ! this needs to be replaced by clubb_liq_deep and clubb_ice+deep accordingly

   cvreffliq(:ncol,top_lev:pver) = 9.0_r8
   cvreffice(:ncol,top_lev:pver) = 37.0_r8

   ! Reassign rate1 if modal aerosols
   if (rate1_cw2pr_st_idx > 0) then
      rate1ord_cw2pr_st(:ncol,top_lev:pver) = rate1cld(:ncol,top_lev:pver)
   end if

   ! Microphysical tendencies for use in the macrophysics at the next time step
   CC_T(:ncol,top_lev:pver)    =  tlat(:ncol,top_lev:pver)/cpair
   CC_qv(:ncol,top_lev:pver)   = qvlat(:ncol,top_lev:pver)
   CC_ql(:ncol,top_lev:pver)   = qcten(:ncol,top_lev:pver)
   CC_qi(:ncol,top_lev:pver)   = qiten(:ncol,top_lev:pver)
   CC_nl(:ncol,top_lev:pver)   = ncten(:ncol,top_lev:pver)
   CC_ni(:ncol,top_lev:pver)   = niten(:ncol,top_lev:pver)
   CC_qlst(:ncol,top_lev:pver) = qcten(:ncol,top_lev:pver)/max(0.01_r8,alst_mic(:ncol,top_lev:pver))

   ! Net micro_mg_cam condensation rate
   qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + cmeiout(:ncol,top_lev:pver)

   !! For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
   !! Other precip output variables are set to 0
   !! Do not subscript by ncol here, because in physpkg we divide the whole
   !! array and need to avoid an FPE due to uninitialized data.
   
   prec_pcw = prect
   prec_sed = 0._r8
   prec_str = prec_pcw + prec_sed
   
   snow_pcw = preci
   snow_sed = 0._r8
   snow_str = snow_pcw + snow_sed


   icecldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)
   liqcldf(:ncol,top_lev:pver) = ast(:ncol,top_lev:pver)

   ! ------------------------------------------------------------ !
   ! Compute in cloud ice and liquid mixing ratios                !
   ! Note that 'iclwp, iciwp' are used for radiation computation. !
   ! ------------------------------------------------------------ !

   icinc = 0._r8
   icwnc = 0._r8
   iciwpst = 0._r8
   iclwpst = 0._r8

   do k = top_lev, pver
      do i = 1, ncol
         ! Limits for in-cloud mixing ratios consistent with P3 microphysics
         ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
         icimrst(i,k)   = min( state_loc%q(i,k,ixcldice) / max(mincld,icecldf(i,k)),0.005_r8 )
         icwmrst(i,k)   = min( state_loc%q(i,k,ixcldliq) / max(mincld,liqcldf(i,k)),0.005_r8 )
         icinc(i,k)     = state_loc%q(i,k,ixnumice) / max(mincld,icecldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         icwnc(i,k)     = state_loc%q(i,k,ixnumliq) / max(mincld,liqcldf(i,k)) * &
              state_loc%pmid(i,k) / (287.15_r8*state_loc%t(i,k))
         ! Calculate micro_p3_acme cloud water paths in each layer
         ! Note: uses stratiform cloud fraction!
         iciwpst(i,k)   = min(state_loc%q(i,k,ixcldice)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
         iclwpst(i,k)   = min(state_loc%q(i,k,ixcldliq)/max(mincld,ast(i,k)),0.005_r8) * state_loc%pdel(i,k) / gravit
      end do
   end do

   ! Calculate cloud fraction for prognostic precip sizes.

      ! Cloud fraction for purposes of precipitation is maximum cloud
      ! fraction out of all the layers that the precipitation may be
      ! falling down from.
      cldmax(:ncol,top_lev:pver) = max(mincld, ast(:ncol,top_lev:pver))
      do k = top_lev+1, pver
         where (state_loc%q(:ncol,k-1,ixrain) >= qsmall .or. &
              state_loc%q(:ncol,k-1,ixcldrim) >= qsmall)
            cldmax(:ncol,k) = max(cldmax(:ncol,k-1), cldmax(:ncol,k))
         end where
      end do


   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 022 -' 
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! average the fields on the grid
   !! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   
   if (use_subcol_microp) then

   !!.......................................................... 
   !!  newly added 
   !!.......................................................... 
   

      call subcol_field_avg(epsc,              ngrdcol, lchnk, epsc_grid)
      call subcol_field_avg(epsi,              ngrdcol, lchnk, epsi_grid)
      call subcol_field_avg(epsr,              ngrdcol, lchnk, epsr_grid)
      call subcol_field_avg(reps,              ngrdcol, lchnk, reps_grid)
      call subcol_field_avg(ac1,               ngrdcol, lchnk, ac1_grid)
      call subcol_field_avg(ac2,               ngrdcol, lchnk, ac2_grid)
      call subcol_field_avg(ac3,               ngrdcol, lchnk, ac3_grid)
      call subcol_field_avg(qccon1,            ngrdcol, lchnk, qccon1_grid)
      call subcol_field_avg(qccon2,            ngrdcol, lchnk, qccon2_grid)
      call subcol_field_avg(qicon1,            ngrdcol, lchnk, qicon1_grid)
      call subcol_field_avg(qicon2,            ngrdcol, lchnk, qicon2_grid)
      call subcol_field_avg(qicon3,            ngrdcol, lchnk, qicon3_grid)
      
      !! liquid 
      
      call subcol_field_avg(qcnuc,             ngrdcol, lchnk, qcnuc_grid)
      call subcol_field_avg(ncnuc,             ngrdcol, lchnk, ncnuc_grid)
      call subcol_field_avg(qccon,             ngrdcol, lchnk, qccon_grid)
      call subcol_field_avg(qrcon,             ngrdcol, lchnk, qrcon_grid)
      call subcol_field_avg(qcaut,             ngrdcol, lchnk, qcaut_grid)
      call subcol_field_avg(ncautc,            ngrdcol, lchnk, ncautc_grid)
      call subcol_field_avg(ncautr,            ngrdcol, lchnk, ncautr_grid)
      call subcol_field_avg(qcacc,             ngrdcol, lchnk, qcacc_grid)
      call subcol_field_avg(ncacc,             ngrdcol, lchnk, ncacc_grid)
      call subcol_field_avg(ncslf,             ngrdcol, lchnk, ncslf_grid)
      call subcol_field_avg(nrslf,             ngrdcol, lchnk, nrslf_grid)
      call subcol_field_avg(prer_evap,         ngrdcol, lchnk, prer_evap_grid)
      call subcol_field_avg(qcevp,             ngrdcol, lchnk, qcevp_grid)
      call subcol_field_avg(qberg,             ngrdcol, lchnk, qberg_grid)
      call subcol_field_avg(qrevp,             ngrdcol, lchnk, qrevp_grid)
      call subcol_field_avg(nrevp,             ngrdcol, lchnk, nrevp_grid)

      !! ice 
      call subcol_field_avg(qccol,             ngrdcol, lchnk, qccol_grid)
      call subcol_field_avg(qidep,             ngrdcol, lchnk, qidep_grid)
      call subcol_field_avg(qrcol,             ngrdcol, lchnk, qrcol_grid)
      call subcol_field_avg(qinuc,             ngrdcol, lchnk, qinuc_grid)
      call subcol_field_avg(nccol,             ngrdcol, lchnk, nccol_grid)
      call subcol_field_avg(nrcol,             ngrdcol, lchnk, nrcol_grid)
      call subcol_field_avg(ninuc,             ngrdcol, lchnk, ninuc_grid)
      call subcol_field_avg(qisub,             ngrdcol, lchnk, qisub_grid)
      call subcol_field_avg(qimlt,             ngrdcol, lchnk, qimlt_grid)
      call subcol_field_avg(nimlt,             ngrdcol, lchnk, nimlt_grid)
      call subcol_field_avg(nisub,             ngrdcol, lchnk, nisub_grid)
      call subcol_field_avg(nislf,             ngrdcol, lchnk, nislf_grid)
      call subcol_field_avg(qchetc,            ngrdcol, lchnk, qchetc_grid)
      call subcol_field_avg(qcheti,            ngrdcol, lchnk, qcheti_grid)
      call subcol_field_avg(qrhetc,            ngrdcol, lchnk, qrhetc_grid)
      call subcol_field_avg(qrheti,            ngrdcol, lchnk, qrheti_grid)
      call subcol_field_avg(nchetc,            ngrdcol, lchnk, nchetc_grid)
      call subcol_field_avg(ncheti,            ngrdcol, lchnk, ncheti_grid)
      call subcol_field_avg(nrhetc,            ngrdcol, lchnk, nrhetc_grid)
      call subcol_field_avg(nrheti,            ngrdcol, lchnk, nrheti_grid)
      call subcol_field_avg(nrshdr,            ngrdcol, lchnk, nrshdr_grid)
      call subcol_field_avg(qcshd,             ngrdcol, lchnk, qcshd_grid)
      call subcol_field_avg(qrmul,             ngrdcol, lchnk, qrmul_grid)
      call subcol_field_avg(nimul,             ngrdcol, lchnk, nimul_grid)
      call subcol_field_avg(ncshdc,            ngrdcol, lchnk, ncshdc_grid)

      call subcol_field_avg(cmeiout,           ngrdcol, lchnk, cmeiout_grid)

   !!.......................................................... 
   !!  others 
   !!.......................................................... 
   
      call subcol_field_avg(prec_str,  ngrdcol, lchnk, prec_str_grid)
      call subcol_field_avg(iclwpst,   ngrdcol, lchnk, iclwpst_grid)
      call subcol_field_avg(cvreffliq, ngrdcol, lchnk, cvreffliq_grid)
      call subcol_field_avg(cvreffice, ngrdcol, lchnk, cvreffice_grid)
      call subcol_field_avg(mgflxprc,  ngrdcol, lchnk, mgflxprc_grid)
      call subcol_field_avg(mgflxsnw,  ngrdcol, lchnk, mgflxsnw_grid)
      call subcol_field_avg(qme,       ngrdcol, lchnk, qme_grid)
      call subcol_field_avg(nevapr,    ngrdcol, lchnk, nevapr_grid)
      call subcol_field_avg(prain,     ngrdcol, lchnk, prain_grid)

      !! average fields which are not in pbuf

      call subcol_field_avg(cld,       ngrdcol, lchnk, cld_grid)
      call subcol_field_avg(icwmrst,   ngrdcol, lchnk, icwmrst_grid)
      call subcol_field_avg(icimrst,   ngrdcol, lchnk, icimrst_grid)
      call subcol_field_avg(liqcldf,   ngrdcol, lchnk, liqcldf_grid)
      call subcol_field_avg(icecldf,   ngrdcol, lchnk, icecldf_grid)
      call subcol_field_avg(icwnc,     ngrdcol, lchnk, icwnc_grid)
      call subcol_field_avg(icinc,     ngrdcol, lchnk, icinc_grid)
      call subcol_field_avg(state_loc%pdel, ngrdcol, lchnk, pdel_grid)

      call subcol_field_avg(state_loc%q(:,:,ixnumliq), ngrdcol, lchnk, nc_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumice), ngrdcol, lchnk, ni_grid)
      call subcol_field_avg(state_loc%q(:,:,ixrain),   ngrdcol, lchnk, qr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixnumrain),ngrdcol, lchnk, nr_grid)
      call subcol_field_avg(state_loc%q(:,:,ixcldrim), ngrdcol, lchnk, qirim_grid)
      call subcol_field_avg(state_loc%q(:,:,ixbvrim),  ngrdcol, lchnk, bvrim_grid)


      call subcol_field_avg(cldmax,    ngrdcol, lchnk, cldmax_grid)

   
   else
      ! These pbuf fields need to be assigned.  There is no corresponding subcol_field_avg
      ! as they are reset before being used, so it would be a needless calculation
      lambdac_grid    => lambdac
      mu_grid         => mu
      rel_grid        => rel
      rei_grid        => rei
      dei_grid        => dei 

      ! fields already on grids, so just assign
      prec_str_grid   => prec_str
      iclwpst_grid    => iclwpst
      cvreffliq_grid  => cvreffliq
      cvreffice_grid  => cvreffice
      mgflxprc_grid   => mgflxprc
      mgflxsnw_grid   => mgflxsnw
      qme_grid        => qme
      nevapr_grid     => nevapr
      prain_grid      => prain
 
      cld_grid        = cld
      
      icwmrst_grid    = icwmrst
      icimrst_grid    = icimrst
      liqcldf_grid    = liqcldf
      icecldf_grid    = icecldf
      icwnc_grid      = icwnc
      icinc_grid      = icinc
      pdel_grid       = state_loc%pdel
      prao_grid       = prao
      prco_grid       = prco

      nc_grid = state_loc%q(:,:,ixnumliq)
      ni_grid = state_loc%q(:,:,ixnumice)
      qr_grid = state_loc%q(:,:,ixrain)
      nr_grid = state_loc%q(:,:,ixnumrain)
      qirim_grid = state_loc%q(:,:,ixcldrim)
      bvrim_grid = state_loc%q(:,:,ixbvrim)

      cldmax_grid = cldmax

     
      epsc_grid              = epsc
      epsi_grid              = epsi
      epsr_grid              = epsr
      reps_grid              = reps
      ac1_grid               = ac1
      ac2_grid               = ac2
      ac3_grid               = ac3
      qccon1_grid              = qccon1
      qccon2_grid              = qccon2
      qicon1_grid              = qicon1
      qicon2_grid              = qicon2
      qicon3_grid              = qicon3
      
      !! liquid 
      
      qcnuc_grid             = qcnuc     
      ncnuc_grid             = ncnuc     
      qccon_grid             = qccon
      qrcon_grid             = qrcon
      qcaut_grid             = qcaut
      ncautc_grid            = ncautc
      ncautr_grid            = ncautr            
      qcacc_grid             = qcacc            
      ncacc_grid             = ncacc
      ncslf_grid             = ncslf
      nrslf_grid             = nrslf
      prer_evap_grid         = prer_evap
      qcevp_grid             = qcevp
      qberg_grid             = qberg
      qrevp_grid             = qrevp
      nrevp_grid             = nrevp

      !! ice 
      qccol_grid             = qccol     
      qidep_grid             = qidep     
      qrcol_grid             = qrcol     
      qinuc_grid             = qinuc     
      nccol_grid             = nccol     
      nrcol_grid             = nrcol     
      ninuc_grid             = ninuc     
      qisub_grid             = qisub     
      qimlt_grid             = qimlt     
      nimlt_grid             = nimlt     
      nisub_grid             = nisub     
      nislf_grid             = nislf     
      qchetc_grid            = qchetc     
      qcheti_grid            = qcheti     
      qrhetc_grid            = qrhetc     
      qrheti_grid            = qrheti     
      nchetc_grid            = nchetc     
      ncheti_grid            = ncheti     
      nrhetc_grid            = nrhetc     
      nrheti_grid            = nrheti     
      nrshdr_grid            = nrshdr     
      qcshd_grid             = qcshd     
      qrmul_grid             = qrmul     
      nimul_grid             = nimul     
      ncshdc_grid            = ncshdc
      
      cmeiout_grid           = cmeiout

   end if

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 023 -' 
   
   !! If on subcolumns, average the rest of the pbuf fields which were modified on 
   !! subcolumns but are not used further in this parameterization  
   !! (no need to assign in the non-subcolumn case -- the else step)
   
   if (use_subcol_microp) then
      call subcol_field_avg(snow_str,    ngrdcol, lchnk, snow_str_grid)
      call subcol_field_avg(prec_pcw,    ngrdcol, lchnk, prec_pcw_grid)
      call subcol_field_avg(snow_pcw,    ngrdcol, lchnk, snow_pcw_grid)
      call subcol_field_avg(prec_sed,    ngrdcol, lchnk, prec_sed_grid)
      call subcol_field_avg(snow_sed,    ngrdcol, lchnk, snow_sed_grid)
      call subcol_field_avg(cldo,        ngrdcol, lchnk, cldo_grid)
      call subcol_field_avg(cc_t,        ngrdcol, lchnk, cc_t_grid)
      call subcol_field_avg(cc_qv,       ngrdcol, lchnk, cc_qv_grid)
      call subcol_field_avg(cc_ql,       ngrdcol, lchnk, cc_ql_grid)
      call subcol_field_avg(cc_qi,       ngrdcol, lchnk, cc_qi_grid)
      call subcol_field_avg(cc_nl,       ngrdcol, lchnk, cc_nl_grid)
      call subcol_field_avg(cc_ni,       ngrdcol, lchnk, cc_ni_grid)
      call subcol_field_avg(cc_qlst,     ngrdcol, lchnk, cc_qlst_grid)
      call subcol_field_avg(iciwpst,     ngrdcol, lchnk, iciwpst_grid)

      if (rate1_cw2pr_st_idx > 0) then
         call subcol_field_avg(rate1ord_cw2pr_st, ngrdcol, lchnk, rate1ord_cw2pr_st_grid)
      end if

   end if
      
   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 024 -' 
   
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! derived fields 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

   ! Calculate rho (on subcolumns if turned on) for size distribution
   ! parameter calculations and average it if needed
   !
   ! State instead of state_loc to preserve answers for P31 (and in any
   ! case, it is unlikely to make much difference).
   
   rho(:ncol,top_lev:) = &
      state%pmid(:ncol,top_lev:) / (rair*state%t(:ncol,top_lev:))

   if (use_subcol_microp) then
      call subcol_field_avg(rho, ngrdcol, lchnk, rho_grid)
   else
      rho_grid = rho
   end if

   !!
   !! Effective radius for cloud liquid, fixed number.
   !!
   
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_fn_grid = 10._r8

   ncic_grid = 1.e8_r8

   !! size distribution 
   
   call size_dist_param_liq( &
                             mg_liq_props, &
                             icwmrst_grid(:ngrdcol,top_lev:), &
                             ncic_grid   (:ngrdcol,top_lev:), &
                             rho_grid    (:ngrdcol,top_lev:), &
                             mu_grid     (:ngrdcol,top_lev:), &
                             lambdac_grid(:ngrdcol,top_lev:))

   where (icwmrst_grid(:ngrdcol,top_lev:) > qsmall)
      rel_fn_grid(:ngrdcol,top_lev:) = &
                    (mu_grid(:ngrdcol,top_lev:) + 3._r8)/ &
                    lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   end where

   !!
   !! Effective radius for cloud liquid, and size parameters
   !! mu_grid and lambdac_grid.
   !!
   
   mu_grid = 0._r8
   lambdac_grid = 0._r8
   rel_grid = 10._r8

   !!
   !! Calculate ncic on the grid
   !!
   
   ncic_grid(:ngrdcol,top_lev:) = nc_grid(:ngrdcol,top_lev:) / &
        max(mincld,liqcldf_grid(:ngrdcol,top_lev:))

   call size_dist_param_liq(&
           mg_liq_props, &
           icwmrst_grid(:ngrdcol,top_lev:), & 
           ncic_grid(:ngrdcol,top_lev:), &
           rho_grid(:ngrdcol,top_lev:), &
           mu_grid(:ngrdcol,top_lev:), &
           lambdac_grid(:ngrdcol,top_lev:))

   where (icwmrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rel_grid(:ngrdcol,top_lev:) = &
           (mu_grid(:ngrdcol,top_lev:) + 3._r8) / &
           lambdac_grid(:ngrdcol,top_lev:)/2._r8 * 1.e6_r8
   elsewhere
      ! Deal with the fact that size_dist_param_liq sets mu_grid to -100
      ! wherever there is no cloud.
      mu_grid(:ngrdcol,top_lev:) = 0._r8
   end where

   !!
   !! Rain/Snow effective diameter
   !!
   
   drout2_grid = 0._r8
   reff_rain_grid = 0._r8

      ! Prognostic precipitation

      where (qr_grid(:ngrdcol,top_lev:) >= 1.e-7_r8)
         drout2_grid(:ngrdcol,top_lev:) = avg_diameter( &
              qr_grid(:ngrdcol,top_lev:), &
              nr_grid(:ngrdcol,top_lev:) * rho_grid(:ngrdcol,top_lev:), &
              rho_grid(:ngrdcol,top_lev:), rhow)

         reff_rain_grid(:ngrdcol,top_lev:) = drout2_grid(:ngrdcol,top_lev:) * &
              1.5_r8 * 1.e6_r8
      end where


   !!
   !! Effective radius and diameter for cloud ice
   !!
   
   rei_grid = 25._r8

   niic_grid(:ngrdcol,top_lev:) = ni_grid(:ngrdcol,top_lev:) / &
        max(mincld,icecldf_grid(:ngrdcol,top_lev:))

   call size_dist_param_basic( &
           mg_ice_props, &
           icimrst_grid(:ngrdcol,top_lev:), &
           niic_grid(:ngrdcol,top_lev:), &
           rei_grid(:ngrdcol,top_lev:))

   where (icimrst_grid(:ngrdcol,top_lev:) >= qsmall)
      rei_grid(:ngrdcol,top_lev:) = &
         1.5_r8/rei_grid(:ngrdcol,top_lev:) * 1.e6_r8
   elsewhere
      rei_grid(:ngrdcol,top_lev:) = 25._r8
   end where

   dei_grid = rei_grid * rhoi/rhows * 2._r8

   !!
   !! Limiters for low cloud fraction
   !!
   
   do k = top_lev, pver
      do i = 1, ngrdcol
         if ( ast_grid(i,k) < 1.e-4_r8 ) then
            mu_grid(i,k) = mucon
            lambdac_grid(i,k) = (mucon + 1._r8)/dcon
            dei_grid(i,k) = deicon
         end if
      end do
   end do

   mgreffrain_grid(:ngrdcol,top_lev:pver) = reff_rain_grid(:ngrdcol,top_lev:pver)
   mgreffsnow_grid(:ngrdcol,top_lev:pver) = 1000._r8 !! dummy value 

   ! ------------------------------------- !
   ! Precipitation efficiency Calculation  !
   ! ------------------------------------- !

   !-----------------------------------------------------------------------
   ! Liquid water path

   ! Compute liquid water paths, and column condensation
   tgliqwp_grid(:ngrdcol) = 0._r8
   tgcmeliq_grid(:ngrdcol) = 0._r8
   do k = top_lev, pver
      do i = 1, ngrdcol
         tgliqwp_grid(i)  = tgliqwp_grid(i) + iclwpst_grid(i,k)*cld_grid(i,k)

         if (cmeliq_grid(i,k) > 1.e-12_r8) then
            !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
            tgcmeliq_grid(i) = tgcmeliq_grid(i) + cmeliq_grid(i,k) * &
                 (pdel_grid(i,k) / gravit) / rhoh2o
         end if
      end do
   end do

   ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
   ! this is 1ppmv of h2o in 10hpa
   ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

   !-----------------------------------------------------------------------
   ! precipitation efficiency calculation  (accumulate cme and precip)

   minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

   ! zero out precip efficiency and total averaged precip
   pe_grid(:ngrdcol)     = 0._r8
   tpr_grid(:ngrdcol)    = 0._r8
   pefrac_grid(:ngrdcol) = 0._r8

   ! accumulate precip and condensation
   do i = 1, ngrdcol

      acgcme_grid(i)  = acgcme_grid(i) + tgcmeliq_grid(i)
      acprecl_grid(i) = acprecl_grid(i) + prec_str_grid(i)
      acnum_grid(i)   = acnum_grid(i) + 1

      ! if LWP is zero, then 'end of cloud': calculate precip efficiency
      if (tgliqwp_grid(i) < minlwp) then
         if (acprecl_grid(i) > 5.e-8_r8) then
            tpr_grid(i) = max(acprecl_grid(i)/acnum_grid(i), 1.e-15_r8)
            if (acgcme_grid(i) > 1.e-10_r8) then
               pe_grid(i) = min(max(acprecl_grid(i)/acgcme_grid(i), 1.e-15_r8), 1.e5_r8)
               pefrac_grid(i) = 1._r8
            end if
         end if

         ! reset counters
!        if (pe_grid(i) /= 0._r8 .and. (pe_grid(i) < 1.e-8_r8 .or. pe_grid(i) > 1.e3_r8)) then
!           write (iulog,*) 'PE_grid:ANOMALY  pe_grid, acprecl_grid, acgcme_grid, tpr_grid, acnum_grid ', &
!                           pe_grid(i),acprecl_grid(i), acgcme_grid(i), tpr_grid(i), acnum_grid(i)
!        endif

         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
         acnum_grid(i)   = 0
      end if               ! end LWP zero conditional

      ! if never find any rain....(after 10^3 timesteps...)
      if (acnum_grid(i) > 1000) then
         acnum_grid(i)   = 0
         acprecl_grid(i) = 0._r8
         acgcme_grid(i)  = 0._r8
      end if

   end do

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 025 -' 
   
   ! --------------------- !
   ! History Output Fields !
   ! --------------------- !

   ! Column droplet concentration
   
   cdnumc_grid(:ngrdcol) = &
       sum( nc_grid(:ngrdcol,top_lev:pver) * &
            pdel_grid(:ngrdcol,top_lev:pver) / gravit, &
            dim=2)


   !! prepare following output 
   !! 
   !! AREL 
   !! AREI 
   !! AWNC 
   !! AWNI 
   !! FREQL 
   !! FREQI
   !!  

   ! Averaging for new output fields
   efcout_grid      = 0._r8
   efiout_grid      = 0._r8
   ncout_grid       = 0._r8
   niout_grid       = 0._r8
   freql_grid       = 0._r8
   freqi_grid       = 0._r8
   icwmrst_grid_out = 0._r8
   icimrst_grid_out = 0._r8

   do k = top_lev, pver
      do i = 1, ngrdcol
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 5.e-5_r8 ) then
            efcout_grid(i,k) = rel_grid(i,k) * liqcldf_grid(i,k)
            ncout_grid(i,k)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            freql_grid(i,k)  = liqcldf_grid(i,k)
            icwmrst_grid_out(i,k) = icwmrst_grid(i,k)
         end if
         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-6_r8 ) then
            efiout_grid(i,k) = rei_grid(i,k) * icecldf_grid(i,k)
            niout_grid(i,k)  = icinc_grid(i,k) * icecldf_grid(i,k)
            freqi_grid(i,k)  = icecldf_grid(i,k)
            icimrst_grid_out(i,k) = icimrst_grid(i,k)
         end if
      end do
   end do

   !!.......................................................... 
   !! cloud top values 
   !!.......................................................... 
   
   fcti_grid  = 0._r8
   fctl_grid  = 0._r8
   ctrel_grid = 0._r8
   ctrei_grid = 0._r8
   ctnl_grid  = 0._r8
   ctni_grid  = 0._r8
   
   do i = 1, ngrdcol
      do k = top_lev, pver
         if ( liqcldf_grid(i,k) > 0.01_r8 .and. icwmrst_grid(i,k) > 1.e-7_r8 ) then
            ctrel_grid(i) = rel_grid(i,k) * liqcldf_grid(i,k)
            ctnl_grid(i)  = icwnc_grid(i,k) * liqcldf_grid(i,k)
            fctl_grid(i)  = liqcldf_grid(i,k)
            exit
         end if
         if ( icecldf_grid(i,k) > 0.01_r8 .and. icimrst_grid(i,k) > 1.e-7_r8 ) then
            ctrei_grid(i) = rei_grid(i,k) * icecldf_grid(i,k)
            ctni_grid(i)  = icinc_grid(i,k) * icecldf_grid(i,k)
            fcti_grid(i)  = icecldf_grid(i,k)
            exit
         end if
      end do
   end do

   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 026 -' 
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! output 
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !!.......................................................... 
   !! averaged already 
   !!.......................................................... 
   

   call outfld('P3_EPSC',   epsc_grid,  pcols, lchnk)
   call outfld('P3_EPSI',   epsi_grid,  pcols, lchnk)
   call outfld('P3_EPSR',   epsr_grid,  pcols, lchnk)
   call outfld('P3_REPS',   reps_grid,  pcols, lchnk)
   call outfld('P3_AC1',    ac1_grid,   pcols, lchnk)
   call outfld('P3_AC2',    ac2_grid,   pcols, lchnk)
   call outfld('P3_AC3',    ac3_grid,   pcols, lchnk)
   call outfld('P3_QCCON1',   qccon1_grid,  pcols, lchnk)
   call outfld('P3_QCCON2',   qccon2_grid,  pcols, lchnk)
   call outfld('P3_QICON1',   qicon1_grid,  pcols, lchnk)
   call outfld('P3_QICON2',   qicon2_grid,  pcols, lchnk)
   call outfld('P3_QICON3',   qicon3_grid,  pcols, lchnk)
   
   !! liquid processes 
   
   call outfld('P3_QCAUT',   qcaut_grid,  pcols, lchnk)
   call outfld('P3_NCAUTC',  ncautc_grid, pcols, lchnk)
   call outfld('P3_NCAUTR',  ncautr_grid, pcols, lchnk)
   call outfld('P3_QCACC',   qcacc_grid,  pcols, lchnk)
   call outfld('P3_NCACC',   ncacc_grid,  pcols, lchnk)
   call outfld('P3_NCSLF',   ncslf_grid,  pcols, lchnk)
   call outfld('P3_NRSLF',   nrslf_grid,  pcols, lchnk)
   call outfld('P3_QCNUC',   qcnuc_grid,  pcols, lchnk)
   call outfld('P3_NCNUC',   ncnuc_grid,  pcols, lchnk)
   call outfld('P3_QREVP',   qrevp_grid,  pcols, lchnk)
   call outfld('P3_QCEVP',   qcevp_grid,  pcols, lchnk)
   call outfld('P3_QBERG',   qberg_grid,  pcols, lchnk)
   call outfld('P3_QRCON',   qrcon_grid,  pcols, lchnk)
   call outfld('P3_QCCON',   qccon_grid,  pcols, lchnk)
   
   !! ice processes 
   
   call outfld('P3_QCCOL',   qccol_grid,  pcols, lchnk)
   call outfld('P3_QIDEP',   qidep_grid,  pcols, lchnk)
   call outfld('P3_QRCOL',   qrcol_grid,  pcols, lchnk)
   call outfld('P3_NCCOL',   nccol_grid,  pcols, lchnk)
   call outfld('P3_NRCOL',   nrcol_grid,  pcols, lchnk)
   call outfld('P3_QINUC',   qinuc_grid,  pcols, lchnk)
   call outfld('P3_NINUC',   ninuc_grid,  pcols, lchnk)
   call outfld('P3_QISUB',   qisub_grid,  pcols, lchnk)
   call outfld('P3_QIMLT',   qimlt_grid,  pcols, lchnk)
   call outfld('P3_NIMLT',   nimlt_grid,  pcols, lchnk)
   call outfld('P3_NISUB',   nisub_grid,  pcols, lchnk)
   call outfld('P3_NISLF',   nislf_grid,  pcols, lchnk)
   call outfld('P3_QCHETI',  qcheti_grid, pcols, lchnk)
   call outfld('P3_NCHETI',  ncheti_grid, pcols, lchnk)
   call outfld('P3_QCHETC',  qchetc_grid, pcols, lchnk)
   call outfld('P3_NCHETC',  nchetc_grid, pcols, lchnk)
   call outfld('P3_QRHETI',  qrheti_grid, pcols, lchnk)
   call outfld('P3_NRHETI',  nrheti_grid, pcols, lchnk)
   call outfld('P3_QRHETC',  qrhetc_grid, pcols, lchnk)
   call outfld('P3_NRHETC',  nrhetc_grid, pcols, lchnk)
   call outfld('P3_QCSHD',   qcshd_grid,  pcols, lchnk)
   call outfld('P3_NCSHDC',  ncshdc_grid, pcols, lchnk)
   call outfld('P3_NRSHDR',  nrshdr_grid, pcols, lchnk)
   call outfld('P3_QRMUL',   qrmul_grid,  pcols, lchnk)
   call outfld('P3_NIMUL',   nimul_grid,  pcols, lchnk)
   
   call outfld('P3_CMEIOUT', cmeiout_grid,  pcols, lchnk)
   call outfld('P3_ICWNC',   icwnc_grid,    pcols, lchnk)
   call outfld('P3_ICINC',   icinc_grid,    pcols, lchnk)
   call outfld('P3_CDNUMC',  cdnumc_grid,   pcols, lchnk)
   call outfld('P3_REL',     rel_grid,         pcols, lchnk)
   call outfld('P3_REI',     rei_grid,         pcols, lchnk)
   call outfld('P3_ICIMRST', icimrst_grid_out, pcols, lchnk)
   call outfld('P3_ICWMRST', icwmrst_grid_out, pcols, lchnk) 
   call outfld('P3_CME',     qme_grid,         pcols, lchnk)
   
   call outfld('P3_ADRAIN',  drout2_grid,      pcols, lchnk)
   
   call outfld('P3_AREL',        efcout_grid,      pcols, lchnk)
   call outfld('P3_AREI',        efiout_grid,      pcols, lchnk)
   call outfld('P3_AWNC' ,       ncout_grid,       pcols, lchnk)
   call outfld('P3_AWNI' ,       niout_grid,       pcols, lchnk)
   call outfld('P3_FREQL',       freql_grid,       pcols, lchnk)
   call outfld('P3_FREQI',       freqi_grid,       pcols, lchnk)
   
   !!.......................................................... 
   !! not averaged already
   !!.......................................................... 
   
   call outfld('P3_UMC', umc, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_UMI', umi, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_UMR', umr, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)

   call outfld('P3_QCSED', qcsedten, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_QISED', qisedten, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_QRSED', qrsedten, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)

   call outfld('P3_NCAL', ncal, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_NCAI', ncai, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)

   call outfld('P3_MPICLWPI', iclwpi, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)
   call outfld('P3_MPICIWPI', iciwpi, psetcols, lchnk, &
                              avg_subcol_field=use_subcol_microp)

!!!   ! Output a handle of variables which are calculated on the fly
!!!   ftem_grid = 0._r8
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) =  qcreso_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDW2V', ftem_grid, pcols, lchnk)
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) =  melto_grid(:ngrdcol,top_lev:pver) - mnuccco_grid(:ngrdcol,top_lev:pver)&
!!!        - mnuccto_grid(:ngrdcol,top_lev:pver) -  bergo_grid(:ngrdcol,top_lev:pver) - homoo_grid(:ngrdcol,top_lev:pver)&
!!!        - msacwio_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDW2I', ftem_grid, pcols, lchnk)
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) = -prao_grid(:ngrdcol,top_lev:pver) - prco_grid(:ngrdcol,top_lev:pver)&
!!!        - psacwso_grid(:ngrdcol,top_lev:pver) - bergso_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDW2P', ftem_grid, pcols, lchnk)
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) =  cmeiout_grid(:ngrdcol,top_lev:pver) + qireso_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDI2V', ftem_grid, pcols, lchnk)
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) = -melto_grid(:ngrdcol,top_lev:pver) + mnuccco_grid(:ngrdcol,top_lev:pver) &
!!!        + mnuccto_grid(:ngrdcol,top_lev:pver) +  bergo_grid(:ngrdcol,top_lev:pver) + homoo_grid(:ngrdcol,top_lev:pver)&
!!!        + msacwio_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDI2W', ftem_grid, pcols, lchnk)
!!!
!!!   ftem_grid(:ngrdcol,top_lev:pver) = - praio_grid(:ngrdcol,top_lev:pver)
!!!   call outfld( 'P3_MPDI2P', ftem_grid, pcols, lchnk)
!!!
!!!   ! Output fields which have not been averaged already, averaging if use_subcol_microp is true
!!!   call outfld('P3_REFL',        refl,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_AREFL',       arefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_AREFLZ',      areflz,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_FREFL',       frefl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_CSRFL',       csrfl,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_ACSRFL',      acsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_FCSRFL',      fcsrfl,      psetcols, lchnk, avg_subcol_field=use_subcol_microp)

!!!   call outfld('P3_FREQR',       freqr,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_MPDT',        tlat,        psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_MPDQ',        qvlat,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_MPDLIQ',      qcten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_MPDICE',      qiten,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!   call outfld('P3_FICE',        nfice,       psetcols, lchnk, avg_subcol_field=use_subcol_microp)
!!!
!!!   ! Example subcolumn outfld call
!!!   if (use_subcol_microp) then
!!!      call outfld('P3_FICE_SCOL',   nfice,       psubcols*pcols, lchnk)
!!!   end if

!!!   call outfld('P3_CV_REFFLIQ',  cvreffliq_grid,   pcols, lchnk)
!!!   call outfld('P3_CV_REFFICE',  cvreffice_grid,   pcols, lchnk)
!!!   call outfld('P3_LS_FLXPRC',   mgflxprc_grid,    pcols, lchnk)
!!!   call outfld('P3_LS_REFFRAIN', mgreffrain_grid,  pcols, lchnk)
!!!   call outfld('P3_LS_REFFSNOW', mgreffsnow_grid,  pcols, lchnk)
!!!
!!!   call outfld('P3_PE',          pe_grid,          pcols, lchnk)
!!!   call outfld('P3_PEFRAC',      pefrac_grid,      pcols, lchnk)
!!!   call outfld('P3_APRL',        tpr_grid,         pcols, lchnk)

!!!   call outfld('P3_AREL',        efcout_grid,      pcols, lchnk)
!!!   call outfld('P3_AREI',        efiout_grid,      pcols, lchnk)
!!!   call outfld('P3_AWNC' ,       ncout_grid,       pcols, lchnk)
!!!   call outfld('P3_AWNI' ,       niout_grid,       pcols, lchnk)
!!!   call outfld('P3_FREQL',       freql_grid,       pcols, lchnk)
!!!   call outfld('P3_FREQI',       freqi_grid,       pcols, lchnk)
!!!   call outfld('P3_ACTREL',      ctrel_grid,       pcols, lchnk)
!!!   call outfld('P3_ACTREI',      ctrei_grid,       pcols, lchnk)
!!!   call outfld('P3_ACTNL',       ctnl_grid,        pcols, lchnk)
!!!   call outfld('P3_ACTNI',       ctni_grid,        pcols, lchnk)
!!!   call outfld('P3_FCTL',        fctl_grid,        pcols, lchnk)
!!!   call outfld('P3_FCTI',        fcti_grid,        pcols, lchnk)
!!!   call outfld('P3_EFFLIQ_IND',  rel_fn_grid,      pcols, lchnk)


   ! ptend_loc is deallocated in physics_update above
   
   call physics_state_dealloc(state_loc)
   
   if(l_summary_debug) write(6,*) 'micro_p3_acme_tend - 027 -' 
   
   call t_stopf('micro_p3_acme_tend_fini')
   
   
   call t_stopf('micro_p3_acme_tend')
   
end subroutine micro_p3_acme_tend



function p1(tin) result(pout)
  real(r8), target, intent(in) :: tin(:)
  real(r8), pointer :: pout(:)
  pout => tin
end function p1



function p2(tin) result(pout)
  real(r8), target, intent(in) :: tin(:,:)
  real(r8), pointer :: pout(:,:)
  pout => tin
end function p2

end module micro_p3_acme



