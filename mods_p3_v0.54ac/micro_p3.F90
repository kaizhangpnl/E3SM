MODULE micro_p3
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Revised P3 microphysics scheme for E3SM 
!!
!! Original P3 code from Hugh Morrison 
!!
!! New development by Kai Zhang
!!
!! New Features 
!!   . Separate in-cloud and grid-mean calculations 
!!   . SILHS 
!!   . Coupled to new ice nucleation schemes 
!!   . Alternative process coupling methods (under development)  
!! 
!! Current Version: 0.54ac
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef HAVE_GAMMA_INTRINSICS
 USE shr_spfn_mod, ONLY: gamma => shr_spfn_gamma
#endif

 USE wv_sat_methods,    ONLY: qsat_water => wv_sat_qsat_water, &
                              qsat_ice => wv_sat_qsat_ice
 USE micro_p3_utils,    ONLY: omsm,   &
                              mincld, &
                              rhosn,  &
                              rhows,  &
                              ac, bc, &
                              ai, bi, &
                              as, bs
 USE cam_logfile,       ONLY: iulog
 USE cam_abortutils,    ONLY: endrun
 USE error_messages,    ONLY: handle_errmsg
 USE spmd_utils,        ONLY: masterproc
 USE cam_debug,         ONLY: l_debug, l_masscon_debug 

 IMPLICIT NONE
 PRIVATE
 SAVE
 
 PUBLIC :: micro_p3_init
 PUBLIC :: micro_p3_get_cols
 PUBLIC :: micro_p3_tend
 PUBLIC :: micro_p3_lookuptable_init 
 
  
 INTEGER, PARAMETER :: r8 = selected_real_kind(12)
 INTEGER, PARAMETER :: i8 = selected_int_kind(18)
 
 INTEGER, PARAMETER :: nCat = 1 

 LOGICAL, PARAMETER :: l_conservation_clipping = .True.  
 LOGICAL, PARAMETER :: l_cons_check = .True.

 LOGICAL, PARAMETER :: l_caut     = .True. 
 LOGICAL, PARAMETER :: l_cacc     = .True. 
 LOGICAL, PARAMETER :: l_cnuc     = .False. 
 LOGICAL, PARAMETER :: l_cnuc2    = .True. 
 LOGICAL, PARAMETER :: l_chetc    = .True. 
 LOGICAL, PARAMETER :: l_cslf     = .True. 
 LOGICAL, PARAMETER :: l_rcol     = .True. 
 LOGICAL, PARAMETER :: l_rmul     = .True. 
 LOGICAL, PARAMETER :: l_rslf     = .True. 
 LOGICAL, PARAMETER :: l_rhetc    = .True. 
 LOGICAL, PARAMETER :: l_rheti    = .True. 
 LOGICAL, PARAMETER :: l_idepsub  = .True. 
 LOGICAL, PARAMETER :: l_islf     = .True. 
 LOGICAL, PARAMETER :: l_icol     = .True. 
 LOGICAL, PARAMETER :: l_imul     = .True.  
 LOGICAL, PARAMETER :: l_csed     = .True. 
 LOGICAL, PARAMETER :: l_rsed     = .True. 
 LOGICAL, PARAMETER :: l_ised     = .True. 
 LOGICAL, PARAMETER :: l_crhomo   = .True. 

 REAL(r8), PARAMETER:: cldm_min   = 1.e-20_r8 !! threshold min value for cloud fraction 
 REAL(r8), PARAMETER:: dum_min    = 1.e-20_r8 !! threshold min value 
 
 REAL(r8), PARAMETER:: reffc_def  = 10.e-6_r8 !! default effective radius for droplet 
 REAL(r8), PARAMETER:: reffr_def  = 25.e-6_r8 !! default effective radius for rain   
 REAL(r8), PARAMETER:: reffi_def  = 25.e-6_r8 !! default effective radius for ice   
 
 REAL(r8), PARAMETER:: zerodegc   = 273.15_r8 
 
 !! Range of cloudsat reflectivities (dBz) for analytic simulator

 REAL(r8), PARAMETER:: csmin   = -30._r8
 REAL(r8), PARAMETER:: csmax   =  26._r8
 REAL(r8), PARAMETER:: mindbz  = -99._r8
 REAL(r8), PARAMETER:: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)
 
 !! ice microphysics lookup table array dimensions
 
 INTEGER, PARAMETER :: isize        = 50
 INTEGER, PARAMETER :: iisize       = 25
 INTEGER, PARAMETER :: zsize        = 20  ! size of mom6 array in lookup_table (for future 3-moment)
 INTEGER, PARAMETER :: densize      =  5
 INTEGER, PARAMETER :: rimsize      =  4
 INTEGER, PARAMETER :: rcollsize    = 30
 
 LOGICAL :: log_predictNc = .True. !! IF true, predict droplet number 
 
 CHARACTER(len=16) :: precip_frac_method = 'max_overlap' 
 
 !! tabsize      : number of quantities 
 !! colltabsize  : number of ice-rain collection quantities 
 !! collitabsize : number of ice-ice collection quantities 
 
 INTEGER, PARAMETER :: tabsize = 12  
 INTEGER, PARAMETER :: colltabsize = 2  
 INTEGER, PARAMETER :: collitabsize = 2  
 
 REAL(r8), PARAMETER:: real_rcollsize = real(rcollsize)
 
 !! itab         : ice lookup table values
 !! itabcoll     : ice lookup table values for ice-rain collision/collection
 
 REAL(r8) :: itab(densize,rimsize,isize,tabsize)
 REAL(r8) :: itabcoll(densize,rimsize,isize,rcollsize,colltabsize)
 
 !! separated into itabcolli1 and itabcolli2, due to max of 7 dimensional arrays on some compilers
 
 REAL(r8) :: itabcolli1(iisize,rimsize,densize,iisize,rimsize,densize)
 REAL(r8) :: itabcolli2(iisize,rimsize,densize,iisize,rimsize,densize)
 
 !! switch for warm-rain autoconversion/accretion parameterization
 !! = 1 Seifert and Beheng 2001
 !! = 2 Beheng 1994
 !! = 3 Khairoutdinov and Kogan 2000

 INTEGER, PARAMETER :: iparam = 3
 
 !! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
 !! warm rain autoconversion/accretion option only (iparam = 1)
 
 REAL(r8) :: dnu(16)
 
 !! lookup table values for rain shape parameter mu_r
 
 REAL(r8) :: mu_r_table(150)
 
 !! lookup table values for rain number- and mass-weighted fallspeeds and ventilation parameters
 
 REAL(r8) :: vn_table(300,10), vm_table(300,10), revap_table(300,10)
 
 !! physical and mathematical constants
 
 REAL(r8) :: rhosur
 REAL(r8) :: rhosui
 REAL(r8) :: ar
 REAL(r8) :: br
 REAL(r8) :: f1r
 REAL(r8) :: f2r
 REAL(r8) :: ecr
 REAL(r8) :: rhow
 REAL(r8) :: kr
 REAL(r8) :: kc
 REAL(r8) :: bimm
 REAL(r8) :: aimm
 REAL(r8) :: rin
 REAL(r8) :: mi0
 REAL(r8) :: nccnst
 REAL(r8) :: eci
 REAL(r8) :: eri
 REAL(r8) :: bcn
 REAL(r8) :: cpw
 REAL(r8) :: e0
 REAL(r8) :: cons1
 REAL(r8) :: cons2
 REAL(r8) :: cons3
 REAL(r8) :: cons4
 REAL(r8) :: cons5
 REAL(r8) :: cons6
 REAL(r8) :: cons7
 REAL(r8) :: inv_rhow
 REAL(r8) :: qsmall
 REAL(r8) :: nsmall
 REAL(r8) :: bsmall
 REAL(r8) :: zsmall
 REAL(r8) :: cp
 REAL(r8) :: rd
 REAL(r8) :: ep_2
 REAL(r8) :: inv_cp
 REAL(r8) :: mw
 REAL(r8) :: osm
 REAL(r8) :: vi
 REAL(r8) :: epsm
 REAL(r8) :: rhoa
 REAL(r8) :: map
 REAL(r8) :: ma
 REAL(r8) :: rr
 REAL(r8) :: bact
 REAL(r8) :: inv_rm1
 REAL(r8) :: inv_rm2
 REAL(r8) :: sig1
 REAL(r8) :: nanew1
 REAL(r8) :: f11
 REAL(r8) :: f21
 REAL(r8) :: sig2
 REAL(r8) :: nanew2
 REAL(r8) :: f12
 REAL(r8) :: f22
 REAL(r8) :: pi
 REAL(r8) :: thrd
 REAL(r8) :: sxth
 REAL(r8) :: piov3
 REAL(r8) :: piov6
 REAL(r8) :: diff_nucthrs
 REAL(r8) :: rho_rimeMin
 REAL(r8) :: rho_rimeMax
 REAL(r8) :: inv_rho_rimeMax
 REAL(r8) :: max_total_Ni
 REAL(r8) :: dbrk
 REAL(r8) :: nmltratio
 REAL(r8) :: mi0l
  
 !! minimum mass of new crystal due to freezing of cloud droplets done externally (kg)
 
 REAL(r8) :: mi0l_min
 
 
!!....................................................................
!! Constants set in initialization
!!....................................................................
 
 ! Set using arguments to micro_p3_init
 REAL(r8) :: g           ! gravity
 REAL(r8) :: r           ! dry air gas constant
 REAL(r8) :: rv          ! water vapor gas constant
 REAL(r8) :: cpp         ! specific heat of dry air
 REAL(r8) :: tmelt       ! freezing point of water (K)
 REAL(r8) :: prc_coef1 = huge(1.0_r8)
 REAL(r8) :: prc_exp   = huge(1.0_r8)
 REAL(r8) :: prc_exp1  = huge(1.0_r8)
 REAL(r8) :: cld_sed   = huge(1.0_r8) !scale fac for cld sedimentation velocity
 
 REAL(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.
 
 ! flags
 LOGICAL :: microp_uniform
 LOGICAL :: do_cldice
 LOGICAL :: do_nccons
 LOGICAL :: do_nicons
 LOGICAL :: use_hetfrz_classnuc
 
 REAL(r8) :: ncnst ! constant droplet concentration
 REAL(r8) :: ninst ! constant ice concentration
 
 REAL(r8) :: rhosu       ! typical 850mn air density
 
 REAL(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C
 
 REAL(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
 REAL(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C
 
 ! additional constants to help speed up code
 REAL(r8) :: gamma_br_plus1
 REAL(r8) :: gamma_br_plus4
 REAL(r8) :: gamma_bs_plus1
 REAL(r8) :: gamma_bs_plus4
 REAL(r8) :: gamma_bi_plus1
 REAL(r8) :: gamma_bi_plus4
 
 REAL(r8)           :: micro_p3_berg_eff_factor     ! berg efficiency factor
 
 LOGICAL  :: allow_sed_supersat ! Allow supersaturated conditions after sedimentation loop
 
 ! switch for specification rather than prediction of droplet and crystal number
 ! note: number will be adjusted as needed to keep mean size within bounds,
 ! even when specified droplet or ice number is used
 
 ! IF constant cloud ice number is set (nicons = .True.),
 ! then all microphysical processes except mass transfer due to ice nucleation
 ! (mnuccd) are based on the fixed cloud ice number. Calculation of
 ! mnuccd follows from the prognosed ice crystal number ni.
 
 ! nccons = .True. to specify constant cloud droplet number
 ! nicons = .True. to specify constant cloud ice number
 
 LOGICAL :: nccons 
 LOGICAL :: nicons
 

!!....................................................................
!! interface 
!!....................................................................
 
 INTERFACE rising_factorial
    MODULE PROCEDURE rising_factorial_r8
    MODULE PROCEDURE rising_factorial_integer
 END INTERFACE rising_factorial
 
 INTERFACE var_coef
    MODULE PROCEDURE var_coef_r8
    MODULE PROCEDURE var_coef_integer
 END INTERFACE var_coef




CONTAINS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                    formula   
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! USE gamma function to implement rising factorial extended to the reals.
PURE FUNCTION rising_factorial_r8(x, n) result(res)

   REAL(r8), INTENT(in) :: x, n
   REAL(r8) :: res

   res = gamma(x+n)/gamma(x)

END FUNCTION rising_factorial_r8


! Rising factorial can be performed much cheaper IF n is a small integer.
PURE FUNCTION rising_factorial_integer(x, n) result(res)

   REAL(r8), INTENT(in) :: x
   INTEGER, INTENT(in) :: n
   REAL(r8) :: res

   INTEGER :: i
   REAL(r8) :: factor

   res = 1._r8
   factor = x

   DO i = 1, n
      res = res * factor
      factor = factor + 1._r8
   END DO

END FUNCTION rising_factorial_integer



ELEMENTAL FUNCTION var_coef_r8(relvar, a) result(res)

  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  REAL(r8), INTENT(in) :: relvar
  REAL(r8), INTENT(in) :: a
  REAL(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

END FUNCTION var_coef_r8

ELEMENTAL FUNCTION var_coef_integer(relvar, a) result(res)

  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  REAL(r8), INTENT(in) :: relvar
  INTEGER, INTENT(in) :: a
  REAL(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

END FUNCTION var_coef_integer




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! P3 init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE micro_p3_init( &
     kind,                &
     gravit,              &
     rair,                &
     rh2o,                &
     cpair,               &
     tmelt_in,            &
     latvap,              &
     latice,              &
     rhmini_in,           &
     microp_uniform_in,   &
     do_cldice_in,        &
     use_hetfrz_classnuc_in, &
     do_nccons_in,        &
     do_nicons_in,        &
     ncnst_in,            &
     ninst_in,            &
     allow_sed_supersat_in, &
     ice_sed_ai,          &
     prc_coef1_in,        &
     prc_exp_in,          &
     prc_exp1_in,         &
     cld_sed_in,          &
     errstring            )

 USE micro_p3_utils,   ONLY: micro_p3_utils_init
 
!!....................................................................
!!
!! initialize constants for P3 microphysics
!!
!!....................................................................

 INTEGER,  INTENT(in)  :: kind 
 
 REAL(r8), INTENT(in)  :: & ! 
    gravit,               & !
    rair,                 & !
    rh2o,                 & !
    cpair,                & !
    tmelt_in,             & ! freezing point of water (K)
    latvap,               & !
    latice,               & ! 
    rhmini_in,            & ! minimum rh for ice cloud fraction > 0. 
    prc_coef1_in,         & !
    prc_exp_in,           & !
    prc_exp1_in,          & !
    cld_sed_in,           & !
    ncnst_in,             & !
    ninst_in,             & !
    ice_sed_ai              ! fall speed parameter for cloud ice

 LOGICAL,  INTENT(in)  ::   & !
    microp_uniform_in,      & ! configure uniform for sub-columns IF true 
    do_cldice_in,           & ! skip all processes affecting IF false 
    do_nccons_in,           & ! set cloud droplet to constant IF true 
    do_nicons_in,           & ! set ice concentration to constant IF true 				    
    use_hetfrz_classnuc_in, & ! USE heterogeneous freezing
    allow_sed_supersat_in     ! allow supersaturated conditions after sedimentation loop

 CHARACTER(128), INTENT(out) :: &
    errstring                 ! output status (non-blank for error return)


!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_init - 001 -' 
  
 prc_coef1 = prc_coef1_in
 prc_exp   = prc_exp_in
 prc_exp1  = prc_exp1_in
 cld_sed   = cld_sed_in

 ! Initialize subordinate utilities module.
 CALL micro_p3_utils_init(kind, rh2o, cpair, tmelt_in, latvap, latice, &
                          ice_sed_ai, errstring)


 IF (trim(errstring) /= "") return

 ! declarations for MG code (transforms variable names)

 g= gravit                 ! gravity
 r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
 rv= rh2o                  ! water vapor gas constant
 cpp = cpair               ! specific heat of dry air
 tmelt = tmelt_in
 rhmini = rhmini_in


 ! flags
 microp_uniform = microp_uniform_in
 do_cldice  = do_cldice_in
 nccons = do_nccons_in
 nicons = do_nicons_in
 ncnst = ncnst_in
 ninst = ninst_in
 use_hetfrz_classnuc = use_hetfrz_classnuc_in

 ! typical air density at 850 mb

 rhosu = 85000._r8/(rair * tmelt)

 !! maximum temperature at which snow is allowed to exist
 
 snowmelt = tmelt + 2._r8
 
 !! minimum temperature at which rain is allowed to exist
 
 rainfrze = tmelt - 40._r8

 !! ice nucleation temperature 
 
 icenuct  = tmelt - 5._r8

 ! Define constants to help speed up code (this limits calls to gamma function)
 
! gamma_br_plus1=gamma(1._r8+br)
! gamma_br_plus4=gamma(4._r8+br)
! gamma_bs_plus1=gamma(1._r8+bs)
! gamma_bs_plus4=gamma(4._r8+bs)
! gamma_bi_plus1=gamma(1._r8+bi)
! gamma_bi_plus4=gamma(4._r8+bi)


!!!!!!!!!!!!!!!!!!!!!!!!!!
!! from P3 
!!!!!!!!!!!!!!!!!!!!!!!!!!

! maximum total ice concentration (sum of all categories)
 max_total_Ni = 500.e3_r8  !(m)



! droplet concentration (m-3)
 nccnst = 400.e+6_r8

! parameters for Seifert and Beheng (2001) autoconversion/accretion
 kc     = 9.44e+9_r8
 kr     = 5.78e+3_r8

! physical constants
 cp     = 1005._r8
 inv_cp = 1./cp
 g      = 9.816_r8
 rd     = 287.15_r8
 rv     = 461.51_r8
 ep_2   = 0.622
 rhosur = 100000./(rd*zerodegc)
 rhosui = 60000./(rd*253.15)
 ar     = 841.99667
 br     = 0.8
 f1r    = 0.78
 f2r    = 0.32
 ecr    = 1.
 rhow   = 997.
 cpw    = 4218.
 inv_rhow = 1.e-3  !inverse of (max.) density of liquid water

! limits for rime density [kg m-3]
 rho_rimeMin     =  50.
 rho_rimeMax     = 900.
 inv_rho_rimeMax = 1./rho_rimeMax

! minium allowable prognostic variables
 qsmall = 1.e-14  !! mg2 REAL(r8), parameter, PUBLIC :: qsmall = 1.e-18_r8
 nsmall = 1.e-16
 bsmall = qsmall*inv_rho_rimeMax
!zsmall = 1.e-35



!! REAL(r8), parameter, PUBLIC :: mi0 = 4._r8/3._r8*pi*rhoi(500._r8)*(10.e-6_r8)**3

 eci    = 0.5
 eri    = 1.
 bcn    = 2.

! mean size for soft lambda_r limiter [microns]
 dbrk   = 600.e-6
! ratio of rain number produced to ice number loss from melting
 nmltratio = 0.2

! saturation pressure at T = 0 C
 e0    = polysvp1(zerodegc,0)



! mathematical/optimization constants
 pi    = 3.14159265_r8
 thrd  = 1._r8/3._r8
 sxth  = 1._r8/6._r8
 piov3 = pi*thrd
 piov6 = pi*sxth
 

! Bigg (1953)
!bimm   = 100.
!aimm   = 0.66
! Barklie and Gokhale (1959)
 bimm   = 2.
 aimm   = 0.65
 rin    = 0.1e-6
!!! mi0    = 4.*piov3*900.*1.e-18
 mi0    = 4.*piov3*900.*(10.e-6)**3 !! set it same as MG2, initial size 10um 
 
 mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3
 
 cons1 = piov6*rhow
 cons2 = 4.*piov3*rhow
 cons3 = 1./(cons2*(25.e-6)**3)
 cons4 = 1./(dbrk**3*pi*rhow)
 cons5 = piov6*bimm
 cons6 = piov6**2*rhow*bimm
 cons7 = 4.*piov3*rhow*(1.e-6)**3

! aerosol/droplet activation parameters
 mw     = 0.018
 osm    = 1.
 vi     = 3.
 epsm   = 0.9
 rhoa   = 1777.
 map    = 0.132
 ma     = 0.0284
 rr     = 8.3187
 bact   = vi*osm*epsm*mw*rhoa/(map*rhow)
! inv_bact = (map*rhow)/(vi*osm*epsm*mw*rhoa)    *** to replace /bact **

! mode 1
 inv_rm1 = 2.e+7           ! inverse aerosol mean size (m-1)
 sig1    = 2.0             ! aerosol standard deviation
 nanew1  = 300.e6          ! aerosol number mixing ratio (kg-1)
 f11     = 0.5*exp(2.5*(log(sig1))**2)
 f21     = 1. + 0.25*log(sig1)

! note: currently only set for a single mode, droplet activation code needs to
!       be modified to include the second mode
! mode 2
 inv_rm2 = 7.6923076e+5    ! inverse aerosol mean size (m-1)
 sig2    = 2.5             ! aerosol standard deviation
 nanew2  = 0._r8              ! aerosol number mixing ratio (kg-1)
 f12     = 0.5*exp(2.5*(log(sig2))**2)
 f22     = 1. + 0.25*log(sig2)

! parameters for droplet mass spectral shape, used by Seifert and Beheng (2001)
! warm rain scheme only (iparam = 1)
 dnu(1)  =  0.
 dnu(2)  = -0.557
 dnu(3)  = -0.430
 dnu(4)  = -0.307
 dnu(5)  = -0.186
 dnu(6)  = -0.067
 dnu(7)  =  0.050
 dnu(8)  =  0.167
 dnu(9)  =  0.282
 dnu(10) =  0.397
 dnu(11) =  0.512
 dnu(12) =  0.626
 dnu(13) =  0.739
 dnu(14) =  0.853
 dnu(15) =  0.966
 dnu(16) =  0.966
 

!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_init - 002 -' 

END SUBROUTINE micro_p3_init



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                 lookup-table  
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE micro_p3_lookuptable_init(lookup_file_1,lookup_file_2)

!!....................................................................
!! This SUBROUTINE must be called at the first model time step 
!!....................................................................

 IMPLICIT NONE

 character*(*), INTENT(in) :: lookup_file_1 !! lookup table for main processes
 character*(*), INTENT(in) :: lookup_file_2 !! lookup table for ice category interactions

 INTEGER  :: i,j,k,ii,jj,kk,jjj,jjj2,jjjj,jjjj2
 REAL(r8) :: lamr,mu_r,lamold,dum,initlamr
 REAL(r8) :: dm,dum1,dum2,dum3,dum4,dum5,dum6
 REAL(r8) :: dd,amg,vt,dia,vn,vm
 
 
!!....................................................................
! read in ice microphysics table
!!....................................................................

 IF(masterproc) WRITE(6,*) ' P3 microphysics, version: 2.3.1 '
 IF(masterproc) WRITE(6,*) ' micro_p3_lookuptable_init ' 
 IF(masterproc) WRITE(6,*) ' open file : ', trim(lookup_file_1)


 open(unit=10,file=trim(lookup_file_1), status='old')

 DO jj = 1,densize
    DO ii = 1,rimsize
       DO i = 1,isize
       
          read(10,*) dum,dum,dum,dum,itab(jj,ii,i,1),itab(jj,ii,i,2),           &
               itab(jj,ii,i,3),itab(jj,ii,i,4),itab(jj,ii,i,5),                 &
               itab(jj,ii,i,6),itab(jj,ii,i,7),itab(jj,ii,i,8),dum,             &
               itab(jj,ii,i,9),itab(jj,ii,i,10),itab(jj,ii,i,11),               &
               itab(jj,ii,i,12)
               
       END DO
       
       !! table for ice-rain collection
       
       DO i = 1,isize
          DO j = 1,rcollsize
             read(10,*) dum,dum,dum,dum,dum,itabcoll(jj,ii,i,j,1),              &
              itabcoll(jj,ii,i,j,2),dum
              itabcoll(jj,ii,i,j,1) = dlog10(itabcoll(jj,ii,i,j,1))
              itabcoll(jj,ii,i,j,2) = dlog10(itabcoll(jj,ii,i,j,2))
          END DO
       END DO
       
    END DO
 END DO

 close(10)

!!....................................................................
!! ice-ice collision lookup table
!!....................................................................

 IF (nCat>1) THEN

    IF(masterproc) WRITE(6,*) ' reading lookup-table for multicategory case' 
    IF(masterproc) WRITE(6,*) ' file : ', trim(lookup_file_2) 

    open(unit=10,file=trim(lookup_file_2),status='old')

    DO i = 1,iisize
       DO jjj = 1,rimsize
          DO jjjj = 1,densize
             DO ii = 1,iisize
                DO jjj2 = 1,rimsize
                   DO jjjj2 = 1,densize
                      read(10,*) dum,dum,dum,dum,dum,dum,dum, &
                      itabcolli1(i,jjj,jjjj,ii,jjj2,jjjj2),   &
                      itabcolli2(i,jjj,jjjj,ii,jjj2,jjjj2) 
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO

    close(unit=10)

 else ! for single cat

    itabcolli1 = 0.
    itabcolli2 = 0.

 END IF

!!....................................................................
!! Generate lookup table for rain shape parameter mu_r
!! this is very fast so it can be generated at the start of each run
!! make a 150x1 1D lookup table, this is done in parameter
!! space of a scaled mean size proportional qr/Nr -- initlamr
!!....................................................................

 !! loop over lookup table values

 DO i = 1,150
 
    initlamr = 1./((real(i)*2.)*1.e-6 + 250.e-6)

    !! iterate to get mu_r
    !! mu_r-lambda relationship is from Cao et al. (2008), eq. (7)
    !! start with first guess, mu_r = 0

    mu_r = 0.

    DO ii=1,50
    
       lamr = initlamr*((mu_r+3.)*(mu_r+2.)*(mu_r+1.)/6.)**thrd

       !! new estimate for mu_r based on lambda
       !! set max lambda in formula for mu_r to 20 mm-1, so Cao et al.
       !! formula is not extrapolated beyond Cao et al. data range

       dum  = min(20.,lamr*1.e-3)
       mu_r = max(0._r8,-0.0201*dum**2+0.902*dum-1.718)

       !! IF lambda is converged within 0.1%, then exit loop

       IF (ii.ge.2) THEN
          IF (abs((lamold-lamr)/lamr).lt.0.001) GOTO 111
       END IF

       lamold = lamr

    END DO

111 CONTINUE

    !! assign lookup table values
    
    mu_r_table(i) = mu_r

 END DO

!!....................................................................
!! Generate lookup table for rain fallspeed and ventilation parameters
!! the lookup table is two dimensional as a function of number-weighted mean size
!! proportional to qr/Nr and shape parameter mu_r
!!....................................................................

 mu_r_loop: DO ii = 1,10   !** change 10 to 9, since range of mu_r is 0-8  CONFIRM
!mu_r_loop: DO ii = 1,9   !** change 10 to 9, since range of mu_r is 0-8

    mu_r = real(ii-1)  ! values of mu

    !! loop over number-weighted mean size

    meansize_loop: DO jj = 1,300

       IF (jj.le.20) THEN
          dm = (real(jj)*10.-5.)*1.e-6      ! mean size [m]
       ELSEIF (jj.gt.20) THEN
          dm = (real(jj-20)*30.+195.)*1.e-6 ! mean size [m]
       END IF

       lamr = (mu_r+1)/dm

       !! DO numerical integration over PSD

       dum1 = 0._r8 ! numerator,   number-weighted fallspeed
       dum2 = 0._r8 ! denominator, number-weighted fallspeed
       dum3 = 0._r8 ! numerator,   mass-weighted fallspeed
       dum4 = 0._r8 ! denominator, mass-weighted fallspeed
       dum5 = 0._r8 ! term for ventilation factor in evap
       dd   = 2.

       !! loop over PSD to numerically integrate number and mass-weighted mean fallspeeds
       
       DO kk = 1,10000

          dia = (real(kk,8)*dd-dd*0.5)*1.e-6  ! size bin [m]
          amg = piov6*997.*dia**3           ! mass [kg]
          amg = amg*1000._r8                   ! convert [kg] to [g]

         !get fallspeed as a function of size [m s-1]
          IF (dia*1.e+6.le.134.43)      then
            vt = 4.5795e+3*amg**(2.*thrd)
          ELSEIF (dia*1.e+6.lt.1511.64) THEN
            vt = 4.962e+1*amg**thrd
          ELSEIF (dia*1.e+6.lt.3477.84) THEN
            vt = 1.732e+1*amg**sxth
          else
            vt = 9.17
          END IF

          dum1 = dum1 + vt*10.**(mu_r*log10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
          dum2 = dum2 + 10.**(mu_r*log10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
          dum3 = dum3 + vt*10.**((mu_r+3.)*log10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
          dum4 = dum4 + 10.**((mu_r+3.)*log10(dia)+4.*mu_r)*exp(-lamr*dia)*dd*1.e-6
          dum5 = dum5 + (vt*dia)**0.5*10.**((mu_r+1.)*log10(dia)+3.*mu_r)*exp(-lamr*dia)*dd*1.e-6

       END DO ! kk-loop (over PSD)

       dum2 = max(dum2, 1.e-30)  !to prevent divide-by-zero below
       dum4 = max(dum4, 1.e-30)  !to prevent divide-by-zero below
       dum5 = max(dum5, 1.e-30)  !to prevent log10-of-zero below

       vn_table(jj,ii)    = dum1/dum2
       vm_table(jj,ii)    = dum3/dum4
       revap_table(jj,ii) = 10.**(log10(dum5)+(mu_r+1.)*log10(lamr)-(3.*mu_r))

    END DO meansize_loop

 END DO mu_r_loop

!!....................................................................

!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_lookuptable_init - 002 -' 

END SUBROUTINE micro_p3_lookuptable_init








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!                                    main  
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
SUBROUTINE micro_p3_tend( &
!!....................................................................
!! in
!!....................................................................
   mncol,         &
   ktop,          &
   kbot,          & 
   is_first_step, &
   dt,            &
   scale_epsi,    & 
   scale_berg,    & 
   l_mg2_qidep,   &
   l_satadj, & 
   l_crconevp,    &
   l_massclip ,   &
   l_limit_qidep_qinuc, &
   l_limit_qisub_qrevp, &
   l_cshd,        &
   l_imlt,        &
   l_ccol,        &
   p3_opt_cheti,  &
   p3_opt_inuc,   & 
   pres,          &
   pdel,          &
   uzpl,          &
   dzq,           &
   cldn,          &
   liqcldf,       &
   icecldf,       &
   ql_relvar,     &
   ql_accfac,     &
   naai,          &
   npccn,         &
   rndst,         &
   nacon,         &
   frzimm,        & 
   frzcnt,        & 
   frzdep,        & 
!!....................................................................
!! in/out 
!!....................................................................
   ttn,           & 
   qvn,           &
   qcn,           &
   ncn,           &
   qrn,           &
   nrn,           &
   qitotn,        &
   nitotn,        &
   qirimn,        &
   birimn,        &
!!....................................................................
!! tendency out 
!!....................................................................
   stend,         &
   qvtend,        &
   qctend,        &
   nctend,        &
   qrtend,        &
   nrtend,        &
   qitend,        &
   nitend,        &
   qirimtend,     &
   birimtend,     &
!!....................................................................
!! diag out 
!!....................................................................
   pcprt_tot,     & 
   pcprt_sol,     & 
   diag_prain,    &
   diag_pevap,    &
   diag_nevapr,   &
   diag_ze,       &
   diag_mu_c,     &
   diag_lamc,     & 
   diag_effc,     &
   diag_effi,     &
   diag_deffi,    &
   diag_rflx,     &
   diag_sflx,     &
   diag_vmc,      &
   diag_vmr,      &
   diag_vmi,      &
   diag_ncal,     &
   diag_ncai,     &
   diag_di,       &
   diag_rhopo,    &
!!....................................................................
!! timescale 
!!....................................................................
   diag_epsc,      & ! 
   diag_epsi,      & ! 
   diag_epsr,      & ! 
   diag_reps,      & ! 
   diag_ac1,       & ! 
   diag_ac2,       & ! 
   diag_ac3,       & ! 
   diag_qccon1,    & ! 
   diag_qccon2,    & ! 
   diag_qicon1,    & ! 
   diag_qicon2,    & ! 
   diag_qicon3,    & ! 
!!....................................................................
!! liquid process 
!!....................................................................
   diag_1stqc2qr, & 
   diag_qccon,    & ! cloud droplet condensation (currently zero)
   diag_qrcon,    & ! rain condensation (currently zero)
   diag_qcaut,    & ! qctend due to autoconversion to rain
   diag_ncautc,   & ! nctend due to autoconversion to rain
   diag_ncautr,   & ! nrtend due to autoconversion of cloud water
   diag_qcacc,    & ! qctend due to accretion by rain
   diag_ncacc,    & ! nctend due to accretion by rain
   diag_ncslf,    & ! nctend due to cloud droplet self-collection
   diag_nrslf,    & ! nrtend due to rain self-collection
   diag_ncnuc,    & ! nctend due to activation of CCN
   diag_qcnuc,    & ! qctend due to activation of CCN
   diag_qcevp,    & ! qctend due to cloud droplet evaporation
   diag_qberg,    & ! qctend/qitend due to Bergeron process 
   diag_qrevp,    & ! qrtend due to rain evaporation
   diag_nrevp,    & ! nrtend due to rain evaporation
!!....................................................................
!! ice process
!!....................................................................
   diag_qccol,    &
   diag_qidep,    &
   diag_qrcol,    &
   diag_qinuc,    &
   diag_nccol,    &
   diag_nrcol,    &
   diag_ninuc,    &
   diag_qisub,    &
   diag_qimlt,    &
   diag_nimlt,    &
   diag_nisub,    &
   diag_nislf,    &
   diag_qchetc,   & ! qctend due to contact freezing droplets
   diag_qcheti,   & ! qctend due to immersion freezing droplets
   diag_qrhetc,   &
   diag_qrheti,   &
   diag_nchetc,   &
   diag_ncheti,   &
   diag_nrhetc,   &
   diag_nrheti,   &
   diag_nrshdr,   &
   diag_qcshd,    &
   diag_qrmul,    &
   diag_nimul,    &
   diag_ncshdc,   &
!!....................................................................
!! sedimentation 
!!....................................................................
   diag_sedqc,    & ! qctend due to sedimentation (1/s)
   diag_sednc,    & ! nctend due to sedimentation (1/s)
   diag_sedqr,    & ! qrtend due to sedimentation (1/s)
   diag_sednr,    & ! nrtend due to sedimentation (1/s)
   diag_sedqi,    & ! qitend due to sedimentation (1/s)
   diag_sedni,    & ! nitend due to sedimentation (1/s)
!!....................................................................
!! other stuff 
!!....................................................................
   diag_cmei,     &
   errstring      ) 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Main SUBROUTINE for the P3 microphysics scheme 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 IMPLICIT NONE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! interface   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 INTEGER, INTENT(in) :: &
    mncol,                              & ! number of columns 
    ktop,                               & ! top level index 
    kbot                                  ! bottom level index 
    
 REAL(r8), INTENT(in) :: &
    dt                                    ! model time step (s) 

 LOGICAL, INTENT(in) :: &
    is_first_step                         ! true IF this is the first model step  

 REAL(r8), INTENT(in) :: &
    scale_epsi           !! scaling parameter for deposition time scale  

 REAL(r8), INTENT(in) :: &
    scale_berg           !! scaling parameter for Bergeron process rate 

 LOGICAL, INTENT(in) :: &
    l_mg2_qidep,        &
    l_satadj, & 
    l_crconevp,         &
    l_massclip ,        &
    l_limit_qidep_qinuc,&
    l_limit_qisub_qrevp,&
    l_cshd,             &
    l_imlt,             &
    l_ccol

 INTEGER, INTENT(in) :: & 
    p3_opt_inuc,        &           ! 1: p3  2: e3sm    
    p3_opt_cheti                    ! 1: p3  2: e3sm 

 REAL(r8), INTENT(in) :: &
    pres(1:mncol,ktop:kbot),            & ! pressure (Pa)
    pdel(1:mncol,ktop:kbot),            & ! pressure difference across level (Pa)
    uzpl(1:mncol,ktop:kbot),            & ! vertical air velocity (m s-1)
    dzq(1:mncol,ktop:kbot),             & ! vertical grid spacing (m) 
    cldn(1:mncol,ktop:kbot),            & ! cloud fraction (/)
    liqcldf(1:mncol,ktop:kbot),         & ! liquid cloud fraction (/)
    icecldf(1:mncol,ktop:kbot),         & ! ice cloud fraction (/)
    ql_relvar(1:mncol,ktop:kbot),       & ! cloud water relative variance (/)
    ql_accfac(1:mncol,ktop:kbot),       & ! optional accretion enhancement factor (/)
    naai(1:mncol,ktop:kbot),            & ! ice nucleation number (from microp_aero_ts) (1/kg)
    npccn(1:mncol,ktop:kbot),           & ! ccn activated number tendency (from microp_aero_ts) (1/kg/s)
    rndst(1:mncol,ktop:kbot,nCat),      & ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
    nacon(1:mncol,ktop:kbot,nCat)         ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

 REAL(r8), INTENT(in), POINTER :: &
    frzimm(:,:),                        & ! Number tendency due to immersion freezing (1/cm3)
    frzcnt(:,:),                        & ! Number tendency due to contact freezing (1/cm3)
    frzdep(:,:)                           ! Number tendency due to deposition nucleation (1/cm3)

 REAL(r8), INTENT(in) :: &
    ttn(1:mncol,ktop:kbot),             & ! temperature (K) 
    qvn(1:mncol,ktop:kbot),             & ! water vapor mixing ratio (kg/kg) 
    qcn(1:mncol,ktop:kbot),             & ! cloud, mass mixing ratio (kg/kg)
    ncn(1:mncol,ktop:kbot),             & ! cloud, number mixing ratio (#/kg)
    qrn(1:mncol,ktop:kbot),             & ! rain, mass mixing ratio (kg/kg)
    nrn(1:mncol,ktop:kbot),             & ! rain, number mixing ratio (#/kg)
    qitotn(1:mncol,ktop:kbot,nCat),     & ! ice, total mass mixing ratio (kg/kg)
    qirimn(1:mncol,ktop:kbot,nCat),     & ! ice, rime mass mixing ratio (kg/kg)
    nitotn(1:mncol,ktop:kbot,nCat),     & ! ice, total number mixing ratio (#/kg)
    birimn(1:mncol,ktop:kbot,nCat)        ! ice, rime volume mixing ratio (m3/kg)

 REAL(r8), INTENT(out) :: &
    pcprt_tot (1:mncol),                & ! precipitation rate, total (m s-1) 
    pcprt_sol (1:mncol),                & ! precipitation rate, solid (m s-1)
    diag_prain(1:mncol,ktop:kbot),      & ! total precipitation (rain + snow)
    diag_pevap(1:mncol,ktop:kbot),      & ! rain evaporation 
    diag_nevapr(1:mncol,ktop:kbot,nCat),& ! evaporation of total precipitation (rain + snow)
    diag_ze   (1:mncol,ktop:kbot),      & ! equivalent reflectivity (dBZ) 
    diag_mu_c (1:mncol,ktop:kbot),      & ! size distribution parameter 
    diag_lamc (1:mncol,ktop:kbot),      & ! size distribution parameter 
    diag_effc (1:mncol,ktop:kbot),      & ! effective radius, cloud (m) 
    diag_effi (1:mncol,ktop:kbot,nCat), & ! effective radius, ice (m) 
    diag_deffi(1:mncol,ktop:kbot,nCat), & ! effective diameter, ice (m) 
    diag_rflx (1:mncol,ktop:kbot+1),    & ! grid-box average rain flux (kg m^-2 s^-1) pverp
    diag_sflx (1:mncol,ktop:kbot+1),    & ! grid-box average ice/snow flux (kg m^-2 s^-1) pverp
    diag_vmc  (1:mncol,ktop:kbot),      & ! mass-weighted vc (m s-1) 
    diag_vmr  (1:mncol,ktop:kbot),      & ! mass-weighted vr (m s-1) 
    diag_vmi  (1:mncol,ktop:kbot,nCat), & ! mass-weighted Vi (m s-1) 
    diag_ncal (1:mncol,ktop:kbot),      & ! nc activated (1/m3)
    diag_ncai (1:mncol,ktop:kbot),      & ! ni nucleated (1/m3) 
    diag_di   (1:mncol,ktop:kbot,nCat), & ! mean size of ice (m) 
    diag_rhopo(1:mncol,ktop:kbot,nCat)    ! bulk ice density (kg m-3) 

 REAL(r8) :: &
    pcprt_liq (1:mncol),                & ! precipitation rate, liquid (m s-1) 
    cldm(1:mncol,ktop:kbot),            & ! cloud fraction 
    lcldm(1:mncol,ktop:kbot),           & ! liquid cloud fraction 
    icldm(1:mncol,ktop:kbot),           & ! ice cloud fraction 
    rcldm(1:mncol,ktop:kbot),           & ! precipitation fraction 
    diag_mu_r (1:mncol,ktop:kbot)         ! shape parameter of rain

!!........................................................................................
!! tendency output 
!!........................................................................................

 REAL(r8), INTENT(out) :: &
     stend(1:mncol,ktop:kbot),          & ! microphysical tendency s (W/kg)
     qvtend(1:mncol,ktop:kbot),         & ! microphysical tendency qv (1/s)
     qctend(1:mncol,ktop:kbot),         & ! microphysical tendency qc (1/s)
     nctend(1:mncol,ktop:kbot),         & ! microphysical tendency nc (1/(kg*s))
     qrtend(1:mncol,ktop:kbot),         & ! microphysical tendency qr (1/s)
     nrtend(1:mncol,ktop:kbot),         & ! microphysical tendency nr (1/(kg*s))
     nitend(1:mncol,ktop:kbot,nCat),    & ! microphysical tendency ni (1/(kg*s))
     qitend(1:mncol,ktop:kbot,nCat),    & ! microphysical tendency qi (1/s)
     qirimtend(1:mncol,ktop:kbot,nCat), & ! microphysical tendency qirim (1/s)
     birimtend(1:mncol,ktop:kbot,nCat)    ! microphysical tendency birim (1/(kg*s))

!!........................................................................................
!! time scales 
!!........................................................................................

 REAL(r8), INTENT(out) :: &
   diag_epsc(1:mncol,ktop:kbot),      & ! 
   diag_epsi(1:mncol,ktop:kbot),      & ! 
   diag_epsr(1:mncol,ktop:kbot),      & ! 
   diag_reps(1:mncol,ktop:kbot),      & ! 
   diag_ac1(1:mncol,ktop:kbot),       & ! 
   diag_ac2(1:mncol,ktop:kbot),       & ! 
   diag_ac3(1:mncol,ktop:kbot),       & ! 
   diag_qccon1(1:mncol,ktop:kbot),    & ! 
   diag_qccon2(1:mncol,ktop:kbot),    & ! 
   diag_qicon1(1:mncol,ktop:kbot),    & ! 
   diag_qicon2(1:mncol,ktop:kbot),    & ! 
   diag_qicon3(1:mncol,ktop:kbot)       ! 
   
!!........................................................................................
!! liquid process rates for output 
!!  (all Q process rates in kg kg-1 s-1)
!!  (all N process rates in # kg-1 s-1)
!!........................................................................................

 REAL(r8), INTENT(out) :: &
     diag_1stqc2qr(1:mncol,ktop:kbot),   & ! 1st order rate for direct cw to precip conversion
     diag_qccon(1:mncol,ktop:kbot),      & ! cloud droplet condensation (currently zero) 
     diag_qrcon(1:mncol,ktop:kbot),      & ! rain condensation 
     diag_qcaut(1:mncol,ktop:kbot),      & ! qctend due to autoconversion to rain
     diag_ncautc(1:mncol,ktop:kbot),     & ! nctend due to autoconversion to rain
     diag_ncautr(1:mncol,ktop:kbot),     & ! nrtend due to autoconversion of cloud water
     diag_qcacc(1:mncol,ktop:kbot),      & ! qctend due to accretion by rain
     diag_ncacc(1:mncol,ktop:kbot),      & ! nctend due to accretion by rain
     diag_ncslf(1:mncol,ktop:kbot),      & ! nctend due to cloud droplet self-collection
     diag_nrslf(1:mncol,ktop:kbot),      & ! nrtend due to rain self-collection
     diag_ncnuc(1:mncol,ktop:kbot),      & ! nctend due to activation of CCN
     diag_qcnuc(1:mncol,ktop:kbot),      & ! qctend due to activation of CCN
     diag_qrevp(1:mncol,ktop:kbot),      & ! qrtend due to rain evaporation
     diag_nrevp(1:mncol,ktop:kbot),      & ! nrtend due to rain evaporation
     diag_sedqc(1:mncol,ktop:kbot),      & ! qctend due to sedimentation (1/s)
     diag_sednc(1:mncol,ktop:kbot),      & ! nctend due to sedimentation (1/s)
     diag_sedqr(1:mncol,ktop:kbot),      & ! qrtend due to sedimentation (1/s)
     diag_sednr(1:mncol,ktop:kbot),      & ! nrtend due to sedimentation (1/s)
     diag_qberg(1:mncol,ktop:kbot),      & ! qctend/qitend due to Bergeron process
     diag_qcevp(1:mncol,ktop:kbot)         ! qctend due to cloud droplet evaporation

!!........................................................................................
!! ice process rates for output 
!!  (all Q process rates in kg kg-1 s-1)
!!  (all N process rates in # kg-1 s-1)
!!........................................................................................

  REAL(r8), INTENT(out) :: &
     diag_qccol(1:mncol,ktop:kbot,nCat),    & ! qc/qitend due to collection cloud water
     diag_qidep(1:mncol,ktop:kbot,nCat),    & ! qv/qitend due to vapor deposition
     diag_qrcol(1:mncol,ktop:kbot,nCat),    & ! qi/qrtend due to collection rain mass by ice
     diag_qinuc(1:mncol,ktop:kbot,nCat),    & ! qitend due to deposition/condensation freezing nuc
     diag_nccol(1:mncol,ktop:kbot,nCat),    & ! change in cloud droplet number from collection by ice
     diag_nrcol(1:mncol,ktop:kbot,nCat),    & ! change in rain number from collection by ice
     diag_ninuc(1:mncol,ktop:kbot,nCat),    & ! change in ice number from deposition/cond-freezing nucleation
     diag_qisub(1:mncol,ktop:kbot,nCat),    & ! sublimation of ice
     diag_cmei(1:mncol,ktop:kbot,nCat),     & ! qitend due to deposition/sublimation 
     diag_qimlt(1:mncol,ktop:kbot,nCat),    & ! qitend due to melting of ice
     diag_nimlt(1:mncol,ktop:kbot,nCat),    & ! nitend due to melting of ice
     diag_nisub(1:mncol,ktop:kbot,nCat),    & ! change in ice number from sublimation
     diag_nislf(1:mncol,ktop:kbot,nCat),    & ! change in ice number from collection within a category
     diag_qchetc(1:mncol,ktop:kbot,nCat),   & ! contact freezing droplets
     diag_qcheti(1:mncol,ktop:kbot,nCat),   & ! immersion freezing droplets
     diag_qrhetc(1:mncol,ktop:kbot,nCat),   & ! contact freezing rain
     diag_qrheti(1:mncol,ktop:kbot,nCat),   & ! immersion freezing rain
     diag_nchetc(1:mncol,ktop:kbot,nCat),   & ! contact freezing droplets
     diag_ncheti(1:mncol,ktop:kbot,nCat),   & ! immersion freezing droplets
     diag_nrhetc(1:mncol,ktop:kbot,nCat),   & ! contact freezing rain
     diag_nrheti(1:mncol,ktop:kbot,nCat),   & ! immersion freezing rain
     diag_sedqi(1:mncol,ktop:kbot,nCat),    & ! qitend due to sedimentation (1/s)
     diag_sedni(1:mncol,ktop:kbot,nCat),    & ! nitend due to sedimentation (1/s)
     diag_nrshdr(1:mncol,ktop:kbot,nCat),   & ! source for rain number from collision of rain/ice above freezing and shedding
     diag_qcshd(1:mncol,ktop:kbot,nCat),    & ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
     diag_qrmul(1:mncol,ktop:kbot,nCat),    & ! change in q, ice multiplication from rime-splitnering of rain (not included in the paper)
     diag_nimul(1:mncol,ktop:kbot,nCat),    & ! change in Ni, ice multiplication from rime-splintering (not included in the paper)
     diag_ncshdc(1:mncol,ktop:kbot,nCat)      ! source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)

  REAL(r8) :: &
     diag_qwgrth(1:mncol,ktop:kbot,nCat) ! wet growth rate 
     
  REAL(r8) :: &
     nimax(1:mncol,ktop:kbot) 
     
  CHARACTER(128), INTENT(out) :: &
     errstring ! output status (non-blank for error return)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! local  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!........................................................................................
!! size distribution and fallspeed parameter
!!........................................................................................
 
 REAL(r8) :: &
    lamr     (1:mncol,ktop:kbot),  &
    logn0r   (1:mncol,ktop:kbot),  & 
    diag_effr(1:mncol,ktop:kbot),  & 
    nu       (1:mncol,ktop:kbot),  & 
    cdist    (1:mncol,ktop:kbot),  & 
    cdist1   (1:mncol,ktop:kbot),  & 
    cdistr   (1:mncol,ktop:kbot),  & 
    Vt_nc    (1:mncol,ktop:kbot),  & 
    Vt_qc    (1:mncol,ktop:kbot),  & 
    Vt_nr    (1:mncol,ktop:kbot),  & 
    Vt_qr    (1:mncol,ktop:kbot),  & 
    Vt_qit   (1:mncol,ktop:kbot),  & 
    Vt_nit   (1:mncol,ktop:kbot),  & 
    Vt_zit   (1:mncol,ktop:kbot)

!!........................................................................................
!! in-cloud values 
!!........................................................................................

 REAL(r8) :: &
    tt(1:mncol,ktop:kbot),         &  ! temperature (K)  
    qv(1:mncol,ktop:kbot),         &  ! water vapor mixing ratio (kg kg-1) 
    qc(1:mncol,ktop:kbot),         &  ! cloud, mass mixing ratio (kg kg-1)
    nc(1:mncol,ktop:kbot),         &  ! cloud, number mixing ratio (# kg-1)
    qr(1:mncol,ktop:kbot),         &  ! rain, mass mixing ratio (kg kg-1)
    nr(1:mncol,ktop:kbot),         &  ! rain, number mixing ratio (# kg-1) 
    qitot(1:mncol,ktop:kbot,nCat), &  ! ice, total mass mixing ratio (kg kg-1) 
    qirim(1:mncol,ktop:kbot,nCat), &  ! ice, rime mass mixing ratio (kg kg-1) 
    nitot(1:mncol,ktop:kbot,nCat), &  ! ice, total number mixing ratio (# kg-1) 
    birim(1:mncol,ktop:kbot,nCat)     ! ice, rime volume mixing ratio (m3 kg-1) 

!!........................................................................................
!! grid box mean values 
!!........................................................................................

 REAL(r8) :: &
    gtt(1:mncol,ktop:kbot),         &  ! temperature (K)  
    gqc(1:mncol,ktop:kbot),         &  ! cloud, mass mixing ratio (kg kg-1)
    gnc(1:mncol,ktop:kbot),         &  ! cloud, number mixing ratio (# kg-1)
    gqr(1:mncol,ktop:kbot),         &  ! rain, mass mixing ratio (kg kg-1)
    gnr(1:mncol,ktop:kbot),         &  ! rain, number mixing ratio (# kg-1) 
    gqitot(1:mncol,ktop:kbot,nCat), &  ! ice, total mass mixing ratio (kg kg-1) 
    gqirim(1:mncol,ktop:kbot,nCat), &  ! ice, rime mass mixing ratio (kg kg-1) 
    gnitot(1:mncol,ktop:kbot,nCat), &  ! ice, total number mixing ratio (# kg-1) 
    gbirim(1:mncol,ktop:kbot,nCat)     ! ice, rime volume mixing ratio (m3 kg-1) 
    
    
!!!  REAL(r8), INTENT(out) :: &  
!!!     diag_qrevp(1:mncol,ktop:kbot),  & ! evaporation rate of rain + snow (1/s)
!!!     evapsnow(1:mncol,ktop:kbot),    & ! sublimation rate of snow (1/s)
!!!     prain(1:mncol,ktop:kbot),       & ! production of rain + snow (1/s)
!!!     prodsnow(1:mncol,ktop:kbot),    & ! production of snow (1/s)
!!!     cmeout(1:mncol,ktop:kbot),      & ! evap/sub of cloud (1/s)
!!!     deffi(1:mncol,ktop:kbot),       & ! ice effective diameter for optics (radiation), (micron)
!!!     pgamrad(1:mncol,ktop:kbot),     & ! ice gamma parameter for optics (radiation), (no units)
!!!     lamcrad(1:mncol,ktop:kbot),     & ! slope of droplet distribution for optics (radiation), (1/m)
!!!     qsout(1:mncol,ktop:kbot),       & ! snow mixing ratio (kg/kg)
!!!     dsout(1:mncol,ktop:kbot),       & ! snow diameter (m)
!!!     qrout(1:mncol,ktop:kbot),       & ! grid-box average rain mixing ratio (kg/kg)
!!!     reff_rain(1:mncol,ktop:kbot),   & ! rain effective radius (micron)
!!!     reff_snow(1:mncol,ktop:kbot),   & ! snow effective radius (micron)
!!!     qcsevap(1:mncol,ktop:kbot),     & ! cloud water evaporation due to sedimentation (1/s)
!!!     qisevap(1:mncol,ktop:kbot),     & ! cloud ice sublimation due to sublimation (1/s)
!!!     qvres(1:mncol,ktop:kbot),       & ! residual condensation term to ensure RH .lt. 100% (1/s)
!!!     cmeitot(1:mncol,ktop:kbot),     & ! grid-mean cloud ice sub/dep (1/s)

 REAL(r8) :: &
    gqc_tend,   & 
    gnc_tend,   & 
    gqr_tend,   & 
    gnr_tend,   & 
    gqitot_tend,   & 
    gnitot_tend,   & 
    gqirim_tend,   & 
    gbirim_tend
 
!!........................................................................................
!! liquid-phase microphysical process rates:
!!  (all Q process rates in kg kg-1 s-1)
!!  (all N process rates in # kg-1)
!!........................................................................................

 REAL(r8) :: &
    qrcon,   & ! rain condensation
    qcacc,   & ! cloud droplet accretion by rain
    qcaut,   & ! cloud droplet autoconversion to rain
    ncacc,   & ! change in cloud droplet number from accretion by rain
    ncautc,  & ! change in cloud droplet number from autoconversion
    ncslf,   & ! change in cloud droplet number from self-collection
    nrslf,   & ! change in rain number from self-collection
    ncnuc,   & ! change in cloud droplet number from activation of CCN
    qccon,   & ! cloud droplet condensation
    qcnuc,   & ! activation of cloud droplets from CCN
    qrevp,   & ! rain evaporation
    qcevp,   & ! cloud droplet evaporation
    qberg,   & ! Bergeron process rate (qc->qi) 
    nrevp,   & ! change in rain number from evaporation
    ncautr     ! change in rain number from autoconversion of cloud water

!!........................................................................................
!! ice-phase microphysical process rates:
!!  (all Q process rates in kg kg-1 s-1)
!!  (all N process rates in # kg-1)
!!........................................................................................

 REAL(r8) :: &
    qccol(nCat),     & ! collection cloud water
    qwgrth(nCat),    & ! wet growth rate
    qidep(nCat),     & ! vapor deposition
    qrcol(nCat),     & ! collection rain mass by ice
    qinuc(nCat),     & ! deposition/condensation freezing nuc
    qinuc_ci(nCat),  & ! provisional ice nucleation number from in-situ ice nucl
    nccol(nCat),     & ! change in cloud droplet number from collection by ice
    nrcol(nCat),     & ! change in rain number from collection by ice
    ninuc(nCat),     & ! change in ice number from deposition/cond-freezing nucleation
    qisub(nCat),     & ! sublimation of ice
    qimlt(nCat),     & ! melting of ice
    nimlt(nCat),     & ! melting of ice
    nisub(nCat),     & ! change in ice number from sublimation
    nislf(nCat),     & ! change in ice number from collection within a category
    qchetc(nCat),    & ! contact freezing droplets
    qcheti(nCat),    & ! contact freezing droplets
    qrhetc(nCat),    & ! contact freezing rain
    qrheti(nCat),    & ! immersion freezing rain
    nchetc(nCat),    & ! contact freezing droplets
    ncheti(nCat),    & ! immersion freezing droplets
    nrhetc(nCat),    & ! contact freezing rain
    nrheti(nCat),    & ! immersion freezing rain
    nrshdr(nCat),    & ! source for rain number from collision of rain/ice above freezing and shedding
    qcshd(nCat),     & ! source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
    qrmul(nCat),     & ! change in q, ice multiplication from rime-splitnering of rain 
    nimul(nCat),     & ! change in Ni, ice multiplication from rime-splintering
    ncshdc(nCat),    & ! source for rain number due to cloud water/ice collision above freezing and shedding
    rhorime_c(nCat), & ! density of rime (from cloud)
    rhorime_r(nCat)    ! density of rime (from rain)

 REAL(r8) :: &
    nicol(nCat,nCat), & ! change of N due to ice-ice collision between categories
    qicol(nCat,nCat)    ! change of q due to ice-ice collision between categories

 LOGICAL, DIMENSION(nCat)   :: log_wetgrowth

 REAL(r8) :: &
    epsi(nCat),         &
    Eii_fact(nCat),     &
    diam_ice(1:mncol,ktop:kbot,nCat) 

 REAL(r8) ::        &
    inv_dzq(1:mncol,ktop:kbot),   &
    inv_rho(1:mncol,ktop:kbot),   &
    ze_ice(1:mncol,ktop:kbot),    &
    ze_rain(1:mncol,ktop:kbot),   &
    prec(1:mncol,ktop:kbot),      &
    rho(1:mncol,ktop:kbot),       &
    rhofacr(1:mncol,ktop:kbot),   &
    rhofaci(1:mncol,ktop:kbot),   &
    acn(1:mncol,ktop:kbot),       &
    xxls(1:mncol,ktop:kbot),      &
    xxlv(1:mncol,ktop:kbot),      &
    xlf(1:mncol,ktop:kbot),       &
    qvs(1:mncol,ktop:kbot),       &
    qvi(1:mncol,ktop:kbot),       &
    sup(1:mncol,ktop:kbot),       &
    supi(1:mncol,ktop:kbot),      &
    ssat(1:mncol,ktop:kbot),      &
    minstsm(1:mncol,ktop:kbot),     & !
    ninstsm(1:mncol,ktop:kbot),     & !
    minstrf(1:mncol,ktop:kbot),     & !
    ninstrf(1:mncol,ktop:kbot),     & !
    ss(1:mncol,ktop:kbot),        &
    vtrmi1(1:mncol,ktop:kbot),    &
    vtrnitot

 REAL(r8) :: & ! 
             ptmp, xqvs, qcadj(1:mncol,ktop:kbot), qvadj(1:mncol,ktop:kbot), qiadj(1:mncol,ktop:kbot)

 REAL(r8) :: & ! 
    dum_qit(ktop:kbot),     & !
    dum_qr(ktop:kbot),      & !
    dum_nit(ktop:kbot),     & !
    dum_qir(ktop:kbot),     & !
    dum_bir(ktop:kbot),     & !
    dum_zit(ktop:kbot),     & !
    dum_nr(ktop:kbot),      & ! 
    dum_qc(ktop:kbot),      & !
    dum_nc(ktop:kbot),      & !
    V_qr(ktop:kbot),        & !
    V_qit(ktop:kbot),       & !
    V_nit(ktop:kbot),       & !
    V_nr(ktop:kbot),        & !
    V_qc(ktop:kbot),        & !
    V_nc(ktop:kbot),        & !
    V_zit(ktop:kbot),       & !
    flux_qr(ktop:kbot),     & !
    flux_qit(ktop:kbot),    & !
    flux_nit(ktop:kbot),    & !
    flux_nr(ktop:kbot),     & !
    flux_qir(ktop:kbot),    & !
    flux_bir(ktop:kbot),    & !
    flux_zit(ktop:kbot),    & !
    flux_qc(ktop:kbot),     & !
    flux_nc(ktop:kbot),     & !
    tend_qc(ktop:kbot),     & !
    tend_qr(ktop:kbot),     & !
    tend_nr(ktop:kbot),     & !
    tend_qit(ktop:kbot),    & !
    tend_qir(ktop:kbot),    & !
    tend_bir(ktop:kbot),    & !
    tend_nit(ktop:kbot),    & !
    tend_nc(ktop:kbot),     & !
    tend_zit(ktop:kbot)

 REAL(r8) ::   & !
    eii,       & ! temperature dependent aggregation efficiency
    lammax,    & ! 
    lammin,    & ! 
    mu,        & ! 
    dv,        & ! 
    sc,        & ! 
    dqsdt,     & ! 
    ab,        & ! 
    kap,       & ! 
    epsr,      & ! 
    epsc,      & ! 
    xx,        & ! 
    aaa,       & ! 
    epsilon,   & ! 
    sigvl,     & ! 
    epsi_tot,  & ! 
    aact,      & ! 
    alpha,     & ! 
    gamm,      & ! 
    gg,        & ! 
    psi,       & ! 
    eta1,      & ! 
    eta2,      & ! 
    sm1,       & ! 
    sm2,       & ! 
    smax,      & ! 
    uu1,       & ! 
    uu2,       & ! 
    dum,       & ! 
    dum0,      & ! 
    dum1,      & ! 
    dum2,      & ! 
    dumqv,     & ! 
    dumqvs,    & ! 
    dums,      & ! 
    dumqc,     & ! 
    ratio,     & ! 
    qsat0,     & ! 
    udiff,     & ! 
    dum3,      & ! 
    dum4,      & ! 
    dum5,      & ! 
    dum6,      & ! 
    lamold,    & ! 
    rdumii,    & !
    rdumjj,    & ! 
    dqsidt,    & ! 
    abi,       & ! 
    dumqvi,    & ! 
    dap,       & ! 
    nacnt,     & ! 
    rhop,      & ! 
    v_impact,  & ! 
    ri,        & ! 
    iTc,       & ! 
    D_c,       & ! 
    D_r,       & ! 
    dumlr,     & ! 
    tmp1,      & ! 
    tmp2,      & ! 
    tmp3,      & ! 
    inv_nstep, & ! 
    inv_dum,   & ! 
    inv_dum3,  & ! 
    odt,       & ! 
    oxx,       & ! 
    oabi,      & ! 
    zero,      & ! 
    test,      & ! 
    test2,     & ! 
    test3,     & ! 
    onstep,    & ! 
    fluxdiv_qr,    & ! 
    fluxdiv_qit,   & ! 
    fluxdiv_nit,   & ! 
    fluxdiv_qir,   & ! 
    fluxdiv_bir,   & !  
    fluxdiv_zit,   & ! 
    fluxdiv_qc,    & ! 
    fluxdiv_nc,    & ! 
    fluxdiv_nr,    & ! 
    rgvm,          & ! 
    D_new,         & ! 
    Q_nuc,         & ! 
    N_nuc,         & ! 
    deltaD_init,   & ! 
    dum1c,         & ! 
    dum4c,         & ! 
    dum5c,         & ! 
    dumt,          & ! 
    drhop,         & ! 
    timeScaleFactor, & !
    qvt,           & ! 
    mtime

 INTEGER :: &
    dumi,          & !
    i,             & !
    k,             & !
    kk,            & !
    ii,            & !
    jj,            & !
    iice,          & !
    iice_dest,     & !
    j,             & !
    dumk,          & !
    dumj,          & !
    dumii,         & !
    dumjj,         & !
    dumzz,         & !
    n,             & !
    nstep,         & !
    tmpint1,       & !
    tmpint2,       & !
    kdir,          & !
    qcindex,       & !
    qrindex,       & !
    qiindex,       & !
    dumic,         & !
    dumiic,        & !
    dumjjc,        & !
    catcoll

 LOGICAL :: &
    log_nucleationPossible, & !
    log_hydrometeorsPresent, & !
    log_predictSsat, & !
    log_tmp1,      & !
    log_exitlevel, & !
    log_hmossopOn, & !
    log_qcpresent, & !
    log_qrpresent, & !
    log_qipresent, & !
    log_ni_add

! quantities related to process rates/parameters, interpolated from lookup tables:

 REAL(r8) :: &
    f1pr01,   & ! number-weighted fallspeed
    f1pr02,   & ! mass-weighted fallspeed
    f1pr03,   & ! ice collection within a category
    f1pr04,   & ! collection of cloud water by ice
    f1pr05,   & ! melting
    f1pr06,   & ! effective radius
    f1pr07,   & ! collection of rain number by ice
    f1pr08,   & ! collection of rain mass by ice
    f1pr09,   & ! minimum ice number (lambda limiter)
    f1pr10,   & ! maximum ice number (lambda limiter)
    f1pr11,   & ! not used
    f1pr12,   & ! not used
    f1pr13,   & ! reflectivity
    f1pr14,   & ! melting (ventilation term)
    f1pr15,   & ! mass-weighted mean diameter
    f1pr16,   & ! mass-weighted mean particle density
    f1pr17,   & ! ice-ice category collection change in number
    f1pr18      ! ice-ice category collection change in mass

 REAL(r8) :: &
    prc_coef, &  
    pra_coef, &
    pra_coef0
    
 REAL(r8) :: &
    xdum1, xdum2, xdum3, xdum4, &
    xdum5, xdum6, xdum7, xdum8, &
    qtmp, ttmp 

 INTEGER :: &
    bad_state(1:mncol,ktop:kbot) 
    
 INTEGER :: &
    n_bad_state 

 LOGICAL :: & 
    l_possible, & 
    l_dum1,     & 
    l_dum2,     & 
    l_dum3,     & 
    l_dum4,     & 
    l_dum5 
 
 !!...............................................
 !! Setup level index 
 !!...............................................
 
 kdir = -1 

 !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 001 -' 

 !!...............................................
 !! default return error message
 !!...............................................
  
 errstring = ' '

 IF (use_hetfrz_classnuc       .and. &
    (.not. (associated(frzimm) .and. &
            associated(frzcnt) .and. &
            associated(frzdep) ) )   ) THEN
     errstring = "External heterogeneous freezing is enabled, but the &
                 &required tendencies were not all passed in."
 END IF


!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 002 -' 


 !!...............................................
 !! grid mean values 
 !!...............................................
 
 tt  = ttn 
 qv  = qvn 
 
 gqc = qcn
 gnc = ncn 
 gqr = qrn 
 gnr = nrn 
 
 gqitot = qitotn 
 gqirim = qirimn 
 gnitot = nitotn 
 gbirim = birimn 

 !!...............................................
 !! init in-cloud values  
 !!...............................................

 qc = 0._r8
 nc = 0._r8
 qr = 0._r8
 nr = 0._r8

 qitot = 0._r8
 nitot = 0._r8
 qirim = 0._r8 
 birim = 0._r8
 
 !!...............................................
 !! cloud fraction 
 !!...............................................

 lcldm = 0._r8 
 icldm = 0._r8 
 
 IF (microp_uniform) THEN
 
    !! subcolumns, set cloud fraction variables to one
    !! IF cloud water or ice is present, IF not present
    !! set to mincld (mincld used instead of zero, to prevent
    !! possible division by zero errors).

    WHERE  (gqc .ge. qsmall)
       lcldm = 1._r8
    ELSEWHERE
       lcldm = mincld
    END WHERE 

    WHERE (gqitot(:,:,1) .ge. qsmall)
       icldm = 1._r8
    ELSEWHERE
       icldm = mincld
    END WHERE 

    cldm = max(icldm, lcldm)

 ELSE
 
    !! get cloud fraction, check for minimum
    
    cldm  = max(cldn,   mincld)
    lcldm = max(liqcldf,mincld)
    icldm = max(icecldf,mincld)
    
 END IF

 !!...............................................
 !! precipitation fraction 
 !!...............................................
 
 rcldm = mincld

 !!DO k = kbot,ktop,-1
 DO k = ktop,kbot
 
    DO i=1,mncol
    
       !! 
       !! precipitation fraction 
       !! 
  
       rcldm(i,k) = cldm(i,k)
  
       IF (trim(precip_frac_method) == 'in_cloud') THEN
  
          IF (k /= ktop) THEN
             IF (gqc(i,k) .lt. qsmall .and. sum(gqitot(i,k,:)) .lt. qsmall) THEN 
                rcldm(i,k) = rcldm(i,k-1)
             END IF
          END IF
  
       ELSE IF (trim(precip_frac_method) == 'max_overlap') THEN
  
          ! calculate precip fraction based on maximum overlap assumption
  
          ! IF rain or snow mix ratios are smaller than threshold,
          ! then leave rcldm as cloud fraction at current level

          IF (k /= ktop) THEN
             IF (gqr(i,k-1) .ge. qsmall .or. sum(gqitot(i,k-1,:)) .ge. qsmall) THEN 
                rcldm(i,k) = max(rcldm(i,k-1),rcldm(i,k))
             END IF
          END IF
  
       END IF
 
    END DO 
 END DO 

!! DO k = ktop,kbot
!!    WRITE(6,*) '# ##99 max, min of rcldm : ', k, maxval(rcldm(:,k)), minval(rcldm(:,k))
!! END DO
 
 diag_epsc   = 0._r8
 diag_epsi   = 0._r8
 diag_epsr   = 0._r8

 !!...............................................
 !! calculate some time-varying atmospheric variables
 !!...............................................
 
 DO k = kbot,ktop,kdir
 
    DO i=1,mncol
    
       rho(i,k)     = pres(i,k)/(rd*tt(i,k))
       inv_rho(i,k) = 1./rho(i,k)
       
       !! xxlv    vaporization
       !! xlf     freezing
       !! xxls    sublimation

       xxlv(i,k)    = 3.1484e6_r8 - 2370._r8 * tt(i,k)
       xxls(i,k)    = xxlv(i,k)+0.3337e6_r8
       xlf(i,k)     = xxls(i,k)-xxlv(i,k)
       
       !! saturation vapor pressure 
       
       dum0 = polysvp1(tt(i,k),0)
       dum1 = polysvp1(tt(i,k),1)

       !! saturation qs and qi 
       
       qvs(i,k) = ep_2*dum0/max(1.e-3,(pres(i,k)-dum0)) !! saturation qs
       qvi(i,k) = ep_2*dum1/max(1.e-3,(pres(i,k)-dum1)) !! saturation qi
       
       !! supersaturation 

       ssat(i,k)    = qv(i,k) - qvs(i,k)
       sup(i,k)     = qv(i,k) / qvs(i,k) - 1._r8
       supi(i,k)    = qv(i,k) / qvi(i,k) - 1._r8

       !! air density adjustment for fallspeed parameters includes air density correction 
       !! factor to the power of 0.54 following Heymsfield and Bansemer 2007

       rhofacr(i,k) = (rhosur*inv_rho(i,k))**0.54_r8
       rhofaci(i,k) = (rhosui*inv_rho(i,k))**0.54_r8
       
       !! mu 
       
       dum = 1.496e-6 * tt(i,k)**1.5_r8 / (tt(i,k)+120._r8)  
       
       !! 'a' parameter for droplet fallspeed (Stokes' law)
       
       acn(i,k) = g*rhow/(18._r8*dum)  

    END DO 
 END DO 


 !!.......................................................................................
 !! droplet activation (as in MG2)  
 !!.......................................................................................

 diag_ncal   = 0._r8
 diag_ncai   = 0._r8
 
 IF(l_cnuc) THEN

  where (gqc >= qsmall)
     gnc = max(gnc + npccn*dt, 0._r8)
     diag_ncal = gnc*rho/lcldm !! sghan minimum in #/cm3
  elsewhere
     diag_ncal = 0._r8
  end where

!!! DO k = kbot,ktop,kdir
!!! 
!!!    DO i=1,mncol
!!!    
!!!     !!...................................
!!!     !! 
!!!     !! get provisional droplet number after activation. This is used for
!!!     !! all microphysical process calculations, for consistency with update of
!!!     !! droplet mass before microphysics
!!!     !!
!!!     !! calculate potential for droplet activation IF cloud water is present
!!!     !! npccn (activation tendency) is from microp_aero 
!!!     !!
!!!     !! npccn is grid-mean value 
!!!     !! 
!!!     !! output activated liquid and ice (convert from #/kg -> #/m3)
!!!     !!...................................
!!!
!!!     l_dum1 = gqc(i,k) .ge. qsmall 
!!!     l_dum2 = lcldm(i,k) .ge. cldm_min
!!!     l_possible = l_dum1 .and. l_dum2
!!!     
!!!     IF (l_possible) THEN 
!!!        gnc(i,k) = max(gnc(i,k) + npccn(i,k)*dt, 0._r8)
!!!        diag_ncal(i,k) = gnc(i,k) * rho(i,k) / lcldm(i,k)
!!!     ELSE
!!!        diag_ncal(i,k) = 0._r8
!!!     END IF
!!!
!!!    END DO 
!!! END DO 

 END IF 


!!! !!.......................................................................................
!!! !! calculate instantaneous precip processes (melting and homogeneous freezing)
!!! !!.......................................................................................
!!!
!!! DO k = kbot,ktop,kdir
!!!    DO i=1,mncol
!!!    
!!!        ! melting of snow at +2 C
!!!
!!!        if (tt(i,k) > snowmelt) then
!!!           if (qvs(i,k) > 0._r8) then
!!!
!!!              ! make sure melting snow doesn't reduce temperature below threshold
!!!              dum = -xlf(i,k)/cpp*gqitot(i,k,1)
!!!              if (tt(i,k)+dum < snowmelt) then
!!!                 dum = (tt(i,k)-snowmelt)*cpp/xlf(i,k)
!!!                 dum = dum/gqitot(i,k,1)
!!!                 dum = max(0._r8,dum)
!!!                 dum = min(1._r8,dum)
!!!              else
!!!                 dum = 1._r8
!!!              end if
!!!
!!!              minstsm(i,k) = dum*gqitot(i,k,1)
!!!              ninstsm(i,k) = dum*gnitot(i,k,1)
!!!
!!!              dum1=-xlf(i,k)*minstsm(i,k)/dt
!!!              !!tlat(i,k)=tlat(i,k)+dum1
!!!              tt(i,k)=tt(i,k)+dum1
!!!              !!meltsdttot(i,k)=meltsdttot(i,k) + dum1
!!!
!!!              gqitot(i,k,1) = max(gqitot(i,k,1) - minstsm(i,k), 0._r8)
!!!              gnitot(i,k,1) = max(gnitot(i,k,1) - ninstsm(i,k), 0._r8)
!!!              gqr(i,k) = max(gqr(i,k) + minstsm(i,k), 0._r8)
!!!              gnr(i,k) = max(gnr(i,k) + ninstsm(i,k), 0._r8)
!!!           end if
!!!        end if
!!!
!!!        ! freezing of rain at -5 C
!!!
!!!        if (tt(i,k) < rainfrze) then
!!!
!!!           if (gqr(i,k) > 0._r8) then
!!!
!!!              ! make sure freezing rain doesn't increase temperature above threshold
!!!              dum = xlf(i,k)/cpp*gqr(i,k)
!!!              if (tt(i,k)+dum > rainfrze) then
!!!                 dum = -(tt(i,k)-rainfrze)*cpp/xlf(i,k)
!!!                 dum = dum/gqr(i,k)
!!!                 dum = max(0._r8,dum)
!!!                 dum = min(1._r8,dum)
!!!              else
!!!                 dum = 1._r8
!!!              end if
!!!
!!!              minstrf(i,k) = dum*gqr(i,k)
!!!              ninstrf(i,k) = dum*gnr(i,k)
!!!
!!!              ! heating tendency
!!!              dum1 = xlf(i,k)*minstrf(i,k)/dt
!!!              !!tlat(i,k)=tlat(i,k)+dum1
!!!              tt(i,k)=tt(i,k)+dum1
!!!              !!!frzrdttot(i,k)=frzrdttot(i,k) + dum1
!!!
!!!              gqr(i,k) = max(gqr(i,k) - minstrf(i,k), 0._r8)
!!!              gnr(i,k) = max(gnr(i,k) - ninstrf(i,k), 0._r8)
!!!              gqitot(i,k,1) = max(gqitot(i,k,1) + minstrf(i,k), 0._r8)
!!!              gnitot(i,k,1) = max(gnitot(i,k,1) + ninstrf(i,k), 0._r8)
!!!
!!!           end if
!!!        end if
!!!
!!!    END DO 
!!! END DO 

 !!.......................................................................................
 !! get in-cloud and in-precipitation fields 
 !! 
 !!   note in-cloud values are initialized to zero above 
 !!   limit in-cloud liquid and ice values to 0.005 kg/kg 
 !!   limit in-precip mixing ratios to 0.01 kg/kg 
 !!.......................................................................................

 DO k = kbot,ktop,kdir
 
    DO i=1,mncol
          
       !!.......................
       !! droplet  
       !!.......................
       
       l_possible = gqc(i,k) .ge. qsmall
       
       IF (l_possible) THEN
           
          qc(i,k) = min(gqc(i,k)/lcldm(i,k),5.e-3_r8)
          nc(i,k) = max(gnc(i,k)/lcldm(i,k),0._r8)

          ! specify droplet concentration
          IF (nccons) THEN
             nc(i,k) = ncnst/rho(i,k)
          END IF
       
       END IF

       DO iice = 1,nCat
        
          !!.......................
          !! ice 
          !!.......................
       
          l_possible = gqitot(i,k,iice) .ge. qsmall
       
          IF (l_possible) THEN
             
             qitot(i,k,iice) = min(gqitot(i,k,iice) / icldm(i,k),5.e-3_r8)
             nitot(i,k,iice) = max(gnitot(i,k,iice) / icldm(i,k),0._r8)

             ! switch for specification of cloud ice number
             IF (nicons .and. iice.eq.1) THEN
                nitot(i,k,iice) = ninst/rho(i,k)
             END IF
             
          END IF

          !!.......................
          !! rime 
          !!.......................

          l_dum1 = gqirim(i,k,iice) .ge. qsmall
          l_dum2 = gqitot(i,k,iice) .ge. qsmall 

          l_possible = l_dum1 .and. l_dum2 
       
          IF (l_possible) THEN
             qirim(i,k,iice) = min(gqirim(i,k,iice) / icldm(i,k),5.e-3_r8)
             birim(i,k,iice) = max(gbirim(i,k,iice) / icldm(i,k),0._r8)
          END IF
          
       END DO 

       !!.......................
       !! rain 
       !!.......................

       l_possible = gqr(i,k).ge.qsmall

       IF(l_possible) THEN
          qr(i,k) = min(gqr(i,k) / rcldm(i,k), 0.01_r8) 
          nr(i,k) = max(gnr(i,k) / rcldm(i,k), 0._r8) 
       END IF 

    END DO 
 END DO 


!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 003 -' 
  
 !!.......................................................................................
 !! Determine threshold size difference [m] as a function of nCat
 !! (used for destination category upon ice initiation) 
 !!.......................................................................................
 
 SELECT CASE (nCat)
    CASE (1)
       deltaD_init = 999._r8
    CASE (2)
       deltaD_init = 500.e-6_r8
    CASE (3)
       deltaD_init = 400.e-6_r8
    CASE (4)
       deltaD_init = 235.e-6_r8
    CASE (5)
       deltaD_init = 175.e-6_r8
    CASE (6:)
       deltaD_init = 150.e-6_r8
 END SELECT 



 log_predictSsat = .false.   ! for prediction of supersaturation

 log_hmossopOn   = (nCat.gt.1)      !default: off for nCat=1, off for nCat>1 

 inv_dzq    = 1./dzq  ! inverse of thickness of layers
 odt        = 1./dt   ! inverse model time step

 !!
 !! Compute time scale factor over which to apply soft rain lambda limiter
 !!
 
 timeScaleFactor = min(1./120., odt)

 !!................................
 !! initialize diag fields 
 !!................................

 pcprt_tot   = 0._r8
 pcprt_liq   = 0._r8
 pcprt_sol   = 0._r8
 diag_prain  = 0._r8
 diag_pevap  = 0._r8
 diag_nevapr = 0._r8
 diag_ze     = 0._r8 
 diag_mu_c   = 1._r8
 diag_lamc   = 1._r8
 diag_effc   = reffc_def
 diag_effr   = reffr_def
 diag_effi   = reffi_def
 diag_deffi  = reffi_def
 diag_rflx   = 0._r8
 diag_sflx   = 0._r8 
 diag_vmc    = 0._r8
 diag_vmr    = 0._r8
 diag_vmi    = 0._r8
 diag_di     = 0._r8
 diag_rhopo  = 0._r8
 
 prec        = 0._r8
 diag_mu_r   = 0._r8
 diam_ice    = 0._r8
 ze_ice      = 1.e-22_r8
 ze_rain     = 1.e-22_r8
 rhorime_c   = 400.
!rhorime_r   = 400.


 !!................................
 !! time scales 
 !!................................

 diag_reps   = 0._r8
 diag_ac1    = 0._r8
 diag_ac2    = 0._r8
 diag_ac3    = 0._r8
 diag_qccon1 = 0._r8
 diag_qccon2 = 0._r8
 diag_qicon1 = 0._r8
 diag_qicon2 = 0._r8
 diag_qicon3 = 0._r8
   
 !!!diag_qicon1(:,:) = ql_relvar(:,:)

 !!................................
 !! warm-phase process rates
 !!................................

 diag_1stqc2qr = 0._r8 !! ### 
 
 diag_qccon   = 0._r8
 diag_qrcon   = 0._r8
 diag_qcaut   = 0._r8
 diag_ncautc  = 0._r8
 diag_ncautr  = 0._r8
 diag_qcacc   = 0._r8
 diag_ncacc   = 0._r8
 diag_ncslf   = 0._r8
 diag_nrslf   = 0._r8
 diag_ncnuc   = 0._r8
 diag_qcnuc   = 0._r8
 diag_qcevp   = 0._r8
 diag_qberg   = 0._r8
 diag_qrevp   = 0._r8
 diag_nrevp   = 0._r8
 diag_sedqc   = 0._r8 
 diag_sednc   = 0._r8 
 diag_sedqr   = 0._r8 
 diag_sednr   = 0._r8 

 !!................................
 !! ice-phase process rates
 !!................................
     
 diag_qccol   = 0._r8
 diag_qidep   = 0._r8
 diag_qrcol   = 0._r8
 diag_qinuc   = 0._r8
 diag_nccol   = 0._r8
 diag_nrcol   = 0._r8
 diag_ninuc   = 0._r8
 diag_qisub   = 0._r8 
 diag_qimlt   = 0._r8
 diag_nimlt   = 0._r8
 diag_nisub   = 0._r8
 diag_nislf   = 0._r8
 diag_qchetc  = 0._r8
 diag_qcheti  = 0._r8
 diag_qrhetc  = 0._r8
 diag_qrheti  = 0._r8
 diag_nchetc  = 0._r8
 diag_ncheti  = 0._r8
 diag_nrhetc  = 0._r8
 diag_nrheti  = 0._r8
 diag_nrshdr  = 0._r8
 diag_sedqi   = 0._r8 
 diag_sedni   = 0._r8 
 diag_qcshd   = 0._r8
 diag_qrmul   = 0._r8
 diag_nimul   = 0._r8
 diag_ncshdc  = 0._r8
 diag_cmei    = 0._r8
 diag_qwgrth  = 0._r8
 nimax        = 0._r8 

 !!................................
 !! tendency 
 !!................................
 
 stend        = 0._r8 
 qvtend       = 0._r8 
 qctend       = 0._r8 
 nctend       = 0._r8 
 qrtend       = 0._r8 
 nrtend       = 0._r8 
 nitend       = 0._r8 
 qitend       = 0._r8 
 qirimtend    = 0._r8 
 birimtend    = 0._r8 


 !!.......................................................................................
 !! main i-loop (around the entire scheme) 
 !!.......................................................................................
 
 
 DO i = 1,mncol
 

    log_hydrometeorsPresent = .false.
    log_nucleationPossible  = .false.
    log_ni_add              = .false. 
 
    !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 005 -', i 

    !!...........................................
    !! initialize some variables 
    !!...........................................

    DO k = kbot,ktop,kdir
    
       l_dum1 = tt(i,k).lt.zerodegc 
       l_dum2 = supi(i,k).ge.-0.05_r8 
       l_dum3 = tt(i,k).ge.zerodegc
       l_dum4 = sup(i,k) .ge.-0.05_r8
       
       l_possible = (l_dum1 .and. l_dum2) .or. (l_dum3 .and. l_dum4) 
       
       IF (l_possible) THEN 
          log_nucleationPossible = .True.
       END IF 

    END DO !! k_loop_1

    DO k = kbot,ktop,kdir
    
       !!.................................................................................
       !! apply mass clipping if dry and mass is sufficiently small
       !!.................................................................................

       !!..................
       !! droplet 
       !!..................
       
       l_dum1 = qc(i,k) .lt. qsmall
       l_dum2 = qc(i,k) .lt. 1.e-8_r8
       l_dum3 = sup(i,k) .lt. -0.1_r8
       
       l_possible = l_dum1 .or. (l_dum2 .and. l_dum3) 
       
       IF (l_massclip .and. l_possible) THEN
          
          qv(i,k) = qv(i,k) + qc(i,k) * lcldm(i,k) 
          
          IF(qv(i,k) .lt. 0._r8) THEN 
             WRITE(6,*) '### L2097 ### qv .lt. 0' 
          END IF 
          
          tt(i,k) = tt(i,k) - qc(i,k) * lcldm(i,k) * xxlv(i,k) * inv_cp
          
          qc(i,k) = 0._r8
          nc(i,k) = 0._r8
          
       ELSE
          log_hydrometeorsPresent = .True. ! updated further down
       END IF

    END DO !! k_loop_1


    DO k = kbot,ktop,kdir
    
       !!..................
       !! rain 
       !!..................
       
       l_dum1 = qr(i,k) .lt. qsmall
       l_dum2 = qr(i,k) .lt. 1.e-8_r8
       l_dum3 = sup(i,k) .lt. -0.1_r8
       
       l_possible = l_dum1 .or. (l_dum2 .and. l_dum3) 
       
       IF (l_massclip .and. l_possible) THEN
          
          qv(i,k) = qv(i,k) + qr(i,k) * rcldm(i,k) 
          
          IF(qv(i,k) .lt. 0._r8) THEN 
             WRITE(6,*) '### L2110 ### qv .lt. 0' 
          END IF 
          
          tt(i,k) = tt(i,k) - qr(i,k) * rcldm(i,k) * xxlv(i,k) * inv_cp
          
          qr(i,k) = 0._r8
          nr(i,k) = 0._r8
          
       ELSE
          log_hydrometeorsPresent = .True. ! updated further down
       END IF

    END DO !! k_loop_1


    DO k = kbot,ktop,kdir
    
       !!..................
       !! ice  
       !!..................

       DO iice = 1,nCat
       
          l_dum1 = qitot(i,k,iice).lt.qsmall
          l_dum2 = qitot(i,k,iice).lt.1.e-8_r8
          l_dum3 = supi(i,k).lt.-0.1_r8
       
          l_possible = l_dum1 .or. (l_dum2 .and. l_dum3) 
       
          IF (l_massclip .and. l_possible) THEN
             
             qv(i,k) = qv(i,k) + qitot(i,k,iice) * icldm(i,k) 
             
             IF(qv(i,k) .lt. 0._r8) THEN 
                WRITE(6,*) '### L2067 ### qv .lt. 0'  
             END IF 
             
             tt(i,k) = tt(i,k) - qitot(i,k,iice) * icldm(i,k) * xxls(i,k) * inv_cp
             
             qitot(i,k,iice) = 0._r8
             nitot(i,k,iice) = 0._r8
             qirim(i,k,iice) = 0._r8
             birim(i,k,iice) = 0._r8
          ELSE
             log_hydrometeorsPresent = .True.    ! final update
          END IF

       END DO !! iice 

    END DO !! k_loop_1


    DO k = kbot,ktop,kdir
    
       DO iice = 1,nCat
       
          l_dum1 = qitot(i,k,iice) .ge. qsmall
          l_dum2 = qitot(i,k,iice) .lt. 1.e-8_r8
          l_dum3 = tt(i,k) .ge. zerodegc
       
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
          
          IF (l_possible) THEN
              
             qr(i,k) = qr(i,k) + qitot(i,k,iice) * icldm(i,k) / rcldm(i,k) 
             tt(i,k) = tt(i,k) - qitot(i,k,iice) * icldm(i,k) * xlf(i,k) * inv_cp
             
             qitot(i,k,iice) = 0._r8
             nitot(i,k,iice) = 0._r8
             qirim(i,k,iice) = 0._r8
             birim(i,k,iice) = 0._r8
             
          END IF

       END DO !! iice 


    END DO !! k_loop_1


    !!.................................................................................
    !! jump to end of i-loop if OK  
    !!.................................................................................
       
    l_possible = .not. (log_nucleationPossible .or. log_hydrometeorsPresent)
          
    IF (l_possible) THEN 
       GOTO 333
    END IF 

    log_hydrometeorsPresent = .false.   ! reset value; used again below


    !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 006 -' 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! main k-loop (for processes)
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    DO k = kbot,ktop,kdir

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 007 - k = ', k 

       !! if relatively dry and no hydrometeors at this level, 
       !! skip to end of k-loop (i.e. skip this level)
     
       log_exitlevel = .True.
       
       IF (qc(i,k).ge.qsmall .or. qr(i,k).ge.qsmall) THEN
          log_exitlevel = .false.
       END IF
       
       DO iice = 1,nCat
          IF (qitot(i,k,iice).ge.qsmall) log_exitlevel = .false.
       END DO

       !! skip all process rate if possible 
       
       l_dum1 = tt(i,k).lt.zerodegc
       
       l_dum2 = supi(i,k).lt.-0.05
       
       l_dum3 = tt(i,k).ge.zerodegc
       
       l_dum4 = sup(i,k) .lt.-0.05
       
       l_possible = (l_dum1 .and. l_dum2) .or. (l_dum3 .and. l_dum4) 
       
       IF (log_exitlevel .and. l_possible) THEN  
          GOTO 555   
       END IF 

       !! warm-phase process rates
       
       qcacc   = 0._r8
       qrevp   = 0._r8
       qccon   = 0._r8
       qcaut   = 0._r8
       qcevp   = 0._r8
       qberg   = 0._r8
       qrcon   = 0._r8
       ncacc   = 0._r8
       ncnuc   = 0._r8
       ncslf   = 0._r8
       ncautc  = 0._r8
       qcnuc   = 0._r8
       nrslf   = 0._r8
       nrevp   = 0._r8
       ncautr  = 0._r8

       !! ice-phase  process rates
       
       qchetc  = 0._r8
       qisub   = 0._r8
       nrshdr  = 0._r8
       qcheti  = 0._r8
       qrcol   = 0._r8
       qcshd   = 0._r8
       qrhetc  = 0._r8
       qimlt   = 0._r8
       qccol   = 0._r8
       qrheti  = 0._r8
       qinuc   = 0._r8
       nimlt   = 0._r8
       nchetc  = 0._r8
       nccol   = 0._r8
       ncshdc  = 0._r8
       ncheti  = 0._r8
       nrcol   = 0._r8
       nislf   = 0._r8
       nrhetc  = 0._r8
       ninuc   = 0._r8
       qidep   = 0._r8
       nrheti  = 0._r8
       nisub   = 0._r8
       qwgrth  = 0._r8
       qrmul   = 0._r8
       nimul   = 0._r8
       qicol   = 0._r8
       nicol   = 0._r8

       qinuc_ci = 0._r8
       
       log_wetgrowth = .false.

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 008 -'


!!                                                                   DO k = kbot,ktop,kdir
    
!!!KZ       IF (log_predictSsat) THEN
!!!KZ
!!!KZ          !! Adjust cloud water and thermodynamics to prognostic supersaturation
!!!KZ          !! following the method in Grabowski and Morrison (2008).
!!!KZ          !! Note that the effects of vertical motion are assumed to dominate the 
!!!KZ          !! production term for supersaturation, and the effects are sub-grid
!!!KZ          !! scale mixing and radiation are not explicitly included.
!!!KZ
!!!KZ          dqsdt   = xxlv(i,k)*qvs(i,k)/(rv*tt(i,k)*tt(i,k))
!!!KZ          ab      = 1. + dqsdt*xxlv(i,k)*inv_cp
!!!KZ          epsilon = (qv(i,k)-qvs(i,k)-ssat(i,k))/ab
!!!KZ          epsilon = max(epsilon,-qc(i,k))   ! limit adjustment to available water
!!!KZ
!!!KZ          !! don't adjust upward IF subsaturated
!!!KZ          !! otherwise this could result in positive adjustment
!!!KZ          !! (spurious generation ofcloud water) in subsaturated conditions
!!!KZ
!!!KZ          IF (ssat(i,k) .lt. 0._r8) epsilon = min(0._r8,epsilon)
!!!KZ
!!!KZ          !! now DO the adjustment
!!!KZ
!!!KZ          IF (abs(epsilon).ge.1.e-15) THEN
!!!KZ
!!!KZ             qc(i,k)   = qc(i,k)+epsilon
!!!KZ             qv(i,k)   = qv(i,k)-epsilon * lcldm(i,k) 
!!!KZ             tt(i,k)   = tt(i,k)+epsilon * lcldm(i,k) * xxlv(i,k) * inv_cp
!!!KZ             
!!!KZ             !! recalculate variables IF there was adjustment
!!!KZ             !! t(i,k)    = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)
!!!KZ
!!!KZ             dum0      = polysvp1(tt(i,k),0)
!!!KZ             dum1      = polysvp1(tt(i,k),1)
!!!KZ             qvs(i,k)  = ep_2*dum0/max(1.e-3,(pres(i,k)-dum0))
!!!KZ             qvi(i,k)  = ep_2*dum1/max(1.e-3,(pres(i,k)-dum1))
!!!KZ             sup(i,k)  = qv(i,k)/qvs(i,k)-1.
!!!KZ             supi(i,k) = qv(i,k)/qvi(i,k)-1.
!!!KZ             ssat(i,k) = qv(i,k)-qvs(i,k)
!!!KZ          END IF
!!!KZ
!!!KZ       END IF !! predict_supersaturation
       

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 009 -'


       !! skip micro process calculations except nucleation/acvtivation IF there no hydrometeors are present

       log_exitlevel = .True.

       IF (qc(i,k).ge.qsmall .or. qr(i,k).ge.qsmall) THEN
          log_exitlevel = .false.
       END IF 

       DO iice = 1,nCat
          IF (qitot(i,k,iice).ge.qsmall) THEN
             log_exitlevel=.false.
          END IF
       END DO

       IF (log_exitlevel) THEN
          GOTO 444   !i.e. skip to nucleation
       END IF 

       !! time/space varying physical variables

       mu     = 1.496e-6*tt(i,k)**1.5/(tt(i,k)+120._r8)
       dv     = 8.794e-5*tt(i,k)**1.81/pres(i,k)
       sc     = mu/(rho(i,k)*dv)
       dum    = 1./(rv*tt(i,k)**2)
       dqsdt  = xxlv(i,k)*qvs(i,k)*dum
       dqsidt = xxls(i,k)*qvi(i,k)*dum
       ab     = 1.+dqsdt*xxlv(i,k)*inv_cp
       abi    = 1.+dqsidt*xxls(i,k)*inv_cp
       kap    = 1.414e+3*mu

       !! very simple temperature dependent aggregation efficiency
       !! linear ramp from 0.1 to 1 between 253.15 and 268.15 K
       
       xdum1 = 253.15_r8
       xdum2 = 268.15_r8
       
       IF (tt(i,k).lt.xdum1) THEN
          eii = 0.1_r8 
       ELSE IF (tt(i,k).ge.xdum1.and.tt(i,k).lt.xdum2) THEN
          eii = 0.1_r8+(tt(i,k)-xdum1)/15._r8*0.9_r8  
       ELSE IF (tt(i,k).ge.xdum2) THEN
          eii = 1._r8 
       END IF

       CALL get_cloud_dsd(qc(i,k),        &
                          nc(i,k),        &
                          diag_mu_c(i,k), &
                          rho(i,k),       &
                          nu(i,k),        &
                          dnu,            &
                          diag_lamc(i,k), &
                          lammin,         &
                          lammax,         &
                          k,              &
                          cdist(i,k),     &
                          cdist1(i,k),    &
                          tmpint1,        &
                          log_tmp1        )

       CALL get_rain_dsd(qr(i,k),         &
                         nr(i,k),         &
                         diag_mu_r(i,k),  &
                         rdumii,          &
                         dumii,           &
                         lamr(i,k),       &
                         mu_r_table,      &
                         cdistr(i,k),     &
                         logn0r(i,k),     &
                         log_tmp1,        &
                         tmpint1,         &
                         tmpint2          )

       !! note: log_tmp1,tmpint1,tmpint2 are not used in this section

       !! initialize inverse supersaturation relaxation timescale for combined ice categories
       epsi_tot = 0._r8

       CALL impose_max_total_Ni(nitot(i,k,:),max_total_Ni,inv_rho(i,k))


       !!......................................................................
       !! loop over ice categories 
       !!......................................................................

       DO iice = 1,nCat

          IF (qitot(i,k,iice).ge.qsmall) THEN

             !! impose lower limits to prevent taking log of # < 0

             nitot(i,k,iice) = max(nitot(i,k,iice),nsmall)
             nr(i,k)         = max(nr(i,k),nsmall)

             !! compute mean-mass ice diameters 
             !! (estimated; rigorous approach to be implemented later)
            
             dum2 = 500._r8 !! ice density
             
             diam_ice(i,k,iice) = ((qitot(i,k,iice)*6.)/(nitot(i,k,iice)*dum2*pi))**thrd

             CALL calc_bulkRhoRime(qitot(i,k,iice), &
                                   qirim(i,k,iice), &
                                   birim(i,k,iice), &
                                   rhop             ) 

             CALL find_lookupTable_indices_1( &
                                   dumi,dumj,dumjj,dumii,dumzz, &
                                   dum1,dum3,dum4,dum5,dum6,    &
                                   isize,rimsize,densize,zsize,rcollsize, &
                                   qr(i,k),         &
                                   nr(i,k),         &
                                   qitot(i,k,iice), &
                                   nitot(i,k,iice), &
                                   qirim(i,k,iice), &
                                   999._r8,rhop,    &
                                   100              )

             !! lookup table interpolation to get process rates 
             
             CALL access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,f1pr02)
             CALL access_lookup_table(dumjj,dumii,dumi, 3,dum1,dum4,dum5,f1pr03)
             CALL access_lookup_table(dumjj,dumii,dumi, 4,dum1,dum4,dum5,f1pr04)
             CALL access_lookup_table(dumjj,dumii,dumi, 5,dum1,dum4,dum5,f1pr05)
             CALL access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
             CALL access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
             CALL access_lookup_table(dumjj,dumii,dumi,10,dum1,dum4,dum5,f1pr14)

             !! ice-rain collection processes
             IF (qr(i,k).ge.qsmall) THEN
                CALL access_lookup_table_coll(dumjj,dumii,dumj,dumi,1,dum1, &
                                              dum3,dum4,dum5,f1pr07)
                CALL access_lookup_table_coll(dumjj,dumii,dumj,dumi,2,dum1, &
                                              dum3,dum4,dum5,f1pr08)
             ELSE
                f1pr07 = 0._r8
                f1pr08 = 0._r8
             END IF

             !! adjust Ni IF needed to make sure mean size is in bounds (i.e. apply lambda limiters)
             !! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
             
             nitot(i,k,iice) = min(nitot(i,k,iice),f1pr09*nitot(i,k,iice))
             nitot(i,k,iice) = max(nitot(i,k,iice),f1pr10*nitot(i,k,iice))

             !! Determine additional collection efficiency factor to be applied to ice-ice collection.
             !! The computed values of qicol and nicol are multipiled by Eii_fact to gradually shut off collection
             !! IF the ice in iice is highly rimed.
             !! linear ramp from 1 to 0 for Fr between 0.6 and 0.9
             
             l_dum1 = qirim(i,k,iice) .gt. 0._r8 
             l_dum2 = qitot(i,k,iice) .gt. 0._r8
             
             l_possible = l_dum1 .and. l_dum2
             
             xdum1 = 0.6_r8
             xdum2 = 0.9_r8 
             
             IF (l_possible) THEN
                tmp1 = qirim(i,k,iice) / qitot(i,k,iice)   !! rime mass fraction
                IF (tmp1.lt.xdum1) THEN
                   Eii_fact(iice) = 1._r8
                ELSE IF (tmp1.ge.xdum1.and.tmp1.lt.xdum2) THEN
                   Eii_fact(iice) = 1.-(tmp1-xdum1)/0.3_r8
                ELSE IF (tmp1.ge.xdum2) THEN
                   Eii_fact(iice) = 0._r8
                END IF
             ELSE
                Eii_fact(iice) = 1._r8
             END IF

          END IF !! qitot > qsmall


          !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 010 -'


!!                                                                   DO k = kbot,ktop,kdir

!!######################################################################
!!#
!!#
!!# ice microphysical processes
!!#
!!#
!!######################################################################


          !!...................................
          !! collection of droplets
          !!...................................

          !! here we multiply rates by air density, air density fallspeed correction
          !! factor, and collection efficiency since these parameters are not
          !! included in lookup table calculations

          IF(l_ccol) THEN 

          !! for T < 273.15_r8, assume collected cloud water is instantly frozen
          !! note 'f1pr' values are normalized, so we need to multiply by N

          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = qc(i,k).ge.qsmall
          l_dum3 = tt(i,k).le.zerodegc
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
       
          IF (l_possible) THEN
             qccol(iice) = rhofaci(i,k)*f1pr04*qc(i,k)*eci*rho(i,k)*nitot(i,k,iice)
             nccol(iice) = rhofaci(i,k)*f1pr04*nc(i,k)*eci*rho(i,k)*nitot(i,k,iice)
          END IF
          
          qccol(iice) = max(qccol(iice),0._r8) !! ### bugfix 
          nccol(iice) = max(nccol(iice),0._r8) !! ### bugfix 

          !! for T > 273.15, assume cloud water is collected and shed as rain drops

          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = qc(i,k).ge.qsmall
          l_dum3 = tt(i,k).gt.zerodegc
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
          
          IF (l_possible) THEN
          
             !! sink for cloud water mass and number, note qcshed is source for rain mass
             qcshd(iice) = rhofaci(i,k)*f1pr04*qc(i,k)*eci*rho(i,k)*nitot(i,k,iice)
             nccol(iice) = rhofaci(i,k)*f1pr04*nc(i,k)*eci*rho(i,k)*nitot(i,k,iice)
             
             !! source for rain number, assume 1mm drops are shed
             ncshdc(iice) = qcshd(iice)*1.923e+6
             
          END IF

          qcshd(iice) = max(qcshd(iice),0._r8) !! ### bugfix 
          nccol(iice) = max(nccol(iice),0._r8) !! ### bugfix 
          ncshdc(iice) = max(ncshdc(iice),0._r8) !! ### bugfix 
          
          END IF !! IF(l_ccol) THEN 


          !!...................................
          !! collection of rain
          !!...................................

          IF(l_rcol) THEN 

          !! here we multiply rates by air density, air density fallspeed correction
          !! factor, collection efficiency, and n0r since these parameters are not
          !! included in lookup table calculations

          !! for T < 273.15, assume all collected rain mass freezes
          !! note this is a sink for rain mass and number and a source
          !! for ice mass
          !! note 'f1pr' values are normalized, so we need to multiply by N

          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = qr(i,k).ge.qsmall
          l_dum3 = tt(i,k).le.zerodegc
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
          
          IF (l_possible) THEN
             qrcol(iice) = 10.**(f1pr08+logn0r(i,k))*rho(i,k)*rhofaci(i,k)*eri*nitot(i,k,iice)
             nrcol(iice) = 10.**(f1pr07+logn0r(i,k))*rho(i,k)*rhofaci(i,k)*eri*nitot(i,k,iice)
          END IF

          !! for T > 273.15, assume collected rain number is shed as
          !! 1 mm drops
          !! note that melting of ice number is scaled to the loss
          !! rate of ice mass due to melting
          !! collection of rain above freezing does not impact total rain mass
     
          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = qr(i,k).ge.qsmall
          l_dum3 = tt(i,k).gt.zerodegc
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
          
          IF (l_possible) THEN
             !! rain number sink due to collection
             nrcol(iice)  = 10.**(f1pr07 + logn0r(i,k))*rho(i,k)*rhofaci(i,k)*eri*nitot(i,k,iice)
             !! rain number source due to shedding = collected rain mass/mass of 1 mm drop
             dum    = 10.**(f1pr08 + logn0r(i,k))*rho(i,k)*rhofaci(i,k)*eri*nitot(i,k,iice)
             !! for now neglect shedding of ice collecting rain above freezing, since snow is
             !! not expected to shed in these conditions (though more hevaily rimed ice would be
             !! expected to lead to shedding)
             !! nrshdr(iice) = dum*1.923e+6   ! 1./5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop
          END IF

          qrcol(iice) = max(qrcol(iice),0._r8) !! ### bugfix 
          nrcol(iice) = max(nrcol(iice),0._r8) !! ### bugfix 
          
          END IF 


!!!KZ          !!...................................
!!!KZ          !! collection between ice categories
!!!KZ          !!...................................
!!!KZ
!!!KZ          IF(l_icol) THEN !! for multi-categories only 
!!!KZ
!!!KZ          IF (iice.ge.2) THEN
!!!KZ             IF (qitot(i,k,iice).ge.qsmall) THEN
!!!KZ                DO catcoll = 1,iice-1
!!!KZ                   IF (qitot(i,k,catcoll).ge.qsmall) THEN
!!!KZ
!!!KZ                      !! first, calculate collection of catcoll category by iice category
!!!KZ                      !! IF (.not. tripleMoment_on) zitot(i,k,iice) = diag_mom6(qitot(i,k,iice),nitot(i,k,iice),rho(i,k))
!!!KZ
!!!KZ                      CALL find_lookupTable_indices_2(dumi,dumii,dumjj,dumic,dumiic,   &
!!!KZ                                 dumjjc,dum1,dum4,dum5,dum1c,dum4c,dum5c,              &
!!!KZ                                 iisize,rimsize,densize,                               &
!!!KZ                                 qitot(i,k,iice),qitot(i,k,catcoll),nitot(i,k,iice),   &
!!!KZ                                 nitot(i,k,catcoll),qirim(i,k,iice),qirim(i,k,catcoll), &
!!!KZ                                 birim(i,k,iice),birim(i,k,catcoll))
!!!KZ
!!!KZ                      CALL access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,  &
!!!KZ                                 dumi,1,dum1c,dum4c,dum5c,dum1,dum4,dum5,f1pr17)
!!!KZ                                 
!!!KZ                      CALL access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,  &
!!!KZ                                 dumi,2,dum1c,dum4c,dum5c,dum1,dum4,dum5,f1pr18)
!!!KZ
!!!KZ                      !! note: need to multiply by air density, air density fallspeed correction factor,
!!!KZ                      !!       and N of the collectee and collector categories for process rates nicol and qicol,
!!!KZ                      !!       first index is the collectee, second is the collector
!!!KZ                    
!!!KZ                      nicol(catcoll,iice) = f1pr17*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
!!!KZ                                            nitot(i,k,catcoll)*nitot(i,k,iice)
!!!KZ                      qicol(catcoll,iice) = f1pr18*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
!!!KZ                                            nitot(i,k,catcoll)*nitot(i,k,iice)
!!!KZ
!!!KZ                      nicol(catcoll,iice) = eii*Eii_fact(iice)*nicol(catcoll,iice)
!!!KZ                      qicol(catcoll,iice) = eii*Eii_fact(iice)*qicol(catcoll,iice)
!!!KZ
!!!KZ                      !! second, calculate collection of iice category by catcoll category
!!!KZ                      !! IF (.not. tripleMoment_on) zitot(i,k,iice) = diag_mom6(qitot(i,k,iice),nitot(i,k,iice),rho(i,k))
!!!KZ                      !! needed to force consistency between qirim(catcoll) and birim(catcoll) (not for rhop)
!!!KZ                      
!!!KZ                      CALL calc_bulkRhoRime(qitot(i,k,catcoll), &
!!!KZ                                            qirim(i,k,catcoll), &
!!!KZ                                            birim(i,k,catcoll), &
!!!KZ                                            rhop)
!!!KZ
!!!KZ                      CALL find_lookupTable_indices_2(dumi,dumii,dumjj,dumic,dumiic,  &
!!!KZ                                 dumjjc,dum1,dum4,dum5,dum1c,dum4c,dum5c,             &
!!!KZ                                 iisize,rimsize,densize,                              &
!!!KZ                                 qitot(i,k,catcoll),qitot(i,k,iice),nitot(i,k,catcoll),    &
!!!KZ                                 nitot(i,k,iice),qirim(i,k,catcoll),qirim(i,k,iice),       &
!!!KZ                                 birim(i,k,catcoll),birim(i,k,iice))
!!!KZ
!!!KZ                      CALL access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj, &
!!!KZ                                 dumi,1,dum1c,dum4c,dum5c,dum1,dum4,dum5,f1pr17)
!!!KZ
!!!KZ                      CALL access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj, &
!!!KZ                                 dumi,2,dum1c,dum4c,dum5c,dum1,dum4,dum5,f1pr18)
!!!KZ
!!!KZ                      nicol(iice,catcoll) = f1pr17*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
!!!KZ                                            nitot(i,k,iice)*nitot(i,k,catcoll)
!!!KZ                      qicol(iice,catcoll) = f1pr18*rhofaci(i,k)*rhofaci(i,k)*rho(i,k)*     &
!!!KZ                                            nitot(i,k,iice)*nitot(i,k,catcoll)
!!!KZ                     ! note: Eii_fact applied to the collector category
!!!KZ                      nicol(iice,catcoll) = eii*Eii_fact(catcoll)*nicol(iice,catcoll)
!!!KZ                      qicol(iice,catcoll) = eii*Eii_fact(catcoll)*qicol(iice,catcoll)
!!!KZ
!!!KZ                   END IF !! qitotcatcoll_notsmall
!!!KZ                END DO !! catcoll_loop
!!!KZ             END IF !! qitot_notsmall
!!!KZ          END IF !! iceice_interaction1
!!!KZ
!!!KZ          END IF 


          !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 011 -'


          !!...................................
          !! self-collection of ice (in a given category) 
          !!...................................

          IF(l_islf) THEN 


          !! here we multiply rates by collection efficiency, air density,
          !! and air density correction factor since these are not included
          !! in the lookup table calculations
          !! note 'f1pr' values are normalized, so we need to multiply by N

          IF (qitot(i,k,iice).ge.qsmall) THEN
             nislf(iice) = f1pr03*rho(i,k)*eii*Eii_fact(iice)*rhofaci(i,k)*nitot(i,k,iice)
          END IF

          END IF


          IF(l_imlt) THEN 

          !!...................................
          !! melting
          !!...................................

          !! nned to add back accelerated melting due to collection of ice mass by rain (pracsw1)
          !! note 'f1pr' values are normalized, so we need to multiply by N

          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = tt(i,k).gt.zerodegc
          l_possible = l_dum1 .and. l_dum2
          
          IF (l_possible) THEN
             qsat0 = 0.622_r8*e0/(pres(i,k)-e0)
             dum   = 0._r8
             !! include RH dependence
             qimlt(iice) = ( ( f1pr05 + &
                               f1pr14*sc**thrd*(rhofaci(i,k)*rho(i,k)/mu)**0.5 ) &
                             * ( (tt(i,k)-zerodegc)*kap - &
                                 rho(i,k)*xxlv(i,k)*dv*(qsat0-qv(i,k)) ) &
                             * 2.*pi/xlf(i,k)+dum &
                           )*nitot(i,k,iice)
             qimlt(iice) = max(qimlt(iice),0._r8)
             qimlt(iice) = min(qimlt(iice),qitot(i,k,iice)/dt) !! ### bugfix 
             dum         = -qimlt(iice)*dt/qitot(i,k,iice)
             dum         = max(-1.,dum)
             nimlt(iice) = dum*nitot(i,k,iice)*odt
          END IF

          END IF 


          !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 012 -'

          !!...................................
          !! calculate wet growth
          !!...................................

          IF(l_cshd) THEN 

          !! similar to Musil (1970), JAS
          !! note 'f1pr' values are normalized, so we need to multiply by N

          l_dum1 = qitot(i,k,iice) .ge. qsmall
          l_dum2 = (qc(i,k)+qr(i,k)) .ge. 1.e-6_r8 
          l_dum3 = tt(i,k) .lt. zerodegc
          
          l_possible = l_dum1 .and. l_dum2 .and. l_dum3
          
          IF (l_possible) THEN

             qsat0  = 0.622_r8*e0/(pres(i,k)-e0)
             
             xdum1 = sc**thrd*(rhofaci(i,k)*rho(i,k)/mu)**0.5_r8  
             xdum2 = rho(i,k)*xxlv(i,k)*dv*(qsat0-qv(i,k)) 
             xdum3 = (tt(i,k)-zerodegc)*kap 
             xdum4 = xlf(i,k)+cpw*(tt(i,k)-zerodegc)
             
             qwgrth(iice) = (f1pr05 + f1pr14*xdum1) * 2._r8 * pi * (xdum2-xdum3)/xdum4 &
                            * nitot(i,k,iice) 
                       
             qwgrth(iice) = max(qwgrth(iice),0._r8)
             
             !! calculate shedding for wet growth
             
             dum = max(0._r8,(qccol(iice)+qrcol(iice))-qwgrth(iice))
             
             IF (dum.ge.1.e-10) THEN
             
                !! 1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop 
                nrshdr(iice) = nrshdr(iice) + dum*1.923e+6_r8   
                
                IF ((qccol(iice)+qrcol(iice)).ge.1.e-10) THEN
                   dum1  = 1._r8/(qccol(iice)+qrcol(iice))
                   qcshd(iice) = qcshd(iice) + dum*qccol(iice)*dum1
                   qccol(iice) = qccol(iice) - dum*qccol(iice)*dum1
                   qrcol(iice) = qrcol(iice) - dum*qrcol(iice)*dum1
               END IF
               
               !! densify due to wet growth
               log_wetgrowth(iice) = .True.
             END IF

          END IF

          END IF 

          !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 013 -'


          !!........................................................................
          !! supersaturation relaxation time scales for qi
          !!........................................................................

          !! note 'f1pr' values are normalized, so we need to multiply by N

          l_dum1 = qitot(i,k,iice).ge.qsmall
          l_dum2 = tt(i,k).lt.zerodegc
          
          l_possible = l_dum1 .and. l_dum2 
          
          IF (l_possible) THEN 
             xdum1 = sc**thrd*(rhofaci(i,k)*rho(i,k)/mu)**0.5
!!!          epsi(iice) = ( (f1pr05+f1pr14*xdum1) * 2.*pi*rho(i,k)*dv ) * nitot(i,k,iice)
             epsi(iice) = scale_epsi * ( (f1pr05+f1pr14*xdum1) * 2.*pi*rho(i,k)*dv ) * nitot(i,k,iice)
             epsi_tot   = epsi_tot + epsi(iice)
          ELSE
             epsi(iice) = 0.
          END IF

          diag_epsi(i,k) = epsi_tot !! ### 

          !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 014 -'

          !!...................................
          !! calculate rime density
          !!...................................

          !! FUTURE:  Add source term for birim (=qccol/rhorime_c) so that all process rates calculations
          !!          are done together, before conservation.
          !! 
          !! NOTE: Tc (ambient) is assumed for the surface temperature.  Technically,
          !! we should diagose graupel surface temperature from heat balance equation.
          !! (but the ambient temperature is a reasonable approximation; tests show
          !! very little sensitivity to different assumed values, Milbrandt and Morrison 2012).
          !! 
          !! Compute rime density: (based on parameterization of Cober and List, 1993 [JAS])
          !! for simplicty USE mass-weighted ice and droplet/rain fallspeeds
          !! 
          !! IF (qitot(i,k,iice).ge.qsmall .and. tt(i,k).lt.273.15_r8) THEN
          !!  NOTE:  condition applicable for cloud only; modify when rain is added back
        
          l_dum1 = qccol(iice).ge.qsmall 
          l_dum2 = tt(i,k).lt.zerodegc
          
          l_possible = l_dum1 .and. l_dum2 
          
          IF (l_possible) THEN

             !! get mass-weighted mean ice fallspeed
           
             vtrmi1(i,k) = f1pr02*rhofaci(i,k)
             iTc         = 1./min(-0.001,tt(i,k)-zerodegc)

             !!............
             !! cloud 
             !!............
          
             IF (qc(i,k).ge.qsmall) THEN
             
                !! droplet fall speed
                !! (use Stokes' formulation (thus USE analytic solution) 
                
                xdum1 = gamma( 4._r8 + bcn + diag_mu_c(i,k) )
                xdum2 = gamma( diag_mu_c(i,k) + 4._r8 ) 
                xdum3 = diag_lamc(i,k)**bcn 
                
                Vt_qc(i,k) = acn(i,k) * xdum1 /(diag_lamc(i,k)**bcn*xdum2)
                
                !! USE mass-weighted mean size
                D_c = (diag_mu_c(i,k)+4.)/diag_lamc(i,k)
                V_impact  = abs(vtrmi1(i,k)-Vt_qc(i,k))
                Ri        = -(0.5e+6*D_c)*V_impact*iTc
                Ri        = max(1.,min(Ri,12.))
                IF (Ri.le.8.) THEN
                   rhorime_c(iice)  = (0.051 + 0.114*Ri - 0.0055*Ri**2)*1000.
                ELSE
                   !! for Ri > 8 assume a linear fit between 8 and 12,
                   !! rhorime = 900 kg m-3 at Ri = 12
                   !! this is somewhat ad-hoc but allows a smoother transition
                   !! in rime density up to wet growth
                   rhorime_c(iice)  = 611.+72.25*(Ri-8.)
                END IF

             END IF    !if qc>qsmall

             !!............
             !! rain 
             !!............
             
             !! assume rime density for rain collecting ice is 900 kg/m3
             
             !! IF (qr(i,k).ge.qsmall) THEN
             !!    D_r = (diag_mu_r(i,k)+1.)/lamr(i,k)
             !!    V_impact  = abs(vtrmi1(i,k)-Vt_qr(i,k))
             !!    Ri        = -(0.5e+6*D_r)*V_impact*iTc
             !!    Ri        = max(1.,min(Ri,8.))
             !!    rhorime_r(iice)  = (0.051 + 0.114*Ri - 0.0055*Ri*Ri)*1000.
             !! ELSE
             !!    rhorime_r(iice) = 400.
             !! END IF

          ELSE
             rhorime_c(iice) = 400._r8
          END IF ! qi .gt. qsmall and T .lt. 273.15_r8


       END DO !! iice_loop1
       
       !!......................................................................
       !! end loop over ice categories 
       !!......................................................................


       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 015 -'
       
       !!...................................
       !! contact and immersion freezing droplets 
       !!...................................

       IF(p3_opt_cheti.eq.1) THEN 

       !! contact freezing currently turned off
       !!         dum=7.37*t(i,k)/(288.*10.*pres(i,k))/100.
       !!         dap=4.*pi*1.38e-23*t(i,k)*(1.+dum/rin)/ &
       !!                (6.*pi*rin*mu)
       !!         nacnt=exp(-2.80+0.262*(273.15_r8-tt(i,k)))*1000.

          l_dum1 = qc(i,k).ge.qsmall
          l_dum2 = tt(i,k).le.269.15
             
          l_possible = l_dum1 .and. l_dum2 
             
          IF (l_possible) THEN

             !! subgrid variability
             
             IF (.not. microp_uniform) THEN
                prc_coef = var_coef(ql_relvar(i,k), 2)
             ELSE
                prc_coef = 1._r8
             END IF
  
             dum = (1._r8/diag_lamc(i,k))**3
             
             xdum1 = gamma(7._r8+diag_mu_c(i,k)) 
             xdum2 = exp( aimm*(zerodegc-tt(i,k)) ) 
             xdum3 = gamma(diag_mu_c(i,k)+4._r8)
             xdum4 = exp( aimm*(zerodegc-tt(i,k)) ) 
             
             Q_nuc = prc_coef * cons6*cdist1(i,k)*xdum1*xdum2*dum**2
             N_nuc = prc_coef * cons5*cdist1(i,k)*xdum3*xdum4*dum 

             !! determine destination ice-phase category

             dum1 = 900._r8 !! density of new ice 

             IF(N_nuc.gt.dum_min) THEN
                D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
             END IF

!!!             CALL icecat_destination(qitot(i,k,:),    &
!!!                                     diam_ice(i,k,:), &
!!!                                     D_new,           &
!!!                                     deltaD_init,     &
!!!                                     log_ni_add,      &
!!!                                     iice_dest        )

             !!IF(l_debug) WRITE(6,*) 'micro_p3_tend - 016 - iice_dest 2584 : ',
             !iice_dest 

             iice_dest = 1

             qcheti(iice_dest) = Q_nuc

!!!          IF (log_ni_add .and. Q_nuc.ge. dum_min) THEN
             IF (Q_nuc.ge. dum_min) THEN
                ncheti(iice_dest) = N_nuc
             END IF

          END IF

       END IF


       IF(p3_opt_cheti.eq.2) THEN 

          !! contact freezing currently turned off
          !!         dum=7.37*t(i,k)/(288.*10.*pres(i,k))/100.
          !!         dap=4.*pi*1.38e-23*t(i,k)*(1.+dum/rin)/ &
          !!                (6.*pi*rin*mu)
          !!         nacnt=exp(-2.80+0.262*(273.15_r8-tt(i,k)))*1000.
   
          l_dum1 = qc(i,k).ge.qsmall
          l_dum2 = tt(i,k).le.269.15
   
          l_possible = l_dum1 .and. l_dum2
   
          IF (l_possible) THEN
             
   !!!write(6,*) 'qc(i,k), nc(i,k), frzimm(i,k), frzcnt(i,k), frzdep(i,k) : ', &
   !!!qc(i,k), nc(i,k), frzimm(i,k), frzcnt(i,k), frzdep(i,k)  
   !!!
            ! Mass of droplets frozen is the average droplet mass, except
            ! with two limiters: concentration must be at least 1/cm^3, and
            ! mass must be at least the minimum defined above.
            mi0l = qc(i,k)/max(nc(i,k), 1.0e6_r8/rho(i,k))
            mi0l = max(mi0l_min, mi0l)
   
            N_nuc = (frzimm(i,k)+frzcnt(i,k)+frzdep(i,k))*1.0e6_r8/rho(i,k) !! BUGFIX
            Q_nuc = (frzimm(i,k)+frzcnt(i,k))*1.0e6_r8/rho(i,k)*mi0l + &    !! mi0l
                     frzdep(i,k)*1.0e6_r8/rho(i,k)*mi0                      !! mi0 
   
            if(N_nuc*dt.gt.nc(i,k)) then
                write(6,*) '#### N_nuc*dt.gt.nc(i,k) ' 
            end if  
            if(Q_nuc*dt.gt.qc(i,k)) then
                write(6,*) '#### Q_nuc*dt.gt.qc(i,k) ' 
            end if  

            !! determine destination ice-phase category
            
            dum1 = 900._r8 !! density of new ice 
            
            IF(N_nuc.gt.dum_min) THEN
               D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd 
            END IF 
            
!!!            CALL icecat_destination(qitot(i,k,:),    &
!!!                                    diam_ice(i,k,:), &
!!!                                    D_new,           &
!!!                                    deltaD_init,     &
!!!                                    log_ni_add,      &
!!!                                    iice_dest        )
   
            !!IF(l_debug) WRITE(6,*) 'micro_p3_tend - 016 - iice_dest 2584 : ', iice_dest 
  
            iice_dest = 1 
 
            qcheti(iice_dest) = Q_nuc
            
!!!         IF (log_ni_add .and. Q_nuc.ge. dum_min) THEN
            IF (Q_nuc.ge. dum_min) THEN
               ncheti(iice_dest) = N_nuc
            END IF
   
          END IF

       END IF 

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 016' 


       !!...................................
       !! immersion freezing of rain
       !!...................................

       IF(l_rheti) THEN 

       !! for future: get rid of log statements below for rain freezing

       l_dum1 = qr(i,k).ge.qsmall
       l_dum2 = tt(i,k).le.269.15
          
       l_possible = l_dum1 .and. l_dum2 
       
       IF (l_possible) THEN
       
          xdum1 = log(cdistr(i,k)) 
          xdum2 = log(gamma(7._r8+diag_mu_r(i,k)))
          xdum3 = log(lamr(i,k))
          xdum4 = aimm*(zerodegc-tt(i,k)) 
          
          Q_nuc = cons6*exp(xdum1+xdum2-6._r8*xdum3)*exp(xdum4)
          
          xdum1 = log(cdistr(i,k)) 
          xdum2 = log(gamma(diag_mu_r(i,k)+4._r8))
          xdum3 = log(lamr(i,k))
          xdum4 = aimm*(zerodegc-tt(i,k))
          
          N_nuc = cons5*exp(xdum1+xdum2-3._r8*xdum3)*exp(xdum4)
                  
          !!--determine destination ice-phase category
          
          dum1 = 900._r8 !! density of new ice
          
          IF(N_nuc.gt.dum_min) THEN
             D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd 
          END IF 
          
!!KZ          CALL icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,        &
!!KZ                                  log_ni_add,iice_dest)


          !!IF(l_debug) WRITE(6,*) 'micro_p3_tend - 017 - iice_dest 2609 : ', iice_dest 

          iice_dest = 1

          qrheti(iice_dest) = Q_nuc
          
!!!       IF (log_ni_add .and. Q_nuc .ge. dum_min) THEN
          IF (Q_nuc .ge. dum_min) THEN
             nrheti(iice_dest) = N_nuc
          END IF 
          
       END IF

       END IF 

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 017 '


       !!...................................
       !! rime splintering (Hallet-Mossop 1974) 
       !!...................................
       
       IF(l_rmul) THEN 

       IF (log_hmossopOn) THEN

          !! assumes ice crystals from rime splintering are tiny
          
          D_new = 10.e-6_r8 
          
          !! determine destination ice-phase category
        
!!          CALL icecat_destination(qitot(i,k,:),    &
!!                                  diam_ice(i,k,:), &
!!                                  D_new,           & 
!!                                  deltaD_init,     &
!!                                  log_ni_add,      &
!!                                  iice_dest        )

          iice_dest = 1 

          !!IF(l_debug) WRITE(6,*) 'micro_p3_tend - 018 - iice_dest 2629 : ', iice_dest 

          DO iice = 1,nCat

             l_dum1 = qitot(i,k,iice).ge.qsmall
             l_dum2 = qccol(iice) .gt. 0._r8
             l_dum3 = qrcol(iice) .gt. 0._r8
          
             l_possible = l_dum1 .and. (l_dum2 .or. l_dum3)  
       
             IF(l_possible) THEN

                IF (tt(i,k).gt.270.15) THEN
                   dum = 0.
                ELSEIF (tt(i,k).le.270.15 .and. tt(i,k).gt.268.15) THEN
                   dum = (270.15-tt(i,k))*0.5
                ELSEIF (tt(i,k).le.268.15 .and. tt(i,k).ge.265.15) THEN
                   dum = (tt(i,k)-265.15)*thrd
                ELSEIF (tt(i,k).lt.265.15) THEN
                   dum = 0.
                END IF

                !! rime splintering from riming of cloud droplets
                !! dum1 = 35.e+4*qccol(iice)*dum*1000._r8 ! 1000 is to convert kg to g
                !! dum2 = dum1*piov6*900.*(10.e-6)**3  ! assume 10 micron splinters
                !! qccol(iice) = qccol(iice)-dum2 ! subtract splintering from rime mass transfer
                !! IF (qccol(iice) .lt. 0._r8) THEN
                !!    dum2 = qccol(iice)
                !!    qccol(iice) = 0.
                !! END IF
                !! qcmul(iice_dest) = qcmul(iice_dest)+dum2
                !! IF (log_ni_add) THEN
                !! nimul(iice_dest) = nimul(iice_dest)+dum2/(piov6*900.*(10.e-6)**3)
                !! END IF

                !! rime splintering from riming of raindrops
               
                dum1 = 35.e+4_r8*qrcol(iice)*dum*1000._r8 ! 1000 is to convert kg to g
                dum2 = dum1*piov6*900._r8*(10.e-6_r8)**3  ! assume 10 micron splinters
                
                qrcol(iice) = qrcol(iice)-dum2      ! subtract splintering from rime mass transfer
                
                IF (qrcol(iice) .lt. 0._r8) THEN
                   dum2 = qrcol(iice)
                   qrcol(iice) = 0._r8
                END IF

                qrmul(iice_dest) = qrmul(iice_dest) + dum2
                
!!!             IF (log_ni_add) THEN
                   nimul(iice_dest) = nimul(iice_dest) + dum2/(piov6*900._r8*(10.e-6_r8)**3)
!!!                END IF 
             END IF

          END DO !! iice_loop_HM

       END IF !! rimesplintering_on

       END IF 

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 018 ' 

       !!...................................
       !! condensation/evaporation/deposition/sublimation 
       !!...................................

       !! (use semi-analytic formulation)
       !! 
       !! IF the bulk rain fall speed is to be used to calculate rhorime_r, this section
       !! could be relocated to above
       !! 
       !! calculate rain evaporation including ventilation
       
       IF (qr(i,k).ge.qsmall) THEN

          CALL find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3, &
                                          diag_mu_r(i,k), &
                                          lamr(i,k)       )
          
          !! interpolate value at diag_mu_r
         
          xdum1 = revap_table(dumii,dumjj) 
          xdum2 = (rdumii-real(dumii)) * inv_dum3 
          xdum3 = revap_table(dumii+1,dumjj)
          xdum4 = revap_table(dumii,dumjj)   !! bugfix 
          
          dum1 = xdum1 + xdum2 * (xdum3-xdum4) 
                 
          !! interoplate value at diag_mu_r+1 
         
          xdum1 = revap_table(dumii,dumjj+1)
          xdum2 = (rdumii-real(dumii)) * inv_dum3
          xdum3 = revap_table(dumii+1,dumjj+1) 
          xdum4 = revap_table(dumii,dumjj+1) !! bugfix 
          
          dum2 = xdum1 + xdum2 * (xdum3-xdum4) 
                 
          !! final interpolation
         
          dum  = dum1+(rdumjj-real(dumjj))*(dum2-dum1)


          !!........................................................................
          !! supersaturation relaxation time scales for qr
          !!........................................................................  
       
          xdum1 = cdistr(i,k)*rho(i,k)*dv 
          xdum2 = gamma(diag_mu_r(i,k)+2._r8)
          xdum3 = f1r*xdum2/lamr(i,k) 
          xdum4 = f2r*(rho(i,k)/mu)**0.5*sc**thrd*dum
          
          epsr = 2._r8*pi*xdum1*( xdum3 + xdum4 ) 

       ELSE
          epsr = 0._r8 
       END IF

       diag_epsr(i,k) = epsr !! ### 
       
       !!........................................................................
       !! supersaturation relaxation time scales for qc
       !!........................................................................
       
       IF (qc(i,k).ge.qsmall) THEN
       !! dv = 8.794e-5*tt(i,k)**1.81/pres(i,k)
          epsc = 2.*pi*rho(i,k)*dv*cdist(i,k)
       ELSE
          epsc = 0.
       END IF

       diag_epsc(i,k) = epsc !! ### 


       !! abi = 1.+dqsidt*xxls(i,k)*inv_cp ;; see (C3) 
       !!
       !! C2 term - M2015 App C. 
       !! 1 / tau
       !!  
       
       IF (tt(i,k).lt.zerodegc) THEN
          oabi = 1./abi
          xx   = epsc + epsr + epsi_tot*(1.+xxls(i,k)*inv_cp*dqsdt)*oabi
       ELSE
          xx   = epsc + epsr
       END IF


       dumqvi = qvi(i,k)   !no modification due to latent heating

       !! !      ! modify due to latent heating from riming rate
       !! !      !   - currently this is done by simple linear interpolation
       !! !      !     between conditions for dry and wet growth --> in wet growth it is assumed
       !! !      !     that particle surface temperature is at 0 C and saturation vapor pressure
       !! !      !     is that with respect to liquid. This simple treatment could be improved in the future.
       !! !        IF (qwgrth(iice).ge.1.e-20) THEN
       !! !           dum = (qccol(iice)+qrcol(iice))/qwgrth(iice)
       !! !        ELSE
       !! !           dum = 0.
       !! !        END IF
       !! !        dumqvi = qvi(i,k) + dum*(qvs(i,k)-qvi(i,k))
       !! !        dumqvi = min(qvs(i,k),dumqvi)

       !! 'A' term including ice (Bergeron process)
       !! note: qv and T tendencies due to mixing and radiation are
       !! currently neglected --> assumed to be much smaller than cooling
       !! due to vertical motion which IS included

       xdum1 = pres(i,k) - polysvp1(tt(i,k),0) 
       xdum2 = max(1.e-3_r8,xdum1) 
          
       dum = qvs(i,k)*rho(i,k)*g*uzpl(i,k)/xdum2

       !! abi    = 1.+dqsidt*xxls(i,k)*inv_cp ;; see (C3) 
       !! oabi   = 1. / abi
       !!
       !! Ac (equation C4) 
       
       xdum1 = -uzpl(i,k)*g*inv_cp
       xdum2 = qvs(i,k)-dumqvi 
       xdum3 = 1._r8 + xxls(i,k)*inv_cp*dqsdt 
          
       IF (tt(i,k).lt.zerodegc) THEN
          aaa = -dum - dqsdt*xdum1 - xdum2*xdum3*oabi*epsi_tot * scale_berg

          diag_ac1(i,k) = - dum 
          diag_ac2(i,k) = - dqsdt*xdum1 
          diag_ac3(i,k) = - xdum2*xdum3*oabi*epsi_tot !! ### 

       ELSE
          aaa = -dum - dqsdt*xdum1 
       END IF

       xx = max(1.e-8_r8,xx)   !! set lower bound on xx to prevent division by zero
       
       !! tau 

       oxx = 1._r8/xx

       diag_reps(i,k) = oxx !! ### 

       xdum1 = ssat(i,k)-aaa*oxx
       
       xdum2 = 1._r8 - dexp( -dble(xx*dt) ) 
       xdum3 = 1._r8 - dexp( -dble(xx*dt) ) 
      

       !!...................................
       !! CLUBB handles cond/evap, here only for testing
       !!...................................
 
       IF(l_crconevp) THEN 

       IF (qc(i,k).ge.qsmall) THEN
         qccon = ( aaa*epsc*oxx + xdum1*odt*epsc*oxx*xdum2 ) / ab 
         diag_qccon1(i,k) = aaa*epsc*oxx / ab 
         diag_qccon2(i,k) = xdum1*odt*epsc*oxx*xdum2 / ab !! ### 
       END IF 
       
       IF (qr(i,k).ge.qsmall) THEN 
         qrcon = ( aaa*epsr*oxx + xdum1*odt*epsr*oxx*xdum3 ) / ab 
       END IF 
       
       !! no condensation; this is handled by macrophysics 
       
       !! before v0.33 qccon = 0._r8 
       !! before v0.33 qrcon = 0._r8 

       !!...................................
       !! condensation/evaporation 
       !!...................................
       
       !! for very small water contents, evaporate instantly
        
       IF (sup(i,k).lt.-0.001 .and. qc(i,k).lt.1.e-12) THEN
          qccon = -qc(i,k)*odt
       END IF
       
       IF (sup(i,k).lt.-0.001 .and. qr(i,k).lt.1.e-12) THEN 
          qrcon = -qr(i,k)*odt
       END IF
       
       IF (qccon .lt. 0._r8) THEN
          qcevp = -qccon
          qccon = 0.
       END IF
       
       IF (qrcon .lt. 0._r8) THEN
          qrevp = -qrcon
          qrcon = 0.
       END IF

       IF (ssat(i,k) .lt. 0._r8) THEN 
          qccon = min(0._r8,qccon)
          qrcon = min(0._r8,qrcon)
       END IF

       !! limit total condensation/evaporation to saturation adjustment
       !! 
       !! xdum1 > 0 supersaturation 
       !! xdum1 < 0 subsaturation 
      
       xdum1  = (qv(i,k)-qvs(i,k))/(1.+xxlv(i,k)**2*qvs(i,k)/(cp*rv*tt(i,k)**2))*odt
       
       if (qccon+qrcon.gt.0._r8) then
          ratio = max(0._r8,xdum1)/(qccon+qrcon)
          ratio = min(1._r8,ratio)
          qccon = qccon*ratio
          qrcon = qrcon*ratio
       elseif (qcevp+qrevp.gt.0._r8) then
          ratio = max(0._r8,-xdum1)/(qcevp+qrevp)
          ratio = min(1._r8,ratio)
          qcevp = qcevp*ratio
          qrevp = qrevp*ratio
       endif
       
       !!qcevp = 0._r8 !! v0.49a
       
       IF (qrevp .gt. 0._r8 .and. qr(i,k).ge.qsmall) THEN
           dum1  = -qrevp*dt/qr(i,k)
           dum1  = max(-1.,dum1)
           !dum2  = exp(-0.2*diag_mu_r(i,k)) ! mu dependence from Seifert (2008), neglecting size dependence
           dum2  = 1.                  ! or, neglect mu dependence
           nrevp = dum2*dum1*nr(i,k)*odt
       END IF

       END IF !! IF(l_crconevp) THEN 


       !!...................................
       !! deposition/sublimation 
       !!...................................

       IF(l_idepsub) THEN 

       !! note that qidep is in-cloud (calculated based in-cloud qi), this is different 
       !! from what is in mg2 
       
       DO iice = 1,nCat


          !!  IF (qitot(i,k,iice).ge.qsmall.and.tt(i,k).lt.273.15_r8) THEN
          !!     qidep(iice) = (aaa*epsi(iice)*oxx+(ssat(i,k)-aaa*oxx)*odt*epsi(iice)*oxx*   &
          !!                   (1.-dexp(-dble(xx*dt))))*oabi+(qvs(i,k)-dumqvi)*epsi(iice)*oabi
          !!  END IF
          !!
          !!  IF (qv(i,k)-qvi(i,k) .lt. 0._r8) qidep(iice) = min(0._r8,qidep(iice))
          !!
          !! !for very small water contents, evaporate instantly
          !!  IF (supi(i,k).lt.-0.001 .and. qitot(i,k,iice).lt.1.e-12) &
          !!     qidep(iice) = -qitot(i,k,iice)*odt
          !!
          !!  IF (qidep(iice) .lt. 0._r8) THEN
          !!     qisub(iice) = -qidep(iice)
          !!     qidep(iice) = 0.
          !!  END IF
          !!
          !!  IF (qisub(iice) .gt. 0._r8 .and. qitot(i,k,iice).ge.qsmall) THEN
          !!      dum         = -qisub(iice)*dt/qitot(i,k,iice)
          !!      dum         = max(-1.,dum)
          !!      nisub(iice) = dum*nitot(i,k,iice)*odt
          !!  END IF

          IF (qitot(i,k,iice).ge.qsmall) THEN 
          
             IF( tt(i,k).lt.zerodegc) THEN
             
                !! 
                !! 
                !! aaa : Ac - is the change in d due to vertical motion, 
                !!            turbulent mixing, radiation, and the Bergeron–Findeisen process
                !! epsi: 1/tau_i 
                !! oabi: 1/gamma_i 
                !! oxx : tau 
                !! ssat: supersaturation 
                !! 
                !!qidep(iice) = (aaa*epsi(iice)*oxx+(ssat(i,k)-aaa*oxx)*odt*epsi(iice)*oxx*   &
                !!              (1.-dexp(-dble(xx*dt))))*oabi+(qvs(i,k)-dumqvi)*epsi(iice)*oabi 

                if (.not. l_mg2_qidep) then 

                   xdum1 = aaa*epsi(iice)*oxx
                   xdum2 = ( ssat(i,k) - aaa*oxx ) * odt * epsi(iice) * oxx * &
                             ( 1.-dexp(-dble(xx*dt)) )
                   xdum3 = (qvs(i,k)-dumqvi)*epsi(iice) * scale_berg
                   
                   qidep(iice) = ( xdum1 + xdum2 + xdum3 )*oabi 
                
!!!                   qberg = max((qvs(i,k) - qvi(i,k))*epsi(iice)*oabi*scale_berg, 0._r8)

!!!                   diag_qicon1(i,k) = xdum1 
!!!                   diag_qicon2(i,k) = xdum2
!!!                   diag_qicon3(i,k) = xdum3 !! ### 
                   
                else

                   qidep(iice) = ( qv(i,k) - qvi(i,k) )*epsi(iice)*oabi 

                   if(qc(i,k).lt.qsmall) then  

                      !!................................................
                      !! only ice exists, no bergeron process  
                      !!................................................
                      qberg = 0._r8

                   else

                      !!................................................
                      !! both ice and liquid exists
                      !!................................................

                      if(qv(i,k).gt.qvi(i,k)) then 
                         !!................................................
                         !! case #2.1 qv > qvi, bergeron process happens
                         !!................................................
                         qberg = max((qvs(i,k) - qvi(i,k))*epsi(iice)*oabi*scale_berg, 0._r8)
                         dum1  = max(0._r8,qc(i,k)/(qberg*dt) ) 
                         if(dum1.lt.1._r8) qberg = max(0._r8,qc(i,k)/dt)
                      else
                         qberg = 0._r8
                      end if

                   end if 

                   qberg = qberg * scale_berg  

                end if 

             ELSE 
                qidep(iice) = 0._r8 
                qberg       = 0._r8
             END IF




!!!          if (qitot(i,k,iice).ge.qsmall.and.tt(i,k).lt.273.15) then
!!!             qidep(iice) = (aaa*epsi(iice)*oxx+(ssat(i,k)-aaa*oxx)*odt*epsi(iice)*oxx*   &
!!!                           (1.-dexp(-dble(xx*dt))))*oabi+(qvs(i,k)-qvi(i,k))*epsi(iice)*oabi
!!!          endif
       
!!!         !for very small ice contents in dry air, sublimate all ice instantly
!!!          if (supi(i,k).lt.-0.001 .and. qitot(i,k,iice).lt.1.e-12) &
!!!             qidep(iice) = -qitot(i,k,iice)*odt
!!!
!!!          if (qidep(iice).lt.0.) then
!!!           !note: limit to saturation adjustment (for dep and subl) is applied
!!!           !later
!!!             qisub(iice) = -qidep(iice)
!!!             qisub(iice) = min(qisub(iice), qitot(i,k,iice)*dt)
!!!             nisub(iice) = qisub(iice)*(nitot(i,k,iice)/qitot(i,k,iice))
!!!             qidep(iice) = 0.
!!!          else
!!!             qidep(iice) = qidep(iice)
!!!          endif

             !! IF sub-saturated, qidep to be negative 
             
             IF ( (qv(i,k)-qvi(i,k)) .lt. 0._r8) THEN 
                qidep(iice) = min(0._r8, qidep(iice))
                !!IF(masterproc) WRITE(6,*) 'qidep: #1  sub-saturated, qidep to be negative ', qidep(iice) 
             END IF
       
             !! for very small water contents, evaporate instantly
             
             IF (supi(i,k)       .lt. -0.001 .and. &
                 qitot(i,k,iice) .lt. 1.e-12       ) THEN 
                 qidep(iice) = -qitot(i,k,iice)*odt 
                 !!IF(masterproc) WRITE(6,*) 'qidep: #2 evaporate very small water contents', qidep(iice) 
             END IF

             !! IF qidep is negative, let it be sublimation  
             
             IF (qidep(iice) .lt. 0._r8) THEN
                qisub(iice) = -qidep(iice)
                qidep(iice) = 0.
                !!IF(masterproc) WRITE(6,*) 'qidep: #3 qidep is negative, let it be sublimation ', qisub(iice) 
             END IF

             IF (qidep(iice) .gt. 0._r8) THEN
                 qidep(iice) = min(qidep(iice),qv(i,k)/dt) !! ### bugfix 
                 qisub(iice) = 0._r8
             END IF

             IF (qisub(iice) .gt. 0._r8) THEN
                 qisub(iice) = min(qisub(iice),qitot(i,k,iice)/dt) !! ### bugfix 
             END IF

             IF (qisub(iice)     .gt. 0._r8 .and. &
                 qitot(i,k,iice) .ge. qsmall      ) THEN
                 dum         = -qisub(iice)*dt/qitot(i,k,iice)
                 dum         = max(-1.,dum)
                 nisub(iice) = dum*nitot(i,k,iice)*odt
             END IF
             

          ELSE 
             qisub(iice) = 0._r8 
             qidep(iice) = 0._r8 
             qberg       = 0._r8
          END IF
          
       END DO !! iice_loop_depsub

       END IF

       !!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 019 -'

       !!in wetgrowth regime, add ice deposition to rain condensation and shed
       !!neglect for now, instead assume deposition onto wet growth hail goes to ice
       !!       IF (log_wetgrowth(iice)) THEN
       !!          IF (qidep(iice) .gt. 0._r8) THEN
       !!             !!1/5.2e-7, 5.2e-7 is the mass of a 1 mm raindrop 
       !!             nrshdr(iice) = nrshdr(iice) + qidep(iice)*1.923e+6          
       !!             qrcon  = qrcon + qidep(iice)
       !!             qidep(iice)  = 0.
       !!          END IF
       !!       END IF


444   CONTINUE


!!                                                                   DO k = kbot,ktop,kdir

     !!...................................
     !! droplet activation
     !!...................................

     IF(l_cnuc2) THEN

     !!...................................
     !! 
     !! get provisional droplet number after activation. This is used for
     !! all microphysical process calculations, for consistency with update of
     !! droplet mass before microphysics
     !!
     !! calculate potential for droplet activation IF cloud water is present
     !! npccn (activation tendency) is from microp_aero 
     !!
     !! npccn is grid-mean value 
     !!...................................

     l_dum1 = qc(i,k) .ge. qsmall 
     
     l_dum2 = .not. is_first_step 
     
     l_dum3 = lcldm(i,k) .ge. cldm_min
     
     l_possible = l_dum1 .and. l_dum2 .and. l_dum3 
     
     IF (l_possible) THEN 
        ncnuc = npccn(i,k) / lcldm(i,k)
        qcnuc = ncnuc*cons7
     ELSE
        qcnuc = 0._r8 
        ncnuc = 0._r8 
     END IF

     !! output activated liquid and ice (convert from #/kg -> #/m3)
     
     l_dum1 = qc(i,k) .ge. qsmall
     
     l_dum2 = lcldm(i,k) .ge. cldm_min 
     
     l_possible = l_dum1 .and. l_dum2 
     
     IF (l_possible) THEN 
        diag_ncal(i,k) = max(nc(i,k) + npccn(i,k)/lcldm(i,k)*dt, 0._r8) *rho(i,k)
     ELSE
        diag_ncal(i,k) = 0._r8
     END IF


     END IF 


     !!...................................
     !! ice nucleation
     !!...................................

     IF(p3_opt_inuc.eq.1) THEN

        IF (tt(i,k) .lt. icenuct) THEN 
           diag_ncai(i,k) = naai(i,k)*rho(i,k) ! already in-cloud 
        ELSE
           diag_ncai(i,k) = 0._r8
        END IF
        
        
        !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 020 -'

        !! nucleation time scale 
        
        mtime = dt 

        !!...................................
        !!
        !! IF activated nuclei exist at t<-5C and rhmini + 5% then ice nucleates 
        !!
        !!...................................
 
        l_dum1 = naai(i,k) .gt. 0._r8 
        
        l_dum2 = tt(i,k) .lt. icenuct
       
        l_dum3 = (supi(i,k) + 1._r8 ).ge. rhmini+0.05_r8 !! bugfix ###  
        l_dum4 = (supi(i,k) ).ge. 0.05_r8 !! bugfix ###  
        !!l_dum3 = supi(i,k) .ge. 0.05_r8 !! p3

!!if(l_dum4) then 
!!write(6,*) 'DUM4 supi(i,k), rhmini, tt(i,k), naai(i,k) : ', supi(i,k), rhmini, tt(i,k), naai(i,k) 
!!end if 
 
        l_possible = l_dum1 .and. l_dum2 .and. l_dum3
        
        IF (l_possible) THEN 

!!write(6,*) 'DUM3 supi(i,k), rhmini, tt(i,k), naai(i,k) : ', supi(i,k), rhmini, tt(i,k), naai(i,k) 

           !if NAAI .gt. 0._r8 then set numice = naai (as before)
           !note: this is gridbox averaged
           
           N_nuc = (naai(i,k) - sum(nitot(i,k,:)) ) / mtime
           Q_nuc = (naai(i,k) - sum(nitot(i,k,:)) ) * mi0 / mtime 
           
           N_nuc = max(0._r8,N_nuc) 
           Q_nuc = max(0._r8,Q_nuc) 
        
!!!           diag_qccon1(i,k) = sum(nitot(i,k,:)) 
!!!           diag_qccon2(i,k) = N_nuc
!!!           diag_qicon1(i,k) = Q_nuc
!!!           diag_qicon2(i,k) = supi(i,k) + 1._r8 

   
           IF (N_nuc.ge.1.e-20) THEN
           
              !--determine destination ice-phase category:
              
              dum1 = 900._r8 !density of new ice
              
              D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
 
!!KZ              CALL icecat_destination(qitot(i,k,:),      &
!!KZ                                      diam_ice(i,k,:),   &
!!KZ                                      D_new,             &
!!KZ                                      deltaD_init,       &
!!KZ                                      log_ni_add,        &
!!KZ                                      iice_dest)
                                      
              !!IF(l_debug) WRITE(6,*) 'micro_p3_tend - 2932 - iice_dest : ', iice_dest 

              iice_dest = 1

              !! bugfix 2017-12-17 
!!!           IF (log_ni_add) THEN  
                 ninuc(iice_dest) = N_nuc 
!!!           END IF 
           
              qinuc(iice_dest) = Q_nuc
                              
           END IF 
           
           nimax(i,k) = naai(i,k)*icldm(i,k)
           
           
        END IF

     !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 021' 

     END IF


     IF(p3_opt_inuc.eq.2) THEN

     !!................................................................
     !! deposition/condensation-freezing nucleation
     !! allow ice nucleation IF .lt. -15 C and .gt. 5% ice supersaturation
     
        IF (tt(i,k).lt.258.15 .and. supi(i,k).ge.0.05) THEN
     
     !!    dum = exp(-0.639+0.1296*100.*supi(i,k))*1000.*inv_rho(i,k)  !Meyers et al. (1992)
           dum = 0.005*exp(0.304*(273.15_r8-tt(i,k)))*1000.*inv_rho(i,k)   !Cooper (1986)
           dum = min(dum,100.e3*inv_rho(i,k))
           N_nuc = max(0._r8,(dum-sum(nitot(i,k,:)))*odt)
     
           IF (N_nuc.ge.1.e-20) THEN
              Q_nuc = max(0._r8,(dum-sum(nitot(i,k,:)))*mi0*odt)
              !--determine destination ice-phase category:
              dum1      = 900._r8     !density of new ice
              D_new     = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
              CALL icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,    &
                                      log_ni_add,iice_dest)
              !==
              qinuc(iice_dest) = Q_nuc
              IF (log_ni_add) ninuc(iice_dest) = N_nuc
           END IF
     
        END IF
     
     END IF

     !!!
     !!!!.................................................................
     !!!! droplet activation
     !!!
     !!!! for specified Nc, make sure droplets are present IF conditions are supersaturated
     !!!! this is not applied at the first time step, since saturation adjustment is applied at the first step
     !!!
     !!!          IF (.not.(log_predictNc).and.sup(i,k).gt.1.e-6.and.it.gt.1) THEN
     !!!             dum   = nccnst*inv_rho(i,k)*cons7-qc(i,k)
     !!!             dum   = max(0._r8,dum)
     !!!             dqsdt = xxlv(i,k)*qvs(i,k)/(rv*tt(i,k)*tt(i,k))
     !!!             ab    = 1. + dqsdt*xxlv(i,k)*inv_cp
     !!!             dum   = min(dum,(qv(i,k)-qvs(i,k))/ab)  ! prevent overdepletion of supersaturation
     !!!             qcnuc = dum*odt
     !!!          END IF
     !!!
     !!!       predict_supersaturation2: IF (log_predictSsat) THEN
     !!!
     !!!! for predicted Nc (and therefore predicted Ssat), calculate activation explicitly from supersaturation
     !!!! note that this is also applied at the first time step.
     !!!
     !!!          IF (sup(i,k).gt.1.e-6) THEN
     !!!             dum1  = 1./bact**0.5
     !!!             sigvl = 0.0761 - 1.55e-4*(tt(i,k)-273.15_r8)
     !!!             aact  = 2.*mw/(rhow*rr*tt(i,k))*sigvl
     !!!             sm1   = 2.*dum1*(aact*thrd*inv_rm1)**1.5
     !!!             sm2   = 2.*dum1*(aact*thrd*inv_rm2)**1.5
     !!!             uu1   = 2.*log(sm1/sup(i,k))/(4.242*log(sig1))
     !!!             uu2   = 2.*log(sm2/sup(i,k))/(4.242*log(sig2))
     !!!             dum1  = nanew1*0.5*(1.-derf(uu1)) ! activated number in kg-1 mode 1
     !!!             dum2  = nanew2*0.5*(1.-derf(uu2)) ! activated number in kg-1 mode 2
     !!!           ! make sure this value is not greater than total number of aerosol
     !!!             dum2  = min((nanew1+nanew2),dum1+dum2)
     !!!             dum2  = (dum2-nc(i,k))*odt
     !!!             dum2  = max(0._r8,dum2)
     !!!             ncnuc = dum2
     !!!           ! don't include mass increase from droplet activation during first time step
     !!!           ! since this is already accounted for by saturation adjustment below
     !!!             IF (is_first_step) THEN
     !!!                qcnuc = 0.
     !!!             ELSE
     !!!                qcnuc = ncnuc*cons7
     !!!             END IF
     !!!          END IF
     !!!
     !!!       END IF predict_supersaturation2


!!                                                                   DO k = kbot,ktop,kdir

!!######################################################################
!!#
!!#
!!# warm cloud microphysical processes 
!!#
!!#
!!######################################################################

       !!...................................
       !! saturation adjustment to get initial cloud water
       !!...................................

       !! This is only called once at the beginning of the simulation
       !! to remove any supersaturation in the intial conditions

       IF (is_first_step) THEN
          dumt   = tt(i,k)
          dumqv  = qv(i,k)
          dumqvs = ep_2*polysvp1(dumt,0)/max(1.e-3,(pres(i,k)-polysvp1(dumt,0)))
          dums   = dumqv-dumqvs
          qccon  = dums/(1.+xxlv(i,k)**2*dumqvs/(cp*rv*dumt**2))*odt
          qccon  = max(0._r8,qccon)
          IF (qccon.le.1.e-7) qccon = 0.
       END IF

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 022 -'


       !!...................................
       !! autoconversion
       !!...................................


       IF(l_caut) THEN 
     
       l_dum1 = qc(i,k).ge.1.e-8
       l_dum2 = nc(i,k).ge.1.e-8
       l_possible = l_dum1 .and. l_dum2
     
       IF (l_possible) THEN

!!!KZ          IF (iparam.eq.1) THEN
!!!KZ
!!!KZ             !! Seifert and Beheng (2001)
!!!KZ            
!!!KZ             dum   = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
!!!KZ             dum1  = 600.*dum**0.68*(1.-dum**0.68)**3
!!!KZ             qcaut =  kc*1.9230769e-5*(nu(i,k)+2.)*(nu(i,k)+4.)/(nu(i,k)+1.)**2* &
!!!KZ                      (rho(i,k)*qc(i,k)*1.e-3)**4/(rho(i,k)*nc(i,k)*1.e-6)**2*(1.+ &
!!!KZ                      dum1/(1.-dum)**2)*1000.*inv_rho(i,k)
!!!KZ                      
!!!KZ             ncautc = qcaut*7.6923076e+9
!!!KZ
!!!KZ          ELSEIF (iparam.eq.2) THEN
!!!KZ
!!!KZ             !! Beheng (1994)
!!!KZ            
!!!KZ             IF (nc(i,k)*rho(i,k)*1.e-6 .lt. 100._r8) THEN
!!!KZ                qcaut = 6.e+28*inv_rho(i,k)*diag_mu_c(i,k)**(-1.7)*(1.e-6*rho(i,k)* &
!!!KZ                        nc(i,k))**(-3.3)*(1.e-3*rho(i,k)*qc(i,k))**4.7
!!!KZ             ELSE
!!!KZ               !2D interpolation of tabled logarithmic values
!!!KZ                dum   = 41.46 + (nc(i,k)*1.e-6*rho(i,k)-100._r8)*(37.53-41.46)*5.e-3
!!!KZ                dum1  = 39.36 + (nc(i,k)*1.e-6*rho(i,k)-100._r8)*(30.72-39.36)*5.e-3
!!!KZ                qcaut = dum+(diag_mu_c(i,k)-5.)*(dum1-dum)*0.1
!!!KZ              ! 1000/rho is for conversion from g cm-3/s to kg/kg
!!!KZ                qcaut = exp(qcaut)*(1.e-3*rho(i,k)*qc(i,k))**4.7*1000.*inv_rho(i,k)
!!!KZ             END IF
!!!KZ             ncautc = 7.7e+9*qcaut
!!!KZ
!!!KZ          ELSEIF (iparam.eq.3) THEN

             !! Khroutdinov and Kogan (2000)

             prc_coef1= 1350._r8   !! KK2000 original value 
             prc_exp  =  2.47_r8   !! KK2000 original value 
             prc_exp1 = -1.79_r8   !! KK2000 original value 
             
             prc_coef1= 30500._r8   !! 04P2 original value ### 
             prc_exp  =  3.19_r8    !! 04P2 original value 
             prc_exp1 = -1.20_r8    !! 04P2 original value 

             !! subgrid variability 

             IF (.not. microp_uniform) THEN !! ### subgrid 
                prc_coef = var_coef(ql_relvar(i,k), prc_exp)
             ELSE
                prc_coef = 1._r8
             END IF
             
             qcaut = prc_coef * prc_coef1 &
                     * qc(i,k)**prc_exp & 
                     * (nc(i,k)*1.e-6*rho(i,k))**prc_exp1
             
             ! note: ncautr is change in Nr; ncautc is change in Nc
             
             ncautr = qcaut*cons3 !! cons3 = 1/droplet_mass_25um 
             
             ncautc = qcaut*nc(i,k)/qc(i,k)

             IF (qcaut .eq.0._r8) ncautc = 0.
             IF (ncautc.eq.0._r8) qcaut  = 0.

       END IF !! qc_not_small

       END IF 

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 023 -'



       IF(l_cslf) THEN 

       !!...................................
       !! self-collection of droplets
       !!...................................

       IF (qc(i,k).ge.qsmall) THEN

          IF (iparam.eq.1) THEN
             !! Seifert and Beheng (2001) 
             xdum1 = 1.e-3_r8*rho(i,k)*qc(i,k) 
             xdum2 = (nu(i,k)+2._r8)/(nu(i,k)+1._r8) 
             ncslf = -kc*xdum1**2 * xdum2 * 1.e+6 * inv_rho(i,k) + ncautc
          ELSEIF (iparam.eq.2) THEN
             !! Beheng (994)
             xdum1 = diag_mu_c(i,k)**(-0.63) 
             xdum2 = 1.e-3*rho(i,k)*qc(i,k) 
             ncslf = -5.5e+16 * inv_rho(i,k) * xdum1 * xdum2**2
          ELSEIF (iparam.eq.3) THEN
             !! Khroutdinov and Kogan (2000)
             ncslf = 0._r8 
          END IF

       END IF

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 024 -'

       END IF 



       !!...................................
       !! accretion of cloud by rain
       !!...................................

       IF(l_cacc) THEN 

 
       l_dum1 = qr(i,k).ge.qsmall
     
       l_dum2 = qc(i,k).ge.qsmall
     
       l_possible = l_dum1 .and. l_dum2 
       
       IF (l_possible) THEN

          IF (iparam.eq.1) THEN 
          
             !! Seifert and Beheng (2001)
           
             dum   = 1.-qc(i,k)/(qc(i,k)+qr(i,k))
             dum1  = (dum/(dum+5.e-4))**4
             
             qcacc = kr*rho(i,k)*0.001*qc(i,k)*qr(i,k)*dum1
             ncacc = qcacc*rho(i,k)*0.001*(nc(i,k)*rho(i,k)*1.e-6)/(qc(i,k)*rho(i,k)*   &
                     0.001)*1.e+6*inv_rho(i,k)
                     
          ELSEIF (iparam.eq.2) THEN
          
             !! Beheng (1994) 
             
             qcacc = 6.*rho(i,k)*(qc(i,k)*qr(i,k))
             
             ncacc = qcacc*rho(i,k)*1.e-3*(nc(i,k)*rho(i,k)*1.e-6)/(qc(i,k)*rho(i,k)*1.e-3)* &
                     1.e+6*inv_rho(i,k)
                     
          ELSEIF (iparam.eq.3) THEN
          
             !! Khroutdinov and Kogan (2000) 
             
             pra_coef0 = 1.0_r8 !! KK2000 default 
             pra_coef0 = 1.5_r8 !! 04P2 default 
             
             IF (.not. microp_uniform) THEN !! ### subgrid 
                pra_coef = pra_coef0 * var_coef(ql_relvar(i,k), 1.15_r8) 
             ELSE
                pra_coef = 1._r8
             END IF
  
             qcacc = pra_coef * 67._r8 * ( qc(i,k)*qr(i,k) ) **1.15_r8 
             
             ncacc = qcacc*nc(i,k)/qc(i,k) 
             
          END IF

          IF (qcacc.eq.0._r8) ncacc = 0._r8
          IF (ncacc.eq.0._r8) qcacc = 0._r8

       END IF

       END IF


       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 025 -'



       IF(l_rslf) THEN 

       !!...................................................................................
       !! self-collection and breakup of rain 
       !! breakup following modified Verlinde and Cotton scheme 
       !!...................................................................................

       l_dum1 = qr(i,k).ge.qsmall
     
       l_dum2 = nr(i,k).ge.qsmall
     
       l_possible = l_dum1 .and. l_dum2 
       
       IF (l_possible) THEN

        ! include breakup
          dum1 = 280.e-6

        ! USE mass-mean diameter (do this by using
        ! the old version of lambda w/o mu dependence)
        ! note there should be a factor of 6^(1/3), but we
        ! want to keep breakup threshold consistent so 'dum'
        ! is expressed in terms of lambda rather than mass-mean D

          dum2 = (qr(i,k)/(pi*rhow*nr(i,k)))**thrd
          
          IF (dum2.lt.dum1) THEN
             dum = 1.
          ELSE IF (dum2.ge.dum1) THEN
             dum = 2.-exp(2300.*(dum2-dum1))
          END IF

          IF (iparam.eq.1.) THEN
             nrslf = -dum*kr*1.e-3*qr(i,k)*nr(i,k)*rho(i,k)
          ELSEIF (iparam.eq.2 .or. iparam.eq.3) THEN
             nrslf = -dum*5.78*nr(i,k)*qr(i,k)*rho(i,k)
          END IF

       END IF

       END IF 


!!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 026 -'


!!######################################################################
!!#
!!#
!!# conservation check, if not, do clipping
!!#
!!#
!!######################################################################

if(l_conservation_clipping) then 


       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 027 -'

       !!.................................................................
       !! limit ice deposition/nucleation
       !!.................................................................

       IF(l_limit_qidep_qinuc) THEN

       DO iice = 1,nCat

          xdum1 = qidep(iice)*icldm(i,k)
          xdum2 = qinuc(iice)*icldm(i,k)

          dum1 = xdum1 + xdum2

          IF (dum1 .gt. dum_min .and. icldm(i,k) .ge. cldm_min ) THEN

             !! available vapor for deposition/nucleation 

             dum = (qv(i,k)-qvi(i,k)) / &
                   (1._r8 + xxls(i,k)**2*qvi(i,k) / (cpp*rv*tt(i,k)**2) ) / dt

             dum = max(dum,0._r8)

             IF (dum1 .gt. dum) THEN
                dum2 = xdum2 / dum1 !! ratio 
                qinuc(iice) = dum * dum2 / icldm(i,k)
                qidep(iice) = dum / icldm(i,k)  - qinuc(iice)
             END IF

          END IF

       END DO !! iice 

       END IF !! 


       !!.................................................................
       !! limit ice sublimation and rain evaporation
       !!.................................................................


       IF(l_limit_qisub_qrevp) THEN

       DO iice = 1,nCat

          xdum2 = qisub(iice)*icldm(i,k) !! sublimation q+ t-
          xdum3 = qidep(iice)*icldm(i,k) !! deposition  q- t+
          xdum4 = qinuc(iice)*icldm(i,k) !! nucleation  q- t+
          xdum5 = qrevp*rcldm(i,k)       !! evaporation q+ t-
          xdum6 = qberg*lcldm(i,k)       !! bergeron    t+

          xdum1 = xdum5 + xdum2

          IF (xdum1 .gt. dum_min) THEN !! all postive terms 

             qtmp = qv(i,k) + (xdum2 + xdum5 - xdum3 - xdum4) * dt

             ttmp = tt(i,k) + ( (xdum3 + xdum4 - xdum2)*xxls(i,k) + xdum6*xlf(i,k) - xdum5*xxlv(i,k)) * dt / cpp

             !! modify ice/precip evaporation rate IF q .gt. qsat

             IF (qtmp .gt. qvs(i,k)) THEN

                dum1 = xdum5 / (xdum5+xdum2) !! contribution of rain evaporation 

                dum2 = 0._r8 !! snow evaporation 

                !! recalculate q and t after vap_dep and mnuccd but without evap
                !or sublim

                qtmp = qv(i,k) - (xdum3+xdum4) * dt

                ttmp = tt(i,k) + ( (xdum3+xdum4) * xxls(i,k) + xdum6*xlf(i,k) ) * dt / cpp

                !! sub-saturation (negative)  

                dum = (qtmp-qvs(i,k)) / &
                      ( 1._r8 + xxlv(i,k)**2 * qvs(i,k) / (cpp*rv*ttmp**2) )

                dum = min(dum,0._r8)

                !! evaporation 

                IF(rcldm(i,k) .ge. cldm_min) THEN
                   qrevp = dum * dum1 / rcldm(i,k) / dt
                END IF

                !! sublimation 

                dum1 = 1._r8 - dum1 - dum2

                IF(icldm(i,k) .ge. cldm_min) THEN
                   qisub(iice) = dum * dum1 / icldm(i,k) / dt
                END IF

             END IF

          END IF

       END DO !! iice 

       END IF !! l_limit_qisub_qrevp 



       qcnuc = 0._r8 

       !!.................................................................
       !! avoid negative water vapor 
       !!.................................................................

       IF(qcnuc .gt. 0._r8) THEN 
          dum  = ( qcnuc      * lcldm(i,k) + &
                   sum(qidep) * icldm(i,k) + & 
                   sum(qinuc) * icldm(i,k)   )*dt
       ELSE
          dum  = (-qcnuc      * lcldm(i,k) + &
                   sum(qidep) * icldm(i,k) + & 
                   sum(qinuc) * icldm(i,k)   )*dt
       END IF

       dum1 = qv(i,k) + ( sum(qisub) * icldm(i,k) + &
                          qcevp      * lcldm(i,k) + &
                          qrevp      * rcldm(i,k)    ) * dt 
       
       l_dum1 = dum.gt.dum1
     
       l_dum2 = dum.ge.dum_min 
     
       l_possible = l_dum1 .and. l_dum2 
      
       ratio = 1._r8 
 
       IF (l_possible) THEN

          ratio  = dum1/dum

          !!IF(l_masscon_debug) THEN 
          !!IF(masterproc) WRITE(6,*) '     qcaut ', qcnuc      * lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qidep ', sum(qidep) * icldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qinuc ', sum(qinuc) * icldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qv    ', qv(i,k) 
          !!IF(masterproc) WRITE(6,*) 'ratio: #0  conservation of vapor ', ratio 
          !!END IF 

          qcnuc = qcnuc*ratio

          DO iice = 1,nCat
             qidep(iice) = qidep(iice)*ratio
             qinuc(iice) = qinuc(iice)*ratio
          END DO 
      
       END IF

!!KZ       diag_ac3(i,k) = ratio 

      
       !!.................................................................
       !! conservation of cloud liquid 
       !!.................................................................

       IF(qcnuc .gt. 0._r8) THEN 
          dum1 = qc(i,k) * lcldm(i,k) + &
                 (qccon+qcnuc)*lcldm(i,k)*dt              
       ELSE
          dum1 = qc(i,k) * lcldm(i,k) + &
                 (qccon-qcnuc)*lcldm(i,k)*dt
       END IF

       !! sink 
       dum  = ( qcaut       + &
                qcacc       + &
                sum(qccol)  + &
                qcevp       + &
                qberg       + & !### 
                sum(qchetc) + &
                sum(qcheti) + &
                sum(qcshd)     )*lcldm(i,k)*dt

       l_dum1 = dum.gt.dum1
     
       l_dum2 = dum.ge.1.e-20
     
       l_possible = l_dum1 .and. l_dum2 
       
       ratio = 1._r8 
 
       IF (l_possible) THEN
          ratio  = dum1/dum

          !!IF(l_masscon_debug) THEN 
          !!IF(masterproc) WRITE(6,*) '     qcaut ', qcaut*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcacc ', qcacc*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcevp ', qcevp*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qccol ', sum(qccol)*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qchetc ', sum(qchetc)*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcheti ', sum(qcheti)*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcshd ', sum(qcshd)*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qc    ', qc(i,k)*lcldm(i,k) 
          !!IF(masterproc) WRITE(6,*) 'ratio: #1  conservation of water ', ratio 
          !!END IF 
       
          qcaut  = qcaut*ratio
          qcacc  = qcacc*ratio
          qcevp  = qcevp*ratio
          qccol  = qccol*ratio
          qchetc = qchetc*ratio
          qcheti = qcheti*ratio
          qcshd  = qcshd*ratio
          qberg  = qberg*ratio
       END IF

!!KZ       diag_ac1(i,k) = ratio 

       !!.................................................................
       !! conservation of rain  
       !!.................................................................

       dum  = (qrevp       + &
               sum(qrcol)  + &
               sum(qrhetc) + &
               sum(qrheti) + &
               sum(qrmul)   )*rcldm(i,k)*dt

       dum1 = qr(i,k)*rcldm(i,k)      + &
              ( qrcon*rcldm(i,k)      + &
                qcaut*lcldm(i,k)      + &
                qcacc*lcldm(i,k)      + &
                sum(qimlt)*icldm(i,k) + &
                sum(qcshd)*lcldm(i,k) )*dt
                
       l_dum1 = dum.gt.dum1
     
       l_dum2 = dum.ge.dum_min 
     
       l_possible = l_dum1 .and. l_dum2 
       
       IF (l_possible) THEN
       
          !!IF(l_masscon_debug) THEN 
          !!IF(masterproc) WRITE(6,*) '     qrevp  ', qrevp*rcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qrcol  ', sum(qrcol)*rcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qrhetc ', sum(qrhetc)*rcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qrheti ', sum(qrheti)*rcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qrmul  ', sum(qrmul)*rcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcaut  ', qcaut*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcacc  ', qcacc*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qimlt  ', sum(qimlt)*icldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qcshd  ', sum(qcshd)*lcldm(i,k) *dt
          !!IF(masterproc) WRITE(6,*) '     qr     ', qr(i,k)*rcldm(i,k) 
          !!IF(masterproc) WRITE(6,*) 'ratio: #2  conservation of rain ', ratio 
          !!END IF
       
          ratio  = dum1/dum
          qrevp  = qrevp*ratio
          qrcol  = qrcol*ratio
          qrhetc = qrhetc*ratio
          qrheti = qrheti*ratio
          qrmul  = qrmul*ratio
       END IF

       !!.................................................................
       !! conservation of cloud ice   
       !!.................................................................

       DO iice = 1,nCat

          !! sink 
          dum  = (qisub(iice) + &
                  qimlt(iice) )*icldm(i,k)*dt

          !! source 
          dum1 = qitot(i,k,iice)*icldm(i,k) + &
                 ( qidep(iice)*icldm(i,k)   + &
                   qberg      *lcldm(i,k)   + & !### 
                   qinuc(iice)*icldm(i,k)   + &
                   qrcol(iice)*rcldm(i,k)   + &
                   qccol(iice)*lcldm(i,k)   + &
                   qrhetc(iice)*rcldm(i,k)  + &
                   qrheti(iice)*rcldm(i,k)  + &
                   qchetc(iice)*lcldm(i,k)  + &
                   qcheti(iice)*lcldm(i,k)  + &
                   qrmul(iice)*rcldm(i,k)  )*dt

          DO catcoll = 1,nCat
            !category interaction leading to source for iice category
             dum1 = dum1 + qicol(catcoll,iice)*icldm(i,k)*dt
            !category interaction leading to sink for iice category
             dum = dum + qicol(iice,catcoll)*icldm(i,k)*dt
          END DO
          
          l_dum1 = dum.gt.dum1
     
          l_dum2 = dum.ge.dum_min 
     
          l_possible = l_dum1 .and. l_dum2 
      
          ratio = 1._r8 
 
          IF (l_possible) THEN
          
             ratio = dum1/dum
             
             !!IF(l_masscon_debug) THEN 
             !!IF(masterproc) WRITE(6,*) '     qisub ', qisub(iice)*icldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qimlt ', qimlt(iice)*icldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qinuc ', qinuc(iice)*icldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qrcol ', qrcol(iice)*rcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qccol ', qccol(iice)*lcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qrhetc ', qrhetc(iice)*rcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qrheti ', qrheti(iice)*rcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qchetc ', qchetc(iice)*lcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qcheti ', qcheti(iice)*lcldm(i,k) *dt
             !!IF(masterproc) WRITE(6,*) '     qrmul  ', qrmul(iice)*rcldm(i,k)  *dt 
             !!IF(masterproc) WRITE(6,*) '     qitot ', qitot(i,k,iice)*icldm(i,k)
             !!IF(masterproc) WRITE(6,*) 'ratio: #3  conservation of ice ', ratio 
             !!END IF
       
             qisub(iice) = qisub(iice)*ratio
             qimlt(iice) = qimlt(iice)*ratio
             
             DO catcoll = 1,nCat
                !!IF(masterproc .and. l_masscon_debug) &
                !!   WRITE(6,*) '     qicol ', qicol(iice,catcoll) *dt
                qicol(iice,catcoll) = qicol(iice,catcoll)*ratio
             END DO
             
!!KZ             diag_ac2(i,k) = ratio 

          END IF
          
       END DO !! iice 

END IF !!if(l_conservation_clipping) then 


!!######################################################################
!!#
!!#
!!# calculate tendency and update state 
!!#
!!#
!!######################################################################

       gqc_tend = 0._r8 
       gnc_tend = 0._r8 
       gqr_tend = 0._r8 
       gnr_tend = 0._r8 

       !!.................................................................
       !! update prognostic microphysics and thermodynamics variables
       !!
       !! grid mean calculation !! 
       !! 
       !! ice-phase dependent processes:
       !!.................................................................

       DO iice = 1,nCat

          gqc_tend = gqc_tend + (-qchetc(iice)*lcldm(i,k) &
                                 -qcheti(iice)*lcldm(i,k) &
                                 -qberg*lcldm(i,k)        & !! added qberg !### lcldm 
                                 -qccol(iice)*lcldm(i,k)  &
                                 -qcshd(iice)*lcldm(i,k)  )*dt
          
          !!gqc(i,k) = gqc(i,k) + (-qchetc(iice)*lcldm(i,k) &
          !!                       -qcheti(iice)*lcldm(i,k) &
          !!                       -qberg*lcldm(i,k)        & !! added qberg !### lcldm 
          !!                       -qccol(iice)*lcldm(i,k)  &
          !!                       -qcshd(iice)*lcldm(i,k)  )*dt
          !!
          IF (log_predictNc) THEN

             gnc_tend = gnc_tend + (-nccol(iice)*lcldm(i,k)  &
                         -nchetc(iice)*lcldm(i,k) &
                         -ncheti(iice)*lcldm(i,k) )*dt

             !!gnc(i,k) = gnc(i,k) + (-nccol(iice)*lcldm(i,k)  &
             !!                       -nchetc(iice)*lcldm(i,k) &
             !!                       -ncheti(iice)*lcldm(i,k) )*dt
          END IF

          gqr_tend = gqr_tend + (-qrcol(iice)*rcldm(i,k)  &
                      +qimlt(iice)*icldm(i,k)  &
                      -qrhetc(iice)*rcldm(i,k) &
                      -qrheti(iice)*rcldm(i,k) &
                      +qcshd(iice)*lcldm(i,k)  &     !! qcshd positive 
                      -qrmul(iice)*rcldm(i,k)  )*dt  !! qrmul zero 

          !!gqr(i,k) = gqr(i,k) + (-qrcol(iice)*rcldm(i,k)  &
          !!                       +qimlt(iice)*rcldm(i,k)  &
          !!                       -qrhetc(iice)*rcldm(i,k) &
          !!                       -qrheti(iice)*rcldm(i,k) &
          !!                       +qcshd(iice)*rcldm(i,k)  &     !! qcshd positive 
          !!                       -qrmul(iice)*rcldm(i,k)  )*dt  !! qrmul zero 

          
          !! apply factor to source for rain number from melting of ice, (ad-hoc
          !! but accounts for rapid evaporation of small melting ice particles)
        
          gnr_tend = gnr_tend + (-nrcol(iice)*rcldm(i,k)  & ! nrcol positive 
                      -nrhetc(iice)*rcldm(i,k) & ! nrhetc positive (now zero) 
                      -nrheti(iice)*rcldm(i,k) & ! nrheti positive 
                      -nmltratio*nimlt(iice)*rcldm(i,k) &
                      +nrshdr(iice)*rcldm(i,k) & ! nrshdr positive 
                      +ncshdc(iice)*rcldm(i,k) )*dt

          !!gnr(i,k) = gnr(i,k) + (-nrcol(iice)*rcldm(i,k)  & ! nrcol positive 
          !!                       -nrhetc(iice)*rcldm(i,k) & ! nrhetc positive (now zero) 
          !!                       -nrheti(iice)*rcldm(i,k) & ! nrheti positive 
          !!                       -nmltratio*nimlt(iice)*rcldm(i,k) &
          !!                       +nrshdr(iice)*rcldm(i,k) & ! nrshdr positive 
          !!                       +ncshdc(iice)*rcldm(i,k) )*dt

          IF (qitot(i,k,iice).ge.qsmall) THEN
             !! add sink terms, assume density stays constant for sink terms
             xdum1 = (qisub(iice)+qimlt(iice))/qitot(i,k,iice) 
             xdum2 = qisub(iice)+qimlt(iice)
             gbirim(i,k,iice) = gbirim(i,k,iice) - gbirim(i,k,iice) * xdum1 * dt 
             gqirim(i,k,iice) = gqirim(i,k,iice) - gqirim(i,k,iice) * xdum1 * dt
             gqitot(i,k,iice) = gqitot(i,k,iice) - xdum2*icldm(i,k)*dt 
          END IF

          !! source for ice
          dum = (qrcol(iice)*rcldm(i,k)  & !! qrcol positive 
                +qccol(iice)*lcldm(i,k)  & !! qccol positive 
                +qrhetc(iice)*rcldm(i,k) & !! qrhetc positive (if not zero) 
                +qrheti(iice)*rcldm(i,k) & !! qrheti positive
                +qchetc(iice)*icldm(i,k) & !! qchetc positive (if not zero) 
                +qcheti(iice)*icldm(i,k) & !! qcheti positive 
                +qberg*lcldm(i,k)        & !! added qberg !### lcldm
                +qrmul(iice)*rcldm(i,k)  )*dt 

          xdum1 = qidep(iice)+qinuc(iice)
          
          gqitot(i,k,iice) = gqitot(i,k,iice) + xdum1*icldm(i,k)*dt + dum
          
          gqirim(i,k,iice) = gqirim(i,k,iice) + dum
          
          xdum1 = qrcol(iice)*rcldm(i,k)*inv_rho_rimeMax
          xdum2 = qccol(iice)*lcldm(i,k)/rhorime_c(iice) 
          xdum3 = (qrhetc(iice)*rcldm(i,k) &
                  +qrheti(iice)*rcldm(i,k) &
                  +qchetc(iice)*icldm(i,k) &
                  +qcheti(iice)*icldm(i,k) &
                  +qrmul(iice)*rcldm(i,k)  )*inv_rho_rimeMax 
          
          gbirim(i,k,iice) = gbirim(i,k,iice) + (xdum1 + xdum2 + xdum3 )*dt

          gnitot(i,k,iice) = gnitot(i,k,iice) + &
                             (ninuc(iice)*icldm(i,k)  & ! ninuc positive 
                             +nimlt(iice)*icldm(i,k)  & ! nimlt negative 
                             +nisub(iice)*icldm(i,k)  & ! nisub negative 
                             -nislf(iice)*icldm(i,k)  & ! nislf postive 
                             +nrhetc(iice)*rcldm(i,k) &
                             +nrheti(iice)*rcldm(i,k) &
                             +nchetc(iice)*lcldm(i,k) &
                             +ncheti(iice)*lcldm(i,k) &
                             +nimul(iice)*icldm(i,k)  )*dt 

          DO catcoll = 1,nCat
          
             !! add ice-ice category interaction collection tendencies
             !! note: nicol is a sink for the collectee category, but NOT a source for collector

             gqitot(i,k,catcoll) = gqitot(i,k,catcoll) - qicol(catcoll,iice)*icldm(i,k)*dt
             gnitot(i,k,catcoll) = gnitot(i,k,catcoll) - nicol(catcoll,iice)*icldm(i,k)*dt
             gqitot(i,k,iice)    = gqitot(i,k,iice)    + qicol(catcoll,iice)*icldm(i,k)*dt
             
             ! now modify rime mass and density, assume collection does not modify rime mass
             ! fraction or density of the collectee, consistent with the assumption that
             ! these are constant over the PSD
             
             IF (gqitot(i,k,catcoll).ge.qsmall) THEN
                
                xdum1 = qicol(catcoll,iice)*icldm(i,k)*dt / gqitot(i,k,catcoll) 
             
                !! source for collector category
                
                gqirim(i,k,iice) = gqirim(i,k,iice)+gqirim(i,k,catcoll)*xdum1
                gbirim(i,k,iice) = gbirim(i,k,iice)+gbirim(i,k,catcoll)*xdum1

                !! sink for collectee category
                
                gqirim(i,k,catcoll) = gqirim(i,k,catcoll)-gqirim(i,k,catcoll)*xdum1
                gbirim(i,k,catcoll) = gbirim(i,k,catcoll)-gbirim(i,k,catcoll)*xdum1
                
             END IF

          END DO !! interactions_loop ! catcoll loop


          IF (gqirim(i,k,iice) .lt. 0._r8) THEN
             gqirim(i,k,iice) = 0._r8
             gbirim(i,k,iice) = 0._r8
          END IF

          !! densify under wet growth
          !! -- to be removed post-v2.1.  Densification automatically happens
          !!    during wet growth due to parameterized rime density --
        
          IF (log_wetgrowth(iice)) THEN
             gqirim(i,k,iice) = gqitot(i,k,iice)
             gbirim(i,k,iice) = gqirim(i,k,iice)*inv_rho_rimeMax
          END IF

          ! densify in above freezing conditions and melting
          ! -- future work --
          !   Ideally, this will be treated with the predicted liquid fraction in ice.
          !   Alternatively, it can be simplified by tending qirim -- qitot
          !   and birim such that rho_rim (qirim/birim) --> rho_liq during melting.
          ! ==


          !!qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
          !!     vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)
          !!     
          !!tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k)) &
          !!     *xxlv+(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
          !!     ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
          !!     pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)
          
          xdum1 = (qidep(iice) - qisub(iice) + qinuc(iice)) * icldm(i,k) 
          xdum2 = qrcol(iice)*rcldm(i,k) + &
                  qccol(iice)*lcldm(i,k) + &
                  qchetc(iice)*lcldm(i,k) + &
                  qcheti(iice)*lcldm(i,k)
          xdum3 = qrhetc(iice)*rcldm(i,k) + qrheti(iice)*rcldm(i,k)
          xdum4 = qimlt(iice)*icldm(i,k) 
          
          !!IF(masterproc) WRITE(6,*) 'qv change: #1 qidep ', k, qidep(iice) 
          !!IF(masterproc) WRITE(6,*) 'qv change: #2 qisub ', k, qisub(iice) 
          !!IF(masterproc) WRITE(6,*) 'qv change: #3 qinuc ', k, qinuc(iice) 
          !!IF(masterproc) WRITE(6,*) 'qv change: #4 total ', k, xdum1 
          !!IF(masterproc) WRITE(6,*) 'Ttend - qidep  : ', k, (qidep*xxls(i,k)*inv_cp)*dt
          !!IF(masterproc) WRITE(6,*) 'Ttend - qisub  : ', k, -(qisub*xxls(i,k)*inv_cp)*dt
          !!IF(masterproc) WRITE(6,*) 'Ttend - qinuc  : ', k, (qinuc*xxls(i,k)*inv_cp)*dt
          
          qv(i,k) = qv(i,k) - xdum1 * dt
          
          IF(qv(i,k).lt.dum_min) THEN 
             IF(masterproc) WRITE(6,*) & 
                               '### L3633 qv is negative : ', i, k, qv(i,k), xdum1, icldm(i,k), dt, &
                               qidep(iice), qisub(iice), qinuc(iice) 
          END IF

          tt(i,k) = tt(i,k) + ( xdum1 * xxls(i,k) & !! nucleation/deposition/sublimation
                              + xdum2 * xlf(i,k)  & !! freezing of droplets 
                              + xdum3 * xlf(i,k)  & !! freezing of rain 
                              + qberg * xlf(i,k)  & !! freezing of droplets through Bergeron process
                              - xdum4 * xlf(i,k)  & !! melting of ice 
                              ) * inv_cp * dt 

#ifdef P3_SLOW
   IF(l_debug) THEN 
       WRITE(6,*) 'Ttend - ice phase : ', k, &
                              ( xdum1 * icldm(i,k) * xxls(i,k) & !! nucleation/deposition/sublimation
                              + xdum2 * lcldm(i,k) * xlf(i,k)  & !! freezing of droplets 
                              + xdum3 * rcldm(i,k) * xlf(i,k) &  !! freezing of rain 
                              - xdum4 * icldm(i,k) * xlf(i,k)      &  !! melting of ice 
                              ) * inv_cp * dt 
       WRITE(6,*) 'Ttend - qidep  : ', k, (qidep*xxls(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qisub  : ', k, -(qisub*xxls(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qinuc  : ', k, (qinuc*xxls(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qrcol  : ', k, (qrcol*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qccol  : ', k, (qccol*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qchetc : ', k, (qchetc*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qcheti : ', k, (qcheti*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qrhetc : ', k, (qrhetc*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qrheti : ', k, (qrheti*xlf(i,k)*inv_cp)*dt
       WRITE(6,*) 'Ttend - qimlt  : ', k, -(qimlt*xlf(i,k)*inv_cp)*dt
   END IF
#endif

       END DO !! iice_loop2



       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 028 -'


       !!.................................................................
       !! warm-phase only processes:
       !!.................................................................

       gqc_tend = gqc_tend + (-qcacc*lcldm(i,k) &
                              -qcaut*lcldm(i,k) &
                              +qcnuc*lcldm(i,k) &
                              +qccon*lcldm(i,k) &
                              -qcevp*lcldm(i,k) )*dt
                            
       !!gqc(i,k) = gqc(i,k) + (-qcacc*lcldm(i,k) &
       !!                       -qcaut*lcldm(i,k) &
       !!                       +qcnuc*lcldm(i,k) &
       !!                       +qccon*lcldm(i,k) &
       !!                       -qcevp*lcldm(i,k) )*dt
                            
       gqr_tend = gqr_tend + ( qcacc*lcldm(i,k) &
                              +qcaut*lcldm(i,k) &
                              +qrcon*rcldm(i,k) &
                              -qrevp*rcldm(i,k) )*dt

       !!gqr(i,k) = gqr(i,k) + ( qcacc*rcldm(i,k) &
       !!                       +qcaut*rcldm(i,k) &
       !!                       +qrcon*rcldm(i,k) &
       !!                       -qrevp*rcldm(i,k) )*dt

       IF (log_predictNc) THEN
          gnc_tend = gnc_tend + (-ncacc *lcldm(i,k)  &
                                 -ncautc*lcldm(i,k)  &
                                 +ncslf *lcldm(i,k)  &
                                 +ncnuc *lcldm(i,k)  )*dt
          !!gnc(i,k) = gnc(i,k) + (-ncacc *lcldm(i,k)  &
          !!                       -ncautc*lcldm(i,k)  &
          !!                       +ncslf *lcldm(i,k)  &
          !!                       +ncnuc *lcldm(i,k)  )*dt
       ELSE
          gnc(i,k) = nccnst*inv_rho(i,k)*lcldm(i,k)
       END IF
       
       IF (iparam.eq.1 .or. iparam.eq.2) THEN
          gnr_tend = gnr_tend + ( 0.5*ncautc*rcldm(i,k) + &
                                  nrslf*rcldm(i,k)  + &
                                  nrevp*rcldm(i,k)      )*dt
          !!gnr(i,k) = gnr(i,k) + ( 0.5*ncautc*rcldm(i,k) + &
          !!                            nrslf*rcldm(i,k)  + &
          !!                            nrevp*rcldm(i,k)      )*dt
       ELSE
          gnr_tend = gnr_tend + (ncautr*rcldm(i,k) + & ! positive 
                                 nrslf*rcldm(i,k)  + & ! negative 
                                 nrevp*rcldm(i,k)      )*dt !
          !!gnr(i,k) = gnr(i,k) + (ncautr*rcldm(i,k) + & ! positive 
          !!                       nrslf*rcldm(i,k)  + & ! negative 
          !!                       nrevp*rcldm(i,k)      )*dt !
       END IF

        !! MG2 
        !!qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
        !!     vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)
        !!     
        !!tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k)) &
        !!     *xxlv+(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
        !!     ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
        !!     pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)
        !! original 
        !!       qv(i,k) = qv(i,k) + (-qcnuc-qccon-qrcon+qcevp+qrevp)*dt
        !!       tt(i,k) = tt(i,k) + ((qcnuc+qccon+qrcon-qcevp-qrevp)*xxlv(i,k)*    &
        !!                 inv_cp)*dt
      
       !!.................................................................
       !! Re-apply droplet activation tendency
       !!.................................................................
     
       if(l_cnuc) then  
          gnc_tend = gnc_tend + npccn(i,k) * dt 
       end if 

       gqc(i,k) = qcn(i,k) + gqc_tend 
       gnc(i,k) = ncn(i,k) + gnc_tend 
       gqr(i,k) = qrn(i,k) + gqr_tend 
       gnr(i,k) = nrn(i,k) + gnr_tend 
 
       !!.................................................................
       !! note that condensational heating is handled by clubb 
       !!.................................................................
       
       xdum1 = qv(i,k) 
 
       qv(i,k) = qv(i,k) + ( -qcnuc * lcldm(i,k) &
                             -qccon * lcldm(i,k) &
                             -qrcon * rcldm(i,k) &
                             +qcevp * lcldm(i,k) &
                             +qrevp * rcldm(i,k) )*dt 
       
       IF(qv(i,k).lt.1.e-22) THEN 
          IF(masterproc) WRITE(6,*) '### L3705 qv is negative : ', &
                         i, k, xdum1, qv(i,k), qcnuc, qcevp, qrevp, lcldm(i,k), rcldm(i,k) 
       END IF

       tt(i,k) = tt(i,k) + ( ( qcnuc * lcldm(i,k) &
                              +qccon * lcldm(i,k) &
                              +qrcon * rcldm(i,k) &
                              -qcevp * lcldm(i,k) &
                              -qrevp * rcldm(i,k) ) * &
                               xxlv(i,k) * inv_cp ) * dt

#ifdef P3_SLOW
       IF(l_debug) THEN
           WRITE(6,*) 'Ttend - warm phase : ', k, ( ( qcnuc * lcldm(i,k) - &
                                   qcevp * lcldm(i,k) - &
                                   qrevp * rcldm(i,k) ) * &
                                   xxlv(i,k) * inv_cp ) * dt
           WRITE(6,*) 'Ttend - qcnuc : ', k, (qcnuc*xxlv(i,k)*inv_cp)*dt
           WRITE(6,*) 'Ttend - qccon : ', k, (qccon*xxlv(i,k)*inv_cp)*dt
           WRITE(6,*) 'Ttend - qrcon : ', k, (qrcon*xxlv(i,k)*inv_cp)*dt
           WRITE(6,*) 'Ttend - qcevp : ', k, -(qcevp*xxlv(i,k)*inv_cp)*dt
           WRITE(6,*) 'Ttend - qrevp : ', k, -(qrevp*xxlv(i,k)*inv_cp)*dt
       END IF 
#endif


       !!.................................................................
       !! clipping for small hydrometeor values
       !!.................................................................
     
       IF (gqc(i,k).lt.qsmall) THEN
          xdum1 = qv(i,k) 
          qv(i,k) = qv(i,k) + gqc(i,k) 
          IF(qv(i,k) .lt. 0._r8) THEN 
             WRITE(6,*) '### L3663 ### qv .lt. 0', i, k, xdum1, qv(i,k), gqc(i,k), lcldm(i,k) 
             WRITE(6,*) '!!!!! L3663 !!!!!' 
          END IF 
          tt(i,k) = tt(i,k) - gqc(i,k) * xxlv(i,k) * inv_cp
          
          
          !!IF(l_debug) WRITE(6,*) 'Ttend - clipping qc : ', k, - qc(i,k)*xxlv(i,k)*inv_cp
          
          gqc(i,k) = 0.
          gnc(i,k) = 0.
       ELSE
          log_hydrometeorsPresent = .True.
       END IF

       IF (gqr(i,k).lt.qsmall) THEN
          xdum1 = qv(i,k) 
          qv(i,k) = qv(i,k) + gqr(i,k)
          IF(qv(i,k) .lt. 0._r8) THEN 
             WRITE(6,*) '### L3680 ### qv .lt. 0', i, k, xdum1, qv(i,k), gqr(i,k), rcldm(i,k) 
             WRITE(6,*) '!!!!! L3680 !!!!!' 
          END IF 
          tt(i,k) = tt(i,k) - gqr(i,k) * xxlv(i,k) * inv_cp
          
          
          !!IF(l_debug) WRITE(6,*) 'Ttend - clipping qr : ', k, &
          !!KZ!!                - qr(i,k) * rcldm(i,k) * xxlv(i,k) * inv_cp
             
             
          gqr(i,k) = 0.
          gnr(i,k) = 0.
       ELSE
          log_hydrometeorsPresent = .True.
       END IF

       DO iice = 1,nCat
          IF (gqitot(i,k,iice).lt.qsmall) THEN
             xdum1 = qv(i,k) 
             qv(i,k) = qv(i,k) + gqitot(i,k,iice) 
             IF(qv(i,k) .lt. 0._r8) THEN 
                WRITE(6,*) '### L3700 ### qv .lt. 0', i, k, xdum1, qv(i,k), gqitot(i,k,iice), icldm(i,k) 
                WRITE(6,*) '!!!!! L3700 !!!!!' 
             END IF 
             tt(i,k) = tt(i,k) - gqitot(i,k,iice) * xxls(i,k) * inv_cp
             
             
             !!IF(l_debug) WRITE(6,*) 'Ttend - clipping qi : ', k, &
             !!KZ!!             - qitot(i,k,iice) * icldm(i,k) * xxls(i,k) * inv_cp
          
             gqitot(i,k,iice) = 0.
             gnitot(i,k,iice) = 0.
             gqirim(i,k,iice) = 0.
             gbirim(i,k,iice) = 0.
          ELSE
             log_hydrometeorsPresent = .True.
          END IF
       END DO !iice-loop

       CALL impose_max_total_Ni(gnitot(i,k,:),max_total_Ni,inv_rho(i,k))

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 029 -'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DO k = kbot,ktop,kdir  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !!.................................................................
       !! output liquid-phase process rates
       !!.................................................................
       
       diag_qcacc (i,k) = qcacc * lcldm(i,k)  
       diag_qrevp (i,k) = qrevp * rcldm(i,k)  
       diag_qccon (i,k) = qccon * lcldm(i,k)  
       diag_qcaut (i,k) = qcaut * lcldm(i,k)  
       diag_qcevp (i,k) = qcevp * lcldm(i,k)  
       diag_qberg (i,k) = qberg * lcldm(i,k)  !### lcldm 
       diag_qrcon (i,k) = qrcon * rcldm(i,k)  
       diag_ncacc (i,k) = ncacc * lcldm(i,k)  
       diag_ncnuc (i,k) = ncnuc * lcldm(i,k)  
       diag_ncslf (i,k) = ncslf * lcldm(i,k)  
       diag_ncautc(i,k) = ncautc * lcldm(i,k)  
       diag_qcnuc (i,k) = qcnuc * lcldm(i,k)  
       diag_nrslf (i,k) = nrslf * rcldm(i,k)  
       diag_nrevp (i,k) = nrevp * rcldm(i,k)  
       diag_ncautr(i,k) = ncautr * lcldm(i,k)  

       !! rain/snow production rate (no evaporation and sublimation considered) 
       !! 
       !! mg2 
       !!
       !! prain(i,k) = (pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
       !!               mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)
       !! prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
       !!                  pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

       diag_prain(i,k) = qcacc * lcldm(i,k) + &
                         qcaut * lcldm(i,k) + &
                         qrcon * rcldm(i,k) 
       
       DO iice = 1,nCat  !! bugfix for v47 !! collection of liquid and ice by ice considered as snow 
          diag_prain(i,k) = diag_prain(i,k) + qcshd(iice) * lcldm(i,k) &
                                            + qicol(iice,iice) * icldm(i,k) & !! only consider self collectio now 
                                            + qccol(iice) * lcldm(i,k) 
       END DO 
       
       diag_pevap(i,k) = qrevp * rcldm(i,k)    ! positive term 

       !!.................................................................
       !! output ice-phase process rates
       !!.................................................................
 
       DO iice = 1,nCat

          diag_qchetc(i,k,iice) = qchetc(iice) * lcldm(i,k)  
          diag_qisub (i,k,iice) = qisub (iice) * icldm(i,k)  
          diag_nrshdr(i,k,iice) = nrshdr(iice) * rcldm(i,k)  
          diag_qcheti(i,k,iice) = qcheti(iice) * icldm(i,k)  
          diag_qrcol (i,k,iice) = qrcol (iice) * rcldm(i,k)  
          diag_qcshd (i,k,iice) = qcshd (iice) * lcldm(i,k)  
          diag_qrhetc(i,k,iice) = qrhetc(iice) * rcldm(i,k)  
          diag_qimlt (i,k,iice) = qimlt (iice) * icldm(i,k)  
          diag_qccol (i,k,iice) = qccol (iice) * lcldm(i,k)  
          diag_qrheti(i,k,iice) = qrheti(iice) * rcldm(i,k)  
          diag_qinuc (i,k,iice) = qinuc (iice) * icldm(i,k)  
          diag_nimlt (i,k,iice) = nimlt (iice) * icldm(i,k)  
          diag_nchetc(i,k,iice) = nchetc(iice) * lcldm(i,k)  
          diag_nccol (i,k,iice) = nccol (iice) * lcldm(i,k)  
          diag_ncshdc(i,k,iice) = ncshdc(iice) * lcldm(i,k)  
          diag_ncheti(i,k,iice) = ncheti(iice) * lcldm(i,k)  
          diag_nrcol (i,k,iice) = nrcol (iice) * rcldm(i,k)  
          diag_nislf (i,k,iice) = nislf (iice) * icldm(i,k)  
          diag_nrhetc(i,k,iice) = nrhetc(iice) * rcldm(i,k)  
          diag_ninuc (i,k,iice) = ninuc (iice) * icldm(i,k)  
          diag_qidep (i,k,iice) = qidep (iice) * icldm(i,k)  
          diag_nrheti(i,k,iice) = nrheti(iice) * rcldm(i,k)  
          diag_nisub (i,k,iice) = nisub (iice) * icldm(i,k)  
          diag_qwgrth(i,k,iice) = qwgrth(iice) * rcldm(i,k)  
          diag_qrmul (i,k,iice) = qrmul (iice) * rcldm(i,k)  
          diag_nimul (i,k,iice) = nimul (iice) * icldm(i,k)  
          diag_nevapr(i,k,iice) = (qisub(iice) * icldm(i,k) + qrevp * rcldm(i,k))    ! positive term 
          diag_cmei  (i,k,iice) = (qidep(iice) - qisub(iice) + qinuc(iice)) * icldm(i,k)  

          if(qidep(iice) .lt. -1.e-20) then 
             write(6,*) 'ZZZ### qidep ', qidep(iice)
          end if 
 
          if(qisub(iice) .lt. -1.e-20) then 
             write(6,*) 'ZZZ### qisub ', qisub(iice)
          end if 
 
          if(qinuc(iice) .lt. -1.e-20) then 
             write(6,*) 'ZZZ### qinuc ', qinuc(iice)
          end if 
 
       END DO 
       
       !! percentage of droplets -> precipitation 
       
       diag_1stqc2qr(i,k) = qcacc * lcldm(i,k) + &
                            qcaut * lcldm(i,k) + &
                            sum(qcshd(:))*lcldm(i,k)  
          
       diag_1stqc2qr(i,k) = diag_1stqc2qr(i,k) / max(gqc(i,k),dum_min) 

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 030 -'


555 CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DO k = kbot,ktop,kdir  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    END DO !! k_loop_main


    IF (.not. log_hydrometeorsPresent) GOTO 333

 
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 031 -'



    !!.................................................................
    !! get in-cloud fields again and clip negative values 
    !!.................................................................

    DO k = kbot,ktop,kdir

       IF (gqc(i,k).ge.qsmall .and. lcldm(i,k).ge. cldm_min ) THEN
          qc(i,k) = gqc(i,k)/lcldm(i,k)
          nc(i,k) = gnc(i,k)/lcldm(i,k)
       ELSE
          qc(i,k) = 0._r8
          nc(i,k) = 0._r8
       END IF
       
       IF (gqr(i,k).ge.qsmall .and. rcldm(i,k).ge. cldm_min ) THEN
          qr(i,k) = gqr(i,k)/rcldm(i,k)
          nr(i,k) = gnr(i,k)/rcldm(i,k)
       ELSE
          qr(i,k) = 0._r8
          nr(i,k) = 0._r8
       END IF
       
       DO iice = 1,nCat
          IF (gqitot(i,k,iice).ge.qsmall .and. icldm(i,k).ge. cldm_min ) THEN
             qitot(i,k,iice) = gqitot(i,k,iice)/icldm(i,k)
             nitot(i,k,iice) = gnitot(i,k,iice)/icldm(i,k)
          ELSE
             qitot(i,k,iice) = 0._r8
             nitot(i,k,iice) = 0._r8
          END IF
       END DO 

       DO iice = 1,nCat
          IF (qirim(i,k,iice).ge.qsmall .and. icldm(i,k).ge. cldm_min ) THEN
             qirim(i,k,iice) = gqirim(i,k,iice)/icldm(i,k)
             birim(i,k,iice) = gbirim(i,k,iice)/icldm(i,k)
          ELSE
             qirim(i,k,iice) = 0._r8
             birim(i,k,iice) = 0._r8
          END IF
       END DO 

    END DO 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Sedimentation 
    !!
    !!
    !!   size distribution calculated using in-cloud values 
    !!   sedimentation calculated using grid-mean values 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !!.................................................................
    !! Cloud sedimentation 
    !!.................................................................


    IF(l_csed) THEN 

    !! initialize logicals for presence of hydrometeor species to .false.
    
    log_qcpresent = .false.

    DO k = ktop,kbot,-kdir

       inv_dzq(i,k) = 1./dzq(i,k)

       !! calculate Q- and N-weighted fallspeeds and find highest k level that hydrometeor is present

       CALL get_cloud_dsd(qc(i,k),        &
                          nc(i,k),        &
                          diag_mu_c(i,k), &
                          rho(i,k),       &
                          nu(i,k),        &
                          dnu,            &
                          diag_lamc(i,k), &
                          lammin,         &
                          lammax,         &
                          k,              &
                          tmp1,           &
                          tmp2,           &
                          qcindex,        &
                          log_qcpresent   )

       !! droplet fall speed
       !! all droplets in smallest category fallspeed; thus, analytic solution can be used

       IF (gqc(i,k).ge.qsmall) THEN
       
          dum = 1._r8/diag_lamc(i,k)**bcn
          
          xdum1 = gamma(1._r8+bcn+diag_mu_c(i,k))  
          xdum2 = gamma(diag_mu_c(i,k)+1._r8)
          xdum3 = gamma(4.+bcn+diag_mu_c(i,k))
          xdum4 = gamma(diag_mu_c(i,k)+4.)
          
          IF (log_predictNc) THEN
             Vt_nc(i,k) = acn(i,k)*xdum1*dum/xdum2 
          END IF
          
          Vt_qc(i,k) = acn(i,k)*xdum3*dum/xdum4 
          
       ELSE
       
          IF (log_predictNc) THEN
             Vt_nc(i,k) = 0._r8
          END IF
          
          Vt_qc(i,k) = 0._r8
          
       END IF

       diag_vmc(i,k) = Vt_qc(i,k) ! output fallspeed
       
    END DO ! k-loop


    !!.................................................................
    !! sedimentation of droplet mass
    !!.................................................................

    IF (log_qcpresent) THEN
    
       nstep = 1
       
       DO k = qcindex+kdir,kbot,-kdir

         !- weighted fall speed arrays used for sedimentation calculations
         !  (assigned below to highest non-zero level value at lower levels with Vt_x=0)
          V_qc(K)  = Vt_qc(i,k)

          IF (kdir.eq.1) THEN
             IF (k.le.qcindex-kdir) THEN
                IF (V_qc(k).lt.1.e-10_r8) THEN
                   V_qc(k) = V_qc(k+kdir)
                END IF
             END IF
          ELSEIF (kdir.eq.-1) THEN
             IF (k.ge.qcindex-kdir) THEN
                IF (V_qc(k).lt.1.e-10_r8) THEN
                   V_qc(k) = V_qc(k+kdir)
                END IF
             END IF
          END IF


          !! calculate number of split time steps

          rgvm       = V_qc(k)
          nstep      = max(int(rgvm*dt*inv_dzq(i,k)+1.),nstep)
          dum_qc(k)  = gqc(i,k)*rho(i,k)
          tend_qc(K) = 0.

       END DO ! k-loop

       inv_nstep = 1./real(nstep)

       IF (l_debug .and. nstep.ge.100) THEN
          print*,'micro_p3 L3496 CLOUD nstep LARGE:',i,nstep
          !!!stop
       END IF

       !! calculate sedimentation using first-order upwind method
       
       tmp1 = 0._r8
       
       DO n = 1,nstep

          DO k = kbot,qcindex,kdir
             flux_qc(k) = V_qc(k)*dum_qc(k)
          END DO
          
          tmp1 = tmp1 + flux_qc(kbot)  !sum flux_ at lowest level for averaging over sub-stepping

          !! top level with hydrometeor present
          
          k = qcindex
          
          fluxdiv_qc = flux_qc(k)*inv_dzq(i,k)
          tend_qc(k) = tend_qc(k)-fluxdiv_qc*inv_nstep*inv_rho(i,k)
          dum_qc(k)  = dum_qc(k)-fluxdiv_qc*dt*inv_nstep

          !! loop from sceond to top level of hydrometeor to surface
          DO k = qcindex-kdir,kbot,-kdir
             fluxdiv_qc = (flux_qc(k+kdir)-flux_qc(K))*inv_dzq(i,k)
             tend_qc(k) = tend_qc(k)+fluxdiv_qc*inv_nstep*inv_rho(i,k)
             dum_qc(k)  = dum_qc(k)+fluxdiv_qc*dt*inv_nstep
          END DO ! k loop

       END DO ! nstep-loop

       DO k = kbot,qcindex,kdir
          gqc(i,k) = gqc(i,k)+tend_qc(k)*dt
          diag_sedqc(i,k) = tend_qc(k) 
       END DO

       !! compute cloud contribution to liquid precipitation rate at surface
       
       tmp1 = tmp1*inv_nstep           !flux_ at surface, averaged over sub-step
       pcprt_liq(i) = tmp1*inv_rhow    !convert flux_ (kg m-2 s-1) to pcp rate (m s-1)
 
       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 032 -'
       
       !!.................................................................
       !! sedimentation of droplet number
       !!.................................................................

       IF (log_predictNc) THEN

       nstep = 1

       DO k = qcindex+kdir,kbot,-kdir

         !- weighted fall speed arrays used for sedimentation calculations
         !  (assigned below to highest non-zero level value at lower levels with Vt_x=0)
          V_nc(K) = Vt_nc(i,k)

          IF (kdir.eq.1) THEN
             IF (k.le.qcindex-kdir) THEN
                IF (V_nc(k).lt.1.E-10) THEN
                   V_nc(k) = V_nc(k+kdir)
                END IF
             END IF
          ELSEIF (kdir.eq.-1) THEN
             IF (k.ge.qcindex-kdir) THEN
                IF (V_nc(k).lt.1.e-10) THEN
                   V_nc(k) = V_nc(k+kdir)
                END IF
             END IF
          END IF

          !! calculate number of split time steps
          rgvm       = V_nc(k)
          nstep      = max(int(rgvm*dt*inv_dzq(i,k)+1.),nstep)
          dum_nc(k)  = gnc(i,k)*rho(i,k)
          tend_nc(K) = 0.

       END DO ! k-loop

       inv_nstep = 1./real(nstep)

       IF (l_debug .and. nstep.ge.100) THEN
          print*,'micro_p3 L3572 CLOUD nstep LARGE:',i,nstep
          !!!stop
       END IF

       !! calculate sedimentation using first-order upwind method
       DO n = 1,nstep

          DO k = kbot,qcindex,kdir
             flux_nc(k) = V_nc(k)*dum_nc(k)
          END DO

          !! top level with hydrometeor present
          k = qcindex
          fluxdiv_nc = flux_nc(k)*inv_dzq(i,k)
          tend_nc(k) = tend_nc(k)-fluxdiv_nc*inv_nstep*inv_rho(i,k)
          dum_nc(k)  = dum_nc(k)-fluxdiv_nc*dt*inv_nstep

          !! loop from sceond to top level of hydrometeor to surface
          DO k = qcindex-kdir,kbot,-kdir

             fluxdiv_nc = (flux_nc(k+kdir)-flux_nc(K))*inv_dzq(i,k)
             tend_nc(k) = tend_nc(k)+fluxdiv_nc*inv_nstep*inv_rho(i,k)
             dum_nc(k)  = dum_nc(k)+fluxdiv_nc*dt*inv_nstep

          END DO ! k loop

       END DO ! nstep-loop

       DO k = kbot,qcindex,kdir
          gnc(i,k) = gnc(i,k)+tend_nc(k)*dt
          diag_sednc(i,k) = tend_nc(k) !! KZP3  
       END DO
       
    END IF ! log_predictNc

    END IF ! log_qcpresent
 
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 033s -'

    END IF 
    

    !!.................................................................
    !! sedimentation of rain mass and number 
    !!.................................................................

    IF(l_rsed) THEN 

    log_qrpresent = .false.

    DO k = ktop,kbot,-kdir
 
       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0331 -', k

       CALL get_rain_dsd(qr(i,k),nr(i,k),diag_mu_r(i,k),rdumii,dumii,lamr(i,k),mu_r_table,    &
                         tmp1,tmp2,log_qrpresent,qrindex,k)

       !! note: tmp1,tmp2 are not used in this section

       !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0332 -', k

       IF (qr(i,k).ge.qsmall) THEN

          !! read in fall mass- and number-weighted fall speeds from lookup table

          CALL find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,       &
                                          diag_mu_r(i,k),lamr(i,k))

          !!.................................
          !! number-weighted fall speed 
          !!.................................
          
          !! at diag_mu_r:
          
          xdum1 = vn_table(dumii,dumjj)
          xdum2 = rdumii-real(dumii) 
          xdum3 = vn_table(dumii+1,dumjj) 
          xdum4 = vn_table(dumii,dumjj)   !! bugfix 
          
          dum1 = xdum1 + xdum2*inv_dum3*(xdum3-xdum4)
                 
          !! at diag_mu_r+1
          
          xdum1 = vn_table(dumii,dumjj+1)
          xdum2 = rdumii-real(dumii)
          xdum3 = vn_table(dumii+1,dumjj+1)
          xdum4 = vn_table(dumii,dumjj+1) !! bugfix 
          
          dum2 = xdum1 + xdum2*inv_dum3*(xdum3-xdum4)
                 
          !! interpolated:
          
          xdum1 = rdumjj-real(dumjj) 
          
          Vt_nr(i,k) = dum1+xdum1*(dum2-dum1)
          Vt_nr(i,k) = Vt_nr(i,k)*rhofacr(i,k)

          !!.................................
          !! mass-weighted fall speed
          !!.................................
          
          !! at diag_mu_r:
          
          xdum1 = vm_table(dumii,dumjj)
          xdum2 = rdumii-real(dumii)
          xdum3 = vm_table(dumii+1,dumjj)
          xdum4 = vm_table(dumii,dumjj)   !! bugfix 
          
          dum1 = xdum1 + xdum2*inv_dum3*(xdum3-xdum4)

          !! at diag_mu_r+1 
          
          xdum1 = vm_table(dumii,dumjj+1)
          xdum2 = rdumii-real(dumii)
          xdum3 = vm_table(dumii+1,dumjj+1)
          xdum4 = vm_table(dumii,dumjj+1) !! bugfix 
          
          dum2 = xdum1 + xdum2*inv_dum3*(xdum3-xdum4) 

          !! interpolated:
          
          xdum1 = rdumjj-real(dumjj)
          
          Vt_qr(i,k) = dum1 + xdum1*(dum2-dum1)
          Vt_qr(i,k) = Vt_qr(i,k)*rhofacr(i,k)

       ELSE

          Vt_nr(i,k) = 0._r8
          Vt_qr(i,k) = 0._r8

       END IF

       diag_vmr(i,k) = Vt_qr(i,k) ! output fallspeed
       
       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0333 -', k, diag_vmr(i,k) 

    END DO ! k-loop


    IF (log_qrpresent) THEN

       nstep = 1

       DO k = qrindex+kdir,kbot,-kdir

          !- weighted fall speed arrays used for sedimentation calculations
          !  (assigned below to highest non-zero level value at lower levels with Vt_x=0)
          V_qr(k) = Vt_qr(i,k)
          V_nr(k) = Vt_nr(i,k)

          IF (kdir.eq.1) THEN
             IF (k.le.qrindex-kdir) THEN
                IF (V_qr(k).lt.1.e-10) THEN
                   V_qr(k) = V_qr(k+kdir)
                END IF
                IF (V_nr(k).lt.1.e-10) THEN
                   V_nr(k) = V_nr(k+kdir)
                END IF
             END IF
          ELSEIF (kdir.eq.-1) THEN
             IF (k.ge.qrindex-kdir) THEN
                IF (V_qr(k).lt.1.e-10) THEN
                   V_qr(k) = V_qr(k+kdir)
                END IF
                IF (V_nr(k).lt.1.e-10) THEN
                   V_nr(k) = V_nr(k+kdir)
                END IF
             END IF
          END IF

          !! calculate number of split time steps
          rgvm       = max(V_qr(k),V_nr(k))
          nstep      = max(int(rgvm*dt*inv_dzq(i,k)+1.),nstep)
          dum_qr(k)  = gqr(i,k)*rho(i,k)
          dum_nr(k)  = gnr(i,k)*rho(i,k)
          tend_qr(k) = 0.
          tend_nr(k) = 0.

       END DO ! k-loop

       inv_nstep = 1./real(nstep)

       IF (l_debug .and. nstep .ge. 100) THEN
          print*,'micro_p3 L3713 RAIN nstep LARGE:',i,nstep
          !!!stop
       END IF

       !--test:  explicitly calculate pcp rate:
       ! pcprt_liq(i) = qr(i,kbot)*rho(i,kbot)*Vt_qr(i,kbot)*1.e-3  !m s-1

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0334 -' 

       !! calculate sedimentation using first-order upwind method
       tmp1 = 0.
       diag_rflx(i,:) = 0._r8
       
       DO n = 1,nstep

          DO k = kbot,qrindex,kdir
             flux_qr(k) = V_qr(k)*dum_qr(k)
             flux_nr(k) = V_nr(k)*dum_nr(k)
             diag_rflx(i,k+1) = diag_rflx(i,k+1) + flux_qr(k)
          END DO
          
          tmp1 = tmp1 + flux_qr(kbot)  !sum flux_ at lowest level for averaging over sub-stepping

          !! top level with hydrometeor present
          k          = qrindex
          fluxdiv_qr = flux_qr(k)*inv_dzq(i,k)
          fluxdiv_nr = flux_nr(k)*inv_dzq(i,k)
          tend_qr(k) = tend_qr(k) - fluxdiv_qr*inv_nstep*inv_rho(i,k)
          tend_nr(k) = tend_nr(k) - fluxdiv_nr*inv_nstep*inv_rho(i,k)
          dum_qr(k)  = dum_qr(k)  - fluxdiv_qr*dt*inv_nstep
          dum_nr(k)  = dum_nr(k)  - fluxdiv_nr*dt*inv_nstep

          !! loop from second to top level of hydrometeor to surface
          DO k = qrindex-kdir,kbot,-kdir
             fluxdiv_qr = (flux_qr(k+kdir) - flux_qr(K))*inv_dzq(i,k)
             fluxdiv_nr = (flux_nr(k+kdir) - flux_nr(K))*inv_dzq(i,k)
             tend_qr(k) = tend_qr(k) + fluxdiv_qr*inv_nstep*inv_rho(i,k)
             tend_nr(k) = tend_nr(k) + fluxdiv_nr*inv_nstep*inv_rho(i,k)
             dum_qr(k)  = dum_qr(k)  + fluxdiv_qr*dt*inv_nstep
             dum_nr(k)  = dum_nr(k)  + fluxdiv_nr*dt*inv_nstep
          END DO ! k loop

       END DO ! nstep loop

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0335 -' 

       !! update prognostic variables with sedimentation tendencies
       DO k = kbot,qrindex,kdir
          gqr(i,k) = gqr(i,k) + tend_qr(k)*dt
          gnr(i,k) = gnr(i,k) + tend_nr(k)*dt
          diag_sedqr(i,k) = tend_qr(k) !! ### 
          diag_sednr(i,k) = tend_nr(k) !! ### 
       END DO

       !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 0335 -' 

       !! add rain component of liquid precipitation rate at surface
       tmp1 = tmp1*inv_nstep               !flux_ at surface, averaged over sub-step
       tmp1 = tmp1*inv_rhow                !convert flux_ (kg m-2 s-1) to pcp rate (m s-1)

       pcprt_liq(i) = pcprt_liq(i) + tmp1  !add pcp rate from cloud and rain

    END IF ! log_qrpresent
 
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 034 -'

    END IF 


    !!.................................................................
    !! sedimentation of ice mass and number 
    !!.................................................................

    IF(l_ised) THEN 
    
    DO iice = 1,nCat

       log_qipresent = .false.  !note: this applies to ice category 'iice' only

       DO k = ktop,kbot,-kdir

          !! get ice fallspeed for updated variables
          IF (qitot(i,k,iice).ge.qsmall) THEN

             !! impose lower limits to prevent taking log of # .lt. 0
             nitot(i,k,iice) = max(nitot(i,k,iice),nsmall)
             nr(i,k)         = max(nr(i,k),nsmall)

             CALL calc_bulkRhoRime(qitot(i,k,iice),qirim(i,k,iice),birim(i,k,iice),rhop)

             !! IF (.not. tripleMoment_on) zitot(i,k,iice) = diag_mom6(qitot(i,k,iice),nitot(i,k,iice),rho(i,k))

             CALL find_lookupTable_indices_1( &
                                       dumi,dumj,dumjj,dumii,dumzz,                      &
                                       dum1,dum3,dum4,dum5,dum6,                         &
                                       isize,rimsize,densize,zsize,rcollsize,            &
                                       qr(i,k),nr(i,k),qitot(i,k,iice),nitot(i,k,iice),  &
                                       qirim(i,k,iice),999._r8,rhop,200)
                                       !!qirim(i,k,iice),zitot(i,k,iice),rhop,200)

             CALL access_lookup_table(dumjj,dumii,dumi, 1,dum1,dum4,dum5,f1pr01)
             CALL access_lookup_table(dumjj,dumii,dumi, 2,dum1,dum4,dum5,f1pr02)
             CALL access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
             CALL access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
             
             !! future (3-moment ice)
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi, 1,dum1,dum4,dum5,dum6,f1pr01)
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi, 2,dum1,dum4,dum5,dum6,f1pr02)
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi, 7,dum1,dum4,dum5,dum6,f1pr09)
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi, 8,dum1,dum4,dum5,dum6,f1pr10)
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi,13,dum1,dum4,dum5,dum6,f1pr19)   !mom6-weighted V
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi,14,dum1,dum4,dum5,dum6,f1pr020)   !z_max
             !! CALL access_lookup_table(dumzz,dumjj,dumii,dumi,15,dum1,dum4,dum5,dum6,f1pr021)   !z_min

             !! impose mean ice size bounds (i.e. apply lambda limiters)
             !! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
             nitot(i,k,iice) = min(nitot(i,k,iice),f1pr09*nitot(i,k,iice))
             nitot(i,k,iice) = max(nitot(i,k,iice),f1pr10*nitot(i,k,iice))

             !! adjust Zi IF needed to make sure mu_i is in bounds
             !! zitot(i,k,iice) = min(zitot(i,k,iice),f1pr020)
             !! zitot(i,k,iice) = max(zitot(i,k,iice),f1pr021)

             IF (.not. log_qipresent) THEN
                qiindex = k
             END IF
             log_qipresent = .True.

             Vt_nit(i,k) = f1pr01*rhofaci(i,k)     !number-weighted    fall speed (with density factor)
             Vt_qit(i,k) = f1pr02*rhofaci(i,k)     !mass-weighted  fall speed (with density factor)
          !  Vt_zit(i,k) = f1pr19*rhofaci(i,k)     !moment6-weighted fall speed (with density factor)
             diag_vmi(i,k,iice) = Vt_qit(i,k)      !output fallspeed

          ELSE

             Vt_nit(i,k) = 0.
             Vt_qit(i,k) = 0.
           ! Vt_zit(i,k) = 0.

          END IF !!qitot_not_small

       END DO !! k-loop

       qipresent: IF (log_qipresent) THEN

          nstep = 1

          DO k = qiindex+kdir,kbot,-kdir

            !- weighted fall speed arrays used for sedimentation calculations
            !  (assigned below to highest non-zero level value at lower levels with Vt_x=0)
             V_qit(k) = Vt_qit(i,k)
             V_nit(k) = Vt_nit(i,k)
             !! V_zit(k) = Vt_zit(i,k)

            !--fill in fall speeds levels below lowest level with hydrometeors
             IF (kdir.eq.1) THEN
                IF (k.le.qiindex-kdir) THEN
                   IF (V_qit(k).lt.1.e-10)  V_qit(k) = V_qit(k+kdir)
                   IF (V_nit(k).lt.1.e-10)  V_nit(k) = V_nit(k+kdir)
                 ! IF (V_zit(k).lt.1.e-10)  V_zit(k) = V_zit(k+kdir)
                END IF
             ELSEIF (kdir.eq.-1) THEN
                IF (k.ge.qiindex-kdir) THEN
                   IF (V_qit(k).lt.1.e-10)  V_qit(k) = V_qit(k+kdir)
                   IF (V_nit(k).lt.1.e-10)  V_nit(k) = V_nit(k+kdir)
                 ! IF (V_zit(k).lt.1.e-10)  V_zit(k) = V_zit(k+kdir)
                END IF
             END IF ! kdir

             !! calculate number of split time steps
             rgvm        = max(V_qit(k),V_nit(k))
             !!rgvm        = max(V_zit(k),max(V_qit(k),V_nit(k)))
             nstep       = max(int(rgvm*dt*inv_dzq(i,k)+1.),nstep)
             dum_qit(k)  = gqitot(i,k,iice)*rho(i,k)
             dum_qir(k)  = gqirim(i,k,iice)*rho(i,k)
             dum_bir(k)  = gbirim(i,k,iice)*rho(i,k)
             dum_nit(k)  = gnitot(i,k,iice)*rho(i,k)
             !!dum_zit(k)  = zitot(i,k,iice)*rho(i,k)
             tend_qit(k) = 0.
             tend_qir(k) = 0.
             tend_bir(k) = 0.
             tend_nit(k) = 0.
             !!tend_zit(k) = 0.

          END DO ! k loop

          inv_nstep = 1./real(nstep)

          IF (l_debug .and. nstep.ge.200) THEN
             print*,'micro_p3 L3897 ICE nstep LARGE:',i,nstep
             !!!if (nstep.ge.500) stop
          END IF

          !! calculate sedimentation using first-order upwind method
          tmp1 = 0.
          diag_sflx(i,:) = 0._r8
       
          DO n = 1,nstep

             DO k = kbot,qiindex,kdir
                flux_qit(k) = V_qit(k)*dum_qit(k)
                flux_nit(k) = V_nit(k)*dum_nit(k)
                flux_qir(k) = V_qit(k)*dum_qir(k)
                flux_bir(k) = V_qit(k)*dum_bir(k)
                !!flux_zit(k) = V_zit(k)*dum_zit(k)
                diag_sflx(i,k+1) = diag_sflx(i,k+1) + flux_qit(k) 
             END DO
             tmp1 = tmp1 + flux_qit(kbot)  !sum flux_ at lowest level for averaging over sub-stepping

             !!! top level with hydrometeor present
             k = qiindex
             fluxdiv_qit = flux_qit(k)*inv_dzq(i,k)
             fluxdiv_qir = flux_qir(k)*inv_dzq(i,k)
             fluxdiv_bir = flux_bir(k)*inv_dzq(i,k)
             fluxdiv_nit = flux_nit(k)*inv_dzq(i,k)
             !!fluxdiv_zit = flux_zit(k)*inv_dzq(i,k)

             tend_qit(k) = tend_qit(k) - fluxdiv_qit*inv_nstep*inv_rho(i,k)
             tend_qir(k) = tend_qir(k) - fluxdiv_qir*inv_nstep*inv_rho(i,k)
             tend_bir(k) = tend_bir(k) - fluxdiv_bir*inv_nstep*inv_rho(i,k)
             tend_nit(k) = tend_nit(k) - fluxdiv_nit*inv_nstep*inv_rho(i,k)
             !!tend_zit(k) = tend_zit(k) - fluxdiv_zit*inv_nstep*inv_rho(i,k)

             dum_qit(k) = dum_qit(k) - fluxdiv_qit*dt*inv_nstep
             dum_qir(k) = dum_qir(k) - fluxdiv_qir*dt*inv_nstep
             dum_bir(k) = dum_bir(k) - fluxdiv_bir*dt*inv_nstep
             dum_nit(k) = dum_nit(k) - fluxdiv_nit*dt*inv_nstep
             !!dum_zit(k) = dum_zit(k) - fluxdiv_zit*dt*inv_nstep

             !! loop from sceond to top level of hydrometeor to surface
             DO k = qiindex-kdir,kbot,-kdir
                fluxdiv_qit = (flux_qit(k+kdir) - flux_qit(k))*inv_dzq(i,k)
                fluxdiv_qir = (flux_qir(k+kdir) - flux_qir(k))*inv_dzq(i,k)
                fluxdiv_bir = (flux_bir(k+kdir) - flux_bir(k))*inv_dzq(i,k)
                fluxdiv_nit = (flux_nit(k+kdir) - flux_nit(k))*inv_dzq(i,k)
                !!fluxdiv_zit = (flux_zit(k+kdir) - flux_zit(k))*inv_dzq(i,k)

                tend_qit(k) = tend_qit(k) + fluxdiv_qit*inv_nstep*inv_rho(i,k)
                tend_qir(k) = tend_qir(k) + fluxdiv_qir*inv_nstep*inv_rho(i,k)
                tend_bir(k) = tend_bir(k) + fluxdiv_bir*inv_nstep*inv_rho(i,k)
                tend_nit(k) = tend_nit(k) + fluxdiv_nit*inv_nstep*inv_rho(i,k)
                !!tend_zit(k) = tend_zit(k) + fluxdiv_zit*inv_nstep*inv_rho(i,k)

                dum_qit(k) = dum_qit(k) + fluxdiv_qit*dt*inv_nstep
                dum_qir(k) = dum_qir(k) + fluxdiv_qir*dt*inv_nstep
                dum_bir(k) = dum_bir(k) + fluxdiv_bir*dt*inv_nstep
                dum_nit(k) = dum_nit(k) + fluxdiv_nit*dt*inv_nstep
                !!dum_zit(k) = dum_nit(k) + fluxdiv_nit*dt*inv_nstep
             END DO ! k loop

          END DO ! nstep loop

          !! update prognostic variables with sedimentation tendencies
          DO k = kbot,qiindex,kdir
             gqitot(i,k,iice) = gqitot(i,k,iice) + tend_qit(k)*dt
             gqirim(i,k,iice) = gqirim(i,k,iice) + tend_qir(k)*dt
             gbirim(i,k,iice) = gbirim(i,k,iice) + tend_bir(k)*dt
             gnitot(i,k,iice) = gnitot(i,k,iice) + tend_nit(k)*dt
             !!zitot(i,k,iice) = zitot(i,k,iice) + tend_zit(k)*dt

             diag_sedqi(i,k,iice) = tend_qit(k) !! KZP3
             diag_sedni(i,k,iice) = tend_nit(k) !! KZP3
          
          END DO

          !! add contirubtion from iice to solid precipitation rate at surface
          tmp1 = tmp1*inv_nstep   !flux_ at surface, averaged over sub-step
          tmp1 = tmp1*inv_rhow    !convert flux_ (kg m-2 s-1) to pcp rate (m s-1), liquid-equivalent
          pcprt_sol(i) = pcprt_sol(i) + tmp1  !add pcp rate from

       END IF qipresent

    END DO !! iice_loop_sedi_ice  !iice-loop
 
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 035 -'

    END IF 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! End of Sedimentation 
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pcprt_tot(i) = pcprt_liq(i) + pcprt_sol(i)


!! if(masterproc) then
!!
!! DO k = ktop,kbot
!!
!!    write(iulog,*) '## KZ after sedimentation ## '
!!
!!    IF(maxval(gqc(:,k)).gt.1. .or. minval(gqc(:,k)).lt.1.e-20) then
!!        WRITE(iulog,*) '# 001 max, min of qc : ', k, maxval(gqc(:,k)), minval(gqc(:,k))
!!    END IF
!!    IF(maxval(gqr(:,k)).gt.1. .or. minval(gqr(:,k)).lt.1.e-20) then
!!        WRITE(iulog,*) '# 001 max, min of qr : ', k, maxval(gqr(:,k)), minval(gqr(:,k))
!!    END IF
!!    IF(maxval(gqitot(:,k,1)).gt.1. .or. minval(gqitot(:,k,1)).lt.1.e-20) then
!!        WRITE(iulog,*) '# 001 max, min of qi : ', k, maxval(gqitot(:,k,1)), minval(gqitot(:,k,1))
!!    END IF
!!
!! END DO
!!
!! end if

    !!.......................................................................................
    !! homogeneous freezing of cloud and rain
    !!.......................................................................................

    IF(l_crhomo) THEN 

    DO k = kbot,ktop,kdir

       !! compute mean-mass ice diameters (estimated; rigorous approach to be implemented later)
    
       diam_ice(i,k,:) = 0._r8
       
       DO iice = 1,nCat
          IF (gqitot(i,k,iice).ge.qsmall) THEN
             dum1 = max(gnitot(i,k,iice),nsmall)
             dum2 = 500._r8 !ice density
             diam_ice(i,k,iice) = ((gqitot(i,k,iice)*6.)/(dum1*dum2*pi))**thrd
          END IF
       END DO  !iice loop


       l_dum1 = gqc(i,k).ge.qsmall
       l_dum2 = tt(i,k).lt.233.15_r8 
     
       l_possible = l_dum1 .and. l_dum2
       
       IF (l_possible) THEN
          Q_nuc = gqc(i,k)
          N_nuc = max(gnc(i,k),nsmall)

          !! mg microphysics
          !! assume 25 micron mean volume radius of homogeneously frozen droplets
          !! consistent with size of detrained ice in stratiform.F90
          !!N_nuc = 3._r8*Q_nuc/(4._r8*3.14_r8*1.563e-14_r8*500._r8) !! ### 
          
         !--determine destination ice-phase category:
          dum1   = 900._r8     !density of new ice
          D_new  = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          CALL icecat_destination(qitot(i,k,:),diam_ice(i,k,:),D_new,deltaD_init,       &
                                  log_ni_add,iice_dest)
                                  
          !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 036 - iice_dest 4280 : ', iice_dest 

          gqirim(i,k,iice_dest) = gqirim(i,k,iice_dest) + Q_nuc 
          gqitot(i,k,iice_dest) = gqitot(i,k,iice_dest) + Q_nuc 
          gbirim(i,k,iice_dest) = gbirim(i,k,iice_dest) + Q_nuc * inv_rho_rimeMax
          gnitot(i,k,iice_dest) = gnitot(i,k,iice_dest) + N_nuc 
          
          tt(i,k) = tt(i,k) + Q_nuc * lcldm(i,k) * xlf(i,k) * inv_cp
          
          
          !!IF(l_debug) WRITE(6,*) 'Ttend - qchomfrz qc : ', k, Q_nuc*xlf(i,k)*inv_cp
          
          
          gqc(i,k) = 0._r8  != qc(i,k) - Q_nuc
          gnc(i,k) = 0._r8  != nc(i,k) - N_nuc
          
       END IF

       l_dum1 = gqr(i,k).ge.qsmall
       l_dum2 = tt(i,k).lt.233.15_r8 
     
       l_possible = l_dum1 .and. l_dum2
       
       IF (l_possible) THEN
       
          Q_nuc = gqr(i,k)
          N_nuc = max(gnr(i,k),nsmall)
          
          !! determine destination ice-phase category
          
          dum1  = 900._r8     !! density of new ice
          
          D_new = ((Q_nuc*6.)/(pi*dum1*N_nuc))**thrd
          
          CALL icecat_destination(qitot(i,k,:),      &
                                  diam_ice(i,k,:),   &
                                  D_new,             &
                                  deltaD_init,       &
                                  log_ni_add,        &
                                  iice_dest          ) 

          gqirim(i,k,iice_dest) = gqirim(i,k,iice_dest) + Q_nuc 
          gqitot(i,k,iice_dest) = gqitot(i,k,iice_dest) + Q_nuc 
          gbirim(i,k,iice_dest) = gbirim(i,k,iice_dest) + Q_nuc * inv_rho_rimeMax
          gnitot(i,k,iice_dest) = gnitot(i,k,iice_dest) + N_nuc 
          
          tt(i,k) = tt(i,k) + Q_nuc * rcldm(i,k) * xlf(i,k) * inv_cp
          
          gqr(i,k) = 0._r8  ! = qr(i,k) - Q_nuc
          gnr(i,k) = 0._r8  ! = nr(i,k) - N_nuc
          
       END IF

    END DO !! k_loop_fz
 
    END IF 


    ! remove any excess over-saturation

    IF(l_satadj) THEN 

    DO k = kbot,ktop,kdir

    qtmp=qv(i,k)
    ttmp=tt(i,k)
    ptmp=pres(i,k)

    dum0 = polysvp1(ttmp,0)
    xqvs = ep_2*dum0/max(1.e-3,(ptmp-dum0))

    xqvs = min(xqvs,1._r8)

    if (qtmp > xqvs .and. xqvs > 0) then

       !!write(iulog,*) 'L6623 #1 : i, k, qv, qvs, rh = ', i, k, qv(i,k), xqvs, qv(i,k)/xqvs 
 
       ! expression below is approximate since there may be ice
       ! deposition
       dum = (qtmp-xqvs)/(1._r8+xxlv(i,k)**2*xqvs/(cpp*rv*ttmp**2))/dt

       diag_cmei(i,k,1) = diag_cmei(i,k,1)+dum

       ! partitioning between liquid and ice 
       if (ttmp > 268.15_r8) then
          dum1=0.0_r8
       else if (ttmp < 238.15_r8) then
          dum1=1.0_r8
       else
          dum1=(268.15_r8-ttmp)/30._r8
       end if

       dum = (qtmp-xqvs)/(1._r8+(xxls(i,k)*dum1+xxlv(i,k)*(1._r8-dum1))**2 &
            *xqvs/(cpp*rv*ttmp**2))/dt

       gqc(i,k)=gqc(i,k)+dum*(1._r8-dum1)
       qcadj(i,k)=dum*(1._r8-dum1)

       !! assuming no riming production 

       gqitot(i,k,1)=gqitot(i,k,1)+dum*dum1
       qiadj(i,k)=dum*dum1

       qv(i,k)=qv(i,k)-dum
       qvadj(i,k)=-dum

       tt(i,k)=tt(i,k)+dum*(1._r8-dum1)*xxlv(i,k)+dum*dum1*xxls(i,k)

       write(iulog,*) 'L6623 #2 : i, k, qv, qvs, rh = ', i, k, qv(i,k), xqvs, qv(i,k)/xqvs 
 
    end if

    END DO !! k_loop_fz
 
    END IF !! k_loop_fz
 
!!!!!!!!!!!!

    !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 036  ' 
    
    !!.......................................................................................
    !! final checks to ensure consistency of mass/number
    !! and compute diagnostic fields for output
    !!.......................................................................................

    DO k = kbot,ktop,kdir

       !!....................................
       !! cloud: 
       !!....................................
       
       IF (gqc(i,k).ge.qsmall) THEN
       
          CALL get_cloud_dsd(gqc(i,k), &
                             gnc(i,k), &
                             diag_mu_c(i,k), &
                             rho(i,k), &
                             nu(i,k),  &
                             dnu,      &
                             diag_lamc(i,k), &
                             lammin,   &
                             lammax,   &
                             k,        &
                             tmp1,     &
                             tmp2,     &
                             tmpint1,  &
                             log_tmp1)
                             
          diag_effc(i,k) = 0.5_r8*(diag_mu_c(i,k)+3._r8)/diag_lamc(i,k)
          
       ELSE
       
          qv(i,k) = qv(i,k)+gqc(i,k) 
          tt(i,k) = tt(i,k)-gqc(i,k) * xxlv(i,k) * inv_cp
          
          gqc(i,k) = 0._r8
          gnc(i,k) = 0._r8
          
          diag_effc(i,k) = 10.e-6_r8
          
       END IF

       !!....................................
       !! rain:
       !!....................................
       
       IF (gqr(i,k).ge.qsmall) THEN
       
          CALL get_rain_dsd( gqr(i,k),   &
                             gnr(i,k),   &
                             diag_mu_r(i,k), &
                             rdumii,     &
                             dumii,      &
                             lamr(i,k),  &
                             mu_r_table, &
                             tmp1,       &
                             tmp2,       &
                             log_tmp1,   &
                             tmpint1,    & 
                             tmpint2)
                             
         ! hm, turn off soft lambda limiter
         ! impose size limits for rain with 'soft' lambda limiter
         ! (adjusts over a set timescale rather than within one timestep)
         ! dum2 = (qr(i,k)/(pi*rhow*nr(i,k)))**thrd
         ! IF (dum2.gt.dbrk) THEN
         !    dum   = qr(i,k)*cons4
         !   !dum1  = (dum-nr(i,k))/max(60.,dt)  !time scale for adjustment is 60 s
         !    dum1  = (dum-nr(i,k))*timeScaleFactor
         !     nr(i,k) = nr(i,k)+dum1*dt
         ! END IF

          diag_effr(i,k) = 0.5*(diag_mu_r(i,k)+3.)/lamr(i,k)
        ! ze_rain(i,k) = n0r(i,k)*720./lamr(i,k)**3/lamr(i,k)**3/lamr(i,k)
          ! non-exponential rain:
          ze_rain(i,k) = nr(i,k)*(diag_mu_r(i,k)+6.)*(diag_mu_r(i,k)+5.)*(diag_mu_r(i,k)+4.)*           &
                        (diag_mu_r(i,k)+3.)*(diag_mu_r(i,k)+2.)*(diag_mu_r(i,k)+1.)/lamr(i,k)**6
          ze_rain(i,k) = max(ze_rain(i,k),1.e-22)
          
       ELSE
       
          qv(i,k) = qv(i,k)+gqr(i,k) 
          tt(i,k) = tt(i,k)-gqr(i,k) * xxlv(i,k) * inv_cp
          gqr(i,k) = 0.
          gnr(i,k) = 0.
          diag_effr(i,k) = 25.e-6
          
       END IF

       !!....................................
       !! ice:
       !!....................................

       CALL impose_max_total_Ni(gnitot(i,k,:), &
                                max_total_Ni,  &
                                inv_rho(i,k)   )

       DO iice = 1,nCat

          qi_not_small:  IF (gqitot(i,k,iice).ge.qsmall) THEN

             !! impose lower limits to prevent taking log of # .lt. 0
            
             gnitot(i,k,iice) = max(gnitot(i,k,iice),nsmall)
             gnr(i,k)         = max(gnr(i,k),nsmall)

             CALL calc_bulkRhoRime(gqitot(i,k,iice), &
                                   gqirim(i,k,iice), &
                                   gbirim(i,k,iice), &
                                   rhop              )

             !! IF (.not. tripleMoment_on) THEN 
             !!   zitot(i,k,iice) = diag_mom6(qitot(i,k,iice),nitot(i,k,iice),rho(i,k))
             !! END IF
           
             CALL find_lookupTable_indices_1(dumi,dumj,dumjj,dumii,dumzz,              &
                                          dum1,dum3,dum4,dum5,dum6,                    &
                                          isize,rimsize,densize,zsize,rcollsize,       &
                                          gqr(i,k),gnr(i,k),gqitot(i,k,iice),gnitot(i,k,iice),  &
                                          gqirim(i,k,iice),999._r8,rhop,300)
!                                         qirim(i,k,iice),zitot(i,k,iice),rhop,100)  !future (3-moment)

             CALL access_lookup_table(dumjj,dumii,dumi, 6,dum1,dum4,dum5,f1pr06)
             CALL access_lookup_table(dumjj,dumii,dumi, 7,dum1,dum4,dum5,f1pr09)
             CALL access_lookup_table(dumjj,dumii,dumi, 8,dum1,dum4,dum5,f1pr10)
             CALL access_lookup_table(dumjj,dumii,dumi, 9,dum1,dum4,dum5,f1pr13)
             CALL access_lookup_table(dumjj,dumii,dumi,11,dum1,dum4,dum5,f1pr15)
             CALL access_lookup_table(dumjj,dumii,dumi,12,dum1,dum4,dum5,f1pr16)

             !! impose mean ice size bounds (i.e. apply lambda limiters)
             !! note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
             
             gnitot(i,k,iice) = min(gnitot(i,k,iice),f1pr09*gnitot(i,k,iice))
             gnitot(i,k,iice) = max(gnitot(i,k,iice),f1pr10*gnitot(i,k,iice))

             !! this should already be done in s/r 'calc_bulkRhoRime'
             IF (gqirim(i,k,iice).lt.qsmall) THEN
                gqirim(i,k,iice) = 0.
                gbirim(i,k,iice) = 0.
             END IF

             !! note that reflectivity from lookup table is normalized, so we need to multiply by N
             diag_effi(i,k,iice)  = f1pr06 ! units are in m
             diag_di(i,k,iice)    = f1pr15
             diag_rhopo(i,k,iice) = f1pr16
             
             !! sum contribution from each ice category (note: 0.1892 = 0.176/0.93)
             ze_ice(i,k) = ze_ice(i,k) + 0.1892*f1pr13*gnitot(i,k,iice)*rho(i,k)  !! bugfix 2018-04-23   
             ze_ice(i,k) = max(ze_ice(i,k),1.e-22)

          ELSE

             qv(i,k) = qv(i,k) + gqitot(i,k,iice) 
             tt(i,k) = tt(i,k) - gqitot(i,k,iice) * xxls(i,k)*inv_cp
             
             !!IF(l_debug) WRITE(6,*) 'Ttend - consistency check qr : ', k, - qitot(i,k,iice)*xxls(i,k)*inv_cp 
             
             gqitot(i,k,iice) = 0.
             gnitot(i,k,iice) = 0.
             gqirim(i,k,iice) = 0.
             gbirim(i,k,iice) = 0.

             diag_effi(i,k,iice) = 25.e-6
             diag_di(i,k,iice)   = 0.

          END IF qi_not_small

          !! ice effective diameter
          
          diag_deffi(i,k,iice) = diag_effi(i,k,iice)*diag_rhopo(i,k,iice)/rhows*2._r8
             
       END DO !! iice_loop_final_diagnostics 

       !! sum ze components and convert to dBZ
     
       diag_ze(i,k) = 10.*log10((ze_rain(i,k) + ze_ice(i,k))*1.d+18)

       !! if qr is very small then set Nr to 0 (needs to be done here after call
       !! to ice lookup table because a minimum Nr of nsmall will be set otherwise even IF qr=0) 
       
       IF (gqr(i,k).lt.qsmall) THEN
          gnr(i,k) = 0.
       END IF

    END DO !! k_loop_final_diagnostics
 
 
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 037 -'
    
    !!.......................................................................................
    !! merge ice categories with similar properties
    !!   note:  this should be relocated to above, such that the diagnostic
    !!          ice properties are computed after merging
    !!.......................................................................................

    IF (nCat.gt.1) THEN

       DO k = kbot,ktop,kdir
          DO iice = nCat,2,-1

           ! simility condition (similar mean sizes; similar bulk densities)
           
             IF (abs(diag_di(i,k,iice)-diag_di(i,k,iice-1)).le.150.e-6   .and.           &
                 abs(diag_rhopo(i,k,iice)-diag_rhopo(i,k,iice-1)).le.100._r8) THEN

                gqitot(i,k,iice-1) = gqitot(i,k,iice-1) + gqitot(i,k,iice)
                gnitot(i,k,iice-1) = gnitot(i,k,iice-1) + gnitot(i,k,iice)
                gqirim(i,k,iice-1) = gqirim(i,k,iice-1) + gqirim(i,k,iice)
                gbirim(i,k,iice-1) = gbirim(i,k,iice-1) + gbirim(i,k,iice)
             !  zitot(i,k,iice-1) = zitot(i,k,iice-1) + zitot(i,k,iice)

                gqitot(i,k,iice) = 0.
                gnitot(i,k,iice) = 0.
                gqirim(i,k,iice) = 0.
                gbirim(i,k,iice) = 0.
             !  zitot(i,k,iice) = 0.

             END IF

          END DO !iice loop
       END DO !k loop

    END IF !! multicat
 
    !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 038 -'



#ifdef P3_SLOW

    IF(l_debug) THEN 
       DO k = kbot,ktop,kdir
          WRITE(6,*) ' sS : ', k, maxval(tt(:,k)),   minval(tt(:,k))
          WRITE(6,*) ' qS : ', k, maxval(qv(:,k)),   minval(qv(:,k))
          WRITE(6,*) 'qcS : ', k, maxval(gqc(:,k)),  minval(gqc(:,k))  
          WRITE(6,*) 'ncS : ', k, maxval(gnc(:,k)),  minval(gnc(:,k))  
          WRITE(6,*) 'qrS : ', k, maxval(gqr(:,k)),  minval(gqr(:,k))  
          WRITE(6,*) 'nrS : ', k, maxval(gnr(:,k)),  minval(gnr(:,k))  
          WRITE(6,*) 'qiS : ', k, maxval(gqitot(:,k,1)), minval(gqitot(:,k,1))  
          WRITE(6,*) 'qiS : ', k, maxval(gnitot(:,k,1)), minval(gnitot(:,k,1))  
          WRITE(6,*) 'rmS : ', k, maxval(gqirim(:,k,1)), minval(gqirim(:,k,1))  
          WRITE(6,*) 'rvS : ', k, maxval(gbirim(:,k,1)), minval(gbirim(:,k,1))  
       END DO 
    END IF 

#endif


    DO k = kbot,ktop,kdir
       IF(any(tt(:,k).gt.500._r8) .or. any(tt(:,k).lt.50)) THEN 
          STOP'temperature wrong in p3' 
       END IF
    END DO 


    !!....................................
    !! tendency update 
    !!....................................

    stend (i,ktop:kbot) = (tt(i,ktop:kbot) - ttn(i,ktop:kbot)) * cp / dt 
    qvtend(i,ktop:kbot) = (qv(i,ktop:kbot) - qvn(i,ktop:kbot)) / dt 
    qctend(i,ktop:kbot) = (gqc(i,ktop:kbot) - qcn(i,ktop:kbot)) / dt
    nctend(i,ktop:kbot) = (gnc(i,ktop:kbot) - ncn(i,ktop:kbot)) / dt
    qrtend(i,ktop:kbot) = (gqr(i,ktop:kbot) - qrn(i,ktop:kbot)) / dt
    nrtend(i,ktop:kbot) = (gnr(i,ktop:kbot) - nrn(i,ktop:kbot)) / dt
    qitend(i,ktop:kbot,nCat) = (gqitot(i,ktop:kbot,nCat) - qitotn(i,ktop:kbot,nCat)) / dt
    nitend(i,ktop:kbot,nCat) = (gnitot(i,ktop:kbot,nCat) - nitotn(i,ktop:kbot,nCat)) / dt
    qirimtend(i,ktop:kbot,nCat) = (gqirim(i,ktop:kbot,nCat) - qirimn(i,ktop:kbot,nCat)) / dt
    birimtend(i,ktop:kbot,nCat) = (gbirim(i,ktop:kbot,nCat) - birimn(i,ktop:kbot,nCat)) / dt
      
    !! IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 039 -'

#ifdef P3_SLOW
    DO k = kbot,ktop,kdir
       WRITE(6,*) ' sT : ', k, maxval(stend(:,k)),  minval(stend(:,k))
       WRITE(6,*) ' qT : ', k, maxval(qvtend(:,k)), minval(qvtend(:,k))
       WRITE(6,*) 'qcT : ', k, maxval(qctend(:,k)), minval(qctend(:,k))  
       WRITE(6,*) 'ncT : ', k, maxval(nctend(:,k)), minval(nctend(:,k))  
       WRITE(6,*) 'qrT : ', k, maxval(qrtend(:,k)), minval(qrtend(:,k))  
       WRITE(6,*) 'nrT : ', k, maxval(nrtend(:,k)), minval(nrtend(:,k))  
       WRITE(6,*) 'qiT : ', k, maxval(qitend(:,k,1)),    minval(qitend(:,k,1))  
       WRITE(6,*) 'qiT : ', k, maxval(nitend(:,k,1)),    minval(nitend(:,k,1))  
       WRITE(6,*) 'rmT : ', k, maxval(qirimtend(:,k,1)), minval(qirimtend(:,k,1))  
       WRITE(6,*) 'rvT : ', k, maxval(birimtend(:,k,1)), minval(birimtend(:,k,1))  
    END DO 
#endif


     
333 CONTINUE

    !! recalculate supersaturation from T and qv
  
    IF (log_predictSsat) THEN
       DO k = kbot,ktop,kdir
          !!t(i,k)    = th(i,k)*(1.e-5*pres(i,k))**(rd*inv_cp)
          dum0      = polysvp1(tt(i,k),0)
          dum       = 0.622*dum0/max(1.e-3,(pres(i,k)-dum0))
          ssat(i,k) = qv(i,k)-dum
       END DO
    END IF



!.....................................................


 END DO !! i_loop_main
 
 
 !!IF(l_debug .and. masterproc) WRITE(6,*) 'micro_p3_tend - 040 -'


 RETURN

 END SUBROUTINE micro_p3_tend




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




 SUBROUTINE access_lookup_table(dumjj,dumii,dumi,index,dum1,dum4,dum5,proc)

 IMPLICIT NONE

 REAL(r8) :: dum1,dum4,dum5,proc,dproc1,dproc2,iproc1,gproc1,tmp1,tmp2
 INTEGER :: dumjj,dumii,dumi,index
 REAL(r8) :: xdum1, xdum2, xdum3, xdum4 

   !! get value at current density index

   !! first interpolate for current rimed fraction index

   xdum1 = itab(dumjj,dumii,dumi,index) 
   xdum2 = dum1-real(dumi)
   xdum3 = itab(dumjj,dumii,dumi+1,index) 
   xdum4 = itab(dumjj,dumii,dumi,index)
   
   iproc1 = xdum1+xdum2*(xdum3-xdum4)

   !! linearly interpolate to get process rates for rimed fraction index + 1

   xdum1 = itab(dumjj,dumii+1,dumi,index)
   xdum2 = dum1-real(dumi)
   xdum3 = itab(dumjj,dumii+1,dumi+1,index)
   xdum4 = itab(dumjj,dumii+1,dumi,index)
   
   gproc1 = xdum1+xdum2*(xdum3-xdum4)

   xdum1 = dum4-real(dumii)
   
   tmp1 = iproc1 + xdum1*(gproc1-iproc1)

   !! get value at density index + 1

   !! first interpolate for current rimed fraction index

   xdum1 = itab(dumjj+1,dumii,dumi,index)
   xdum2 = dum1-real(dumi)
   xdum3 = itab(dumjj+1,dumii,dumi+1,index)
   xdum4 = itab(dumjj+1,dumii,dumi,index)
   
   iproc1 = xdum1+xdum2*(xdum3-xdum4)

   !! linearly interpolate to get process rates for rimed fraction index + 1

   xdum1 = itab(dumjj+1,dumii+1,dumi,index)
   xdum2 = dum1-real(dumi)
   xdum3 = itab(dumjj+1,dumii+1,dumi+1,index)
   xdum4 = itab(dumjj+1,dumii+1,dumi,index)
   
   gproc1 = xdum1+xdum2*(xdum3-xdum4)

   xdum1 = dum4-real(dumii)
   
   tmp2 = iproc1 + xdum1*(gproc1-iproc1) 

   !! get final process rate
   
   xdum1 = dum5-real(dumjj)
   
   proc   = tmp1+xdum1*(tmp2-tmp1)

END SUBROUTINE access_lookup_table



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE access_lookup_table_coll(dumjj,dumii,dumj,dumi,index,dum1,dum3,     &
                                    dum4,dum5,proc)

 IMPLICIT NONE

 REAL(r8)    :: dum1,dum3,dum4,dum5,proc,dproc1,dproc2,iproc1,gproc1,tmp1,tmp2,dproc11, &
            dproc12,dproc21,dproc22
 INTEGER :: dumjj,dumii,dumj,dumi,index


! This SUBROUTINE interpolates lookup table values for rain/ice collection processes

! current density index

   dproc1  = itabcoll(dumjj,dumii,dumi,dumj,index)+(dum1-real(dumi))*               &
             (itabcoll(dumjj,dumii,dumi+1,dumj,index)-itabcoll(dumjj,dumii,dumi,    &
             dumj,index))

   dproc2  = itabcoll(dumjj,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj,dumii,dumi+1,dumj+1,index)-itabcoll(dumjj,dumii,dumi,  &
             dumj+1,index))

   iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

   dproc1  = itabcoll(dumjj,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj,dumii+1,dumi+1,dumj,index)-itabcoll(dumjj,dumii+1,     &
                 dumi,dumj,index))

   dproc2  = itabcoll(dumjj,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj,dumii+1,dumi+1,dumj+1,index)-itabcoll(dumjj,dumii+1,   &
             dumi,dumj+1,index))

   gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
   tmp1    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! density index + 1

! current rime fraction index

   dproc1  = itabcoll(dumjj+1,dumii,dumi,dumj,index)+(dum1-real(dumi))*             &
             (itabcoll(dumjj+1,dumii,dumi+1,dumj,index)-itabcoll(dumjj+1,dumii,     &
                 dumi,dumj,index))

   dproc2  = itabcoll(dumjj+1,dumii,dumi,dumj+1,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj+1,dumii,dumi+1,dumj+1,index)-itabcoll(dumjj+1,dumii,   &
             dumi,dumj+1,index))

   iproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)

! rime fraction index + 1

   dproc1  = itabcoll(dumjj+1,dumii+1,dumi,dumj,index)+(dum1-real(dumi))*           &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumj,index)-itabcoll(dumjj+1,dumii+1, &
             dumi,dumj,index))

   dproc2  = itabcoll(dumjj+1,dumii+1,dumi,dumj+1,index)+(dum1-real(dumi))*         &
             (itabcoll(dumjj+1,dumii+1,dumi+1,dumj+1,index)-itabcoll(dumjj+1,       &
                 dumii+1,dumi,dumj+1,index))

   gproc1  = dproc1+(dum3-real(dumj))*(dproc2-dproc1)
   tmp2    = iproc1+(dum4-real(dumii))*(gproc1-iproc1)

! interpolate over density to get final values
   proc    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

 END SUBROUTINE access_lookup_table_coll



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 SUBROUTINE access_lookup_table_colli(dumjjc,dumiic,dumic,dumjj,dumii,dumj,dumi,     &
                                     index,dum1c,dum4c,dum5c,dum1,dum4,dum5,proc)

 IMPLICIT NONE

 REAL(r8)    :: dum1,dum3,dum4,dum5,dum1c,dum4c,dum5c,proc,dproc1,dproc2,iproc1,iproc2, &
            gproc1,gproc2,rproc1,rproc2,tmp1,tmp2,dproc11,dproc12
 INTEGER :: dumjj,dumii,dumj,dumi,index,dumjjc,dumiic,dumic


! This SUBROUTINE interpolates lookup table values for rain/ice collection processes

! current density index collectee category

! current rime fraction index for collectee category

! current density index collector category

! current rime fraction index for collector category

  IF (index.eq.1) THEN

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumii,dumjj)+(dum1c-real(dumic))*    &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumii,dumjj)-                     &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumii,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumii+1,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumii+1,dumjj)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumii,dumjj+1)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumii,dumjj+1)-                   &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumii,dumjj+1)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))*&
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi,dumii+1,dumjj+1)-                 &
             itabcolli1(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj)+(dum1c-real(dumic))*   &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumii,dumjj)-                    &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumii+1,dumjj)-                   &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc1  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!.................................
! collectee density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumii,dumjj)-                   &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj))

!   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)
   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj+1)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumii,dumjj+1)-                   &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumii,dumjj)-                   &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli1(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli1(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc2  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!..........................................................................................
! final process rate interpolation over collectee density

   proc    = rproc1+(dum5c-real(dumjjc))*(rproc2-rproc1)

 ELSE IF (index.eq.2) THEN

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumii,dumjj)+(dum1c-real(dumic))*    &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumii,dumjj)-                     &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumii,dumjj)-                   &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumii+1,dumjj)-                   &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj))

!   dproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)
   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

!   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
!             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj)-                  &
!             itabcolli2(dumic,dumiic,dumjjc,dumi,dumk+1,dumii+1,dumjj))

!   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
!             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj)-                  &
!             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumk+1,dumii+1,dumjj))

!   dproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)


!   iproc2  = dproc1+(dum2-real(dumk))*(dproc2-dproc1)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumii,dumjj+1)+(dum1c-real(dumic))*  &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumii,dumjj+1)-                   &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumii+1,dumjj)-                   &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc1  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!......................................... 
! collectee density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj)+(dum1c-real(dumic))*  &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumii,dumjj)-                   &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic,dumjjc+1,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc1    = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

!.......................................................................................................
! collectee rime fraction + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp1    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

! collector density index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii,dumjj+1))

   iproc1  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

! collector rime fraction index + 1

   dproc11 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi,dumii+1,dumjj+1))

   dproc12 = itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1)+(dum1c-real(dumic))* &
             (itabcolli2(dumic+1,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1)-                  &
             itabcolli2(dumic,dumiic+1,dumjjc+1,dumi+1,dumii+1,dumjj+1))

   iproc2  = dproc11+(dum1-real(dumi))*(dproc12-dproc11)

   tmp2    = iproc1+(dum4-real(dumii))*(iproc2-iproc1)

   gproc2  = tmp1+(dum5-real(dumjj))*(tmp2-tmp1)

   rproc2  = gproc1+(dum4c-real(dumiic))*(gproc2-gproc1)

!..........................................................................................
! final process rate interpolation over collectee density

   proc    = rproc1+(dum5c-real(dumjjc))*(rproc2-rproc1)

 END IF ! index =1 or 2

 END SUBROUTINE access_lookup_table_colli



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 REAL(r8) FUNCTION polysvp1(T,i_type)

!-------------------------------------------
!  COMPUTE SATURATION VAPOR PRESSURE
!  POLYSVP1 RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  i_type REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
!-------------------------------------------

      IMPLICIT NONE

       REAL(r8)    :: DUM,T
      INTEGER :: i_type

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
       REAL(r8) a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
        6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/

! liquid
       REAL(r8) a0,a1,a2,a3,a4,a5,a6,a7,a8

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
        6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
       REAL(r8) dt

!-------------------------------------------

      IF (i_type.EQ.1 .and. T.lt.zerodegc) THEN
! ICE

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+       &
                    a8i*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-9.09718*(273.16/T-1.)-3.56654*                 &
!          log10(273.16/T)+0.876793*(1.-T/273.16)+                        &
!          log10(6.1071))*100.


      ELSEIF (i_type.EQ.0 .or. T.ge.zerodegc) THEN
! LIQUID

!       Flatau formulation:
         dt       = max(-80.,t-273.16)
         polysvp1 = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
         polysvp1 = polysvp1*100.

!       Goff-Gratch formulation:
!        POLYSVP1 = 10.**(-7.90298*(373.16/T-1.)+                         &
!             5.02808*log10(373.16/T)-                                    &
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+                  &
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+                &
!             log10(1013.246))*100.

         END IF


 END FUNCTION polysvp1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  REAL(r8) FUNCTION p3_gamma(X)
!----------------------------------------------------------------------
! THIS ROUTINE CALCULATES THE gamma FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE gamma
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!----------------------------------------------------------------------
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH gamma(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  gamma(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!----------------------------------------------------------------------
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, log, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: I,N
      LOGICAL :: l_parity
       REAL(r8) ::                                                       &
          CONV,EPS,FACT,HALF,ONE,res,sum,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL(r8), DIMENSION(7) :: C
      REAL(r8), DIMENSION(8) :: P
      REAL(r8), DIMENSION(8) :: Q
      REAL(r8), PARAMETER   :: constant1 = 0.9189385332046727417803297

!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      data ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      data XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      data P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      data Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      data C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,                 &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,    &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      l_parity=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF (Y.LE.ZERO) THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        res=Y-Y1
        IF (res.NE.ZERO) THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)l_parity=.True.
          FACT=-PI/SIN(PI*res)
          Y=Y+ONE
        ELSE
          res=XINF
          GOTO 900
        END IF
      END IF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF (Y.LT.EPS) THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF (Y.GE.XMININ) THEN
          res=ONE/Y
        ELSE
          res=XINF
          GOTO 900
        END IF
      ELSEIF (Y.LT.TWELVE) THEN
        Y1=Y
        IF (Y.LT.ONE) THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        END IF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        res=XNUM/XDEN+ONE
        IF (Y1.LT.Y) THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          res=res/Y1
        ELSEIF (Y1.GT.Y) THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            res=res*Y
            Y=Y+ONE
          END DO
        END IF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF (Y.LE.XBIG) THEN
          YSQ=Y*Y
          sum=C(7)
          DO I=1,6
            sum=sum/YSQ+C(I)
          END DO
          sum=sum/Y-Y+constant1
          sum=sum+(Y-HALF)*log(Y)
          res=exp(sum)
        ELSE
          res=XINF
          GOTO 900
        END IF
      END IF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF (l_parity)res=-res
      IF (FACT.NE.ONE)res=FACT/res
  900 p3_gamma=res
      return
! ---------- LAST LINE OF gamma ----------

 END FUNCTION p3_gamma


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


REAL(r8) FUNCTION DERF(X)

 IMPLICIT NONE

 REAL(r8) :: X
 REAL(r8), DIMENSION(0 : 64) :: A, B
 REAL(r8) :: W,T,Y
 INTEGER :: K,I
      data A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,                            &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,                            &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0,                            &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,                            &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      data (B(I), I = 0, 12) /                                 &
         -0.00000000029734388465E0,  0.00000000269776334046E0, &
         -0.00000000640788827665E0, -0.00000001667820132100E0, &
         -0.00000021854388148686E0,  0.00000266246030457984E0, &
          0.00001612722157047886E0, -0.00025616361025506629E0, &
          0.00015380842432375365E0,  0.00815533022524927908E0, &
         -0.01402283663896319337E0, -0.19746892495383021487E0, &
          0.71511720328842845913E0 /
      data (B(I), I = 13, 25) /                                &
         -0.00000000001951073787E0, -0.00000000032302692214E0, &
          0.00000000522461866919E0,  0.00000000342940918551E0, &
         -0.00000035772874310272E0,  0.00000019999935792654E0, &
          0.00002687044575042908E0, -0.00011843240273775776E0, &
         -0.00080991728956032271E0,  0.00661062970502241174E0, &
          0.00909530922354827295E0, -0.20160072778491013140E0, &
          0.51169696718727644908E0 /
      data (B(I), I = 26, 38) /                                &
         0.00000000003147682272E0, -0.00000000048465972408E0,  &
         0.00000000063675740242E0,  0.00000003377623323271E0,  &
        -0.00000015451139637086E0, -0.00000203340624738438E0,  &
         0.00001947204525295057E0,  0.00002854147231653228E0,  &
        -0.00101565063152200272E0,  0.00271187003520095655E0,  &
         0.02328095035422810727E0, -0.16725021123116877197E0,  &
         0.32490054966649436974E0 /
      data (B(I), I = 39, 51) /                                &
         0.00000000002319363370E0, -0.00000000006303206648E0,  &
        -0.00000000264888267434E0,  0.00000002050708040581E0,  &
         0.00000011371857327578E0, -0.00000211211337219663E0,  &
         0.00000368797328322935E0,  0.00009823686253424796E0,  &
        -0.00065860243990455368E0, -0.00075285814895230877E0,  &
         0.02585434424202960464E0, -0.11637092784486193258E0,  &
         0.18267336775296612024E0 /
      data (B(I), I = 52, 64) /                                &
        -0.00000000000367789363E0,  0.00000000020876046746E0,  &
        -0.00000000193319027226E0, -0.00000000435953392472E0,  &
         0.00000018006992266137E0, -0.00000078441223763969E0,  &
        -0.00000675407647949153E0,  0.00008428418334440096E0,  &
        -0.00017604388937031815E0, -0.00239729611435071610E0,  &
         0.02064129023876022970E0, -0.06905562880005864105E0,  &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T +    &
              A(K + 11)) * T + A(K + 12)) * W
      ELSEIF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T +     &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T +     &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T +    &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF = Y

 END FUNCTION DERF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 LOGICAL FUNCTION isnan(arg1)
       REAL(r8),INTENT(in) :: arg1
       isnan=( arg1  .ne. arg1 )
       return
 END FUNCTION isnan



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE icecat_destination(Qi,Di,D_nuc,deltaD_init,log_ni_add,iice_dest)

 !--------------------------------------------------------------------------------------!
 ! Returns the index of the destination ice category into which new ice is nucleated.
 !
 ! New ice will be nucleated into the category in which the existing ice is
 ! closest in size to the ice being nucleated.  The exception is that IF the
 ! size difference between the nucleated ice and existing ice exceeds a threshold
 ! value for all categories, then ice is initiated into a new category.
 !
 ! D_nuc        = mean diameter of new particles being added to a category
 ! D(i)         = mean diameter of particles in category i
 ! diff(i)      = |D(i) - D_nuc|
 ! deltaD_init  = threshold size difference to consider a new (empty) category
 ! mindiff      = minimum of all diff(i) (for non-empty categories)
 !
 ! POSSIBLE CASES                      DESTINATION CATEGORY
 !---------------                      --------------------
 ! case 1:  all empty                  category 1
 ! case 2:  all full                   category with smallest diff
 ! case 3:  partly full
 !  case 3a:  mindiff .lt.  diff_thrs     category with smallest diff
 !  case 3b:  mindiff .ge. diff_thrs     first empty category
 !--------------------------------------------------------------------------------------!

 IMPLICIT NONE

! arguments:
 REAL(r8), INTENT(in), DIMENSION(:) :: Qi,Di
 REAL(r8), INTENT(in)               :: D_nuc,deltaD_init
 INTEGER, INTENT(out)           :: iice_dest
 LOGICAL, INTENT(out)           :: log_ni_add

! local variables:
 LOGICAL                        :: all_full,all_empty
 INTEGER                        :: i_firstEmptyCategory,iice,i_mindiff,n_cat
 REAL(r8)                           :: mindiff,diff
 REAL(r8), PARAMETER               :: qsmall_loc = 1.e-14

 !--------------------------------------------------------------------------------------!

 n_cat      = size(Qi)
 log_ni_add = .True.
 iice_dest  = -99

!-- test:
! iice_dest = 1
! return
!==

 IF (sum(Qi(:))<qsmall_loc) THEN

 !case 1:
    iice_dest = 1
    return

 ELSE

    all_full  = .True.
    all_empty = .false.
    mindiff   = 9.e+9
    i_firstEmptyCategory = 0

    DO iice = 1,n_cat
       IF (Qi(iice) .ge. qsmall_loc) THEN
          all_empty = .false.
          diff      = abs(Di(iice)-D_nuc)
          IF (diff .lt. mindiff) THEN
             mindiff   = diff
             i_mindiff = iice
          END IF
       ELSE
          all_full = .false.
          IF (i_firstEmptyCategory.eq.0) i_firstEmptyCategory = iice
       END IF
    END DO

    IF (all_full) THEN
 !case 2:
       iice_dest = i_mindiff
       IF (mindiff .ge. 100.e-6) log_ni_add=.false.
       return
    ELSE
       IF (mindiff .lt. deltaD_init) THEN
 !case 3a:
          iice_dest = i_mindiff
          return
       ELSE
 !case 3b:
          iice_dest = i_firstEmptyCategory
          return
       END IF
    END IF

 END IF

 print*, 'ERROR in s/r icecat_destination -- made it to end'
 stop


 END SUBROUTINE icecat_destination


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE find_lookupTable_indices_1(dumi,dumj,dumjj,dumii,dumzz,                 &
                                       dum1,dum3,dum4,dum5,dum6,                    &
                                       isize,rimsize,densize,zsize,rcollsize,       &
                                       qr,nr,qitot,nitot,qirim,zitot_in,rhop,src_ind)

!------------------------------------------------------------------------------------------!
! Finds indices in 3D ice lookup table
!------------------------------------------------------------------------------------------!

 IMPLICIT NONE

! arguments:
 INTEGER, INTENT(out) :: dumi,dumj,dumjj,dumii,dumzz
 REAL(r8),    INTENT(out) :: dum1,dum3,dum4,dum5,dum6
 INTEGER, INTENT(in)  :: isize,rimsize,densize,zsize,rcollsize,src_ind
 REAL(r8),    INTENT(in)  :: qr,nr,qitot,nitot,qirim,zitot_in,rhop

! local variables:
 REAL(r8)                 :: dumlr,zitot

!------------------------------------------------------------------------------------------!

           ! find index for qi (normalized ice mass mixing ratio = qitot/nitot)
!             dum1 = (log10(qitot)+16.)/0.70757  !orig
!             dum1 = (log10(qitot)+16.)*1.41328
! we are inverting this equation from the lookup table to solve for i:
! qitot/nitot=261.7**((i+10)*0.1)*1.e-18
             dum1 = (log10(qitot/nitot)+18.)/(0.1*log10(261.7))-10.
             dumi = int(dum1)
             ! set limits (to make sure the calculated index doesn't exceed range of lookup table)
             dum1 = min(dum1,real(isize))
             dum1 = max(dum1,1.)
             dumi = max(1,dumi)
             dumi = min(isize-1,dumi)

           ! find index for scaled mean rain size
           ! IF no rain, then just choose dumj = 1 and DO not calculate rain-ice collection processes
             IF (qr.ge.qsmall) THEN
              ! calculate scaled mean size for consistency with ice lookup table
                dumlr = (qr/(pi*rhow*nr))**thrd
                dum3  = (log10(1.*dumlr)+5.)*10.70415
                dumj  = int(dum3)
              ! set limits
                dum3  = min(dum3,real_rcollsize)
                dum3  = max(dum3,1.)
                dumj  = max(1,dumj)
                dumj  = min(rcollsize-1,dumj)
             ELSE
                dumj  = 1
                dum3  = 1.
             END IF

           ! find index for rime mass fraction
             dum4  = (qirim/qitot)*3. + 1.
             dumii = int(dum4)
             ! set limits
             dum4  = min(dum4,real(rimsize))
             dum4  = max(dum4,1.)
             dumii = max(1,dumii)
             dumii = min(rimsize-1,dumii)

           ! find index for bulk rime density
           ! (account for uneven spacing in lookup table for density)
             IF (rhop.le.650._r8) THEN
                dum5 = (rhop-50._r8)*0.005 + 1.
             ELSE
                dum5 =(rhop-650._r8)*0.004 + 4.
             END IF
             dumjj = int(dum5)
             ! set limits
             dum5  = min(dum5,real(densize))
             dum5  = max(dum5,1.)
             dumjj = max(1,dumjj)
             dumjj = min(densize-1,dumjj)

! ! ! find index for moment6
! !             !invert equation in lookupTable1 that assigns mom6 values
! !             !to index values:  Z_value = 9.**i_Z*1.e-30
! !              dum6  = (log10(zitot)+30._r8)*1.04795
! !              dumzz = int(dum6)
! !              ! set limits
! !              dum6  = min(dum6,real(zsize))
! !              dum6  = max(dum6,1.)
! !              dumzz = max(1,dumzz)
! !              dumzz = min(zsize-1,dumzz)
             dum6  = -99
             dumzz = -99

 END SUBROUTINE find_lookupTable_indices_1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



 SUBROUTINE find_lookupTable_indices_2(dumi,   dumii,   dumjj,  dumic, dumiic, dumjjc,  &
                                       dum1,   dum4,    dum5,   dum1c, dum4c,  dum5c,   &
                                       iisize, rimsize, densize,                        &
                                       qitot_1, qitot_2, nitot_1, nitot_2,                      &
                                       qirim_1, qirim_2, birim_1, birim_2)


!------------------------------------------------------------------------------------------!
! Finds indices in ice-ice interaction lookup table (2)
!------------------------------------------------------------------------------------------!

 IMPLICIT NONE

! arguments:
 INTEGER, INTENT(out) :: dumi,   dumii,   dumjj,  dumic, dumiic, dumjjc
 REAL(r8),INTENT(out) :: dum1,   dum4,    dum5,   dum1c, dum4c,  dum5c
 INTEGER, INTENT(in)  :: iisize, rimsize, densize
 REAL(r8),INTENT(in)  :: qitot_1,qitot_2,nitot_1,nitot_2,qirim_1,qirim_2,birim_1,birim_2

! local variables:
 REAL(r8)                 :: drhop

!------------------------------------------------------------------------------------------!

                    ! find index in lookup table for collector category

                    ! find index for qi (total ice mass mixing ratio)
! replace with new inversion for new lookup table 2 w/ reduced dimensionality
                      dum1 = (log10(qitot_1/nitot_1)+18.)/(0.2*log10(261.7))-5.
                      dumi = int(dum1)
                      dum1 = min(dum1,real(iisize))
                      dum1 = max(dum1,1.)
                      dumi = max(1,dumi)
                      dumi = min(iisize-1,dumi)

   ! note that the code below for finding rime mass fraction and density index is
   ! redundant with code for main ice lookup table and can probably be omitted
   ! for efficiency; for now it is left in

                    ! find index for rime mass fraction
                      dum4  = qirim_1/qitot_1*3. + 1.
                      dumii = int(dum4)
                      dum4  = min(dum4,real(rimsize))
                      dum4  = max(dum4,1.)
                      dumii = max(1,dumii)
                      dumii = min(rimsize-1,dumii)


                    ! find index for bulk rime density
                    ! (account for uneven spacing in lookup table for density)
                    ! bulk rime density
                      IF (birim_1.ge.bsmall) THEN
                         drhop = qirim_1/birim_1
                      ELSE
                         drhop = 0.
                      END IF

                      IF (drhop.le.650._r8) THEN
                         dum5 = (drhop-50._r8)*0.005 + 1.
                      ELSE
                         dum5 =(drhop-650._r8)*0.004 + 4.
                      END IF
                      dumjj = int(dum5)
                      dum5  = min(dum5,real(densize))
                      dum5  = max(dum5,1.)
                      dumjj = max(1,dumjj)
                      dumjj = min(densize-1,dumjj)


                    ! find index in lookup table for collectee category, here 'q' is a scaled q/N
                    ! find index for qi (total ice mass mixing ratio)
                      dum1c = (log10(qitot_2/nitot_2)+18.)/(0.2*log10(261.7))-5.
                      dumic = int(dum1c)
                      dum1c = min(dum1c,real(iisize))
                      dum1c = max(dum1c,1.)
                      dumic = max(1,dumic)
                      dumic = min(iisize-1,dumic)


                    ! find index for rime mass fraction
                      dum4c  = qirim_2/qitot_2*3. + 1.
                      dumiic = int(dum4c)
                      dum4c  = min(dum4c,real(rimsize))
                      dum4c  = max(dum4c,1.)
                      dumiic = max(1,dumiic)
                      dumiic = min(rimsize-1,dumiic)
                    ! calculate predicted bulk rime density
                      IF (birim_2.ge.1.e-15) THEN            !*** NOTE:  change to 'bsmall'
                         drhop = qirim_2/birim_2
                      ELSE
                         drhop = 0.
                      END IF

                    ! find index for bulk rime density
                    ! (account for uneven spacing in lookup table for density)
                      IF (drhop.le.650._r8) THEN
                         dum5c = (drhop-50._r8)*0.005 + 1.
                      ELSE
                         dum5c =(drhop-650._r8)*0.004 + 4.
                      END IF
                      dumjjc = int(dum5c)
                      dum5c  = min(dum5c,real(densize))
                      dum5c  = max(dum5c,1.)
                      dumjjc = max(1,dumjjc)
                      dumjjc = min(densize-1,dumjjc)

 END SUBROUTINE find_lookupTable_indices_2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE find_lookupTable_indices_3(dumii,dumjj,dum1,rdumii,rdumjj,inv_dum3,mu_r,lamr)

!------------------------------------------------------------------------------------------!
! Finds indices in rain lookup table (3)
!------------------------------------------------------------------------------------------!

 IMPLICIT NONE

! arguments:
 INTEGER, INTENT(out) :: dumii,dumjj
 REAL(r8),    INTENT(out) :: dum1,rdumii,rdumjj,inv_dum3
 REAL(r8),    INTENT(in)  :: mu_r,lamr

!------------------------------------------------------------------------------------------!

        ! find location in scaled mean size space
          dum1 = (mu_r+1.)/lamr
          IF (dum1.le.195.e-6) THEN
             inv_dum3  = 0.1
             rdumii = (dum1*1.e6+5.)*inv_dum3
             rdumii = max(rdumii, 1.)
             rdumii = min(rdumii,20._r8)
             dumii  = int(rdumii)
             dumii  = max(dumii, 1)
             dumii  = min(dumii,20)
          ELSEIF (dum1.gt.195.e-6) THEN
             inv_dum3  = thrd*0.1            !i.e. 1/30
             rdumii = (dum1*1.e+6-195.)*inv_dum3 + 20.
             rdumii = max(rdumii, 20._r8)
             rdumii = min(rdumii,300._r8)
             dumii  = int(rdumii)
             dumii  = max(dumii, 20)
             dumii  = min(dumii,299)
          END IF

        ! find location in mu_r space
          rdumjj = mu_r+1.
          rdumjj = max(rdumjj,1.)
          rdumjj = min(rdumjj,10._r8)
          dumjj  = int(rdumjj)
          dumjj  = max(dumjj,1)
          dumjj  = min(dumjj,9)

 END SUBROUTINE find_lookupTable_indices_3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE get_cloud_dsd(qc,nc,mu_c,rho,nu,dnu,lamc,lammin,lammax,k,cdist, &
                          cdist1,qcindex,log_qcpresent)

 IMPLICIT NONE

!arguments:
 REAL(r8), DIMENSION(:), INTENT(in)  :: dnu
 REAL(r8), INTENT(in)            :: qc,rho
 REAL(r8), INTENT(inout)         :: nc
 REAL(r8), INTENT(out)           :: mu_c,nu,lamc,cdist,cdist1
 INTEGER,  INTENT(in)            :: k
 INTEGER, INTENT(out)            :: qcindex
 LOGICAL, INTENT(inout)          :: log_qcpresent

!local variables
 REAL(r8)                        :: lammin,lammax
 INTEGER                         :: dumi

!--------------------------------------------------------------------------

       IF (qc.ge.qsmall) THEN

        ! set minimum nc to prevent floating point error
          nc   = max(nc,nsmall)
          mu_c = 0.0005714*(nc*1.e-6*rho)+0.2714
          mu_c = 1./(mu_c**2)-1.
          mu_c = max(mu_c,2.)
          mu_c = min(mu_c,15.)

        ! interpolate for mass distribution spectral shape parameter (for SB warm processes)
          IF (iparam.eq.1) THEN
             dumi = int(mu_c)
             nu   = dnu(dumi)+(dnu(dumi+1)-dnu(dumi))*(mu_c-dumi)
          END IF

        ! calculate lamc
          lamc = (cons1*nc*(mu_c+3.)*(mu_c+2.)*(mu_c+1.)/qc)**thrd

        ! apply lambda limiters
          lammin = (mu_c+1.)*2.5e+4   ! min: 40 micron mean diameter
          lammax = (mu_c+1.)*1.e+6    ! max:  1 micron mean diameter

          IF (lamc.lt.lammin) THEN
             lamc = lammin
             nc   = 6.*lamc**3*qc/(pi*rhow*(mu_c+3.)*(mu_c+2.)*(mu_c+1.))
          ELSEIF (lamc.gt.lammax) THEN
             lamc = lammax
             nc   = 6.*lamc**3*qc/(pi*rhow*(mu_c+3.)*(mu_c+2.)*(mu_c+1.))
             IF (.not. log_qcpresent) THEN
                qcindex = k
             END IF
             log_qcpresent = .True.

          END IF

          cdist  = nc*(mu_c+1.)/lamc
          cdist1 = nc/gamma(mu_c+1.)

       ELSE

          lamc   = 0.
          cdist  = 0.
          cdist1 = 0.

       END IF

 END SUBROUTINE get_cloud_dsd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE get_rain_dsd(qr,nr,mu_r,rdumii,dumii,lamr,mu_r_table,cdistr,logn0r, &
                         log_qrpresent,qrindex,k)

! Computes and returns rain size distribution parameters

 IMPLICIT NONE

!arguments:
 REAL(r8), DIMENSION(:), INTENT(in)  :: mu_r_table
 REAL(r8), INTENT(in)            :: qr
 REAL(r8), INTENT(inout)         :: nr
 REAL(r8), INTENT(out)           :: rdumii,lamr,mu_r,cdistr,logn0r
 INTEGER,  INTENT(in)            :: k
 INTEGER,  INTENT(out)           :: dumii,qrindex
 LOGICAL,  INTENT(inout)         :: log_qrpresent

!local variables:
 REAL(r8)                            :: inv_dum,lammax,lammin

!--------------------------------------------------------------------------

       IF (qr.ge.qsmall) THEN

       ! USE lookup table to get mu
       ! mu-lambda relationship is from Cao et al. (2008), eq. (7)

       ! find spot in lookup table
       ! (scaled N/q for lookup table parameter space_
          nr      = max(nr,nsmall)
          inv_dum = (qr/(cons1*nr*6.))**thrd

          IF (inv_dum.lt.282.e-6) THEN
             mu_r = 8.282
          ELSEIF (inv_dum.ge.282.e-6 .and. inv_dum.lt.502.e-6) THEN
           ! interpolate
             rdumii   = (inv_dum-250.e-6)*1.e+6*0.5
             rdumii   = max(rdumii,1.)
             rdumii   = min(rdumii,150._r8)
             dumii    = int(rdumii)
             dumii    = min(149,dumii)
             mu_r     = mu_r_table(dumii)+(mu_r_table(dumii+1)-mu_r_table(dumii))*(rdumii-  &
                        real(dumii))
          ELSEIF (inv_dum.ge.502.e-6) THEN
             mu_r = 0.
          END IF

          lamr = (cons1*nr*(mu_r+3.)*(mu_r+2)*(mu_r+1.)/(qr))**thrd  ! recalculate slope based on mu_r
          lammax = (mu_r+1.)*1.e+5   ! check for slope
          lammin = (mu_r+1.)*1250._r8   ! set to small value since breakup is explicitly included (mean size 0.8 mm)

        ! apply lambda limiters for rain
          IF (lamr.lt.lammin) THEN
             lamr = lammin
             nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(cons1)
          ELSEIF (lamr.gt.lammax) THEN
             lamr = lammax
             nr   = exp(3.*log(lamr)+log(qr)+log(gamma(mu_r+1.))-log(gamma(mu_r+4.)))/(cons1)
          END IF

          IF (.not. log_qrpresent) THEN
             qrindex = k
          END IF
          log_qrpresent = .True.

          cdistr = nr/gamma(mu_r+1.)
          logn0r    = log10(nr)+(mu_r+1.)*log10(lamr)-log10(gamma(mu_r+1)) !note: logn0r is calculated as log10(n0r)

       ELSE

          lamr = 0.
          cdistr = 0.
          logn0r = 0.

       END IF


 END SUBROUTINE get_rain_dsd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE calc_bulkRhoRime(qi_tot,qi_rim,bi_rim,rho_rime)

!--------------------------------------------------------------------------------
!  Calculates and returns the bulk rime density from the prognostic ice variables
!  and adjusts qirim and birim appropriately.
!--------------------------------------------------------------------------------

 IMPLICIT NONE

!arguments:
 REAL(r8), INTENT(in)    :: qi_tot
 REAL(r8), INTENT(inout) :: qi_rim,bi_rim
 REAL(r8), INTENT(out)   :: rho_rime

 !--------------------------------------------------------------------------

 IF (bi_rim.ge.1.e-15) THEN
!if (bi_rim.ge.bsmall) THEN
    rho_rime = qi_rim/bi_rim
    !impose limits on rho_rime;  adjust bi_rim IF needed
    IF (rho_rime.lt.rho_rimeMin) THEN
       rho_rime = rho_rimeMin
       bi_rim   = qi_rim/rho_rime
    ELSEIF (rho_rime.gt.rho_rimeMax) THEN
       rho_rime = rho_rimeMax
       bi_rim   = qi_rim/rho_rime
    END IF
 ELSE
    qi_rim   = 0.
    bi_rim   = 0.
    rho_rime = 0.
 END IF

 !set upper constraint qi_rim <= qi_tot
 IF (qi_rim.gt.qi_tot .and. rho_rime .gt. 0._r8) THEN
    qi_rim = qi_tot
    bi_rim = qi_rim/rho_rime
 END IF

 !impose consistency
 IF (qi_rim.lt.qsmall) THEN
    qi_rim = 0.
    bi_rim = 0.
 END IF


 END SUBROUTINE calc_bulkRhoRime


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE impose_max_total_Ni(nitot_local,max_total_Ni,inv_rho_local)

!--------------------------------------------------------------------------------
! Impose maximum total ice number concentration (total of all ice categories).
! IF the sum of all nitot(:) exceeds maximum allowable, each category to preserve
! ratio of number between categories.
!--------------------------------------------------------------------------------

 IMPLICIT NONE

!arguments:
 REAL(r8), INTENT(inout), DIMENSION(:) :: nitot_local           !note: dimension (nCat)
 REAL(r8), INTENT(in)                  :: max_total_Ni,inv_rho_local

!local variables:
 REAL(r8)                              :: dum

 IF (sum(nitot_local(:)).ge.1.e-20) THEN
    dum = max_total_Ni*inv_rho_local/sum(nitot_local(:))
    nitot_local(:) = nitot_local(:)*min(dum,1.)
 END IF

 END SUBROUTINE impose_max_total_Ni



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE micro_p3_get_cols(ncol, nlev, top_lev, qcn, qin, &
                             qrn, mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  INTEGER, INTENT(in) :: ncol      ! Number of columns with meaningful data
  INTEGER, INTENT(in) :: nlev      ! Number of levels to use
  INTEGER, INTENT(in) :: top_lev   ! Top level for microphysics

  REAL(r8), INTENT(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  REAL(r8), INTENT(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)
  REAL(r8), INTENT(in) :: qrn(:,:) ! rain mixing ratio (kg/kg)

  INTEGER, INTENT(out) :: mgncol   ! Number of columns MG will use
  INTEGER, ALLOCATABLE, INTENT(out) :: mgcols(:) ! column indices

  INTEGER :: lev_offset  ! top_lev - 1 (defined here for consistency)
  LOGICAL :: ltrue(ncol) ! store tests for each column

  INTEGER :: i, ii ! column indices

  IF (ALLOCATEd(mgcols)) DEALLOCATE(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know IF water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) .ge. qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) .ge. qsmall, 2)
  ltrue = ltrue .or. any(qrn(:ncol,top_lev:(nlev+lev_offset)) .ge. qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  ALLOCATE(mgcols(mgncol))
  i = 0
  DO ii = 1,ncol
     IF (ltrue(ii)) THEN
        i = i + 1
        mgcols(i) = ii
     END IF
  END DO

END SUBROUTINE micro_p3_get_cols


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE micro_p3
       
