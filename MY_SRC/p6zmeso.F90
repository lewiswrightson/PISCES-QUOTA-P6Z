MODULE p6zmeso
   !!======================================================================
   !!                         ***  MODULE p6zmeso  ***
   !! TOP :   PISCES-QUOTA Compute the sources/sinks for mesozooplankton
   !!         Including Explicit Diazotroph PFT
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_meso       : Compute the sources/sinks for mesozooplankton
   !!   p6z_meso_init  : Initialization of the parameters for mesozooplankton
   !!   p6z_meso_alloc : Allocate variables for mesozooplankton 
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p6z_meso              ! called in p4zbio.F90
   PUBLIC   p6z_meso_init         ! called in trcsms_pisces.F90
   PUBLIC   p6z_meso_alloc        ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for zooplankton
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for Diatoms 
   REAL(wp), PUBLIC ::  xpref2m      !: mesozoo preference for mesozoo
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2mes  !: mesozoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::  unass2c      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2n      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  unass2p      !: Non-assimilated fraction of food
   REAL(wp), PUBLIC ::  epsher2      !: Growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  epsher2min   !: Minimum growth efficiency of mesozoo
   REAL(wp), PUBLIC ::  ssigma2      !: Fraction excreted as semi-labile DOM
   REAL(wp), PUBLIC ::  srespir2     !: Active respiration
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate
   REAL(wp), PUBLIC ::  xfracmig     !: Fractional biomass of meso that performs DVM
   REAL(wp), PUBLIC ::  xsigma2      !: Width of the predation window
   REAL(wp), PUBLIC ::  xsigma2del   !: Maximum width of the predation window at low food density
   REAL(wp), PUBLIC ::  xpref2dz     !: mesozoo preference diazotrophs
   REAL(wp), PUBLIC ::  xthresh2dz   !: diazotroph feeding threshold for mesozooplankton
   LOGICAL,  PUBLIC ::  bmetexc2     !: Use of excess carbon for respiration
   LOGICAL , PUBLIC ::  ln_dvm_meso  !: Boolean to activate DVM of mesozooplankton
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) :: depmig  !: DVM of mesozooplankton : migration depth
   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:,:) :: kmig    !: Vertical indice of the the migration depth

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zmeso.F90 13234 2020-07-02 18:40:58Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!                This includes ingestion and assimilation, flux feeding
      !!                and mortality. We use an active prey switching  
      !!                parameterization Morozov and Petrovskii (2013). 
      !!                All living compartments and mesozooplankton
      !!                are potential preys of mesozooplankton as well as small
      !!                sinking particles 
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  :: ji, jj, jk, jkt
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam, zcompames
      REAL(wp) :: zgraze2, zdenom, zfact, zfood, zfoodlim, zproport
      REAL(wp) :: zmortzgoc, zfracc, zfracn, zfracp, zfracfe, zratio, zratio2
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zrespirc, zrespirn, zrespirp, zbasresb, zbasresi
      REAL(wp) :: zgraztotc, zgraztotn, zgraztotp, zgraztotf, zbasresn, zbasresp, zbasresf
      REAL(wp) :: zgratmp, zgradoct, zgradont, zgrareft, zgradopt
      REAL(wp) :: zprcaca, zmortz, zexcess
      REAL(wp) :: zbeta, zrespz, ztortz, zgrasratp, zgrasratn, zgrasratf
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot
      REAL(wp) :: zgrazdc, zgrazz, zgrazm, zgrazpof, zgrazcal, zfracal
      REAL(wp) :: zgraznc, zgrazpoc, zgrazpon, zgrazpop, zgraznf, zgrazdf
      REAL(wp) :: zgraznp, zgraznn, zgrazdn, zgrazdp
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      REAL(wp) :: zgrazffnp, zgrazffng, zgrazffpp, zgrazffpg
      REAL(wp) :: zmigreltime, zrum, zcodel, zargu, zva
      ! Explicit Diazotroph PFT
      REAL(wp) :: zcompadz, ztmp6, zgrazdzc, zgrazdzn, zgrazdzp, zgrazdzf, zdiffdzn, zdiffdzd
      CHARACTER (len=25) :: charout
      REAL(wp) :: zrfact2, zmetexcess, zsigma, zdiffdn
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing, zfezoo2, znitzoo2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrarem, zgraref, zgrapoc, zgrapof
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrarep, zgraren, zgrapon, zgrapop
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgradoc, zgradon, zgradop
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zgramigrem, zgramigref, zgramigpoc, zgramigpof, zstrn
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zgramigrep, zgramigren, zgramigpop, zgramigpon
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zgramigdoc, zgramigdop, zgramigdon
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d, zz2ligprod

      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_meso')
      !
      ! Initialization of local arrays
      zgrazing(:,:,:) = 0._wp  ;  zfezoo2(:,:,:) = 0._wp
      zgrarem (:,:,:) = 0._wp  ;  zgraren(:,:,:) = 0._wp
      zgrarep (:,:,:) = 0._wp  ;  zgraref(:,:,:) = 0._wp
      zgrapoc (:,:,:) = 0._wp  ;  zgrapon(:,:,:) = 0._wp
      zgrapop (:,:,:) = 0._wp  ;  zgrapof(:,:,:) = 0._wp
      zgradoc (:,:,:) = 0._wp  ;  zgradon(:,:,:) = 0._wp
      zgradop (:,:,:) = 0._wp  ;  znitzoo2(:,:,:) = 0._wp
      !
      IF (ln_ligand) THEN
         ALLOCATE( zz2ligprod(jpi,jpj,jpk) )
         zz2ligprod(:,:,:) = 0._wp
      ENDIF

      !
      ! Diurnal vertical migration of mesozooplankton
      ! Computation of the migration depth
      ! ---------------------------------------------
      IF (ln_dvm_meso) CALL p6z_meso_depmig

      ! Use of excess carbon for metabolism
      zmetexcess = 0.0
      IF ( bmetexc2 ) zmetexcess = 1.0

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompam   = MAX( ( trb(ji,jj,jk,jpmes) - 1.e-9 ), 0.e0 )
               zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

               !  linear mortality of mesozooplankton
               !  A michaelis menten modulation term is used to avoid extinction of 
               !  mesozooplankton at very low food concentrations
               !  -----------------------------------------------------------------
               zrespz   = resrat2 * zfact * ( trb(ji,jj,jk,jpmes) / ( xkmort + trb(ji,jj,jk,jpmes) )  &
               &          + 3. * nitrfac(ji,jj,jk) )

               !  Zooplankton quadratic mortality. A square function has been selected with
               !  to mimic predation and disease (density dependent mortality). It also tends
               !  to stabilise the model
               !  -------------------------------------------------------------------------
               ztortz   = mzrat2 * 1.e6 * zfact * trb(ji,jj,jk,jpmes) * (1. - nitrfac(ji,jj,jk))

               !   Computation of the abundance of the preys
               !   A threshold can be specified in the namelist
               !   --------------------------------------------
               zcompadi  = MAX( ( trb(ji,jj,jk,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,jk,jpzoo) - xthresh2zoo ), 0.e0 )
               zcompaph  = MAX( ( trb(ji,jj,jk,jpphy) - xthresh2phy ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,jk,jppoc) - xthresh2poc ), 0.e0 )
               zcompames = MAX( ( trb(ji,jj,jk,jpmes) - xthresh2mes ), 0.e0 )
               zcompadz  = MAX( ( trb(ji,jj,jk,jpcdz) - xthresh2dz  ), 0.e0 )
               !   Mesozooplankton grazing
               ! The total amount of food is the sum of all preys accessible to mesozooplankton 
               ! multiplied by their food preference
               ! A threshold can be specified in the namelist (xthresh2). However, when food 
               ! concentration is close to this threshold, it is decreased to avoid the 
               ! accumulation of food in the mesozoopelagic domain
               ! -------------------------------------------------------------------------------
               zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc   &
               &           + xpref2m * zcompames + xpref2dz * zcompadz
               zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpmes) * (1. - nitrfac(ji,jj,jk)) 

               ! An active switching parameterization is used here.
               ! We don't use the KTW parameterization proposed by 
               ! Vallina et al. because it tends to produce too steady biomass
               ! composition and the variance of Chl is too low as it grazes
               ! too strongly on winning organisms. We use a generalized
               ! switching parameterization proposed by Morozov and 
               ! Petrovskii (2013)
               ! ------------------------------------------------------------  
               ! The width of the selection window is increased when preys
               ! have low abundance, .i.e. zooplankton become less specific 
               ! to avoid starvation.
               ! ----------------------------------------------------------
               zsigma = 1.0 - zdenom**2/(0.05**2+zdenom**2)
               zsigma = xsigma2 + xsigma2del * zsigma
               ! Nanophytoplankton, diatoms and diazotrophs are considered
               ! to be close enough to have potential interference
               ! -----------------------------------------------------------
               zdiffdn = exp( -ABS(log(3.0 * sizen(ji,jj,jk) / (5.0 * sized(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               zdiffdzn = exp( -ABS(log(3.0 * sizen(ji,jj,jk) / (3.0 * sizedz(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               zdiffdzd = exp( -ABS(log(3.0 * sizedz(ji,jj,jk) / (5.0 * sized(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               ztmp1 = xpref2n * zcompaph * ( zcompaph + zdiffdn * zcompadi + zdiffdzn * zcompadz ) / (1.0 + zdiffdn + zdiffdzn )
               ztmp2 = xpref2m * zcompames**2
               ztmp3 = xpref2c * zcompapoc**2
               ztmp4 = xpref2d * zcompadi * ( zdiffdn * zcompadi + zcompaph + zdiffdzd * zcompadz) / (1.0 + zdiffdn + zdiffdzd )
               ztmp5 = xpref2z * zcompaz**2
               ztmp6 = xpref2dz * zcompadz * ( zcompadz + zdiffdzn * zcompaph + zdiffdzd * zcompadi ) / ( 1.0 + zdiffdzd + zdiffdzn )
               ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + ztmp6 + rtrn
               ztmp1 = ztmp1 / ztmptot
               ztmp2 = ztmp2 / ztmptot
               ztmp3 = ztmp3 / ztmptot
               ztmp4 = ztmp4 / ztmptot
               ztmp5 = ztmp5 / ztmptot
               ztmp6 = ztmp6 / ztmptot

               !   Mesozooplankton regular grazing on the different preys
               !   ------------------------------------------------------
               zgrazdc   = zgraze2 * ztmp4 * zdenom
               zgrazdn   = zgrazdc * trb(ji,jj,jk,jpndi) / ( trb(ji,jj,jk,jpdia) + rtrn)
               zgrazdp   = zgrazdc * trb(ji,jj,jk,jppdi) / ( trb(ji,jj,jk,jpdia) + rtrn)
               zgrazdf   = zgrazdc * trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn)
               zgrazz    = zgraze2 * ztmp5 * zdenom
               zgrazm    = zgraze2 * ztmp2 * zdenom
               zgraznc   = zgraze2 * ztmp1 * zdenom
               zgraznn   = zgraznc * trb(ji,jj,jk,jpnph) / ( trb(ji,jj,jk,jpphy) + rtrn)
               zgraznp   = zgraznc * trb(ji,jj,jk,jppph) / ( trb(ji,jj,jk,jpphy) + rtrn)
               zgraznf   = zgraznc * trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) + rtrn)
               zgrazpoc  = zgraze2 * ztmp3 * zdenom
               zgrazpon  = zgrazpoc * trb(ji,jj,jk,jppon) / ( trb(ji,jj,jk,jppoc) + rtrn)
               zgrazpop  = zgrazpoc * trb(ji,jj,jk,jppop) / ( trb(ji,jj,jk,jppoc) + rtrn)
               zgrazpof  = zgrazpoc * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn)
               ! Explicit Diazotroph
               zgrazdzc  = zgraze2 * ztmp6 * zdenom
               zgrazdzn  = zgrazdzc * trb(ji,jj,jk,jpndz) / ( trb(ji,jj,jk,jpcdz) + rtrn)
               zgrazdzp  = zgrazdzc * trb(ji,jj,jk,jppdz) / ( trb(ji,jj,jk,jpcdz) + rtrn)
               zgrazdzf  = zgrazdzc * trb(ji,jj,jk,jpfed) / ( trb(ji,jj,jk,jpcdz) + rtrn)
               !  Mesozooplankton flux feeding on GOC and POC. The feeding pressure
               ! is proportional to the flux
               !  ------------------------------------------------------------------
               zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes)  &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffg = zgrazffeg * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)
               zgrazffng = zgrazffeg * trb(ji,jj,jk,jpgon) / (trb(ji,jj,jk,jpgoc) + rtrn)
               zgrazffpg = zgrazffeg * trb(ji,jj,jk,jpgop) / (trb(ji,jj,jk,jpgoc) + rtrn)
               zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)     &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpmes)   &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazfffp = zgrazffep * trb(ji,jj,jk,jpsfe) / (trb(ji,jj,jk,jppoc) + rtrn)
               zgrazffnp = zgrazffep * trb(ji,jj,jk,jppon) / (trb(ji,jj,jk,jppoc) + rtrn)
               zgrazffpp = zgrazffep * trb(ji,jj,jk,jppop) / (trb(ji,jj,jk,jppoc) + rtrn)
               !
               zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg &
               &            + zgrazdzc

               ! Compute the proportion of filter feeders. It is assumed steady state.
               ! ---------------------------------------------------------------------  
               zproport  = (zgrazffep + zgrazffeg)/(rtrn + zgraztotc)

               !   Compute fractionation of aggregates. It is assumed that 
               !   diatoms based aggregates are more prone to fractionation
               !   since they are more porous (marine snow instead of fecal pellets)
               !   ----------------------------------------------------------------
               zratio    = trb(ji,jj,jk,jpgsi) / ( trb(ji,jj,jk,jpgoc) + rtrn )
               zratio2   = zratio * zratio
               zfracc    = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &          * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes)          &
               &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
               zfracfe   = zfracc * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)
               zfracn    = zfracc * trb(ji,jj,jk,jpgon) / (trb(ji,jj,jk,jpgoc) + rtrn)
               zfracp    = zfracc * trb(ji,jj,jk,jpgop) / (trb(ji,jj,jk,jpgoc) + rtrn)

               ! Flux feeding is multiplied by the fractional biomass of flux feeders
               zgrazffep = zproport * zgrazffep   ;   zgrazffeg = zproport * zgrazffeg
               zgrazfffp = zproport * zgrazfffp   ;   zgrazfffg = zproport * zgrazfffg
               zgrazffnp = zproport * zgrazffnp   ;   zgrazffng = zproport * zgrazffng
               zgrazffpp = zproport * zgrazffpp   ;   zgrazffpg = zproport * zgrazffpg

               zgraztotc  = zgrazdc + zgrazz + zgraznc + zgrazm + zgrazpoc + zgrazffep + zgrazffeg + zgrazdzc
               zgraztotf  = zgrazdf + zgraznf + ( zgrazz + zgrazm ) * ferat3 + zgrazpof &
               &            + zgrazfffp + zgrazfffg + zgrazdzf
               zgraztotn  = zgrazdn + (zgrazm + zgrazz) * no3rat3 + zgraznn + zgrazpon  &
               &            + zgrazffnp + zgrazffng + zgrazdzn
               zgraztotp  = zgrazdp + (zgrazz + zgrazm) * po4rat3 + zgraznp + zgrazpop  &
               &            + zgrazffpp + zgrazffpg + zgrazdzp

               ! Total grazing ( grazing by microzoo is already computed in p6zmicro )
               zgrazing(ji,jj,jk) = zgraztotc
               !   Stoichiometruc ratios of the food ingested by zooplanton 
               !   --------------------------------------------------------
               zgrasratf  =  (zgraztotf + rtrn) / ( zgraztotc + rtrn )
               zgrasratn  =  (zgraztotn + rtrn) / ( zgraztotc + rtrn )
               zgrasratp  =  (zgraztotp + rtrn) / ( zgraztotc + rtrn )

               ! Mesozooplankton efficiency. 
               ! We adopt a formulation proposed by Mitra et al. (2007)
               ! The gross growth efficiency is controled by the most limiting nutrient.
               ! Growth is also further decreased when the food quality is poor. This is currently
               ! hard coded : it can be decreased by up to 50% (zepsherq)
               ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
               ! Fulton, 2012)
               ! -----------------------------------------------------------------------------------
               zepshert  = MIN( 1., zgrasratn/ no3rat3, zgrasratp/ po4rat3, zgrasratf / ferat3)
               zbeta     = MAX(0., (epsher2 - epsher2min) )
               zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
               zepsherv  = zepsherf * zepshert * zepsherq
               !   Respiration of mesozooplankton
               !   Excess carbon in the food is used preferentially
               !   when bmetexc2 is set to .true.
               !   -----------------------------------------------
               zexcess  = zgraztotc * zepsherf * (1.0 - zepshert) * zmetexcess 
               zbasresb = MAX(0., zrespz - zexcess)
               zbasresi = zexcess + MIN(0., zrespz - zexcess)
               zrespirc = srespir2 * zepsherv * zgraztotc + zbasresb

               !   When excess carbon is used, the other elements in excess
               !   are also used proportionally to their abundance
               !   --------------------------------------------------------
               zexcess  = ( zgrasratn/ no3rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresn = zbasresi * zexcess * zgrasratn
               zexcess  = ( zgrasratp/ po4rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresp = zbasresi * zexcess * zgrasratp
               zexcess  = ( zgrasratf/ ferat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresf = zbasresi * zexcess * zgrasratf

               !   Voiding of the excessive elements as organic matter
               !   --------------------------------------------------------
               zgradoct = (1. - unass2c - zepsherv) * zgraztotc - zbasresi
               zgradont = (1. - unass2n) * zgraztotn - zepsherv * no3rat3 * zgraztotc - zbasresn
               zgradopt = (1. - unass2p) * zgraztotp - zepsherv * po4rat3 * zgraztotc - zbasresp
               zgrareft = (1. - unass2c) * zgraztotf - zepsherv * ferat3 * zgraztotc - zbasresf
               ztmp1   = ( 1. - epsher2 - unass2c ) /( 1. - 0.8 * epsher2 ) * ztortz
               zgradoc(ji,jj,jk) = (zgradoct + ztmp1) * ssigma2
               zgradon(ji,jj,jk) = (zgradont + no3rat3 * ztmp1) * ssigma2
               zgradop(ji,jj,jk) = (zgradopt + po4rat3 * ztmp1) * ssigma2
               zgratmp = 0.2 * epsher2 /( 1. - 0.8 * epsher2 ) * ztortz

               !  Since only semilabile DOM is represented in PISCES
               !  part of DOM is in fact labile and is then released
               !  as dissolved inorganic compounds (ssigma2)
               !  --------------------------------------------------
               zgrarem(ji,jj,jk) = zgratmp + ( zgradoct + ztmp1 ) * (1.0 - ssigma2)
               zgraren(ji,jj,jk) = no3rat3 * zgratmp + ( zgradont + no3rat3 * ztmp1 ) * (1.0 - ssigma2)
               zgrarep(ji,jj,jk) = po4rat3 * zgratmp + ( zgradopt + po4rat3 * ztmp1 ) * (1.0 - ssigma2)
               zgraref(ji,jj,jk) = zgrareft + ferat3 * ( ztmp1 + zgratmp )

               !   Defecation as a result of non assimilated products
               !   --------------------------------------------------
               zgrapoc(ji,jj,jk)  = zgraztotc * unass2c + unass2c / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapon(ji,jj,jk)  = zgraztotn * unass2n + no3rat3 * unass2n / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapop(ji,jj,jk)  = zgraztotp * unass2p + po4rat3 * unass2p / ( 1. - 0.8 * epsher2 ) * ztortz
               zgrapof(ji,jj,jk)  = zgraztotf * unass2c + ferat3  * unass2c / ( 1. - 0.8 * epsher2 ) * ztortz

               !  Addition of respiration to the release of inorganic nutrients
               !  -------------------------------------------------------------
               zgrarem(ji,jj,jk) = zgrarem(ji,jj,jk) + zbasresi + zrespirc
               zgraren(ji,jj,jk) = zgraren(ji,jj,jk) + zbasresn + zrespirc * no3rat3
               zgrarep(ji,jj,jk) = zgrarep(ji,jj,jk) + zbasresp + zrespirc * po4rat3
               zgraref(ji,jj,jk) = zgraref(ji,jj,jk) + zbasresf + zrespirc * ferat3

               !   Update the arrays TRA which contain the biological sources and
               !   sinks
               !   --------------------------------------------------------------
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) + zepsherv * zgraztotc - zrespirc   &
               &                     - ztortz - zgrazm
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazdc
               tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) - zgrazdn
               tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) - zgrazdp
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazdf
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgraznc
               tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) - zgraznn
               tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) - zgraznp
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgraznc * trb(ji,jj,jk,jpnch) / ( trb(ji,jj,jk,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazdc * trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazdc * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazdc * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               ! Diazotroph PFT
               tra(ji,jj,jk,jpcdz) = tra(ji,jj,jk,jpcdz) - zgrazdzc
               tra(ji,jj,jk,jpndz) = tra(ji,jj,jk,jpndz) - zgrazdzn
               tra(ji,jj,jk,jppdz) = tra(ji,jj,jk,jppdz) - zgrazdzp
               tra(ji,jj,jk,jpfed) = tra(ji,jj,jk,jpfed) - zgrazdzf
               tra(ji,jj,jk,jpchd) = tra(ji,jj,jk,jpchd) - zgrazdzc * trb(ji,jj,jk,jpchd) / ( trb(ji,jj,jk,jpcdz) + rtrn )
               !
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc - zgrazffep + zfracc
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zfracc
               conspoc(ji,jj,jk)   = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) - zgrazpon - zgrazffnp + zfracn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) - zgrazpop - zgrazffpp + zfracp
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) - zgrazffeg - zfracc
               consgoc(ji,jj,jk)   = consgoc(ji,jj,jk) - zgrazffeg - zfracc
               tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) - zgrazffng - zfracn
               tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) - zgrazffpg - zfracp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof - zgrazfffp + zfracfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) - zgrazfffg - zfracfe
               zfracal = trb(ji,jj,jk,jpcal) / ( trb(ji,jj,jk,jpgoc) + rtrn )
               zgrazcal = zgrazffeg * (1. - part2) * zfracal
               ! Calcite production
               ! Calcite remineralization due to zooplankton activity
               ! part2 of the ingested calcite is dissolving in the acidic gut
               ! -------------------------------------------------------------
               zprcaca = xfracal(ji,jj,jk) * zgraznc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrazcal - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + 2. * ( zgrazcal - zprcaca )
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) - zgrazcal + zprcaca
            END DO
         END DO
      END DO

      ! Computation of the effect of DVM by mesozooplankton
      ! This part is only activated if ln_dvm_meso is set to true
      ! The parameterization has been published in Gorgues et al. (2019).
      ! -----------------------------------------------------------------
      IF (ln_dvm_meso) THEN
         ALLOCATE( zgramigrem(jpi,jpj), zgramigref(jpi,jpj), zgramigpoc(jpi,jpj), zgramigpof(jpi,jpj) )
         ALLOCATE( zgramigrep(jpi,jpj), zgramigren(jpi,jpj), zgramigpop(jpi,jpj), zgramigpon(jpi,jpj) )
         ALLOCATE( zgramigdoc(jpi,jpj), zgramigdon(jpi,jpj), zgramigdop(jpi,jpj) )
         ALLOCATE( zstrn(jpi,jpj) )
         zgramigrem(:,:)  = 0.0   ;   zgramigref(:,:) = 0.0
         zgramigrep(:,:)  = 0.0   ;   zgramigren(:,:) = 0.0
         zgramigpoc(:,:)  = 0.0   ;   zgramigpof(:,:) = 0.0
         zgramigpop(:,:)  = 0.0   ;   zgramigpon(:,:) = 0.0
         zgramigdoc(:,:)  = 0.0   ;   zgramigdon(:,:) = 0.0
         zgramigdop(:,:)  = 0.0   

         ! compute the day length depending on latitude and the day
         zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
         zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )

         ! day length in hours
         zstrn(:,:) = 0.
         DO jj = 1, jpj
            DO ji = 1, jpi
               zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
               zargu = MAX( -1., MIN(  1., zargu ) )
               zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
               zstrn(ji,jj) = MIN(0.75, MAX( 0.25, zstrn(ji,jj) / 24.) )
            END DO
         END DO


        ! Compute the amount of materials that will go into vertical migration
        ! This fraction is sumed over the euphotic zone and is removed from 
        ! the fluxes driven by mesozooplankton in the euphotic zone.
        ! --------------------------------------------------------------------
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zmigreltime = (1. - zstrn(ji,jj))
                  IF ( gdept_n(ji,jj,jk) <= heup(ji,jj) ) THEN
                     zgramigrem(ji,jj) = zgramigrem(ji,jj) + xfracmig * zgrarem(ji,jj,jk) * (1. - zmigreltime )    &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigrep(ji,jj) = zgramigrep(ji,jj) + xfracmig * zgrarep(ji,jj,jk) * (1. - zmigreltime )    &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigrep(ji,jj) = zgramigren(ji,jj) + xfracmig * zgrarep(ji,jj,jk) * (1. - zmigreltime )    &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigref(ji,jj) = zgramigref(ji,jj) + xfracmig * zgraref(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigpoc(ji,jj) = zgramigpoc(ji,jj) + xfracmig * zgrapoc(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigpop(ji,jj) = zgramigpop(ji,jj) + xfracmig * zgrapop(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigpon(ji,jj) = zgramigpon(ji,jj) + xfracmig * zgrapon(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigpof(ji,jj) = zgramigpof(ji,jj) + xfracmig * zgrapof(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigdoc(ji,jj) = zgramigdoc(ji,jj) + xfracmig * zgradoc(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigdop(ji,jj) = zgramigdop(ji,jj) + xfracmig * zgradop(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)
                     zgramigdon(ji,jj) = zgramigdon(ji,jj) + xfracmig * zgradon(ji,jj,jk) * (1. - zmigreltime )   &
                     &                   * e3t_n(ji,jj,jk) * tmask(ji,jj,jk)

                     zgrarem(ji,jj,jk) = zgrarem(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgrarep(ji,jj,jk) = zgrarep(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgraren(ji,jj,jk) = zgraren(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgraref(ji,jj,jk) = zgraref(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgrapoc(ji,jj,jk) = zgrapoc(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgrapop(ji,jj,jk) = zgrapop(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgrapon(ji,jj,jk) = zgrapon(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgrapof(ji,jj,jk) = zgrapof(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgradoc(ji,jj,jk) = zgradoc(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgradop(ji,jj,jk) = zgradop(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                     zgradon(ji,jj,jk) = zgradon(ji,jj,jk) * ( (1.0 - xfracmig) + xfracmig * zmigreltime )
                  ENDIF
               END DO
            END DO
         END DO

         ! The inorganic and organic fluxes induced by migrating organisms are added at the 
         ! the migration depth (corresponding indice is set by kmig)
         ! --------------------------------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF (tmask(ji,jj,1) == 1.) THEN
                  jkt = kmig(ji,jj)
                  zgrarem(ji,jj,jkt) = zgrarem(ji,jj,jkt) + zgramigrem(ji,jj) / e3t_n(ji,jj,jkt)
                  zgrarep(ji,jj,jkt) = zgrarep(ji,jj,jkt) + zgramigrep(ji,jj) / e3t_n(ji,jj,jkt)
                  zgraren(ji,jj,jkt) = zgraren(ji,jj,jkt) + zgramigren(ji,jj) / e3t_n(ji,jj,jkt)
                  zgraref(ji,jj,jkt) = zgraref(ji,jj,jkt) + zgramigref(ji,jj) / e3t_n(ji,jj,jkt)
                  zgrapoc(ji,jj,jkt) = zgrapoc(ji,jj,jkt) + zgramigpoc(ji,jj) / e3t_n(ji,jj,jkt)
                  zgrapop(ji,jj,jkt) = zgrapop(ji,jj,jkt) + zgramigpop(ji,jj) / e3t_n(ji,jj,jkt)
                  zgrapon(ji,jj,jkt) = zgrapon(ji,jj,jkt) + zgramigpon(ji,jj) / e3t_n(ji,jj,jkt)
                  zgrapof(ji,jj,jkt) = zgrapof(ji,jj,jkt) + zgramigpof(ji,jj) / e3t_n(ji,jj,jkt)
                  zgradoc(ji,jj,jkt) = zgradoc(ji,jj,jkt) + zgramigdoc(ji,jj) / e3t_n(ji,jj,jkt)
                  zgradop(ji,jj,jkt) = zgradop(ji,jj,jkt) + zgramigdop(ji,jj) / e3t_n(ji,jj,jkt)
                  zgradon(ji,jj,jkt) = zgradon(ji,jj,jkt) + zgramigdon(ji,jj) / e3t_n(ji,jj,jkt)
               ENDIF
            END DO
         END DO
         !
         ! Deallocate temporary variables
         ! ------------------------------
         DEALLOCATE( zgramigrem, zgramigref, zgramigpoc, zgramigpof )
         DEALLOCATE( zgramigrep, zgramigren, zgramigpop, zgramigpon )
         DEALLOCATE( zgramigdoc, zgramigdon, zgramigdop )
         DEALLOCATE( zstrn )
      ! End of the ln_dvm_meso part
      ENDIF

      !   Update the arrays TRA which contain the biological sources and sinks
      !   This only concerns the variables which are affected by DVM (inorganic 
      !   nutrients, DOC agands, and particulate organic carbon). 
      !   ---------------------------------------------------------------------
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarep(ji,jj,jk) 
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgraren(ji,jj,jk)
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgradoc(ji,jj,jk)
               !
               IF( ln_ligand ) THEN
                  tra(ji,jj,jk,jplgw)  = tra(ji,jj,jk,jplgw) + zgradoc(ji,jj,jk) * ldocz
                  zz2ligprod(ji,jj,jk) = zgradoc(ji,jj,jk) * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zgradon(ji,jj,jk)
               tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + zgradop(ji,jj,jk)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarem(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgraref(ji,jj,jk)
               zfezoo2(ji,jj,jk)   = zgraref(ji,jj,jk)
               znitzoo2(ji,jj,jk)    = zgraren(ji,jj,jk)
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarem(ji,jj,jk)
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgraren(ji,jj,jk)
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zgrapoc(ji,jj,jk)
               prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zgrapoc(ji,jj,jk)
               tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zgrapon(ji,jj,jk)
               tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + zgrapop(ji,jj,jk)
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zgrapof(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(jpi,jpj,jpk) )
         IF( iom_use( "GRAZ2" ) ) THEN
            zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !   Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ2", zw3d )
         ENDIF
         IF( iom_use( "PCAL" ) ) THEN
            zw3d(:,:,:) = prodcal(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)   !  Calcite production
            CALL iom_put( "PCAL", zw3d )
         ENDIF
         IF( iom_use( "FEZOO2" ) ) THEN
            zw3d(:,:,:) = zfezoo2(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
            CALL iom_put( "FEZOO2", zw3d )
         ENDIF
         IF( iom_use( "LPRODZ2" ) .AND. ln_ligand )  THEN
            zw3d(:,:,:) = zz2ligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
            CALL iom_put( "LPRODZ2"  , zw3d )
         ENDIF
         DEALLOCATE( zw3d )
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_meso')
      !
   END SUBROUTINE p6z_meso


   SUBROUTINE p6z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the namp6zmes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp6zmes
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios    ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zmes/part2, bmetexc2, grazrat2, resrat2, mzrat2, xpref2c, xpref2n, xpref2z, &
         &                xpref2m, xpref2d, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2mes, xthresh2, xkgraz2, epsher2, epsher2min, ssigma2, unass2c, &
         &                unass2n, unass2p, srespir2, xsigma2, xsigma2del, grazflux, ln_dvm_meso, xfracmig, &
         &                xpref2dz, xthresh2dz
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp_ref )              ! Namelist namp6zmes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, namp6zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zmes in reference namelist' )
      !
      REWIND( numnatp_cfg )              ! Namelist namp6zmes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, namp6zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zmes in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' ' 
         WRITE(numout,*) ' Namelist parameters for mesozooplankton, namp6zmes'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in mesozoo guts  part2       = ', part2
         WRITE(numout,*) '    mesozoo preference for nano.                   xpref2n     = ', xpref2n
         WRITE(numout,*) '    mesozoo preference for diatoms                 xpref2d     = ', xpref2d
         WRITE(numout,*) '    mesozoo preference for zoo                     xpref2z     = ', xpref2z
         WRITE(numout,*) '    mesozoo preference for mesozoo                 xpref2m     = ', xpref2m
         WRITE(numout,*) '    mesozoo preference for poc                     xpref2c     = ', xpref2c
         WRITE(numout,*) '    microzoo feeding threshold  for mesozoo        xthresh2zoo = ', xthresh2zoo
         WRITE(numout,*) '    diatoms feeding threshold  for mesozoo         xthresh2dia = ', xthresh2dia
         WRITE(numout,*) '    nanophyto feeding threshold for mesozoo        xthresh2phy = ', xthresh2phy
         WRITE(numout,*) '    poc feeding threshold for mesozoo              xthresh2poc = ', xthresh2poc
         WRITE(numout,*) '    mesozoo feeding threshold for mesozoo          xthresh2mes = ', xthresh2mes
         WRITE(numout,*) '    feeding threshold for mesozooplankton          xthresh2    = ', xthresh2
         WRITE(numout,*) '    exsudation rate of mesozooplankton             resrat2     = ', resrat2
         WRITE(numout,*) '    mesozooplankton mortality rate                 mzrat2      = ', mzrat2
         WRITE(numout,*) '    maximal mesozoo grazing rate                   grazrat2    = ', grazrat2
         WRITE(numout,*) '    mesozoo flux feeding rate                      grazflux    = ', grazflux
         WRITE(numout,*) '    C egested fraction of food by mesozoo          unass2c     = ', unass2c
         WRITE(numout,*) '    N egested fraction of food by mesozoo          unass2n     = ', unass2n
         WRITE(numout,*) '    P egested fraction of food by mesozoo          unass2p     = ', unass2p
         WRITE(numout,*) '    Efficicency of Mesozoo growth                  epsher2     = ', epsher2
         WRITE(numout,*) '    Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '    Fraction excreted as semi-labile DOM           ssigma2     = ', ssigma2
         WRITE(numout,*) '    Active respiration                             srespir2    = ', srespir2
         WRITE(numout,*) '    half sturation constant for grazing 2          xkgraz2     = ', xkgraz2
         WRITE(numout,*) '    Use excess carbon for respiration              bmetexc2    = ', bmetexc2
         WRITE(numout,*) '      Width of the grazing window                     xsigma2     =', xsigma2
         WRITE(numout,*) '      Maximum additional width of the grazing window  xsigma2del  =', xsigma2del
         WRITE(numout,*) '      Diurnal vertical migration of mesozoo.         ln_dvm_meso  =', ln_dvm_meso
         WRITE(numout,*) '      Fractional biomass of meso  that performs DVM  xfracmig     =', xfracmig
         WRITE(numout,*) '    mesozoo preference for diazotrophs                xpref2dz    =', xpref2dz
         WRITE(numout,*) '    diazotroph feeding threshold for mesozoo          xthresh2dz  =', xthresh2dz
      ENDIF
      !
   END SUBROUTINE p6z_meso_init

   SUBROUTINE p6z_meso_depmig
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_meso_depmig  ***
      !!
      !! ** Purpose :   Computation the migration depth of mesozooplankton
      !!
      !! ** Method  :   Computes the DVM depth of mesozooplankton from oxygen
      !!      temperature and chlorophylle following the parameterization 
      !!      proposed by Bianchi et al. (2013)
      !!
      !! ** input   :   
      !!----------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      !
      REAL(wp) :: totchl
      REAL(wp), DIMENSION(jpi,jpj) :: oxymoy, tempmoy, zdepmoy

      !!---------------------------------------------------------------------
      !
      IF( ln_timing == 1 )  CALL timing_start('p6z_meso_zdepmig')
      !
      oxymoy(:,:)  = 0.
      tempmoy(:,:) = 0.
      zdepmoy(:,:) = 0.
      depmig (:,:) = 5.
      kmig   (:,:) = 1
      !
      ! Compute the averaged values of oxygen, temperature over the domain 
      ! 150m to 500 m depth.
      !
      DO jk =1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF (tmask(ji,jj,jk) == 1.) THEN
                  IF (gdept_n(ji,jj,jk) >= 150. .AND. gdept_n(ji,jj,jk) <= 500.) THEN
                     oxymoy(ji,jj)  = oxymoy(ji,jj)  + trb(ji,jj,jk,jpoxy)*e3t_n(ji,jj,jk)*1E6
                     tempmoy(ji,jj) = tempmoy(ji,jj) + tsn(ji,jj,jk,jp_tem)*e3t_n(ji,jj,jk)
                     zdepmoy(ji,jj) = zdepmoy(ji,jj) + e3t_n(ji,jj,jk)
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

      ! Compute the difference between surface values and the mean values in the mesopelagic
      ! domain
      ! ------------------------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            oxymoy(ji,jj) = trb(ji,jj,1,jpoxy)*1E6 - oxymoy(ji,jj) / (zdepmoy(ji,jj) + rtrn)
            tempmoy(ji,jj) = tsn(ji,jj,1,jp_tem)-tempmoy(ji,jj) / (zdepmoy(ji,jj) + rtrn)
         END DO
      END DO

      ! Computation of the migration depth based on the parameterization of 
      ! Bianchi et al. (2013)
      ! -------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF (tmask(ji,jj,1) == 1.) THEN
               totchl = (trb(ji,jj,1,jppch)+trb(ji,jj,1,jpnch)+trb(ji,jj,1,jpdch))*1E6
               depmig(ji,jj) = 398. - 0.56 * oxymoy(ji,jj) -115. * log10(totchl) + 0.36 * hmld(ji,jj) -2.4 * tempmoy(ji,jj)
            ENDIF
         END DO
      END DO
      ! 
      ! Computation of the corresponding jk indice 
      ! ------------------------------------------
      DO jk = 1, jpk-1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF (depmig(ji,jj) .GE. gdepw_n(ji,jj,jk) .AND. depmig(ji,jj) .LT. gdepw_n(ji,jj,jk+1) ) THEN
                  kmig(ji,jj) = jk
               ENDIF
            END DO
         END DO
      END DO
      !
      ! Correction of the migration depth and indice based on O2 levels
      ! If O2 is too low, imposing a migration depth at this low O2 levels
      ! would lead to negative O2 concentrations (respiration while O2 is close
      ! to 0. Thus, to avoid that problem, the migration depth is adjusted so
      ! that it falls above the OMZ
      ! -----------------------------------------------------------------------
      DO ji =1, jpi
         DO jj = 1, jpj
            IF (trb(ji,jj,kmig(ji,jj),jpoxy) < 5E-6) THEN
               DO jk = kmig(ji,jj),1,-1
                  IF (trb(ji,jj,jk,jpoxy) >= 5E-6 .AND. trb(ji,jj,jk+1,jpoxy)  < 5E-6) THEN
                     kmig(ji,jj) = jk
                     depmig(ji,jj) = gdept_n(ji,jj,jk)
                  ENDIF
               END DO
            ENDIF
         END DO
      END DO
      !
      IF( ln_timing )   CALL timing_stop('p6z_meso_depmig')
      !
   END SUBROUTINE p6z_meso_depmig

   INTEGER FUNCTION p6z_meso_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_meso_alloc  ***
      !!----------------------------------------------------------------------
      !
      ALLOCATE( depmig(jpi,jpj), kmig(jpi,jpj), STAT= p6z_meso_alloc  )
      !
      IF( p6z_meso_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p6z_meso_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p6z_meso_alloc

   !!======================================================================
END MODULE p6zmeso
