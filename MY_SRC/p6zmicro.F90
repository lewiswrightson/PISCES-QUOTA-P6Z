MODULE p6zmicro
   !!======================================================================
   !!                         ***  MODULE p6zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!         Including Explicit Diazotroph PFT
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_micro       :   Compute the sources/sinks for microzooplankton
   !!   p6z_micro_init  :   Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  PISCES nutrient limitation term of PISCES std
   USE p6zlim          !  Phytoplankton limitation terms
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p6z_micro         ! called in p6zbio.F90
   PUBLIC   p6z_micro_init    ! called in trcsms_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::  xprefc     !: microzoo preference for POC 
   REAL(wp), PUBLIC ::  xprefn     !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::  xprefp     !: microzoo preference for picophyto
   REAL(wp), PUBLIC ::  xprefd     !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::  xprefz     !: microzoo preference for microzoo
   REAL(wp), PUBLIC ::  xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshpic  !: picophyto feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshzoo  !: microzoo threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::  xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::  resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::  mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::  unassc      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  unassn      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  unassp      !: Non-assimilated part of food
   REAL(wp), PUBLIC ::  epsher      !: Growth efficiency for microzoo
   REAL(wp), PUBLIC ::  epshermin   !: Minimum growth efficiency for microzoo
   REAL(wp), PUBLIC ::  srespir     !: half sturation constant for grazing 1 
   REAL(wp), PUBLIC ::  ssigma      !: Fraction excreted as semi-labile DOM
   REAL(wp), PUBLIC ::  xsigma      !: Width of the grazing window
   REAL(wp), PUBLIC ::  xsigmadel   !: Maximum additional width of the grazing window at low food density
   REAL(wp), PUBLIC ::  xprefdz     !: Microzoo preference for diazotrophs
   REAL(wp), PUBLIC ::  xthreshdz   !: Diazotroph feeding threshold for microzooplankton
   LOGICAL,  PUBLIC ::  bmetexc     !: Use of excess carbon for respiration

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zmicro.F90 13233 2020-07-02 18:34:16Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::  kt  ! ocean time step
      INTEGER, INTENT(in) ::  knt 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc, zcompapon, zcompapop
      REAL(wp) :: zcompapi, zgraze  , zdenom, zfact, zfood, zfoodlim
      REAL(wp) :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5, ztmptot
      REAL(wp) :: zepsherf, zepshert, zepsherq, zepsherv, zrespirc, zrespirn, zrespirp, zbasresb, zbasresi
      REAL(wp) :: zgraztotc, zgraztotn, zgraztotp, zgraztotf, zbasresn, zbasresp, zbasresf
      REAL(wp) :: zgradoc, zgradon, zgradop, zgraref, zgradoct, zgradont, zgradopt, zgrareft
      REAL(wp) :: zexcess, zgraren, zgrarep, zgrarem
      REAL(wp) :: zgrapoc, zgrapon, zgrapop, zgrapof, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasratf, zgrasratn, zgrasratp
      REAL(wp) :: zgraznc, zgraznn, zgraznp, zgrazpoc, zgrazpon, zgrazpop, zgrazpof
      REAL(wp) :: zgrazdc, zgrazdn, zgrazdp, zgrazdf, zgraznf, zgrazz
      REAL(wp) :: zgrazpc, zgrazpn, zgrazpp, zgrazpf, zbeta, zrfact2, zmetexcess
      REAL(wp) :: zsigma, zdiffdn, zdiffpn, zdiffdp, zproport, zproport2
      ! Explicit Diazotroph PFT
      REAL(wp) :: zcompadz, ztmp6, zgrazdzc, zgrazdzn, zgrazdzp, zgrazdzf, zdiffdzn, zdiffdzp, zdiffdzd, zproport3
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing, zfezoo, znitzoo
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d, zzligprod
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_micro')
      !
      IF (ln_ligand) THEN
         ALLOCATE( zzligprod(jpi,jpj,jpk) )
         zzligprod(:,:,:) = 0._wp
      ENDIF
      !
      ! Use of excess carbon for metabolism
      zmetexcess = 0.0
      IF ( bmetexc ) zmetexcess = 1.0
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaz = MAX( ( trb(ji,jj,jk,jpzoo) - 1.e-9 ), 0.e0 )
               zfact   = xstep * tgfunc2(ji,jj,jk) * zcompaz
               ! Proportion of nano and diatoms that are within the size range
               ! accessible to microzooplankton. 
               zproport  = min(1.0, exp(-1.1 * MAX(0., ( sized(ji,jj,jk) - 1.8 ))**0.8 ))
               zproport2 = sizen(ji,jj,jk)**(-0.54)
               zproport3 = sizedz(ji,jj,jk)**(-0.54)
               !!! TO fix works for now
            !  zproport2  = min(1.0, exp(-1.1 * MAX(0., ( sizen(ji,jj,jk) - 1.8 ))**0.8 ))
             ! zproport3  = min(1.0, exp(-1.1 * MAX(0., ( sizedz(ji,jj,jk) - 1.8 ))**0.8 ))

               !  linear mortality of mesozooplankton
               !  A michaelis menten modulation term is used to avoid extinction of 
               !  microzooplankton at very low food concentrations. Mortality is 
               !  enhanced in low O2 waters
               !  -----------------------------------------------------------------
               zrespz = resrat * zfact * ( trb(ji,jj,jk,jpzoo) / ( xkmort + trb(ji,jj,jk,jpzoo) )  &
               &        + 3. * nitrfac(ji,jj,jk) )

               !  Zooplankton quadratic mortality. A square function has been selected with
               !  to mimic predation and disease (density dependent mortality). It also tends
               !  to stabilise the model
               !  -------------------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trb(ji,jj,jk,jpzoo) * (1. - nitrfac(ji,jj,jk))

               !   Computation of the abundance of the preys
               !   A threshold can be specified in the namelist
               !   Nanophyto and diatoms have a specific treatment with 
               !   teir preference decreasing with size.
               !   --------------------------------------------------------
               zcompadi  = zproport * MAX( ( trb(ji,jj,jk,jpdia) - xthreshdia ), 0.e0 )
               zcompaph  = zproport2 * MAX( ( trb(ji,jj,jk,jpphy) - xthreshphy ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,jk,jpzoo) - xthreshzoo ), 0.e0 )
               zcompapi  = MAX( ( trb(ji,jj,jk,jppic) - xthreshpic ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,jk,jppoc) - xthreshpoc ), 0.e0 )
               zcompadz  = zproport3 * MAX( ( trb(ji,jj,jk,jpcdz) - xthreshdz ), 0.e0 )               
               ! Microzooplankton grazing
               ! The total amount of food is the sum of all preys accessible to mesozooplankton 
               ! multiplied by their food preference
               ! A threshold can be specified in the namelist (xthresh). However, when food 
               ! concentration is close to this threshold, it is decreased to avoid the 
               ! accumulation of food in the mesozoopelagic domain
               ! -------------------------------------------------------------------------------
               zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi   &
               &           + xprefz * zcompaz + xprefp * zcompapi + xprefdz * zcompadz
               zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
               zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               zgraze    = grazrat * xstep * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpzoo) * (1. - nitrfac(ji,jj,jk)) 

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
               zsigma = xsigma + xsigmadel * zsigma
               zdiffpn = exp( -ABS(log(0.5 * sizep(ji,jj,jk) / (3.0 * sizen(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               zdiffdn = exp( -ABS(log(3.0 * sizen(ji,jj,jk) / (5.0 * sized(ji,jj,jk) + rtrn )) )**2 / zsigma**2)
               zdiffdp = exp( -ABS(log(0.5 * sizep(ji,jj,jk) / (5.0 * sized(ji,jj,jk) + rtrn )) )**2 / zsigma**2)
               zdiffdzp = exp( -ABS(log(0.5 * sizep(ji,jj,jk) / (3.0 * sizedz(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               zdiffdzd = exp( -ABS(log(3.0 * sizedz(ji,jj,jk) / (5.0 * sized(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )
               zdiffdzn = exp( -ABS(log(3.0 * sizen(ji,jj,jk) / (3.0 * sizedz(ji,jj,jk) + rtrn )) )**2 / zsigma**2 )               
               ztmp1 = xprefn * zcompaph * ( zcompaph + zdiffdn * zcompadi + zdiffpn * zcompapi ) / ( 1.0 + zdiffdn + zdiffpn )
               ztmp2 = xprefp * zcompapi * ( zcompapi + zdiffpn * zcompaph + zdiffdp * zcompadi ) / ( 1.0 + zdiffpn + zdiffdp )
               ztmp3 = xprefc * zcompapoc**2
               ztmp4 = xprefd * zcompadi * ( zdiffdp * zcompapi + zdiffdn * zcompaph + zcompadi ) / ( 1.0 + zdiffdn + zdiffdp )
               ztmp5 = xprefz * zcompaz**2
               ztmp6 = xprefdz * zcompadz * ( zcompadz + zdiffdzn * zcompaph + zdiffdzp * zcompapi + zdiffdzd &
               &       * zcompadi ) / ( 1.0 + zdiffdzn + zdiffdzp + zdiffdzd )
               ztmptot = ztmp1 + ztmp2 + ztmp3 + ztmp4 + ztmp5 + ztmp6 + rtrn
               ztmp1 = ztmp1 / ztmptot
               ztmp2 = ztmp2 / ztmptot
               ztmp3 = ztmp3 / ztmptot
               ztmp4 = ztmp4 / ztmptot
               ztmp5 = ztmp5 / ztmptot
               ztmp6 = ztmp6 / ztmptot

               !   Microzooplankton regular grazing on the different preys
               !   -------------------------------------------------------
               !   Nanophytoplankton
               zgraznc   = zgraze  * ztmp1  * zdenom  
               zgraznn   = zgraznc * trb(ji,jj,jk,jpnph) / (trb(ji,jj,jk,jpphy) + rtrn)
               zgraznp   = zgraznc * trb(ji,jj,jk,jppph) / (trb(ji,jj,jk,jpphy) + rtrn)
               zgraznf   = zgraznc * trb(ji,jj,jk,jpnfe) / (trb(ji,jj,jk,jpphy) + rtrn)

               ! Picophytoplankton
               zgrazpc   = zgraze  * ztmp2  * zdenom
               zgrazpn   = zgrazpc * trb(ji,jj,jk,jpnpi) / (trb(ji,jj,jk,jppic) + rtrn)
               zgrazpp   = zgrazpc * trb(ji,jj,jk,jpppi) / (trb(ji,jj,jk,jppic) + rtrn)
               zgrazpf   = zgrazpc * trb(ji,jj,jk,jppfe) / (trb(ji,jj,jk,jppic) + rtrn)
               ! Microzooplankton
               zgrazz    = zgraze  * ztmp5   * zdenom

               ! small POC
               zgrazpoc  = zgraze  * ztmp3   * zdenom
               zgrazpon  = zgrazpoc * trb(ji,jj,jk,jppon) / ( trb(ji,jj,jk,jppoc) + rtrn )
               zgrazpop  = zgrazpoc * trb(ji,jj,jk,jppop) / ( trb(ji,jj,jk,jppoc) + rtrn )
               zgrazpof  = zgrazpoc* trb(ji,jj,jk,jpsfe) / (trb(ji,jj,jk,jppoc) + rtrn)

               ! Diatoms
               zgrazdc   = zgraze  * ztmp4  * zdenom
               zgrazdn   = zgrazdc * trb(ji,jj,jk,jpndi) / (trb(ji,jj,jk,jpdia) + rtrn)
               zgrazdp   = zgrazdc * trb(ji,jj,jk,jppdi) / (trb(ji,jj,jk,jpdia) + rtrn)
               zgrazdf   = zgrazdc * trb(ji,jj,jk,jpdfe) / (trb(ji,jj,jk,jpdia) + rtrn)
               !
               ! Diazotrophs
               zgrazdzc  = zgraze  * ztmp6  * zdenom
               zgrazdzn  = zgrazdzc * trb(ji,jj,jk,jpndz) / (trb(ji,jj,jk,jpcdz) + rtrn)
               zgrazdzp  = zgrazdzc * trb(ji,jj,jk,jppdz) / (trb(ji,jj,jk,jpcdz) + rtrn)
               zgrazdzf  = zgrazdzc * trb(ji,jj,jk,jpfed) / (trb(ji,jj,jk,jpcdz) + rtrn)  
               ! Total ingestion rates in C, P, Fe, N
               zgraztotc = zgraznc + zgrazpoc + zgrazdc + zgrazz + zgrazpc + zgrazdzc
               zgraztotn = zgraznn + zgrazpn + zgrazpon + zgrazdn + zgrazdzn + zgrazz * no3rat3
               zgraztotp = zgraznp + zgrazpp + zgrazpop + zgrazdp + zgrazdzp + zgrazz * po4rat3
               zgraztotf = zgraznf + zgrazpf + zgrazpof + zgrazdf + zgrazdzf + zgrazz * ferat3
               !
               ! Grazing by microzooplankton
               zgrazing(ji,jj,jk) = zgraztotc
               ! Stoichiometruc ratios of the food ingested by zooplanton 
               ! --------------------------------------------------------
               zgrasratf =  (zgraztotf + rtrn) / ( zgraztotc + rtrn )
               zgrasratn =  (zgraztotn + rtrn) / ( zgraztotc + rtrn )
               zgrasratp =  (zgraztotp + rtrn) / ( zgraztotc + rtrn )

               ! Mesozooplankton efficiency. 
               ! We adopt a formulation proposed by Mitra et al. (2007)
               ! The gross growth efficiency is controled by the most limiting nutrient.
               ! Growth is also further decreased when the food quality is poor. This is currently
               ! hard coded : it can be decreased by up to 50% (zepsherq)
               ! GGE can also be decreased when food quantity is high, zepsherf (Montagnes and 
               ! Fulton, 2012)
               ! -----------------------------------------------------------------------------------
               zepshert  = MIN( 1., zgrasratn/ no3rat3, zgrasratp/ po4rat3, zgrasratf / ferat3)
               zbeta     = MAX( 0., (epsher - epshermin) )
               ! Food density deprivation of GGE
               zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               ! Food quality deprivation of GGE
               zepsherq  = 0.5 + (1.0 - 0.5) * zepshert * ( 1.0 + 1.0 ) / ( zepshert + 1.0 )
               ! Actual GGE
               zepsherv  = zepsherf * zepshert * zepsherq

               ! Respiration of microzooplankton
               ! Excess carbon in the food is used preferentially
               ! when activated by zmetexcess
               ! ------------------------------------------------
               zexcess  = zgraztotc * zepsherf * (1.0 - zepshert) * zmetexcess
               zbasresb = MAX(0., zrespz - zexcess)
               zbasresi = zexcess + MIN(0., zrespz - zexcess)  
               zrespirc = srespir * zepsherv * zgraztotc + zbasresb
               
               ! When excess carbon is used, the other elements in excess
               ! are also used proportionally to their abundance
               ! --------------------------------------------------------
               zexcess  = ( zgrasratn/ no3rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresn = zbasresi * zexcess * zgrasratn 
               zexcess  = ( zgrasratp/ po4rat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresp = zbasresi * zexcess * zgrasratp
               zexcess  = ( zgrasratf/ ferat3 - zepshert ) / ( 1.0 - zepshert + rtrn)
               zbasresf = zbasresi * zexcess * zgrasratf

               ! Voiding of the excessive elements as DOM
               ! ----------------------------------------
               zgradoct   = (1. - unassc - zepsherv) * zgraztotc - zbasresi  
               zgradont   = (1. - unassn) * zgraztotn - zepsherv * no3rat3 * zgraztotc - zbasresn
               zgradopt   = (1. - unassp) * zgraztotp - zepsherv * po4rat3 * zgraztotc - zbasresp
               zgrareft   = (1. - unassc) * zgraztotf - zepsherv * ferat3 * zgraztotc - zbasresf

               ! Since only semilabile DOM is represented in PISCES
               ! part of DOM is in fact labile and is then released
               ! as dissolved inorganic compounds (ssigma)
               ! --------------------------------------------------
               zgradoc =  zgradoct * ssigma
               zgradon =  zgradont * ssigma
               zgradop =  zgradopt * ssigma
               zgrarem = (1.0 - ssigma) * zgradoct
               zgraren = (1.0 - ssigma) * zgradont
               zgrarep = (1.0 - ssigma) * zgradopt
               zgraref = zgrareft

               ! Defecation as a result of non assimilated products
               ! --------------------------------------------------
               zgrapoc   = zgraztotc * unassc
               zgrapon   = zgraztotn * unassn
               zgrapop   = zgraztotp * unassp
               zgrapof   = zgraztotf * unassc

               ! Addition of respiration to the release of inorganic nutrients
               ! -------------------------------------------------------------
               zgrarem = zgrarem + zbasresi + zrespirc
               zgraren = zgraren + zbasresn + zrespirc * no3rat3
               zgrarep = zgrarep + zbasresp + zrespirc * po4rat3
               zgraref = zgraref + zbasresf + zrespirc * ferat3

               ! Update of the TRA arrays
               ! ------------------------
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarep
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgraren
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgradoc
               !
               IF( ln_ligand ) THEN 
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zgradoc * ldocz
                  zzligprod(ji,jj,jk) = zgradoc * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zgradon
               tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + zgradop
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarem 
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgraref
               zfezoo(ji,jj,jk)    = zgraref
               znitzoo(ji,jj,jk)     = zgraren
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) + zepsherv * zgraztotc - zrespirc - ztortz - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgraznc
               tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) - zgraznn
               tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) - zgraznp
               tra(ji,jj,jk,jppic) = tra(ji,jj,jk,jppic) - zgrazpc
               tra(ji,jj,jk,jpnpi) = tra(ji,jj,jk,jpnpi) - zgrazpn
               tra(ji,jj,jk,jpppi) = tra(ji,jj,jk,jpppi) - zgrazpp
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazdc
               tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) - zgrazdn
               tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) - zgrazdp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgraznc * trb(ji,jj,jk,jpnch)/(trb(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jppch) = tra(ji,jj,jk,jppch) - zgrazpc * trb(ji,jj,jk,jppch)/(trb(ji,jj,jk,jppic)+rtrn)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazdc * trb(ji,jj,jk,jpdch)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazdc * trb(ji,jj,jk,jpdsi)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazdc * trb(ji,jj,jk,jpdsi)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jppfe) = tra(ji,jj,jk,jppfe) - zgrazpf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazdf
               ! Diazotroph PFT
               tra(ji,jj,jk,jpcdz) = tra(ji,jj,jk,jpcdz) - zgrazdzc
               tra(ji,jj,jk,jpndz) = tra(ji,jj,jk,jpndz) - zgrazdzn
               tra(ji,jj,jk,jppdz) = tra(ji,jj,jk,jppdz) - zgrazdzp
               tra(ji,jj,jk,jpfed) = tra(ji,jj,jk,jpfed) - zgrazdzf
               tra(ji,jj,jk,jpchd) = tra(ji,jj,jk,jpchd) - zgrazdzc * trb(ji,jj,jk,jpchd)/(trb(ji,jj,jk,jpcdz)+rtrn)
               !
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + ztortz + zgrapoc - zgrazpoc 
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + ztortz + zgrapoc
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + no3rat3 * ztortz + zgrapon - zgrazpon
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + po4rat3 * ztortz + zgrapop - zgrazpop
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * ztortz  + zgrapof - zgrazpof
               !
               ! Calcite production
               ! Calcite remineralization due to zooplankton activity
               ! part of the ingested calcite is dissolving in the acidic gut
               ! -------------------------------------------------------------
               zprcaca = xfracal(ji,jj,jk) * zgraznc
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarem - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca     &
               &                     + rno3 * zgraren
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
            END DO
         END DO
      END DO
      !
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            ALLOCATE( zw3d(jpi,jpj,jpk) )
            IF( iom_use( "GRAZ1" ) ) THEN
               zw3d(:,:,:) = zgrazing(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !  Total grazing of phyto by zooplankton
               CALL iom_put( "GRAZ1", zw3d )
            ENDIF
            IF( iom_use( "FEZOO" ) ) THEN
               zw3d(:,:,:) = zfezoo(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
               CALL iom_put( "FEZOO", zw3d )
            ENDIF
            IF( iom_use( "LPRODZ" ) .AND. ln_ligand )  THEN
               zw3d(:,:,:) = zzligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
               CALL iom_put( "LPRODZ"  , zw3d )
            ENDIF
            DEALLOCATE( zw3d )
         ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_micro')
      !
   END SUBROUTINE p6z_micro


   SUBROUTINE p6z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the namp6zzoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp6zzoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namp6zzoo/ part, grazrat, bmetexc, resrat, mzrat, xprefc, xprefn, &
         &                xprefp, xprefd, xprefz, xthreshdia, xthreshphy, &
         &                xthreshpic, xthreshpoc, xthreshzoo, xthresh, xkgraz, &
         &                epsher, epshermin, ssigma, srespir, unassc, unassn, unassp,   &
         &                xsigma, xsigmadel, xprefdz, xthreshdz   
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp_ref )              ! Namelist namp6zzoo in reference namelist : Pisces microzooplankton
      READ  ( numnatp_ref, namp6zzoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zzoo in reference namelist' )
      !
      REWIND( numnatp_cfg )              ! Namelist namp6zzoo in configuration namelist : Pisces microzooplankton
      READ  ( numnatp_cfg, namp6zzoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zzoo in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zzoo )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for microzooplankton, nampiszooq'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '    microzoo preference for POC                     xprefc     =', xprefc
         WRITE(numout,*) '    microzoo preference for nano                    xprefn     =', xprefn
         WRITE(numout,*) '    microzoo preference for pico                    xprefp     =', xprefp
         WRITE(numout,*) '    microzoo preference for diatoms                 xprefd     =', xprefd
         WRITE(numout,*) '    microzoo preference for microzoo                xprefz     =', xprefz
         WRITE(numout,*) '    diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(numout,*) '    nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '    picophyto feeding threshold for microzoo        xthreshpic  =', xthreshpic
         WRITE(numout,*) '    poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '    microzoo feeding threshold for microzoo         xthreshzoo  =', xthreshzoo
         WRITE(numout,*) '    feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '    exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '    microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '    maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '    C egested fraction of fodd by microzoo          unassc      =', unassc
         WRITE(numout,*) '    N egested fraction of fodd by microzoo          unassn      =', unassn
         WRITE(numout,*) '    P egested fraction of fodd by microzoo          unassp      =', unassp
         WRITE(numout,*) '    Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '    Minimum Efficiency of Microzoo growth           epshermin   =', epshermin
         WRITE(numout,*) '    Fraction excreted as semi-labile DOM            ssigma      =', ssigma
         WRITE(numout,*) '    Active respiration                              srespir     =', srespir
         WRITE(numout,*) '    half sturation constant for grazing 1           xkgraz      =', xkgraz
         WRITE(numout,*) '    Use of excess carbon for respiration            bmetexc     =', bmetexc
         WRITE(numout,*) '      Width of the grazing window                     xsigma      =', xsigma
         WRITE(numout,*) '      Maximum additional width of the grazing window  xsigmadel   =', xsigmadel
         WRITE(numout,*) '    microzoo preference for diazotrophs               xprefdz     =', xprefdz
         WRITE(numout,*) '    diazotroph feeding threshold for microzoo         xthreshdz   =', xthreshdz
      ENDIF
      !
   END SUBROUTINE p6z_micro_init

   !!======================================================================
END MODULE p6zmicro
