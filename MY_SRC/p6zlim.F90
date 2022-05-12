MODULE p6zlim
   !!======================================================================
   !!                         ***  MODULE p6zlim  ***
   !! TOP :   PISCES-QUOTA : Computes the various nutrient limitation terms
   !!                        of phytoplankton, with explicit diazotrophy
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-04  (O. Aumont, C. Ethe) Limitation for iron modelled in quota 
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_lim        :   Compute the nutrients limitation terms 
   !!   p6z_lim_init   :   Read the namelist 
   !!----------------------------------------------------------------------
   USE oce_trc         ! Shared ocean-passive tracers variables
   USE trc             ! Tracers defined
   USE p4zlim          ! Nutrient limitation 
   USE sms_pisces      ! PISCES variables
   USE iom             !  I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC p6z_lim           ! called in p4zbio.F90  
   PUBLIC p6z_lim_init      ! called in trcsms_pisces.F90 
   PUBLIC p6z_lim_alloc     ! called in trcini_pisces.F90

   !! * Shared module variables
   REAL(wp), PUBLIC ::  concpno3    !:  NO3 half saturation for picophyto  
   REAL(wp), PUBLIC ::  concpnh4    !:  NH4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concnpo4    !:  PO4 half saturation for nanophyto
   REAL(wp), PUBLIC ::  concppo4    !:  PO4 half saturation for picophyto
   REAL(wp), PUBLIC ::  concdpo4    !:  PO4 half saturation for diatoms
   REAL(wp), PUBLIC ::  concpfer    !:  Iron half saturation for picophyto
   REAL(wp), PUBLIC ::  concbpo4    !:  PO4 half saturation for bacteria
   REAL(wp), PUBLIC ::  xsizepic    !:  Minimum size criteria for picophyto
   REAL(wp), PUBLIC ::  xsizerp     !:  Size ratio for picophytoplankton
   REAL(wp), PUBLIC ::  qfnopt      !:  optimal Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpopt      !:  optimal Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdopt      !:  optimal Fe quota for diatoms
   REAL(wp), PUBLIC ::  qnnmin      !:  minimum N  quota for nanophyto
   REAL(wp), PUBLIC ::  qnnmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qpnmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmin      !:  minimum N quota for nanophyto
   REAL(wp), PUBLIC ::  qnpmax      !:  maximum N quota for nanophyto
   REAL(wp), PUBLIC ::  qppmin      !:  minimum P quota for nanophyto
   REAL(wp), PUBLIC ::  qppmax      !:  maximum P quota for nanophyto
   REAL(wp), PUBLIC ::  qndmin      !:  minimum N quota for diatoms
   REAL(wp), PUBLIC ::  qndmax      !:  maximum N quota for diatoms
   REAL(wp), PUBLIC ::  qpdmin      !:  minimum P quota for diatoms
   REAL(wp), PUBLIC ::  qpdmax      !:  maximum P quota for diatoms
   REAL(wp), PUBLIC ::  qfnmax      !:  maximum Fe quota for nanophyto
   REAL(wp), PUBLIC ::  qfpmax      !:  maximum Fe quota for picophyto
   REAL(wp), PUBLIC ::  qfdmax      !:  maximum Fe quota for diatoms
   REAL(wp), PUBLIC ::  zpsinh4     !:  respiration cost of NH4 assimilation
   REAL(wp), PUBLIC ::  zpsino3     !:  respiration cost of NO3 assimilation
   REAL(wp), PUBLIC ::  zpsiuptk    !:  Mean respiration cost
   REAL(wp), PUBLIC ::  concdzno3   !:  Nitrate half saturation for diazo
   REAL(wp), PUBLIC ::  concdznh4   !:  NH4 half saturation for diazo
   REAL(wp), PUBLIC ::  concdzpo4   !:  PO4 half saturation for diazo
   REAL(wp), PUBLIC ::  concdzfer   !:  Fe half saturation for diazo
   REAL(wp), PUBLIC ::  xsizedz     !:  Minimum size criteria for phyto
   REAL(wp), PUBLIC ::  xsizerdz    !:  size ratio of diazotrophs
   REAL(wp), PUBLIC ::  qfdzopt     !:  optimal Fe quota of diazotrophs 
   REAL(wp), PUBLIC ::  qndzmin     !:  Minimal N quota of diazotrophs
   REAL(wp), PUBLIC ::  qndzmax     !:  Maximal N quota of diazotrophs
   REAL(wp), PUBLIC ::  qpdzmin     !:  Minimal P quota of diazotrophs
   REAL(wp), PUBLIC ::  qpdzmax     !:  Maximal P quota of diazotrophs
   REAL(wp), PUBLIC ::  qfdzmax     !:  Maximal Fe quota of diazotrophs
   REAL(wp), PUBLIC ::  zpsinfix    !:  Cost of biosynthesis associated with Nfix
   REAL(wp), PUBLIC ::  xkdop       !:  half saturation of DOP uptake phytos
   REAL(wp), PUBLIC ::  xkdopdz     !:  half saturation of DOP uptake diazos
   REAL(wp), PUBLIC ::  Facul_lim   !:  Diazos Facultative sensitivity
   REAL(wp), PUBLIC ::  kustkaFe    !: Diazo Fe limitation based on kustka 2003
   REAL(wp), PUBLIC ::  maxFescale  !: MAximum Fe cost of nfix can be scaled
   REAL(wp), PUBLIC ::  maxPminscale!: Maximum Pmin scaling 
   LOGICAL , PUBLIC ::  ln_tiue     !: Boolean to activate temperature dependence on Fe cost of nfix (IUE)
   LOGICAL , PUBLIC ::  ln_tpue     !: Boolean to activate temperature dependence on QPmin for diazo (PUE)
   !!*  Allometric variations of the quotas
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmin    !: Minimum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnnmax    !: Maximum N quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmin    !: Minimum P quota of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpnmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmin    !: Minimum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqnpmax    !: Maximum N quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmin    !: Minimum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqppmax    !: Maximum P quota of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmin    !: Minimum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndmax    !: Maximum N quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmin    !: Minimum P quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdmax    !: Maximum P quota of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndzmin   !: Minimum N quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqndzmax   !: Maximum N quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdzmin   !: Minimum P quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   xqpdzmax   !: Maximum P quota of diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   qfeminn    !: QFe min of nanophytoplankton 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   qfeminp    !: QFe min of picophytoplankton
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   qfemind    !: QFe min of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE,   DIMENSION(:,:,:)  ::   qfemindz   !: QFe min of diazotrophs
   !!* Phytoplankton nutrient limitation terms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicono3   !: Limitation of NO3 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpiconh4   !: Limitation of NH4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicopo4   !: Limitation of PO4 uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xnanodop   !: Limitation of DOP uptake by nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicodop   !: Limitation of DOP uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiatdop   !: Limitation of DOP uptake by diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xpicofer   !: Limitation of Fe uptake by picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpic    !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpics   !: Limitation of picophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimphys   !: Limitation of nanophyto PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdias   !: Limitation of diatoms PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimpfe    !: Limitation of picophyto PP by Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvnuptk    !: Maximum potential uptake rate of nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvpuptk    !: Maximum potential uptake rate of picophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvduptk    !: Maximum potential uptake rate of diatoms
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazno3   !: Limitation of NO3 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiaznh4   !: Limitation of NH4 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazpo4   !: Limitation of PO4 uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazdop   !: Limitation of DOP uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xdiazfer   !: Limitation of Fe uptake by diazos
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdiaz   !: Limitation of diazo PP by C and N (P by proxy)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdzfe   !: Limitation of diazo PP by Fe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   fvdzuptk   !: 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdzp    !: Limitation of diazo Nfix by P
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimdiazo  !: Limitation of diazo PP by nutrients
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   qfenfixdz  !: Fe cost of nitrogen fixation
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   xlimnpdz    !: Limitation of Nfix by N
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   Fe_scale_nfix !: Temp dependant Fe cost scaling
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  ::   Pmin_scale_nfix !: Temp dependent QPmin scaling
   ! Coefficient for iron limitation following Flynn and Hipkin (1999)
   REAL(wp) ::  xcoef1   = 0.00167  / 55.85
   REAL(wp) ::  xcoef2   = 1.21E-5 * 14. / 55.85 / 7.625 * 0.5 * 1.5
   REAL(wp) ::  xcoef3   = 1.15E-4 * 14. / 55.85 / 7.625 * 0.5
   REAL(wp) ::  xcoef4   = 7.5E-4  * 14. / 55.85 / 7.625 * 0.5 
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p5zlim.F90 10070 2018-08-28 14:30:54Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_lim( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_lim  ***
      !!
      !! ** Purpose :   Compute the co-limitations by the various nutrients
      !!                for the various phytoplankton species. Quota based
      !!                approach. The quota model is derived from theoretical
      !!                models proposed by Pahlow and Oschlies (2009) and 
      !!                Flynn (2001). Various adaptations from several 
      !!                publications by these authors have been also adopted.
      !!                Explicit Diazotroph Added 
      !!
      !! ** Method  : Quota based approach. The quota model is derived from 
      !!              theoretical models by Pahlow and Oschlies (2009) and 
      !!              Flynn (2001). Various adaptations from several publications
      !!              by these authors have been also adopted.
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in)  :: kt, knt
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zlim1, zlim2, zlim3, zlim4, zno3, zferlim
      REAL(wp) ::   z1_trndia, z1_trnpic, z1_trnphy, ztem1, ztem2, zetot1
      REAL(wp) ::   zratio, zration, zratiof, znutlim, zfalim, zzpsiuptk
      REAL(wp) ::   zconc1d, zconc1dnh4, zconc0n, zconc0nnh4, zconc0npo4, zconc0dpo4
      REAL(wp) ::   zconc0p, zconc0pnh4, zconc0ppo4, zconcpfe, zconcnfe, zconcdfe
      REAL(wp) ::   fanano, fananop, fananof, fadiat, fadiatp, fadiatf
      REAL(wp) ::   fapico, fapicop, fapicof
      REAL(wp) ::   zrpho, zrass, zcoef, zfuptk, zratchl
      REAL(wp) ::   zfvn, zfvp, zfvf, zsizen, zsizep, zsized, znanochl, zpicochl, zdiatchl
      REAL(wp) ::   zqfemn, zqfemp, zqfemd, zbactno3, zbactnh4
      REAL(wp) ::   zlim1f, zsizetmp
      !! Explicit diazotroph PFT
      REAL(wp) ::   z1_trndiaz, zconc0dz, zconc0dznh4, zconc0dzpo4, zconcdzfe, zzpsinfix
      REAL(wp) ::   fadiaz, fadiazp, fadiazf, zsizedz, zdiazchl, zqfemdz, zratiop, qfenfix
      REAL(wp) ::   facul, IUE_scale, PUE_scale
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: xlimnpn, xlimnpp, xlimnpd
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_lim')
      !
      zratchl = 6.0
      sizena(:,:,:) = 0.0  ;  sizepa(:,:,:) = 0.0  ;  sizeda(:,:,:) = 0.0
      sizedza(:,:,:) = 0.0;
      !! Ensuring size terms never become 0
      sizen(:,:,:)=max(1.,sizen(:,:,:))
      sized(:,:,:)=max(1.,sized(:,:,:))
      sizep(:,:,:)=max(1.,sizep(:,:,:))
      sizedz(:,:,:)=max(1.,sizedz(:,:,:))!
      
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! 
               ! Tuning of the iron concentration to a minimum level that
               ! is set to the detection limit
               ! --------------------------------------------------------
               zno3    = trb(ji,jj,jk,jpno3) / 40.e-6
               zferlim = MAX( 3e-11 * zno3 * zno3, 5e-12 )
               zferlim = MIN( zferlim, 7e-11 )
               trb(ji,jj,jk,jpfer) = MAX( trb(ji,jj,jk,jpfer), zferlim )

               ! Computation of the Chl/C ratio of each phytoplankton group
               ! ----------------------------------------------------------
               z1_trnphy   = 1. / ( trb(ji,jj,jk,jpphy) + rtrn )
               z1_trnpic   = 1. / ( trb(ji,jj,jk,jppic) + rtrn )
               z1_trndia   = 1. / ( trb(ji,jj,jk,jpdia) + rtrn )
               z1_trndiaz  = 1. / ( trb(ji,jj,jk,jpcdz) + rtrn )
               znanochl = trb(ji,jj,jk,jpnch) * z1_trnphy
               zpicochl = trb(ji,jj,jk,jppch) * z1_trnpic
               zdiatchl = trb(ji,jj,jk,jpdch) * z1_trndia
               zdiazchl = trb(ji,jj,jk,jpchd) * z1_trndiaz

               ! Computation of a variable Ks for the different phytoplankton
               ! group as a function of their relative size. Allometry
               ! from Edwards et al. (2012)
               ! ------------------------------------------------------------

               ! diatoms
               zsized            = sized(ji,jj,jk)**0.81
               zconcdfe          = concdfer * zsized
               zconc1d           = concdno3 * zsized
               zconc1dnh4        = concdnh4 * zsized
               zconc0dpo4        = concdpo4 * zsized

               ! picophytoplankton
               zsizep            = sizep(ji,jj,jk)**0.81
               zconcpfe          = concpfer * zsizep
               zconc0p           = concpno3 * zsizep
               zconc0pnh4        = concpnh4 * zsizep
               zconc0ppo4        = concppo4 * zsizep

               ! nanophytoplankton
               zsizen            = sizen(ji,jj,jk)**0.81
               zconcnfe          = concnfer * zsizen
               zconc0n           = concnno3 * zsizen
               zconc0nnh4        = concnnh4 * zsizen
               zconc0npo4        = concnpo4 * zsizen

               ! Explicit Diazotroph PFT
               zsizedz           = sizedz(ji,jj,jk)**0.81
               zconcdzfe         = concdzfer * zsizedz
               zconc0dz          = concdzno3 * zsizedz
               zconc0dznh4       = concdznh4 * zsizedz
               zconc0dzpo4       = concdzpo4 * zsizedz

               ! Allometric variations of the minimum and maximum quotas
               ! From Talmy et al. (2014) and Maranon et al. (2013)
               ! -------------------------------------------------------
               xqnnmin(ji,jj,jk) = qnnmin  * sizen(ji,jj,jk)**(-0.3)
               xqnnmax(ji,jj,jk) = qnnmax
               xqndmin(ji,jj,jk) = qndmin * sized(ji,jj,jk)**(-0.3)
               xqndmax(ji,jj,jk) = qndmax
               xqnpmin(ji,jj,jk) = qnpmin * sizep(ji,jj,jk)**(-0.48)
               xqnpmax(ji,jj,jk) = qnpmax * sizep(ji,jj,jk)**(-0.21)
               xqndzmin(ji,jj,jk) = qndzmin * sizedz(ji,jj,jk)**(-0.3)
               xqndzmax(ji,jj,jk) = qndzmax

               ! Computation of the optimal allocation parameters
               ! Based on the different papers by Pahlow et al., and 
               ! Smith et al.
               ! ---------------------------------------------------

               ! Nanophytoplankton
               znutlim = MAX( trb(ji,jj,jk,jpnh4) / zconc0nnh4,    &
                 &         trb(ji,jj,jk,jpno3) / zconc0n)
               fanano = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,jk,jppo4) / zconc0npo4
               fananop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcnfe
               fananof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

               ! Picophytoplankton
               znutlim = MAX( trb(ji,jj,jk,jpnh4) / zconc0pnh4,    &
                 &         trb(ji,jj,jk,jpno3) / zconc0p)
               fapico = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,jk,jppo4) / zconc0ppo4
               fapicop = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcpfe
               fapicof = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

               ! Diatoms
               znutlim = MAX( trb(ji,jj,jk,jpnh4) / zconc1dnh4,    &
                 &         trb(ji,jj,jk,jpno3) / zconc1d )
               fadiat = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,jk,jppo4) / zconc0dpo4
               fadiatp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcdfe
               fadiatf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
              
               ! Explicit Diazotrophs
               znutlim = MAX( trb(ji,jj,jk,jpnh4) / zconc0dznh4,    &
                 &         trb(ji,jj,jk,jpno3) / zconc0dz)
               fadiaz  = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = trb(ji,jj,jk,jppo4) / zconc0dzpo4
               fadiazp = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )
               znutlim = biron(ji,jj,jk) / zconcdzfe
               fadiazf = MAX(0.01, MIN(0.99, 1. / ( SQRT(znutlim) + 1.) ) )

               ! Michaelis-Menten Limitation term by nutrients of
               ! heterotrophic bacteria
               ! -------------------------------------------------
               zbactnh4 = trb(ji,jj,jk,jpnh4) / ( concbnh4 + trb(ji,jj,jk,jpnh4) )
               zbactno3 = trb(ji,jj,jk,jpno3) / ( concbno3 + trb(ji,jj,jk,jpno3) ) * (1. - zbactnh4)
               !
               zlim1    = zbactno3 + zbactnh4
               zlim2    = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concbpo4)
               zlim3    = biron(ji,jj,jk) / ( concbfe + biron(ji,jj,jk) )
               zlim4    = trb(ji,jj,jk,jpdoc) / ( xkdoc   + trb(ji,jj,jk,jpdoc) )
               ! Xlimbac is used for DOC solubilization whereas xlimbacl
               ! is used for all the other bacterial-dependent terms
               ! -------------------------------------------------------
               xlimbacl(ji,jj,jk) = MIN( zlim1, zlim2, zlim3 )
               xlimbac (ji,jj,jk) = xlimbacl(ji,jj,jk) * zlim4
               
               ! Michaelis-Menten Limitation term by nutrients: Nanophyto
               ! --------------------------------------------------------
               !
               ! Limitation of N based nutrients uptake (NO3 and NH4)
               zfalim = (1.-fanano) / fanano
               xnanonh4(ji,jj,jk) = (1. - fanano) * trb(ji,jj,jk,jpnh4) / ( zfalim * zconc0nnh4 + trb(ji,jj,jk,jpnh4) )
               xnanono3(ji,jj,jk) = (1. - fanano) * trb(ji,jj,jk,jpno3) / ( zfalim * zconc0n + trb(ji,jj,jk,jpno3) )  &
               &                    * (1. - xnanonh4(ji,jj,jk))
               !
               ! Limitation of P based nutrients (PO4 and DOP)
               zfalim = (1.-fananop) / fananop
               xnanopo4(ji,jj,jk) = (1. - fananop) * trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zfalim * zconc0npo4 )
               !xnanodop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop) + xkdoc )   &
               !&                    * ( 1.0 - xnanopo4(ji,jj,jk) )
               xnanodop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / (trb(ji,jj,jk,jpdop) + xkdop )   &
               &                    * ( 1.0 - xnanopo4(ji,jj,jk) )
               !xnanodop(ji,jj,jk) = 0.
               !
               ! Limitation of Fe uptake
               zfalim = (1.-fananof) / fananof
               xnanofer(ji,jj,jk) = (1. - fananof) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcnfe )
               !
               ! The minimum iron quota depends on the size of PSU, respiration
               ! and the reduction of nitrate following the parameterization 
               ! proposed by Flynn and Hipkin (1999)
               zratiof   = trb(ji,jj,jk,jpnfe) * z1_trnphy
               zqfemn = xcoef1 * znanochl + xcoef2 + xcoef3 * xnanono3(ji,jj,jk)
               qfeminn(ji,jj,jk) = zqfemn
               !
               zration = trb(ji,jj,jk,jpnph) * z1_trnphy
               zration = MIN(xqnnmax(ji,jj,jk), MAX( xqnnmin(ji,jj,jk), zration ))
               zzpsiuptk = xqnnmin(ji,jj,jk) * rno3 / zpsiuptk**2
               fvnuptk(ji,jj,jk) = 1. / zzpsiuptk * xqnnmin(ji,jj,jk) / (zration + rtrn)  &
               &                   * MAX(0., (1. - zratchl * znanochl / 12. ) )
               !
               zlim1  = max(0., (zration - xqnnmin(ji,jj,jk) )  &
               &          / (xqnnmax(ji,jj,jk) - xqnnmin(ji,jj,jk) ) ) * xqnnmax(ji,jj,jk)  &
               &          / (zration + rtrn)
               ! The value of the optimal quota in the formulation below
               ! has been found by solving a non linear equation
               zlim1f = max(0., ( 1.086 - xqnnmin(ji,jj,jk) )  &
               &          / (xqnnmax(ji,jj,jk) - xqnnmin(ji,jj,jk) ) ) * xqnnmax(ji,jj,jk)
               zlim3  = MAX( 0.,( zratiof - zqfemn ) / qfnopt )
               ! computation of the various limitation terms of nanophyto
               ! growth and PP
               xlimnfe (ji,jj,jk) = MIN( 1., zlim3 )
               xlimphy (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
               xlimphys(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )
               xlimnpn (ji,jj,jk) = MIN( 1., zlim1)

               ! Michaelis-Menten Limitation term by nutrients: Picophyto
               ! --------------------------------------------------------
               !
               ! Limitation of N based nutrients uptake (NO3 and NH4) 
               zfalim = (1.-fapico) / fapico 
               xpiconh4(ji,jj,jk) = (1. - fapico) * trb(ji,jj,jk,jpnh4) / ( zfalim * zconc0pnh4 + trb(ji,jj,jk,jpnh4) )
               xpicono3(ji,jj,jk) = (1. - fapico) * trb(ji,jj,jk,jpno3) / ( zfalim * zconc0p + trb(ji,jj,jk,jpno3) )  &
               &                    * (1. - xpiconh4(ji,jj,jk))
               !
               ! Limitation of P based nutrients uptake (PO4 and DOP)
               zfalim = (1.-fapicop) / fapicop 
               xpicopo4(ji,jj,jk) = (1. - fapicop) * trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zfalim * zconc0ppo4 )
               !xpicodop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop) + xkdoc )   &
               !&                    * ( 1.0 - xpicopo4(ji,jj,jk) )
               xpicodop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop) + xkdop )   &
               &                    * ( 1.0 - xpicopo4(ji,jj,jk) )               
               !xpicodop(ji,jj,jk) = 0.
               !
               zfalim = (1.-fapicof) / fapicof
               xpicofer(ji,jj,jk) = (1. - fapicof) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcpfe )
               !
               ! The minimum iron quota depends on the size of PSU, respiration
               ! and the reduction of nitrate following the parameterization 
               ! proposed by Flynn and Hipkin (1999)
               zratiof   = trb(ji,jj,jk,jppfe) * z1_trnpic
               zqfemp = xcoef1 * zpicochl + xcoef2 + xcoef3 * xpicono3(ji,jj,jk)
               qfeminp(ji,jj,jk) = zqfemp
               !
               zration   = trb(ji,jj,jk,jpnpi) * z1_trnpic
               zration = MIN(xqnpmax(ji,jj,jk), MAX( xqnpmin(ji,jj,jk), zration ))
               zzpsiuptk = xqnpmin(ji,jj,jk) * rno3 / zpsiuptk**2
               fvpuptk(ji,jj,jk) = 1. / zzpsiuptk * xqnpmin(ji,jj,jk) / (zration + rtrn)  &
               &                   * MAX(0., (1. - zratchl * zpicochl / 12. ) ) 
               !
               zlim1    = max(0., (zration - xqnpmin(ji,jj,jk) )  &
               &          / (xqnpmax(ji,jj,jk) - xqnpmin(ji,jj,jk) ) ) * xqnpmax(ji,jj,jk)  &
               &          / (zration + rtrn)
               ! The value of the optimal quota in the formulation below
               ! has been found by solving a non linear equation
               zlim1f   = max(0., (1.367 - xqnpmin(ji,jj,jk) )  &
               &          / (xqnpmax(ji,jj,jk) - xqnpmin(ji,jj,jk) ) ) * xqnpmax(ji,jj,jk)
               zlim3    = MAX( 0.,( zratiof - zqfemp ) / qfpopt )

               ! computation of the various limitation terms of picophyto
               ! growth and PP
               xlimpfe (ji,jj,jk) = MIN( 1., zlim3 )
               xlimpic (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
               xlimnpp (ji,jj,jk) = MIN( 1., zlim1 )
               xlimpics(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )

               ! Michaelis-Menten Limitation term by nutrients : Diatoms
               ! -------------------------------------------------------
               !
               ! Limitation of N based nutrients uptake (NO3 and NH4)
               zfalim = (1.-fadiat) / fadiat 
               xdiatnh4(ji,jj,jk) = (1. - fadiat) * trb(ji,jj,jk,jpnh4) / ( zfalim * zconc1dnh4 + trb(ji,jj,jk,jpnh4) )
               xdiatno3(ji,jj,jk) = (1. - fadiat) * trb(ji,jj,jk,jpno3) / ( zfalim * zconc1d + trb(ji,jj,jk,jpno3) )  &
               &                    * (1. - xdiatnh4(ji,jj,jk))
               !
               ! Limitation of P based nutrients uptake (PO4 and DOP)
               zfalim = (1.-fadiatp) / fadiatp
               xdiatpo4(ji,jj,jk) = (1. - fadiatp) * trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zfalim * zconc0dpo4 )
               !xdiatdop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop) + xkdoc )  &
               !&                    * ( 1.0 - xdiatpo4(ji,jj,jk) )
               xdiatdop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop) + xkdop )  &
               &                    * ( 1.0 - xdiatpo4(ji,jj,jk) )
               !xdiatdop(ji,jj,jk) = 0.
               !
               ! Limitation of Fe uptake
               zfalim = (1.-fadiatf) / fadiatf
               xdiatfer(ji,jj,jk) = (1. - fadiatf) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcdfe )
               !
               ! The minimum iron quota depends on the size of PSU, respiration
               ! and the reduction of nitrate following the parameterization 
               ! proposed by Flynn and Hipkin (1999)
               zratiof   = trb(ji,jj,jk,jpdfe) * z1_trndia
               zqfemd = xcoef1 * zdiatchl + xcoef2 + xcoef3 * xdiatno3(ji,jj,jk)
               qfemind(ji,jj,jk) = zqfemd               
               !
               zration   = trb(ji,jj,jk,jpndi) * z1_trndia
               zration   = MIN(xqndmax(ji,jj,jk), MAX( xqndmin(ji,jj,jk), zration ))
               zzpsiuptk = xqndmin(ji,jj,jk) * rno3 / zpsiuptk**2
               fvduptk(ji,jj,jk) = 1. / zzpsiuptk * xqndmin(ji,jj,jk) / (zration + rtrn)   &
               &                   * MAX(0., (1. - zratchl * zdiatchl / 12. ) ) 
               !
               zlim1    = max(0., (zration - xqndmin(ji,jj,jk) )    &
               &          / (xqndmax(ji,jj,jk) - xqndmin(ji,jj,jk) ) )   &
               &          * xqndmax(ji,jj,jk) / (zration + rtrn)
               ! The value of the optimal quota in the formulation below
               ! has been found by solving a non linear equation
               zlim1f   = max(0., (1.077 - xqndmin(ji,jj,jk) )    &
               &          / (xqndmax(ji,jj,jk) - xqndmin(ji,jj,jk) ) )   &
               &          * xqndmax(ji,jj,jk)
               zlim3    = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi(ji,jj) )
               zlim4    = MAX( 0., ( zratiof - zqfemd ) / qfdopt )
               ! computation of the various limitation terms of diatoms
               ! growth and PP
               xlimdfe(ji,jj,jk) = MIN( 1., zlim4 )
               xlimdia(ji,jj,jk) = MIN( 1., zlim1, zlim3, zlim4 )
               xlimdias(ji,jj,jk) = MIN (1.0, zlim1 / (zlim1f + rtrn ), zlim3, zlim4 )
               xlimsi(ji,jj,jk)  = MIN( zlim1, zlim4 )
               xlimnpd(ji,jj,jk) = MIN( 1., zlim1 )
               !
               ! Michaelis-Menten Limitation term for nutrients diazotrophs
               ! -----------------------------------------------
               zfalim = (1.-fadiaz) / fadiaz
               xdiaznh4(ji,jj,jk) = (1. - fadiaz) * trb(ji,jj,jk,jpnh4) / ( zfalim * zconc0dznh4 + trb(ji,jj,jk,jpnh4) )
               xdiazno3(ji,jj,jk) = (1. - fadiaz) * trb(ji,jj,jk,jpno3) / ( zfalim * zconc0dz + trb(ji,jj,jk,jpno3) )  &
               &                    * (1. - xdiaznh4(ji,jj,jk))
               !
               zfalim = (1.-fadiazp) / fadiazp
               xdiazpo4(ji,jj,jk) = (1. - fadiazp) * trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + zfalim * zconc0dzpo4 )
               xdiazdop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdop)+ xkdopdz )   &
               &                    * ( 1.0 - xdiazpo4(ji,jj,jk) )
               !xdiazdop(ji,jj,jk) = 0. !Remove to allow DOP uptake
               !
               ! Limitation of Fe uptake
               zfalim = (1.-fadiazf) / fadiazf
               xdiazfer(ji,jj,jk) = (1. - fadiazf) * biron(ji,jj,jk) / ( biron(ji,jj,jk) + zfalim * zconcdzfe )
               !
               ! The minimum iron quota depends on the size of PSU, respiration
               ! and the reduction of nitrate following the parameterization 
               ! proposed by Flynn and Hipkin (1999)
               zratiof   = trb(ji,jj,jk,jpfed) * z1_trndiaz
               
               ! Temperature dependence of diazotroph IUE of nitrogen fixation
               !
               IF ( ln_tiue ) THEN
               ! Trichodesmium values based upon Jiang et al. (2018)
                 IF ( ln_tricho ) THEN
               !Tricho Temp driven Fe cost scaling set to 1 where Trico growth = 0.1 d-1 ~20.5oC
                    IUE_scale = (-(0.001392*(tsn(ji,jj,jk,jp_tem)**5)) + (0.1559*(tsn(ji,jj,jk,jp_tem)**4)) &
               &                -(6.7685*(tsn(ji,jj,jk,jp_tem)**3)) + (141.81*(tsn(ji,jj,jk,jp_tem)**2))   &
               &                -(1421.1*tsn(ji,jj,jk,jp_tem)) + 5388.1)/33.49
                    Fe_scale_nfix(ji,jj,jk) = 1/(IUE_scale+rtrn)
                    IF(tsn(ji,jj,jk,jp_tem) .LT. 17 .OR. tsn(ji,jj,jk,jp_tem).GT. 35 ) THEN
                      Fe_scale_nfix(ji,jj,jk) = 0.
                    ENDIF
                    IF (Fe_scale_nfix(ji,jj,jk) .GT. maxFescale ) THEN
                      Fe_scale_nfix(ji,jj,jk) = maxFescale
                    ENDIF
               ! Crocosphaera values based upon Yang et al. (2021)
                 ELSE
               !Croco Temp driven Fe cost scaling set to 1 where Croco growth =0.1 d-1 ~21.3oC
                    IUE_scale = ((0.02092*(tsn(ji,jj,jk,jp_tem)**4)) - (2.302*(tsn(ji,jj,jk,jp_tem)**3)) &
               &                +(92.08*(tsn(ji,jj,jk,jp_tem)**2)) -(1582*tsn(ji,jj,jk,jp_tem))   &
               &                 + 9881)/20.64
                    Fe_scale_nfix(ji,jj,jk) = 1/(IUE_scale+rtrn)
                    IF(tsn(ji,jj,jk,jp_tem) .LT. 20 .OR. tsn(ji,jj,jk,jp_tem).GT. 35 ) THEN
                      Fe_scale_nfix(ji,jj,jk) = 0.
                    ENDIF
                    IF ( Fe_scale_nfix(ji,jj,jk) .GT. maxFescale ) THEN
                      Fe_scale_nfix(ji,jj,jk) = maxFescale
                    ENDIF
                 ENDIF

               ELSE
               !No temperature dependance on Fe cost
                 Fe_scale_nfix(ji,jj,jk) = 1

               ENDIF
              
               facul = (1-(xdiazno3(ji,jj,jk)+xdiaznh4(ji,jj,jk)))  ! Diazotroph Facultative term (1 = all nfix, 0 = no nfix)
               qfenfix = kustkaFe * facul * Fe_scale_nfix(ji,jj,jk) ! IUE scales the Fe cost of Nfix based upon Kustka et al. (2003)
                                                                    ! kustkaFe = 13 umol/mol 
               zqfemdz = xcoef1 * zdiazchl + xcoef2 + xcoef3 * xdiazno3(ji,jj,jk) + qfenfix
               qfemindz(ji,jj,jk) = zqfemdz
               !
               zration = trb(ji,jj,jk,jpndz) * z1_trndiaz
               zration = MIN(xqndzmax(ji,jj,jk), MAX( xqndzmin(ji,jj,jk), zration ))
               zzpsiuptk = xqndzmin(ji,jj,jk) * rno3 / zpsiuptk**2
               fvdzuptk(ji,jj,jk) = 1. / zzpsiuptk * xqndzmin(ji,jj,jk) / (zration + rtrn)  &
               &                   * MAX(0., (1. - zratchl * zdiazchl / 12. ) )
               !
               zlim1  = max(0., (zration - xqndzmin(ji,jj,jk) )  &
               &          / (xqndzmax(ji,jj,jk) - xqndzmin(ji,jj,jk) ) ) * xqndzmax(ji,jj,jk)  &
               &          / (zration + rtrn)
               ! The value of the optimal quota in the formulation below
               ! has been found by solving a non linear equation
               zlim1f = max(0., ( 1.086 - xqndzmin(ji,jj,jk) )  &
               &          / (xqndzmax(ji,jj,jk) - xqndzmin(ji,jj,jk) ) ) * xqndzmax(ji,jj,jk)
               zlim3  = MAX( 0.,( zratiof - zqfemdz ) / qfdzopt )
               ! computation of the various limitation terms of nanophyto
               ! growth and PP
               xlimdzfe (ji,jj,jk) = MIN( 1., zlim3 )
               xlimdiaz (ji,jj,jk) = MIN( 1., zlim1, zlim3 )
               xlimdiazo(ji,jj,jk) = MIN( 1., zlim1/( zlim1f + rtrn ), zlim3 )
               xlimnpdz (ji,jj,jk) = MIN( 1., zlim1)
            END DO
         END DO
      END DO
      !
      ! Compute the phosphorus quota values. It is based on Litchmann et al., 2004 and Daines et al, 2013.
      ! The relative contribution of three fonctional pools are computed: light harvesting apparatus, 
      ! nutrient uptake pool and assembly machinery. DNA is assumed to represent 1% of the dry mass of 
      ! phytoplankton (see Daines et al., 2013). 
      ! --------------------------------------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Size estimation of nanophytoplankton based on total biomass
               ! Assumes that larger biomass implies addition of larger cells
               ! ------------------------------------------------------------
               zcoef = trb(ji,jj,jk,jpphy) - MIN(xsizephy, trb(ji,jj,jk,jpphy) )
               sizena(ji,jj,jk) = 1. + ( xsizern -1.0 ) * zcoef / ( xsizephy + zcoef+rtrn )
               ! N/P ratio of nanophytoplankton
               ! ------------------------------
               zfuptk = 0.2 + 0.12 / ( 3.0 * sizen(ji,jj,jk) + rtrn )
               zrpho  = 1.54 * trb(ji,jj,jk,jpnch) / ( trb(ji,jj,jk,jpnph) * rno3 * 14. + rtrn )
               zrass = MAX(0.62/4., ( 1. - zrpho - zfuptk ) * xlimnpn(ji,jj,jk) )
               xqpnmin(ji,jj,jk) = ( 0.0 + 0.0078 + 0.62/4. * 0.0783 * xqnnmin(ji,jj,jk) ) * 16.
               xqpnmax(ji,jj,jk) = ( zrpho * 0.0128 + zrass * 0.0783 ) * 16.
               xqpnmax(ji,jj,jk) = xqpnmax(ji,jj,jk) * trb(ji,jj,jk,jpnph) / ( trb(ji,jj,jk,jpphy) + rtrn )  &
               &      + (0.033 + 0.0078 ) * 16.
               xqpnmax(ji,jj,jk) = MIN( qpnmax, xqpnmax(ji,jj,jk) )


               ! Size estimation of picophytoplankton based on total biomass
               ! Assumes that larger biomass implies addition of larger cells
               ! ------------------------------------------------------------
               zcoef = trb(ji,jj,jk,jppic) - MIN(xsizepic, trb(ji,jj,jk,jppic) )
               sizepa(ji,jj,jk) = 1. + ( xsizerp -1.0 ) * zcoef / ( xsizepic + zcoef )

               ! N/P ratio of picophytoplankton
               ! ------------------------------
               zfuptk = 0.2 + 0.12 / ( 0.5 * sizep(ji,jj,jk) + rtrn )
               zrpho = 1.54 * trb(ji,jj,jk,jppch) / ( trb(ji,jj,jk,jpnpi) * rno3 * 14. + rtrn )
               zrass = MAX(0.4/4., ( 1. - zrpho - zfuptk ) * xlimnpp(ji,jj,jk) )
               xqppmin(ji,jj,jk) = ( (0.0 + 0.0078 ) + 0.4/4. * 0.0517 * xqnpmin(ji,jj,jk) ) * 16.
               xqppmax(ji,jj,jk) = ( zrpho * 0.0128 + zrass * 0.0517 ) * 16.
               xqppmax(ji,jj,jk) = xqppmax(ji,jj,jk) * trb(ji,jj,jk,jpnpi) / ( trb(ji,jj,jk,jppic) + rtrn ) &
               &      +  (0.033 + 0.0078 ) * 16
               xqppmax(ji,jj,jk) = MIN( qppmax, xqppmax(ji,jj,jk) )

               ! Size estimation of diatoms based on total biomass
               ! Assumes that larger biomass implies addition of larger cells
               ! ------------------------------------------------------------
               zcoef = trb(ji,jj,jk,jpdia) - MIN(xsizedia, trb(ji,jj,jk,jpdia) )
               sizeda(ji,jj,jk) = 1. + ( xsizerd - 1.0 ) * zcoef / ( xsizedia + zcoef )

               ! N/P ratio of diatoms
               ! --------------------
               zfuptk = 0.2 + 0.12 / ( 5.0 * sized(ji,jj,jk) + rtrn )
               zrpho = 1.54 * trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpndi) * rno3 * 14. + rtrn )
               zrass = MAX(0.66/4., ( 1. - zrpho - zfuptk ) * xlimnpd(ji,jj,jk) )

               xqpdmin(ji,jj,jk) = ( ( 0.0 + 0.0078 ) + 0.66/4. * 0.0783 *  xqndmin(ji,jj,jk) ) * 16.
               xqpdmax(ji,jj,jk) = ( zrpho * 0.0128 + zrass * 0.0783 ) * 16.
               xqpdmax(ji,jj,jk) = xqpdmax(ji,jj,jk) * trb(ji,jj,jk,jpndi) / ( trb(ji,jj,jk,jpdia) + rtrn ) &
               &      + ( 0.0078 + 0.033 ) * 16.
               xqpdmax(ji,jj,jk) = MIN(qpdmax, xqpdmax(ji,jj,jk) )

               ! Size estimation of diazotrophs
               ! ------------------------------------
               zcoef = trb(ji,jj,jk,jpcdz) - MIN(xsizedz, trb(ji,jj,jk,jpcdz) )
               sizedza(ji,jj,jk) = 1. + ( xsizerdz -1.0 ) * zcoef / ( xsizedz + zcoef )
               ! N/P ratio of diazotrophs
               ! --------------------
               zfuptk = 0.2 + 0.12 / ( 3.0 * sizedz(ji,jj,jk) + rtrn )
               zrpho  = 1.54 * trb(ji,jj,jk,jpchd) / ( trb(ji,jj,jk,jpndz) *rno3 * 14. + rtrn )
               zrass = MAX(0.62/4., ( 1. - zrpho - zfuptk ) *xlimnpdz(ji,jj,jk))
               xqpdzmin(ji,jj,jk) = ( 0.0 + 0.0078 + 0.62/4. * 0.0783 *xqndzmin(ji,jj,jk) ) * 16.
               xqpdzmax(ji,jj,jk) = ( zrpho * 0.0128 + zrass * 0.0783 ) * 16.
               xqpdzmax(ji,jj,jk) = xqpdzmax(ji,jj,jk) * trb(ji,jj,jk,jpndz) / (trb(ji,jj,jk,jpcdz) + rtrn )  &
               &      + (0.033 + 0.0078 ) * 16.
               xqpdzmax(ji,jj,jk) = MIN( qpdzmax, xqpdzmax(ji,jj,jk) )

               ! Temperature dependence of diazotroph PUE of nitrogen fixation
               IF ( ln_tpue ) THEN
               ! Trichodesmium PUE values based upon Jiang et al. (2018)
                 IF (ln_tricho ) THEN
                 ! Tricho Temp driven QPmin scaling set to 1 where Tricho growth = 0.1 d-1 20.5oC
                    PUE_scale = ((-0.000133*(tsn(ji,jj,jk,jp_tem)**4)) + (0.012452*(tsn(ji,jj,jk,jp_tem)**3)) &
                    &            -(0.4294*(tsn(ji,jj,jk,jp_tem)**2)) + (6.538*tsn(ji,jj,jk,jp_tem))       &
                    &            - 37.11)/0.25
                    Pmin_scale_nfix(ji,jj,jk) = 1/(PUE_scale+rtrn)
                    IF (tsn(ji,jj,jk,jp_tem) .LT. 17 .OR. tsn(ji,jj,jk,jp_tem).GT. 35 ) THEN
                       Pmin_scale_nfix(ji,jj,jk) = 0.
                    ENDIF
                    IF ( Pmin_scale_nfix(ji,jj,jk) .GT. maxPminscale ) THEN
                       Pmin_scale_nfix(ji,jj,jk) = maxPminscale
                    ENDIF
                ! Crocosphaera PUE values based upon Yang et al. (2021)
                 ELSE
                 !Croco Temp driven QPmin scaling set to 1 where Croco growth =0.1 d-1 ~21.3oC       
                    PUE_scale = ((-0.0004429*(tsn(ji,jj,jk,jp_tem)**4)) + (0.04684*(tsn(ji,jj,jk,jp_tem)**3)) &
                    &            -(1.83905*(tsn(ji,jj,jk,jp_tem)**2)) + (31.829*tsn(ji,jj,jk,jp_tem))       &
                    &            - 204.815)/0.2628
                    Pmin_scale_nfix(ji,jj,jk) = 1/(PUE_scale+rtrn)
                    IF (tsn(ji,jj,jk,jp_tem) .LT. 20 .OR. tsn(ji,jj,jk,jp_tem).GT. 35 ) THEN
                       Pmin_scale_nfix(ji,jj,jk) = 0.
                    ENDIF
                    IF ( Pmin_scale_nfix(ji,jj,jk) .GT. maxPminscale ) THEN
                       Pmin_scale_nfix(ji,jj,jk) = maxPminscale
                    ENDIF
                  ENDIF

                ELSE
                ! No temperature dependance on QPmin
                  Pmin_scale_nfix(ji,jj,jk) = 1

                ENDIF
               xqpdzmin(ji,jj,jk) = xqpdzmin(ji,jj,jk) * Pmin_scale_nfix(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Compute the fraction of nanophytoplankton that is made of calcifiers
      ! This is a purely adhoc formulation described in Aumont et al. (2015)
      ! This fraction depends on nutrient limitation, light, temperature
      ! --------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zlim1 =  trb(ji,jj,jk,jpnh4) / ( trb(ji,jj,jk,jpnh4) + concnnh4 ) + trb(ji,jj,jk,jpno3)    &
               &        / ( trb(ji,jj,jk,jpno3) + concnno3 ) * ( 1.0 - trb(ji,jj,jk,jpnh4)   &
               &        / ( trb(ji,jj,jk,jpnh4) + concnnh4 ) )
               zlim2  = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concnpo4 )
               zlim3  = trb(ji,jj,jk,jpfer) / ( trb(ji,jj,jk,jpfer) +  5.E-11 ) 
               ztem1  = MAX( 0., tsn(ji,jj,jk,jp_tem) )
               ztem2  = tsn(ji,jj,jk,jp_tem) - 10.
               zetot1 = MAX( 0., etot_ndcy(ji,jj,jk) - 1.) / ( 4. + etot_ndcy(ji,jj,jk) ) * 20. / ( 20. + etot_ndcy(ji,jj,jk) ) 

               xfracal(ji,jj,jk) = caco3r * MIN( zlim1, zlim2, zlim3 )    &
               &                   * ztem1 / ( 1. + ztem1 ) * MAX( 1., trb(ji,jj,jk,jpphy)*1E6 )   &
                  &                * ( 1. + EXP(-ztem2 * ztem2 / 25. ) )         &
                  &                * zetot1 * MIN( 1., 50. / ( hmld(ji,jj) + rtrn ) )
               xfracal(ji,jj,jk) = MAX( 0.02, MIN( 0.8 , xfracal(ji,jj,jk) ) )
            END DO
         END DO
      END DO
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! denitrification factor computed from O2 levels
               nitrfac(ji,jj,jk) = MAX(  0.e0, 0.4 * ( 6.e-6  - trb(ji,jj,jk,jpoxy) )    &
                  &                                / ( oxymin + trb(ji,jj,jk,jpoxy) )  )
               nitrfac(ji,jj,jk) = MIN( 1., nitrfac(ji,jj,jk) )
               !
               ! redox factor computed from NO3 levels
               nitrfac2(ji,jj,jk) = MAX( 0.e0,       ( 1.E-6 - trb(ji,jj,jk,jpno3) )  &
                  &                                / ( 1.E-6 + trb(ji,jj,jk,jpno3) ) )
               nitrfac2(ji,jj,jk) = MIN( 1., nitrfac2(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN        ! save output diagnostics
        IF( iom_use( "xfracal" ) ) CALL iom_put( "xfracal", xfracal(:,:,:) * tmask(:,:,:) )  ! euphotic layer deptht
        IF( iom_use( "LNnut"   ) ) CALL iom_put( "LNnut"  , xlimphy(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LPnut"   ) ) CALL iom_put( "LPnut"  , xlimpic(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LDnut"   ) ) CALL iom_put( "LDnut"  , xlimdia(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LNFe"    ) ) CALL iom_put( "LNFe"   , xlimnfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LPFe"    ) ) CALL iom_put( "LPFe"   , xlimpfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDFe"    ) ) CALL iom_put( "LDFe"   , xlimdfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZEN"   ) ) CALL iom_put( "SIZEN"  , sizen(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZEP"   ) ) CALL iom_put( "SIZEP"  , sizep(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZED"   ) ) CALL iom_put( "SIZED"  , sized(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "SIZEDZ"   ) ) CALL iom_put( "SIZEDZ"  , sizedz(:,:,:) * tmask(:,:,:) )  ! Iron limitation term
        IF( iom_use( "LDZnut"   ) ) CALL iom_put( "LDZnut"  , xlimdiaz(:,:,:) * tmask(:,:,:) )  ! Nutrient limitation term
        IF( iom_use( "LDZFe"    ) ) CALL iom_put( "LDZFe"   , xlimdzfe(:,:,:) * tmask(:,:,:) )  ! Iron limitation term 
        IF( iom_use( "LDZNP"     ) ) CALL iom_put( "LDZNP"     , xlimnpdz(:,:,:) * tmask(:,:,:) )  ! N/P limitation diazo
        IF( iom_use( "QFEMINN"  ) ) CALL iom_put( "QFEMINN"  , qfeminn(:,:,:) * tmask(:,:,:) ) ! min Qfe nano
        IF( iom_use( "QFEMINP"  ) ) CALL iom_put( "QFEMINP"  , qfeminp(:,:,:) * tmask(:,:,:) ) ! min Qfe pico
        IF( iom_use( "QFEMIND"  ) ) CALL iom_put( "QFEMIND"  , qfemind(:,:,:) * tmask(:,:,:) ) ! min Qfe diatoms
        IF( iom_use( "QFEMINDZ"  ) ) CALL iom_put( "QFEMINDZ"  , qfemindz(:,:,:) * tmask(:,:,:) ) ! min Qfe diazos
        IF( iom_use( "Fecost_scale") ) CALL iom_put( "Fecost_scale", Fe_scale_nfix(:,:,:) * tmask(:,:,:) ) !Fe cost of nfix Temp dependant scaling
        IF( iom_use( "Pmin_scale") ) CALL iom_put( "Pmin_scale", Pmin_scale_nfix(:,:,:) * tmask(:,:,:) ) !QPmin temperature dependance scaling
      ENDIF
      !
      IF( ln_timing )  CALL timing_stop('p6z_lim')
      !
   END SUBROUTINE p6z_lim


   SUBROUTINE p6z_lim_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_lim_init  ***
      !!
      !! ** Purpose :   Initialization of nutrient limitation parameters
      !!
      !! ** Method  :   Read the namp6zlim and nampisquota namelists and check
      !!      the parameters called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp6zlim
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zlim/ concnno3, concpno3, concdno3, concnnh4, concpnh4, concdnh4,  &
         &                concnfer, concpfer, concdfer, concbfe, concnpo4, concppo4,   &
         &                concdpo4, concbno3, concbnh4, concbpo4, xsizedia, xsizepic,  &
         &                xsizephy, xsizern, xsizerp, xsizerd, xksi1, xksi2, xkdoc,    &
         &                caco3r, oxymin,  concdzno3, concdznh4, concdzpo4,concdzfer,  &
         &                xsizerdz, xsizedz, xkdop, xkdopdz, Facul_lim, kustkaFe,      &
         &                maxFescale, maxPminscale, ln_tiue, ln_tpue         
         !
      NAMELIST/namp6zquota/ qnnmin, qnnmax, qpnmin, qpnmax, qnpmin, qnpmax, qppmin,      &
         &                  qppmax, qndmin, qndmax, qpdmin, qpdmax, qfnmax, qfpmax, qfdmax,  &
         &                  qfnopt, qfpopt, qfdopt, qfdzopt, qndzmin, qndzmax, qpdzmin, qpdzmax, &
         &                  qfdzmax
      !!----------------------------------------------------------------------
      !
      REWIND( numnatp_ref )              ! Namelist namp6zlim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp6zlim, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zlim in reference namelist' )
      !
      REWIND( numnatp_cfg )              ! Namelist namp6zlim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp6zlim, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zlim in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zlim )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp6zlim'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    mean rainratio                           caco3r    = ', caco3r
         WRITE(numout,*) '    NO3 half saturation of nanophyto         concnno3  = ', concnno3
         WRITE(numout,*) '    NO3 half saturation of picophyto         concpno3  = ', concpno3
         WRITE(numout,*) '    NO3 half saturation of diatoms           concdno3  = ', concdno3
         WRITE(numout,*) '    NH4 half saturation for phyto            concnnh4  = ', concnnh4
         WRITE(numout,*) '    NH4 half saturation for pico             concpnh4  = ', concpnh4
         WRITE(numout,*) '    NH4 half saturation for diatoms          concdnh4  = ', concdnh4
         WRITE(numout,*) '    PO4 half saturation for phyto            concnpo4  = ', concnpo4
         WRITE(numout,*) '    PO4 half saturation for pico             concppo4  = ', concppo4
         WRITE(numout,*) '    PO4 half saturation for diatoms          concdpo4  = ', concdpo4
         WRITE(numout,*) '    half saturation constant for Si uptake   xksi1     = ', xksi1
         WRITE(numout,*) '    half saturation constant for Si/C        xksi2     = ', xksi2
         WRITE(numout,*) '    half-sat. of DOC remineralization        xkdoc     = ', xkdoc
         WRITE(numout,*) '    Iron half saturation for nanophyto       concnfer  = ', concnfer
         WRITE(numout,*) '    Iron half saturation for picophyto       concpfer  = ', concpfer
         WRITE(numout,*) '    Iron half saturation for diatoms         concdfer  = ', concdfer
         WRITE(numout,*) '    size ratio for nanophytoplankton         xsizern   = ', xsizern
         WRITE(numout,*) '    size ratio for picophytoplankton         xsizerp   = ', xsizerp
         WRITE(numout,*) '    size ratio for diatoms                   xsizerd   = ', xsizerd
         WRITE(numout,*) '    NO3 half saturation of bacteria          concbno3  = ', concbno3
         WRITE(numout,*) '    NH4 half saturation for bacteria         concbnh4  = ', concbnh4
         WRITE(numout,*) '    Minimum size criteria for diatoms        xsizedia  = ', xsizedia
         WRITE(numout,*) '    Minimum size criteria for picophyto      xsizepic  = ', xsizepic
         WRITE(numout,*) '    Minimum size criteria for nanophyto      xsizephy  = ', xsizephy
         WRITE(numout,*) '    Fe half saturation for bacteria          concbfe   = ', concbfe
         WRITE(numout,*) '    halk saturation constant for anoxia       oxymin   = ', oxymin
         WRITE(numout,*) '    NO3 half saturation of diazotrophs       concdzno3 = ', concdzno3
         WRITE(numout,*) '    NH4 half saturation of diazotrophs       concdznh4 = ', concdznh4
         WRITE(numout,*) '    PO4 half saturation of diazotrophs       concdzpo4 = ', concdzpo4
         WRITE(numout,*) '    Fe half saturation of diazotrophs        concdzfer = ', concdzfer
         WRITE(numout,*) '    Size ratio for diazotrophs               xsizerdz  = ', xsizerdz
         WRITE(numout,*) '    Minimum size criteria for diazotrophs    xsizedz   = ', xsizedz
         WRITE(numout,*) '    half-sat. of DOP uptake                  xkdop     = ', xkdop
         WRITE(numout,*) '    half-sat. of DOP uptake diazos           xkdopdz   = ', xkdopdz
         WRITE(numout,*) '    Diazo Facultative limit                  Facul_lim = ', Facul_lim
         WRITE(numout,*) '    Kustka Fe limitation of diazos           kustkaFe  = ', kustkaFe
         WRITE(numout,*) '    Max Fe cost of nfix can be scaled        maxFescale  = ', maxFescale
         WRITE(numout,*) '    Max QPmin can be scaled                  maxPminscale  = ', maxPminscale
         WRITE(numout,*) '    turn temp dependence of nfix fe cost     ln_tiue   = ', ln_tiue
         WRITE(numout,*) '    turn temp dependence of nfix QPmin       ln_tpue   = ', ln_tpue
      ENDIF

      REWIND( numnatp_ref )              ! Namelist nampislim in reference namelist : Pisces nutrient limitation parameters
      READ  ( numnatp_ref, namp6zquota, IOSTAT = ios, ERR = 903)
903   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nampisquota in reference namelist' )
      !
      REWIND( numnatp_cfg )              ! Namelist nampislim in configuration namelist : Pisces nutrient limitation parameters 
      READ  ( numnatp_cfg, namp6zquota, IOSTAT = ios, ERR = 904 )
904   IF( ios >  0 ) CALL ctl_nam ( ios , 'nampisquota in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zquota )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for nutrient limitations, namp6zquota'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    optimal Fe quota for nano.               qfnopt    = ', qfnopt
         WRITE(numout,*) '    optimal Fe quota for pico.               qfpopt    = ', qfpopt
         WRITE(numout,*) '    Optimal Fe quota for diatoms             qfdopt    = ', qfdopt
         WRITE(numout,*) '    Minimal N quota for nano                 qnnmin    = ', qnnmin
         WRITE(numout,*) '    Maximal N quota for nano                 qnnmax    = ', qnnmax
         WRITE(numout,*) '    Minimal P quota for nano                 qpnmin    = ', qpnmin
         WRITE(numout,*) '    Maximal P quota for nano                 qpnmax    = ', qpnmax
         WRITE(numout,*) '    Minimal N quota for pico                 qnpmin    = ', qnpmin
         WRITE(numout,*) '    Maximal N quota for pico                 qnpmax    = ', qnpmax
         WRITE(numout,*) '    Minimal P quota for pico                 qppmin    = ', qppmin
         WRITE(numout,*) '    Maximal P quota for pico                 qppmax    = ', qppmax
         WRITE(numout,*) '    Minimal N quota for diatoms              qndmin    = ', qndmin
         WRITE(numout,*) '    Maximal N quota for diatoms              qndmax    = ', qndmax
         WRITE(numout,*) '    Minimal P quota for diatoms              qpdmin    = ', qpdmin
         WRITE(numout,*) '    Maximal P quota for diatoms              qpdmax    = ', qpdmax
         WRITE(numout,*) '    Maximal Fe quota for nanophyto.          qfnmax    = ', qfnmax
         WRITE(numout,*) '    Maximal Fe quota for picophyto.          qfpmax    = ', qfpmax
         WRITE(numout,*) '    Maximal Fe quota for diatoms             qfdmax    = ', qfdmax
         WRITE(numout,*) '    optimal Fe quota for diazotrophs         qfdzopt   = ', qfdzopt
         WRITE(numout,*) '    Minimal N quota for diazotrophs          qndzmin   = ', qndzmin
         WRITE(numout,*) '    Maximal N quota for diazotrophs          qndzmax   = ', qndzmax
         WRITE(numout,*) '    Minimal P quota for diazotrophs          qpdzmin   = ', qpdzmin
         WRITE(numout,*) '    Maximal P quota for diazotrophs          qpdzmax   = ', qpdzmax
         WRITE(numout,*) '    Maximal Fe quota for diazotrophs         qfdzmax   = ', qfdzmax
      ENDIF
      !
      ! Metabolic cost of nitrate and ammonium utilisation
      zpsino3  = 2.3 * rno3
      zpsinh4  = 1.8 * rno3
      zpsiuptk = 1.0 / 6.625
      zpsinfix = (2.0 / 0.7 + rtrn) * 2.3 * rno3
      !
      nitrfac (:,:,:) = 0._wp
      !
   END SUBROUTINE p6z_lim_init


   INTEGER FUNCTION p6z_lim_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_lim_alloc  ***
      !!----------------------------------------------------------------------
      USE lib_mpp , ONLY: ctl_stop
      INTEGER ::   ierr(2)        ! Local variables
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      !*  Biological arrays for phytoplankton growth
      ALLOCATE( xpicono3(jpi,jpj,jpk), xpiconh4(jpi,jpj,jpk),       &
         &      xpicopo4(jpi,jpj,jpk), xpicodop(jpi,jpj,jpk),       &
         &      xnanodop(jpi,jpj,jpk), xdiatdop(jpi,jpj,jpk),       &
         &      xpicofer(jpi,jpj,jpk), xlimpfe (jpi,jpj,jpk),       &
         &      fvnuptk (jpi,jpj,jpk), fvduptk (jpi,jpj,jpk),       &
         &      xlimphys(jpi,jpj,jpk), xlimdias(jpi,jpj,jpk),       &
         &      xlimpics(jpi,jpj,jpk), xlimdiazo(jpi,jpj,jpk),      &
         &      fvpuptk (jpi,jpj,jpk), xlimpic (jpi,jpj,jpk),       &
         &      xdiazno3(jpi,jpj,jpk), xdiaznh4(jpi,jpj,jpk),       &
         &      xdiazpo4(jpi,jpj,jpk), xdiazdop(jpi,jpj,jpk),       &
         &      xdiazfer(jpi,jpj,jpk), xlimdiaz(jpi,jpj,jpk),       &
         &      xlimdzfe(jpi,jpj,jpk), fvdzuptk(jpi,jpj,jpk),       &
         &      xlimdzp(jpi,jpj,jpk),  qfeminn(jpi,jpj,jpk),        &
         &      qfeminp(jpi,jpj,jpk),  qfemind(jpi,jpj,jpk),        &
         &      qfemindz(jpi,jpj,jpk), qfenfixdz(jpi,jpj,jpk),     &
         &      xlimnpdz(jpi,jpj,jpk), Fe_scale_nfix(jpi,jpj,jpk),  &
         &      Pmin_scale_nfix(jpi,jpj,jpk),                     STAT=ierr(1) )

         !
      !*  Minimum/maximum quotas of phytoplankton
      ALLOCATE( xqnnmin (jpi,jpj,jpk), xqnnmax(jpi,jpj,jpk),       &
         &      xqpnmin (jpi,jpj,jpk), xqpnmax(jpi,jpj,jpk),       &
         &      xqnpmin (jpi,jpj,jpk), xqnpmax(jpi,jpj,jpk),       &
         &      xqppmin (jpi,jpj,jpk), xqppmax(jpi,jpj,jpk),       &
         &      xqndmin (jpi,jpj,jpk), xqndmax(jpi,jpj,jpk),       &
         &      xqpdmin (jpi,jpj,jpk), xqpdmax(jpi,jpj,jpk),       &
         &      xqndzmin (jpi,jpj,jpk), xqndzmax(jpi,jpj,jpk),     &
         &      xqpdzmin (jpi,jpj,jpk), xqpdzmax(jpi,jpj,jpk),    STAT=ierr(2) )
         !
      p6z_lim_alloc = MAXVAL( ierr )
      !
      IF( p6z_lim_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p6z_lim_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p6z_lim_alloc
   !!======================================================================
END MODULE p6zlim
