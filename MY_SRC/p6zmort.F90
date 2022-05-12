MODULE p6zmort
   !!======================================================================
   !!                         ***  MODULE p6zmort  ***
   !! TOP :   PISCES-QUOTA Compute the mortality terms for phytoplankton
   !!         including an explicit diazotroph PFT
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont)  Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p6z_mort       :   Compute the mortality terms for phytoplankton
   !!   p6z_mort_init  :   Initialize the mortality params for phytoplankton
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  Phytoplankton limitation terms (p4z)
   USE p6zlim          !  Phytoplankton limitation terms (p6z)
   USE prtctl_trc      !  print control for debugging
   USE iom
   IMPLICIT NONE
   PRIVATE

   PUBLIC   p6z_mort           ! Called from p4zbio.F90 
   PUBLIC   p6z_mort_init      ! Called from trcini_pisces.F90 

   !! * Shared module variables
   REAL(wp), PUBLIC :: wchln   !! Quadratic mortality rate of nanophytoplankton
   REAL(wp), PUBLIC :: wchlp   !: Quadratic mortality rate of picophytoplankton
   REAL(wp), PUBLIC :: wchld   !: Quadratic mortality rate of diatoms
   REAL(wp), PUBLIC :: wchldm  !: Maximum quadratic mortality rate of diatoms
   REAL(wp), PUBLIC :: mpratn  !: Linear mortality rate of nanophytoplankton
   REAL(wp), PUBLIC :: mpratp  !: Linear mortality rate of picophytoplankton
   REAL(wp), PUBLIC :: mpratd  !: Linear mortality rate of diatoms
   REAL(wp), PUBLIC :: wchldz  !: Quadratic mortality rate of diazotrophs
   REAL(wp), PUBLIC :: mpratdz !: Linear mortality rate of diazotrophs
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p6zmort.F90 13233 2020-07-02 18:34:16Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE p6z_mort( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort  ***
      !!
      !! ** Purpose :   Calls the different subroutine to compute
      !!                the different phytoplankton mortality terms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt ! ocean time step
      !!---------------------------------------------------------------------

      CALL p6z_mort_nano            ! nanophytoplankton
      CALL p6z_mort_pico            ! picophytoplankton
      CALL p6z_mort_diat            ! diatoms
      CALL p6z_mort_diazo           ! diazotrophs

   END SUBROUTINE p6z_mort


   SUBROUTINE p6z_mort_nano
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_nano  ***
      !!
      !! ** Purpose :   Compute the mortality terms for nanophytoplankton
      !!
      !! ** Method  : - Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp, zprcaca
      REAL(wp) :: ztortp , zrespp , zmortp
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_nano')
      !
      prodcal(:,:,:) = 0.  !: calcite production variable set to zero
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaph = MAX( ( trb(ji,jj,jk,jpphy) - 1e-9 ), 0.e0 )

               ! Quadratic mortality of nano due to aggregation during
               ! blooms (Doney et al. 1996)
               ! -----------------------------------------------------
               zrespp = wchln * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * trb(ji,jj,jk,jpphy)

               ! Phytoplankton linear mortality
               ! A michaelis-menten like term is introduced to avoid 
               ! extinction of nanophyto in highly limited areas
               ! ----------------------------------------------------
               ztortp = mpratn * xstep * zcompaph * trb(ji,jj,jk,jpphy) / ( xkmort + trb(ji,jj,jk,jpphy) )
               zmortp = zrespp + ztortp

               ! Update the arrays TRA which contains the biological sources and sinks

               zfactn  = trb(ji,jj,jk,jpnph)/(trb(ji,jj,jk,jpphy)+rtrn)
               zfactp  = trb(ji,jj,jk,jppph)/(trb(ji,jj,jk,jpphy)+rtrn)
               zfactfe = trb(ji,jj,jk,jpnfe)/(trb(ji,jj,jk,jpphy)+rtrn)
               zfactch = trb(ji,jj,jk,jpnch)/(trb(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zmortp
               tra(ji,jj,jk,jpnph) = tra(ji,jj,jk,jpnph) - zmortp * zfactn
               tra(ji,jj,jk,jppph) = tra(ji,jj,jk,jppph) - zmortp * zfactp
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zmortp * zfactch
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zmortp * zfactfe
               ! Production PIC particles due to mortality
               zprcaca = xfracal(ji,jj,jk) * zmortp
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zmortp * zfactn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + zmortp * zfactp
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp * zfactfe
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nano')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_nano')
      !
   END SUBROUTINE p6z_mort_nano


   SUBROUTINE p6z_mort_pico
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_pico  ***
      !!
      !! ** Purpose :   Compute the mortality terms for picophytoplankton
      !!
      !! ** Method  : - Both quadratic and semilininear terms are used
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp
      REAL(wp) :: ztortp , zrespp , zmortp 
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_pico')
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaph = MAX( ( trb(ji,jj,jk,jppic) - 1e-9 ), 0.e0 )

               ! Quadratic mortality of pico due to aggregation during
               ! blooms (Doney et al. 1996)
               ! -----------------------------------------------------
               zrespp = wchlp * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * trb(ji,jj,jk,jppic)

               ! Phytoplankton linear mortality
               ! A michaelis-menten like term is introduced to avoid 
               ! extinction of picophyto in highly limited areas
               ! ----------------------------------------------------
               ztortp = mpratp * xstep  * zcompaph * trb(ji,jj,jk,jppic) /  ( xkmort + trb(ji,jj,jk,jppic) )
               zmortp = zrespp + ztortp

               !   Update the arrays TRA which contains the biological sources and sinks
               zfactn = trb(ji,jj,jk,jpnpi)/(trb(ji,jj,jk,jppic)+rtrn)
               zfactp = trb(ji,jj,jk,jpppi)/(trb(ji,jj,jk,jppic)+rtrn)
               zfactfe = trb(ji,jj,jk,jppfe)/(trb(ji,jj,jk,jppic)+rtrn)
               zfactch = trb(ji,jj,jk,jppch)/(trb(ji,jj,jk,jppic)+rtrn)
               tra(ji,jj,jk,jppic) = tra(ji,jj,jk,jppic) - zmortp
               tra(ji,jj,jk,jpnpi) = tra(ji,jj,jk,jpnpi) - zmortp * zfactn
               tra(ji,jj,jk,jpppi) = tra(ji,jj,jk,jpppi) - zmortp * zfactp
               tra(ji,jj,jk,jppch) = tra(ji,jj,jk,jppch) - zmortp * zfactch
               tra(ji,jj,jk,jppfe) = tra(ji,jj,jk,jppfe) - zmortp * zfactfe
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortp
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zmortp * zfactn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + zmortp * zfactp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zmortp * zfactfe
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortp
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('pico')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_pico')
      !
   END SUBROUTINE p6z_mort_pico


   SUBROUTINE p6z_mort_diat
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_diat  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diatoms
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER  ::  ji, jj, jk
      REAL(wp) ::  zfactfe,zfactsi,zfactch, zfactn, zfactp, zcompadi
      REAL(wp) ::  zrespp2, ztortp2, zmortp2
      REAL(wp) ::  zlim2, zlim1
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_diat')
      !

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               zcompadi = MAX( ( trb(ji,jj,jk,jpdia) - 1E-9), 0. )

               !   Aggregation term for diatoms is increased in case of nutrient
               !   stress as observed in reality. The stressed cells become more
               !   sticky and coagulate to sink quickly out of the euphotic zone
               !   -------------------------------------------------------------
               !  Phytoplankton squared mortality
               !  -------------------------------
               zlim2   = xlimdia(ji,jj,jk) * xlimdia(ji,jj,jk)
               zlim1   = 0.25 * ( 1. - zlim2 ) / ( 0.25 + zlim2 ) 
               zrespp2 = 1.e6 * xstep * (  wchld + wchldm * zlim1 ) * xdiss(ji,jj,jk) * zcompadi * trb(ji,jj,jk,jpdia)

               ! Phytoplankton linear mortality
               ! A michaelis-menten like term is introduced to avoid 
               ! extinction of diatoms in highly limited areas
               !  ---------------------------------------------------
               ztortp2 = mpratd * xstep  * zcompadi * trb(ji,jj,jk,jpdia) /  ( xkmort + trb(ji,jj,jk,jpdia) )
               zmortp2 = zrespp2 + ztortp2

               !   Update the arrays tra which contains the biological sources and sinks
               !   ---------------------------------------------------------------------
               zfactn  = trb(ji,jj,jk,jpndi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactp  = trb(ji,jj,jk,jppdi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactch = trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactfe = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn )
               zfactsi = trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zmortp2 
               tra(ji,jj,jk,jpndi) = tra(ji,jj,jk,jpndi) - zmortp2 * zfactn
               tra(ji,jj,jk,jppdi) = tra(ji,jj,jk,jppdi) - zmortp2 * zfactp
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zmortp2 * zfactch
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zmortp2 * zfactfe
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zmortp2 * zfactsi
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zmortp2 * zfactsi
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zrespp2 
               tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zrespp2 * zfactn
               tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + zrespp2 * zfactp
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zrespp2 * zfactfe
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + ztortp2
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + ztortp2 * zfactn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + ztortp2 * zfactp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ztortp2 * zfactfe
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + ztortp2
               prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zrespp2
            END DO
         END DO
      END DO
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diat')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_diat')
      !
   END SUBROUTINE p6z_mort_diat

   SUBROUTINE p6z_mort_diazo
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p6z_mort_diazo  ***
      !!
      !! ** Purpose :   Compute the mortality terms for diazotrophs
      !!
      !! ** Method  : - Both quadratic and simili linear mortality terms
      !!---------------------------------------------------------------------
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompaph
      REAL(wp) :: zfactfe, zfactch, zfactn, zfactp, zprcaca
      REAL(wp) :: ztortp , zrespp , zmortp
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: mortdz
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p6z_mort_diazo')
      !
      mortdz(:,:,:) = 0.
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaph = MAX( ( trb(ji,jj,jk,jpcdz) - 1e-9 ), 0.e0 )

               ! Quadratic mortality of nano due to aggregation during
               ! blooms (Doney et al. 1996)
               ! -----------------------------------------------------
               zrespp = wchldz * 1.e6 * xstep * xdiss(ji,jj,jk) * zcompaph * trb(ji,jj,jk,jpcdz)

               ! Phytoplankton linear mortality
               ! A michaelis-menten like term is introduced to avoid 
               ! extinction of nanophyto in highly limited areas
               ! ----------------------------------------------------
               ztortp = mpratdz * xstep * zcompaph * trb(ji,jj,jk,jpcdz) / ( xkmort + trb(ji,jj,jk,jpcdz) )
               zmortp = zrespp + ztortp

               ! Update the arrays TRA which contains the biological sources and
               ! sinks

               zfactn  = trb(ji,jj,jk,jpndz)/(trb(ji,jj,jk,jpcdz)+rtrn)
               zfactp  = trb(ji,jj,jk,jppdz)/(trb(ji,jj,jk,jpcdz)+rtrn)
               zfactfe = trb(ji,jj,jk,jpfed)/(trb(ji,jj,jk,jpcdz)+rtrn)
               zfactch = trb(ji,jj,jk,jpchd)/(trb(ji,jj,jk,jpcdz)+rtrn)
               tra(ji,jj,jk,jpcdz) = tra(ji,jj,jk,jpcdz) - zmortp
               tra(ji,jj,jk,jpndz) = tra(ji,jj,jk,jpndz) - zmortp * zfactn
               tra(ji,jj,jk,jppdz) = tra(ji,jj,jk,jppdz) - zmortp * zfactp
               tra(ji,jj,jk,jpchd) = tra(ji,jj,jk,jpchd) - zmortp * zfactch
               tra(ji,jj,jk,jpfed) = tra(ji,jj,jk,jpfed) - zmortp * zfactfe
               !
               mortdz(ji,jj,jk) = zmortp * zfactn
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zrespp
               tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zrespp * zfactn
               tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + zrespp * zfactp
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zrespp * zfactfe
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + ztortp
               tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + ztortp * zfactn
               tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + ztortp * zfactp
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ztortp * zfactfe
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + ztortp
               prodgoc(ji,jj,jk)   = prodgoc(ji,jj,jk) + zrespp
            END DO
         END DO
      END DO
      !
       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('diazo')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p6z_mort_diazo')
      !
   END SUBROUTINE p6z_mort_diazo

   SUBROUTINE p6z_mort_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p6z_mort_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton mortality parameters
      !!
      !! ** Method  :   Read the namp6zmort namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist namp6zmort
      !!
      !!----------------------------------------------------------------------
      INTEGER :: ios   ! Local integer output status for namelist read
      !!
      NAMELIST/namp6zmort/ wchln, wchlp, wchld, wchldm, mpratn, mpratp, mpratd, wchldz, mpratdz
      !!----------------------------------------------------------------------

      REWIND( numnatp_ref )              ! Namelist namp6zmort in reference namelist : Pisces phytoplankton
      READ  ( numnatp_ref, namp6zmort, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namp6zmort in reference namelist' )

      REWIND( numnatp_cfg )              ! Namelist namp6zmort in configuration namelist : Pisces phytoplankton
      READ  ( numnatp_cfg, namp6zmort, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 ) CALL ctl_nam ( ios , 'namp6zmort in configuration namelist' )
      IF(lwm) WRITE ( numonp, namp6zmort )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) ' '
         WRITE(numout,*) ' Namelist parameters for phytoplankton mortality, namp6zmort'
         WRITE(numout,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '    quadratic mortality of phytoplankton      wchln     =', wchln
         WRITE(numout,*) '    quadratic mortality of picophyto.         wchlp     =', wchlp
         WRITE(numout,*) '    quadratic mortality of diatoms            wchld     =', wchld
         WRITE(numout,*) '    Additional quadratic mortality of diatoms wchldm    =', wchldm
         WRITE(numout,*) '    nanophyto. mortality rate                 mpratn    =', mpratn
         WRITE(numout,*) '    picophyto. mortality rate                 mpratp    =', mpratp
         WRITE(numout,*) '    Diatoms mortality rate                    mpratd    =', mpratd
         WRITE(numout,*) '    quadratic mortality of diazotrophs        wchldz    =', wchldz
         WRITE(numout,*) '    diazotrophs mortality rate                mpratdz   =', mpratdz
      ENDIF

   END SUBROUTINE p6z_mort_init

   !!======================================================================
END MODULE p6zmort
