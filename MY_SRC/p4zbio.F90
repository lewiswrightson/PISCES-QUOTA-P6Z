MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES biogeochemical model
   !!         This module is for both PISCES and PISCES-QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.6  ! 2015 (O. Aumont) PISCES-QUOTA
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p5zlim          !  Co-limitations of differents nutrients
   USE p5zprod         !  Growth rate of the 2 phyto groups
   USE p5zmort         !  Mortality terms for phytoplankton
   USE p5zmicro        !  Sources and sinks of microzooplankton
   USE p5zmeso         !  Sources and sinks of mesozooplankton
   !! Added compartments for Explicit diazotroph model
   USE p6zlim          !  Co-limitations of differents nutrients
   USE p6zprod         !  Growth rate of the 4 phyto groups
   USE p6zmort         !  Mortality terms for phytoplankton
   USE p6zmicro        !  Sources and sinks of microzooplankton
   USE p6zmeso         !  Sources and sinks of mesozooplankton
   !!
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zpoc          !  Remineralization of organic particles
   USE p4zagg          !  Aggregation of particles
   USE p4zfechem
   USE p4zligand       !  Prognostic ligand model
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zbio.F90 13233 2020-07-02 18:34:16Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!                different interactions between the different compartments
      !!                of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER             :: ji, jj, jk, jn
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_bio')

      ! ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      ! OF PHYTOPLANKTON AND DETRITUS. Shear rate is supposed to equal 1
      ! in the mixed layer and 0.1 below the mixed layer.
      xdiss(:,:,:) = 1.
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gdepw_n(ji,jj,jk+1) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

      CALL p4z_opt     ( kt, knt )     ! Optics: PAR in the water column
      CALL p4z_sink    ( kt, knt )     ! vertical flux of particulate organic matter
      CALL p4z_fechem  ( kt, knt )     ! Iron chemistry/scavenging
      !
      IF( ln_p4z ) THEN  ! PISCES standard
         ! Phytoplankton only sources/sinks terms
         CALL p4z_lim  ( kt, knt )     ! co-limitations by the various nutrients
         CALL p4z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
         !                             ! (for each element : C, Si, Fe, Chl )
         CALL p4z_mort ( kt      )     ! phytoplankton mortality
         ! zooplankton sources/sinks routines 
         CALL p4z_micro( kt, knt )     ! microzooplankton
         CALL p4z_meso ( kt, knt )     ! mesozooplankton
      ENDIF
      IF( ln_p5z ) THEN  ! PISCES-QUOTA implicit diazotrophy
         ! Phytoplankton only sources/sinks terms
         CALL p5z_lim  ( kt, knt )     ! co-limitations by the various nutrients
         CALL p5z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
         !                             ! (for each element : C, N, P, Si, Fe, Chl )
         CALL p5z_mort ( kt      )     ! phytoplankton mortality
         ! zooplankton sources/sinks routines 
         CALL p5z_micro( kt, knt )     ! microzooplankton
         CALL p5z_meso ( kt, knt )     ! mesozooplankton
      ENDIF
      IF( ln_p6z ) THEN  ! PISCES-QUOTA Explicit diazotrophy
         ! Phytoplankton only sources/sinks terms
         CALL p6z_lim  ( kt, knt )     ! co-limitations by the various nutrients
         CALL p6z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
         !                             ! (for each element : C, N, P, Si, Fe, Chl )
         CALL p6z_mort ( kt      )     ! phytoplankton mortality
         ! zooplankton sources/sinks routines 
         CALL p6z_micro( kt, knt )     ! microzooplankton
         CALL p6z_meso ( kt, knt )     ! mesozooplankton
      ENDIF 
      !
      CALL p4z_agg     ( kt, knt )     ! Aggregation of particles
      CALL p4z_rem     ( kt, knt )     ! remineralization terms of organic matter
      CALL p4z_poc     ( kt, knt )     ! Remineralization of organic particles
      !
      ! Ligand production. ln_ligand should be set .true. to activate
      IF( ln_ligand )  &
      & CALL p4z_ligand( kt, knt )

      ! Update of the size of the different phytoplankton groups
      sized(:,:,:) = sizeda(:,:,:)
      sizen(:,:,:) = sizena(:,:,:)
      ! Update the size of Picophytoplankton for PISCES QUOTA
      IF (ln_p5z .OR. ln_p6z) THEN
         sizep(:,:,:) = sizepa(:,:,:)
      ENDIF
      ! Update the size of the explicit diazotroph PFT in PISCES QUOTA
      IF ( ln_p6z ) THEN
         sizedz(:,:,:) = sizedza(:,:,:)
      ENDIF
      !                                                             !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_bio')
      !
   END SUBROUTINE p4z_bio

   !!======================================================================
END MODULE p4zbio
