! $Id: wetscav_adj_mod.f,v 1.4 2012/04/20 19:19:31 nicolas Exp $
      MODULE WETSCAV_ADJ_MOD
!
!******************************************************************************
!  Module WETSCAV_MOD contains arrays for used in the wet scavenging of
!  tracer in cloud updrafts, rainout, and washout. (bmy, 2/28/00, 3/5/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NSOLMAX (INTEGER) : Max # of soluble tracers       [unitless]
!  (2 ) NSOL    (INTEGER) : Actual # of soluble tracers    [unitless]
!  (3 ) IDWETD  (INTEGER) : Index array for WETDEP routine [unitless]
!  (4 ) Vud     (REAL*8 ) : Array for updraft velocity     [m/s]
!  (5 ) CLDLIQ  (REAL*8 ) : Array for cloud liquid water   [cm3 H2O/cm3 air]
!  (6 ) CLDICE  (REAL*8 ) : Array for cloud ice content    [cm3 ice/cm3 air]
!  (7 ) C_H2O   (REAL*8 ) : Array for Mixing ratio of ,
!                            water, computed from Eice(T)  [v/v]
!  (8 ) PDOWN   (REAL*8 ) : Precip thru bottom of grid box [cm3 H2O/cm2 area/s]
!  (9 ) QQ      (REAL*8 ) : Rate of new precip formation   [cm3 H2O/cm3 air/s]
!  (10) EPSILON (REAL*8 ) : A very small positive number   [unitless]
!  (11) H2O2s   (REAL*8 ) : Array to save H2O2 for wetdep  [v/v]
!  (12) SO2s    (REAL*8 ) : Array to save SO2 for wetdep   [v/v]
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_QQ           : Constructs the QQ field (precipitable water)
!  (2 ) E_ICE             : Computes saturation vapor pressure for ice
!  (3 ) COMPUTE_L2G       : Computes the ratio [v/v liquid] / [v/v gas]
!  (4 ) COMPUTE_F         : Computes fraction of tracer lost in cloud updrafts
!  (5 ) F_AEROSOL         : Computes fraction of tracer scavenged in updrafts
!  (6 ) GET_ISOL          : Returns correct index for ND38 diagnostic
!  (7 ) RAINOUT           : Computes fraction of soluble tracer lost to rainout
!  (8 ) GET_RAINFRAC      : Computes rainout fraction -- called by RAINOUT
!  (9 ) WASHOUT           : Computes fraction of soluble tracer lost to washout
!  (10) WASHFRAC_AEROSOL  : Computes fraction of aerosol lost to washout
!  (11) WASHFRAC_LIQ_GAS  : Computes fraction of soluble gases lost to washout
!  (12) WETDEP            : Driver routine for computing wet deposition losses
!  (13) LS_K_RAIN         : Computes K_RAIN (for LS precipitation)
!  (14) LS_F_PRIME        : Computes F_PRIME (for LS precipitation)
!  (15) CONV_F_PRIME      : Computes F_PRIME (for convective precipitation)
!  (16) SAFETY            : Stops WETDEP w/ error msg if negative tracer found
!  (17) WETDEPID          : Initalizes the IDWETD array for routine WETDEP
!  (18) GET_WETDEP_NMAX   : Returns max # of soluble tracers per simulation
!  (19) GET_WETDEP_NSOL   : Returns actual # of soluble tracers per simulation
!  (20) GET_WETDEP_IDWETD : Returns CTM tracer # of for a given wetdep species
!  (21) INIT_WETSCAV      : Initializes fields used for computing wetdep losses
!  (22) CLEANUP_WETSCAV   : Deallocates all allocatable module arrays
!
!  GEOS-CHEM modules referenced by wetscav_mod.f
!  ============================================================================
!  (1 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) logical_mod.f  : Module containing GEOS-CHEM logical switches
!  (5 ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!  (6 ) tracer_mod.f   : Module containing GEOS-CHEM tracer array STT etc.
!  (7 ) tracerid_mod.f : Module containing pointers to tracers and emissions
!
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.J. Jacob, I. Bey and R.M. Yantosca, "Constraints from 210Pb
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields",
!        JGR, Vol 106, pp 12109-12128, 2001.
!  (2 ) D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet
!        deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling
!        Group, March 2000.
!  (3 ) Chin, M., D.J. Jacob, G.M. Gardner, M.S. Foreman-Fowler, and P.A.
!        Spiro, "A global three-dimensional model of tropospheric sulfate",
!        J. Geophys. Res., 101, 18667-18690, 1996.
!  (4 ) Balkanski, Y  D.J. Jacob, G.M. Gardner, W.C. Graustein, and K.K.
!        Turekian, "Transport and Residence Times of Tropospheric Aerosols
!        from a Global Three-Dimensional Simulation of 210Pb", JGR, Vol 98,
!        (D11) pp 20573-20586, 1993.
!  (5 ) Giorgi, F, & W.L. Chaimedes, "Rainout Lifetimes of Highly Soluble
!        Aerosols and Gases as Inferred from Simulations With a General
!        Circulation Model", JGR, Vol 86 (D13) pp 14367-14376, 1986.
!
!  NOTES:
!  (1 ) Now trap allocation errors with routine ALLOC_ERR. (bmy, 7/11/00)
!  (2 ) Moved routine MAKE_QQ here from "dao_mod.f" (bmy, 10/12/00)
!  (3 ) Reordered arguments in INIT_PRECIP (bmy, 10/12/00)
!  (4 ) Updated comments (bmy, 9/4/01)
!  (5 ) Bug fix in MAKE_QQ: BXHEIGHT is sized IIPAR,JJPAR,LLPAR (bmy, 10/4/01)
!  (6 ) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!  (7 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (8 ) Now zero allocatable arrays (bmy, 8/5/02)
!  (9 ) Bug fix: ND39 diagnostic now closes the budget.  Also bundled several
!        standalone routines into this module.  Now references F90 module
!        "tracerid_mod.f".  Also set NSOLMAX=10 since we now have sulfate
!        tracers for wetdep.   Now prevent out-of-bounds errors in routine
!        WETDEP.  Added GET_WETDEP_NMAX function to return max # of soluble
!        tracers for allocating diagnostic arrays.  Added functions
!        GET_WETDEP_NSOL and GET_WETDEP_IDWETD.  Now init H2O2s and SO2s
!        to the initial H2O2 and SO2 from STT.  Updated comments.
!        (qli, bmy, 1/14/03)
!  (10) Improvements for SO2/SO4 scavenging (rjp, bmy, 3/23/03)
!  (11) Now references "time_mod.f".  Added driver routine DO_WETDEP to
!        remove cumbersome calling sequence from MAIN program.  Also declared
!        WETDEP and MAKE_QQ PRIVATE to this module. (bmy, 3/27/03)
!  (11) Add parallelization to routine WETDEP (bmy, 3/17/04)
!  (12) Added carbon and dust aerosol tracers (rjp, tdf, bmy, 4/5/04)
!  (13) Added seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (14) Added secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (15) Now references "logical_mod.f" and "tracer_mod.f".  Now move all
!        internal routines to the module and pass arguments explicitly in
!        order to facilitate parallelization on the Altix. (bmy, 7/20/04)
!  (16) Updated for mercury aerosol tracers (eck, bmy, 12/9/04)
!  (17) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also now pass Hg2 wetdep loss
!        to "ocean_mercury_mod.f". (cas, sas, bmy, 1/20/05)
!  (18) Bug fix to avoid numerical blowup in WETDEP.  Now use analytical
!        function for E_ICE(T). (bmy, 3/7/05)
!  (19) Added SO4s, NITs.  Increased NSOLMAX to 31.  Also block out
!        parallel loop in WETDEP for SGI MIPS compiler. (bec, bmy, 5/5/05)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Bug fixes: do not over-deplete H2O2s.  Also include updates for
!        tagged Hg simulation. (dkh, rjp, eck, cdh, bmy, 1/6/06)
!  (22) Now wet deposit SOG4, SOA4. Remove unnecessary variables in WETDEP.
!        (dkh, bmy, 5/18/06)
!  (23) Bug fixes in COMPUTE_F (bmy, 7/26/06)
!  (24) Resize DSTT array in WETDEP to save memory.  Added fixes for GEOS-5
!        wet deposition per Hongyu Liu's suggestions. (bmy, 3/5/08)
!  (25) Add wet scavenging of GLYX, MGLY, GLYC, SOAG, SOAM (tmf, 1/7/09)
!  (26) Effective Henry's law constant and coefficient from
!       Sander, R, 1999, Compilation of Henry's Law Constants for
!          Inorganic and Organic Species of Potential Importance in
!          Environmental Chemistry.
!          http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!       (tmf, 1/7/09)
!  26 Oct 2015 - M. Sulprizio- Add support for MERRA2 meteorology
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      !PUBLIC :: H2O2s
      !PUBLIC :: SO2s

      PUBLIC :: H2O2s_ADJ
      PUBLIC :: SO2s_ADJ

      ! ... and these routines
      PUBLIC :: CLEANUP_WETSCAV_ADJ
      !PUBLIC :: COMPUTE_F
      PUBLIC :: DO_WETDEP_ADJ
      !PUBLIC :: GET_WETDEP_IDWETD
      !PUBLIC :: GET_WETDEP_NMAX
      !PUBLIC :: GET_WETDEP_NSOL
      PUBLIC :: INIT_WETSCAV_ADJ
      PUBLIC :: ADJ_INIT_WETSCAV
      PUBLIC :: WETSCAV_ADJ_FORCE
      !PUBLIC :: WETDEPID

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
!      INTEGER, PARAMETER   :: NSOLMAX = 38
!      REAL*8,  PARAMETER   :: EPSILON = 1d-32
!
!      ! Scalars
!      INTEGER              :: NSOL
!
!      ! Arrays
!      INTEGER              :: IDWETD(NSOLMAX)
!      REAL*8,  ALLOCATABLE :: Vud(:,:)
!      REAL*8,  ALLOCATABLE :: C_H2O(:,:,:)
!      REAL*8,  ALLOCATABLE :: CLDLIQ(:,:,:)
!      REAL*8,  ALLOCATABLE :: CLDICE(:,:,:)
!      REAL*8,  ALLOCATABLE :: PDOWN(:,:,:)
!      REAL*8,  ALLOCATABLE :: QQ(:,:,:)
!      REAL*8,  ALLOCATABLE :: H2O2s(:,:,:)
!      REAL*8,  ALLOCATABLE :: SO2s(:,:,:)

      ! adjoint variables
      REAL*8,  ALLOCATABLE :: H2O2s_ADJ(:,:,:)
      REAL*8,  ALLOCATABLE :: SO2s_ADJ(:,:,:)

      !>>>
      ! Now include adjoint of F (dkh, 10/03/08)
      REAL*8,  ALLOCATABLE :: F_ADJ(:,:,:)
      !<<<

      ! wetdep adj (fp, dkh, 03/04/13)
      REAL*8, ALLOCATABLE :: BOX_DEP(:,:,:,:)
      REAL*8, ALLOCATABLE :: LOWER_DEP(:,:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_WETDEP_ADJ
!
!******************************************************************************
!  Subroutine DO_WETDEP_ADJ is a driver for the adjoint of the wet deposition
!  code. (dkh, 02/??/05) !
!
!  It is based on the subroutine DO_WETDEP, a driver for the wet deposition code,
!  called from the MAIN program. (bmy, 3/27/03, 3/5/08)
!
!  NOTES:
!  (1 ) BUG FIX.  Do adjoint wetdep of convective precip BEFORE adjoint wetdep
!        of stratiform precip.  (dkh, 10/23/05)
!  (2 )  Added support for full chemistry.  Now call ADJ_SO2_WETDEP to do the
!         adjoint for SO2.  Now call forward wetdep routines first
!         in order to recompute variables (H2O2s, SO2s,STT(SO4)) needed for
!         ADJ_SO2_WETDEP. (dkh, 10/23/05)
!  (3 ) Updated to GCv8 adjoint (dkh, 09/28/09)
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE LOGICAL_MOD,        ONLY : LPRT
      ! adj_group debug (dkh, 06/08/09)
      USE ADJ_ARRAYS_MOD,     ONLY : IFD, JFD, LFD, NFD
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : NSPAN
      USE CHECKPT_MOD,        ONLY : CHK_STT_BEFCHEM
      USE LOGICAL_ADJ_MOD,    ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD,    ONLY : LMAX_OBS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_WDEP_LS
      USE WETSCAV_MOD,        ONLY : H2O2s, SO2s
      USE WETSCAV_MOD,        ONLY : MAKE_QQ
      USE WETSCAV_MOD,        ONLY : RESTORE
      USE WETSCAV_MOD,        ONLY : RECALC_SOX_WETDEP
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTSO4
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TRACER_MOD,         ONLY : STT
      USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM
#if defined(EANET_OBS) || defined(EMEP_OBS) || defined(NADP_OBS) || defined(IDAF_OBS)
      USE N_DEPOSITION_OBS_MOD,  ONLY : NDEP_FORCING
#endif


#     include "CMN_SIZE"           ! Size parameters
#     include "define_adj.h"       ! Obs operators

      INTEGER                     :: I, J, L
      REAL*8,  SAVE               :: OBS_COUNT  = 0


      !==================================================================
      ! DO_WETDEP_ADJ begins here!
      !==================================================================

      !=================================================================
      ! Adjoint of wetdep for all species other than SO2
      !  ( need to do adjoint of wetdep of SO4 before adjoint of SO2 )
      !=================================================================

#if   !defined( GEOS_5 ) && !defined( GEOS_FP ) && !defined( MERRA2 )

      ! Adjoint of wetdep by convective precip
      CALL MAKE_QQ( .FALSE. )
      CALL WETDEP_ADJ(  .FALSE. )

      IF ( ITS_A_FULLCHEM_SIM() ) THEN
         ! Adjoint of wetdep by convective precip for SO2
         ! Restore initial values of H2O2s, SO2s, STT(SO4), STT(SO2)
         CALL RESTORE

         IF ( LPRINTFD ) THEN
            WRITE(6,*) ' WETD CHK variables before RECALC_SOX '
            print*, ' H2O2s(FD) = ', H2O2s(IFD,JFD,LFD)
            print*, ' SO2s(FD) = ', SO2s(IFD,JFD,LFD)
            print*, ' SO4(FD) = ', STT(IFD,JFD,LFD,IDTSO4)
         ENDIF

         ! Recompute intemediate values
         CALL MAKE_QQ( .TRUE. )
         ! Just do fwd LS precip for SO2 and SO4
         CALL RECALC_SOX_WETDEP(  .TRUE. )

         IF ( LPRINTFD ) THEN
            WRITE(6,*) ' WETD CHK variables before ADJ_SO2_WETDEP(F) '
            print*, ' H2O2s(FD) = ', H2O2s(IFD,JFD,LFD)
            print*, ' SO2s(FD) = ', SO2s(IFD,JFD,LFD)
            print*, ' SO4(FD) = ', STT(IFD,JFD,LFD,IDTSO4)
            print*, ' STT_ADJ(FD) = ', STT_ADJ(IFD,JFD,LFD,NFD)
            print*, ' STT(FD) = ', STT(IFD,JFD,LFD,NFD)
         ENDIF

         CALL MAKE_QQ( .FALSE. )
         CALL ADJ_SO2_WETDEP(  .FALSE. )

      ENDIF

      IF ( LPRT ) CALL DEBUG_MSG('### DO_WETDEP_ADJ: after conv wetdep')

#endif

      ! Wetdep by large-scale (stratiform) precip
      CALL MAKE_QQ( .TRUE. )
      CALL WETDEP_ADJ(  .TRUE. )

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Restore initial values of H2O2s, SO2s, STT(SO4)
         CALL RESTORE

         IF ( LPRINTFD ) THEN
            WRITE(6,*) ' WETD CHK variables before ADJ_SO2_WETDEP(T) '
            print*, ' H2O2s(FD) = ', H2O2s(IFD,JFD,LFD)
            print*, ' SO2s(FD) = ', SO2s(IFD,JFD,LFD)
            print*, ' SO4(FD) = ', STT(IFD,JFD,LFD,IDTSO4)
            print*, ' SO2(FD) = ', STT(IFD,JFD,LFD,IDTSO2)
         ENDIF

         CALL ADJ_SO2_WETDEP(  .TRUE. )

         ! Reset SO2 and SO4 in STT to values before chemistry
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            STT(I,J,L,IDTSO4) = CHK_STT_BEFCHEM(I,J,L,IDTSO4)
            STT(I,J,L,IDTSO2) = CHK_STT_BEFCHEM(I,J,L,IDTSO2)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_WETDEP_ADJ: after LS wetdep' )

      ! Apply adjoint forcing now because WETSCAV_ADJ_FORCE and NDEP_FORCING
      ! use arrays calculated in WETDEP_ADJ currently stored in memory
      IF ( LMAX_OBS ) THEN
         OBS_COUNT = OBS_COUNT
     &             + REAL(GET_TS_DYN(),8) / REAL(GET_TS_CHEM(),8)

         IF ( OBS_COUNT > NSPAN ) RETURN
      ENDIF

#if defined(EANET_OBS) || defined(EMEP_OBS) || defined(NADP_OBS) || defined(IDAF_OBS)
      CALL NDEP_FORCING
      CALL WETSCAV_ADJ_FORCE
#else
      IF ( LADJ_WDEP_LS ) THEN
         CALL WETSCAV_ADJ_FORCE
      ENDIF
#endif

      ! Return to calling program
      END SUBROUTINE DO_WETDEP_ADJ

!------------------------------------------------------------------------------
!
!      SUBROUTINE MAKE_QQ( LS )
!!
!!*****************************************************************************
!!  Subroutine MAKE_QQ computes the large-scale or convective precipitation
!!  fields for use with wetdep.f. (hyl, bmy, 2/29/00, 11/8/02)
!!
!!  Arguments as Input:
!!  ===========================================================================
!!  (1 ) LS       : = T for Large-scale precip, =F otherwise
!!
!!  DAO met fields from "dao_mod.f:"
!!  ===========================================================================
!!  (1 ) AIRDEN   : Density of air in grid box (I,J,L) [kg air/m^3]
!!  (2 ) BXHEIGHT : Height of grid box (I,J,L) in [m]
!!  (3 ) MOISTQ   : DAO field for change in specific
!!                  humidity due to moist processes    [kg H2O/kg air/s]
!!  (4 ) PREACC   : DAO total accumulated precipitaton [mm/day]
!!  (5 ) PRECON   : DAO convective precipitation       [mm/day]
!!
!!  References (see above for full citations):
!!  ===========================================================================
!!  (1 ) Liu et al, 2000
!!  (2 ) Jacob et al, 2000
!!
!!  NOTES:
!!  (1 ) Now we partition MOISTQ into large-scale and convective parts, using
!!        total precipitation PREACC and convective precipitation PRECON (both
!!        are vertical integral amounts). The precipitation field at altitudes
!!        (PDOWN) is also made (hyl, djj, 10/17/98).
!!  (2 ) MAKE_QQ is written in Fixed-Form Fortran 90. (bmy, 4/2/99)!
!!  (3 ) AIRDEN, MOISTQ, QQ, and PDOWN are dimensioned (LLPAR,IIPAR,JJPAR)
!!       in order to maximize loop efficiency when processing an (I,J)
!!       column layer by layer. (bmy, 3/14/00)
!!  (4 ) MOISTQ is originally [g H2O/kg air/day], and is converted in
!!        READ_A6 to [kg H2O/kg air/s]. (bmy, 3/14/00)
!!  (5 ) Now reference PREACC, PRECON from "dao_mod.f" instead of from
!!        common block header file "CMN_PRECIP" (bmy, 6/26/00)
!!  (6 ) Now pass BXHEIGHT as an argument.  Also added to "dao_mod.f".
!!        (bmy, 6/26/00)
!!  (7 ) Moved from "dao_mod.f" to "wetscav_mod.f".  Also made PREACC
!!        and PRECON into arguments. (bmy, 10/12/00)
!!  (8 ) Updated comments (bmy, 9/4/01)
!!  (9 ) BXHEIGHT is now sized (IIPAR,JJPAR,LLPAR) (bmy, 10/4/01)
!!  (10) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!!  (11) Now reference met field arrays directly from "dao_mod.f" (bmy, 11/8/02)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,   ONLY : AIRDEN, BXHEIGHT, MOISTQ, PREACC, PRECON
!      USE ERROR_MOD, ONLY : ALLOC_ERR
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      LOGICAL, INTENT(IN)  :: LS
!
!      ! Local variables
!      INTEGER              :: I, J, L, AS
!      REAL*8               :: PTEMP, FRAC
!      LOGICAL              :: FIRST = .TRUE.
!
!      !=================================================================
!      ! MAKE_QQ begins here!
!      !=================================================================
!      IF ( FIRST ) THEN
!
!         ! Allocate PDOWN on first call
!         ALLOCATE( PDOWN( LLPAR, IIPAR, JJPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDOWN' )
!         PDOWN = 0d0
!
!         ! Allocate QQ on first call
!         ALLOCATE( QQ( LLPAR, IIPAR, JJPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'QQ' )
!         QQ = 0d0
!
!         ! Reset flag
!         FIRST = .FALSE.
!      ENDIF
!
!      !=================================================================
!      ! Loop over surface grid boxes
!      !=================================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, FRAC, L, PTEMP )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         !==============================================================
!         ! If there is total precipitation in the (I,J) column, then:
!         !
!         ! (1) Compute FRAC, the large scale fraction (if LS = .TRUE.)
!         !     or convective fraction (if LS = .FALSE.) total
!         !     precipitation.  FRAC is computed from PREACC and PRECON.
!         !
!         ! (2) Compute QQ, the rate of formation of precipitation
!         !     [cm3 H2O/cm3 air/s].  From MOISTQ [kg H2O/kg air/s],
!         !     the unit conversion is:
!         !
!         !     kg H2O   |   m^3 H2O   | AIRDEN kg air         m^3 H2O
!         !  ------------+-------------+--------------- ==> -------------
!         !   kg air * s | 1000 kg H2O |    m^3 air          m^3 air * s
!         !
!         ! and
!         !
!         !         m^3 H2O                         cm^3 H2O
!         !      -------------  is equivalent to  --------------
!         !       m^3 air * s                      cm^3 air * s!
!         !
!         ! since the same conversion factor (10^6 cm^3/m^3) is in both
!         ! the numerator and the denominator.
!         !
!         ! Therefore, the equation for QQ is:
!         !
!         !   QQ(L,I,J) = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1000.0
!         !
!         ! (3) Compute PDOWN, the column precipitation
!         !     [cm3 H2O/cm2 air/s], by multiplying QQ(L,I,J) by
!         !     BXHEIGHT(I,J,L) * 100 cm.
!         !
!         ! (4) The reason why we do not force PTEMP to be positive is
!         !     that PREACC is the integral of the MOISTQ field.  MOISTQ
!         !     contains both negative (evap) and positive (precip)
!         !     values.  If we forced PTEMP to be positive, then we would
!         !     be adding extra precipitation to PDOWN (hyl, bmy, 3/6/99).
!         !==============================================================
!         IF ( PREACC(I,J) > 0d0 ) THEN
!
!            ! Large scale or convective fraction of precipitation
!            IF ( LS ) THEN
!               FRAC = ( PREACC(I,J) - PRECON(I,J) ) / PREACC(I,J)
!            ELSE
!               FRAC = PRECON(I,J) / PREACC(I,J)
!            ENDIF
!
!            ! Start at the top of the atmosphere
!            L = LLPAR
!
!            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
!            QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
!            PTEMP        = QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
!            PDOWN(L,I,J) = PTEMP
!
!            ! PDOWN cannot be negative
!            IF ( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0
!
!            ! Loop down from LLPAR to the surface
!            DO L = LLPAR-1, 1, -1
!
!               ! Compute QQ and PDOWN.  Keep PTEMP for the next level.
!               QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
!               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
!               PTEMP        = PDOWN(L,I,J)
!
!               ! PDOWN cannot be negative
!               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
!            ENDDO
!
!         !==============================================================
!         ! If there is no precipitation reaching the surface in the
!         ! (I,J) column, then assume any precipitation at altitude to
!         ! be large-scale.
!         !
!         ! (1) Assume the large scale fraction = 1d0,
!         !                convective fraction  = 0d0
!         ! (2) Compute QQ as described above
!         ! (3) Compute PDOWN as described above
!         !==============================================================
!         ELSE
!
!            ! Assume large-scale precipitation!
!            IF ( LS ) THEN
!               FRAC = 1d0
!            ELSE
!               FRAC = 0d0
!            ENDIF
!
!            ! Start at the top of the atmosphere
!            L = LLPAR
!
!            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
!            QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
!            PTEMP        = QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
!            PDOWN(L,I,J) = PTEMP
!
!            ! PDOWN cannot be negative
!            IF( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0
!
!            ! Loop down from LLPAR to the surface
!            DO L = LLPAR-1, 1, -1
!
!               ! Compute QQ and PDOWN.  Keep PTEMP for the next level
!               QQ(L,I,J)    = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1d3
!               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) * BXHEIGHT(I,J,L) * 1d2
!               PTEMP        = PDOWN(L,I,J)
!
!               ! PDOWN cannot be negative
!               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
!            ENDDO
!         ENDIF
!      ENDDO  ! J
!      ENDDO  ! I
!!$OMP END PARALLEL DO
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_QQ
!
!!------------------------------------------------------------------------------
!
!      FUNCTION E_ICE( TK ) RESULT( VALUE )
!!
!!******************************************************************************
!!  Subroutine E_ICE computes Eice(T), the saturation vapor pressure of ice
!!  at a given Celsius temperature. (bmy, 2/8/05)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) TK (REAL*8) : Ambient temperature [K]
!!
!!  References:
!!  ============================================================================
!!  (1 ) Marti & Mauersberber (GRL '93) formulation of saturation
!!        vapor pressure of ice [Pa] is: log P = A/TK + B
!!
!!  NOTES:
!!  (1 ) Now use the same analytic function as the Goddard CTM (bmy, 2/8/05)
!!******************************************************************************
!!
!      ! Arguments as Input
!      REAL*8, INTENT(IN) :: TK
!
!      ! Return value
!      REAL*8             :: VALUE
!
!      ! Parameters
!      REAL*8, PARAMETER  :: A = -2663.5d0
!      REAL*8, PARAMETER  :: B =  12.537d0
!
!      !=================================================================
!      ! E_ICE begins here!
!      !=================================================================
!
!      ! Saturation vap press of Ice [Pa] -- divide by 100 for [hPa]
!      VALUE = ( 10d0**( A/TK + B ) ) / 100d0
!
!      ! Return to calling program
!      END FUNCTION E_ICE
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE COMPUTE_L2G( Kstar298, H298_R, TK, H2OLIQ, L2G )
!!
!!******************************************************************************
!!  Subroutine COMPUTE_L2G computes the ratio L2G = Cliq / Cgas, which is
!!  the mixing ratio of tracer in the liquid phase, divided by the mixing
!!  ratio of tracer in the gas phase.  (bmy, 2/23/00, 11/8/02)
!!
!!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate
!!  values of Kstar298 and H298_R must be supplied for each tracer.
!!  (cf Jacob et al 2000, p. 3)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Kstar298 (REAL*8) : Eff. Henry's law constant @ 298 K   [moles/atm]
!!  (2 ) H298_R   (REAL*8) : Molar heat of formation @ 298 K / R [K]
!!  (3 ) TK       (REAL*8) : Temperature at grid box (I,J,L)     [K]
!!  (4 ) H2OLIQ   (REAL*8) : Liquid water content at (I,J,L)     [cm3 H2O/cm3 air]
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (5 ) L2G      (REAL*8) : Cliq/Cgas ratio for given tracer  [unitless]
!!
!!  References (see above for full citations):
!!  ===========================================================================
!!  (1 ) Jacob et al, 2000
!!
!!  NOTES:
!!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN)  :: KStar298, H298_R, TK, H2OLIQ
!      REAL*8, INTENT(OUT) :: L2G
!
!      ! Local variables
!      REAL*8              :: Kstar
!
!      ! R = universal gas constant [atm/moles/K]
!      REAL*8, PARAMETER   :: R = 8.32d-2
!
!      ! INV_T0 = 1/298 K
!      REAL*8, PARAMETER   :: INV_T0 = 1d0 / 298d0
!
!      !=================================================================
!      ! COMPUTE_L2G begins here!
!      !=================================================================
!
!      ! Get Kstar, the effective Henry's law constant for temperature TK
!      Kstar = Kstar298 * EXP( -H298_R * ( ( 1d0 / TK ) - INV_T0 ) )
!
!      ! Use Henry's Law to get the ratio:
!      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
!      L2G   = Kstar * H2OLIQ * R * TK
!
!      ! Return to calling program
!      END SUBROUTINE COMPUTE_L2G
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE COMPUTE_F( N, F, ISOL )
!!
!!******************************************************************************
!!  Subroutine COMPUTE_F computes F, the fraction of soluble tracer lost by
!!  scavenging in convective cloud updrafts. (hyl, bmy, djj, 2/23/00, 7/26/06)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) N    (INTEGER) : Tracer number
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (2 ) F    (REAL*8)  : Fraction of tracer scavenged in cloud updraft [0-1]
!!  (3 ) ISOL (INTEGER) : Index number for ND38 diagnostic
!!
!!  References (see above for full citations):
!!  ===========================================================================
!!  (1 ) Jacob et al, 2000
!!  (2 ) Chin et al, 1996
!!
!!  NOTES:
!!  (1 ) Currently works computes scavenging fractions for either full
!!        chemistry simulation (NSRCX == 3) or Rn-Pb-Be chemistry simulation
!!        (NSRCX == 1).  Set the scavenging fraction to zero for other
!!        simulations which do not carry soluble tracers. (bmy, 3/2/00)
!!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ,
!!        and CLDICE fields once per timestep. (bmy, 2/23/00)
!!  (3 ) For aerosols only: now apply Eq. 2 for all temperatures.  Also
!!        use the distance between the grid box centers in Eq. 2.  Updated
!!        comments and made some cosmetic changes (hyl, bmy, 6/18/01)
!!  (4 ) Remove IREF, JREF -- these are obsolete.  T is now dimensioned
!!        (IIPAR,JJPAR,LLPAR).  T(IREF,JREF,L) is now T(I,J,L). (bmy, 9/27/01)
!!  (5 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!!  (6 ) Fix 2 bugs for aerosol scavenging in Rn-Pb-Be simulation:
!!        (a) set F(:,:,1) = 0 since we don't do any scavenging there.
!!        (b) DO L = 2, LLPAR to avoid any subscript range out of bounds
!!        errors (rjp, hyl, bmy, 1/10/02)
!!  (7 ) Now set F=0 in the first level for all tracers.  Also now
!!        compute the distance between grid box centers and use that in
!!        in Eq. 10 from Jacob et al, 2000 to compute F. (hyl, bmy, 1/24/02)
!!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!!  (9 ) Now reference T from "dao_mod.f" instead of from "CMN".  Also reference
!!        BXHEIGHT from "dao_mod.f" instead of from "CMN_NOX".  Now bundled
!!        into "wetscav_mod.f".  Now references IDTHNO3, IDTH2O2, etc, from
!!        F90 module "tracerid_mod.f".  Added internal routines F_AEROSOL
!!        and GET_ISOL.  Rewritten so that we don't duplicate code for
!!        different chemistry simulations. (bmy, 1/17/03)
!!  (10) Now compute F for SO2 in the same way for both fullchem and offline
!!        simulations (rjp, bmy, 3/23/03)
!!  (11) Added slots for carbon aerosol & dust tracers.  Now modified internal
!!        routine GET_ISOL so it's not hardwired anymore. (rjp, bmy, 4/5/04)
!!  (12) Added slots for sea salt aerosol tracers (rjp, bec, bmy, 4/20/04)
!!  (13) Added slots for secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!!  (14) Remove reference to CMN, it's not needed.  Made internal routine
!!        F_AEROSOL a module procedure rather than an internal routine to
!!        COMPUTE_F in order to facilitate parallelization on the Altix.  Also
!!        now pass all arguments explicitly to F_AEROSOL. (bmy, 7/20/04)
!!  (15) Now wet scavenge mercury aerosol tracers (eck, bmy, 12/9/04)
!!  (16) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!!        statement by combining branches for aerosols. (cas, bmy, 12/20/04)
!!  (17) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (19) Bug fix: Now do not over-deplete H2O2s.  Also change Henry's law
!!        constant for Hg2 to 1.0d+14. Now use functions IS_Hg2 and IS_HgP to
!!        determine if a tracer is an Hg2 or HgP tagged tracer.
!!        (dkh, rjp, eck, cdh, bmy, 1/6/06)
!!  (20) Updated for SOG4 and SOA4 (dkh, bmy, 5/18/06)
!!  (21) Bug fix: now use separate conversion factors for H2O2 and NH3.
!!        (havala, bmy, 7/26/06)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,      ONLY : BXHEIGHT, T
!      USE TRACERID_MOD, ONLY : IDTPB,    IDTBE7,   IDTHNO3, IDTH2O2
!      USE TRACERID_MOD, ONLY : IDTCH2O,  IDTMP,    IDTSO2,  IDTSO4
!      USE TRACERID_MOD, ONLY : IDTSO4s,  IDTSO4aq, IDTMSA,  IDTNH3
!      USE TRACERID_MOD, ONLY : IDTNH4,   IDTNH4aq, IDTNIT,  IDTNITs
!      USE TRACERID_MOD, ONLY : IDTAS,    IDTAHS,   IDTLET,  IDTBCPI
!      USE TRACERID_MOD, ONLY : IDTOCPI,  IDTBCPO,  IDTOCPO, IDTDST1
!      USE TRACERID_MOD, ONLY : IDTDST2,  IDTDST3,  IDTDST4, IDTSALA
!      USE TRACERID_MOD, ONLY : IDTSALC,  IDTALPH,  IDTLIMO, IDTALCO
!      USE TRACERID_MOD, ONLY : IDTSOG1,  IDTSOG2,  IDTSOG3, IDTSOG4
!      USE TRACERID_MOD, ONLY : IDTSOA1,  IDTSOA2,  IDTSOA3, IDTSOA4
!      USE TRACERID_MOD, ONLY : IS_Hg2,   IS_HgP
!      USE TRACERID_MOD, ONLY : IDTGLYX,  IDTMGLY,  IDTGLYC
!      USE TRACERID_MOD, ONLY : IDTSOAG,  IDTSOAM
!
!#     include "CMN_SIZE"    ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: N
!      REAL*8,  INTENT(OUT) :: F(IIPAR,JJPAR,LLPAR)
!      INTEGER, INTENT(OUT) :: ISOL
!
!      ! Local variables
!      INTEGER              :: I, J, L, NN
!      REAL*8               :: L2G, I2G, C_TOT, F_L, F_I, K, TMP, SO2LOSS
!
!      ! Kc is the conversion rate from cloud condensate to precip [s^-1]
!      REAL*8, PARAMETER    :: KC        = 5d-3
!
!      ! CONV_H2O2 = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
!      ! 0.6 is ( sticking  coeff H2O2  / sticking  coeff  water )
!      ! 1.9 is ( molecular weight H2O2 / molecular weight water )
!      REAL*8, PARAMETER    :: CONV_H2O2 = 8.27042925126d-1
!
!      ! CONV_NH3 = 0.6 * SQRT( 0.9 ), used for the ice to gas ratio for NH3
!      ! 0.6 is ( sticking  coeff  NH3 / sticking  coeff  water )
!      ! 0.9 is ( molecular weight NH3 / molecular weight water )
!      REAL*8, PARAMETER    :: CONV_NH3  = 5.69209978831d-1
!
!      !=================================================================
!      ! COMPUTE_F begins here!
!      !
!      ! For aerosol tracers, compute F with internal routine F_AEROSOL.
!      ! ISOL = tracer index for the ND38 diagnostic.
!      !=================================================================
!
!      !-------------------------------
!      ! 210Pb and 7Be (aerosols)
!      !-------------------------------
!      IF ( N == IDTPb .or. N == IDTBe7 ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! HNO3 (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTHNO3 ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! H2O2 (liquid & ice phases)
!      !-------------------------------
!      ELSE IF ( N == IDTH2O2 ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute ice to gas ratio for H2O2 by co-condensation
!            ! (Eq. 9, Jacob et al, 2000)
!            IF ( C_H2O(I,J,L) > 0d0 ) THEN
!               I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV_H2O2
!            ELSE
!               I2G = 0d0
!            ENDIF
!
!            ! Compute liquid to gas ratio for H2O2, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 8.3d4,    -7.4d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of H2O2 in liquid & ice phases
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G + I2G
!            F_L   = L2G / C_TOT
!            F_I   = I2G / C_TOT
!
!            ! Compute the rate constant K.  The retention factor for
!            ! liquid H2O2 is 0.05 for 248 K < T < 268 K and 1.0 for
!            ! T >= 268 K. (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * ( F_L + F_I )
!
!            ELSE IF ( T(I,J,L) > 248d0  .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( ( 5d-2 * F_L ) + F_I )
!
!            ELSE
!               K = KC * F_I
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! Compute F, the fraction of scavenged H2O2.
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! CH2O (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTCH2O ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Apply scavenging in levels 2 and higher
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for CH2O, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 3.0d3,    -7.2d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of CH2O in liquid phase
!            ! NOTE: CH2O does not exist in the ice phase!
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for
!            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of CH2O scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      ! Update GLYX and MGLY Henry's Law Const calculations (tmf, 9/13/06)
!      !-------------------------------
!      ! GLYX (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTGLYX ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Apply scavenging in levels 2 and higher
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for GLYX, using
!            ! (1) Zhou and Mopper (1990): Kstar298 = 3.6e5 M/atm
!            ! (2) Schweitzer et al. (1998) showed that the temperature dependence
!            ! for CH2O works well for glyoxal, so we use the same H298_R as CH2O
!            CALL COMPUTE_L2G( 3.6d5,   -7.2d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of GLYX in liquid phase
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! assume same retention factor as CH2O
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for
!            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of GLYX scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! MGLY (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTMGLY ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Apply scavenging in levels 2 and higher
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for MGLY, using
!            ! the appropriate parameters for Henry's law
!            ! from Betterton and Hoffman 1988): Kstar298 = 3.71d3 M/atm;
!            ! H298_R = -7.5d3 K
!            CALL COMPUTE_L2G( 3.7d3,    -7.5d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of MGLY in liquid phase
!            ! NOTE: CH2O does not exist in the ice phase!
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! assume same retention factor as CH2O
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for
!            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of MGLY scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! GLYC (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTGLYC ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Apply scavenging in levels 2 and higher
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for GLYC, using
!            ! the appropriate parameters for Henry's law
!            ! from Betterton and Hoffman 1988): Kstar298 = 4.1d4 M/atm;  H298_R = -4600 K
!            CALL COMPUTE_L2G( 4.1d4,    -4.6d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of MGLY in liquid phase
!            ! NOTE: CH2O does not exist in the ice phase!
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! assume same retention factor as CH2O
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid CH2O is 0.0 for T <= 248K and 0.02 for
!            ! 248 K < T < 268 K. (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of MGLY scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! CH3OOH (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTMP ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for CH3OOH, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 3.1d2,    -5.2d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of CH3OOH in liquid phase
!            ! NOTE: CH3OOH does not exist in the ice phase!
!            ! (Eq. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid CH3OOH is 0.0 for T <= 248 K and 0.02 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of CH3OOH scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !------------------------------
!      ! SO2 (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSO2 ) THEN
!
!         ! Compute fraction of SO2 scavenged
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!         !==============================================================
!         ! Coupled full chemistry/aerosol simulation:
!         ! Use the wet scavenging formula of Chin et al [1996],
!         ! such that a soluble fraction of SO2 is limited by the
!         ! availability of H2O2 in the precipitating grid box.
!         ! Scavenge the soluble SO2 at the same rate as the sulfate.
!         ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
!         !==============================================================
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Make sure to deplete H2O2s the same as SO2s.
!            ! (dkh, rjp, bmy, 11/17/05)
!            IF ( SO2s(I,J,L) > EPSILON ) THEN
!
!               ! Limit F
!               SO2LOSS      = MIN( H2O2s(I,J,L), SO2s(I,J,L) )
!               F(I,J,L)     = F(I,J,L) * SO2LOSS / SO2s(I,J,L)
!               F(I,J,L)     = MAX(F(I,J,L), 0d0)
!
!               ! Update saved H2O2 concentration
!               H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * F(I,J,L) )
!               H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )
!
!            ELSE
!
!               ! Set F = 0 if SO2s < EPSILON (dkh, rjp, bmy, 11/17/05)
!               F(I,J,L)     = 0d0
!
!            ENDIF
!
!            ! Update SO2
!            SO2s(I,J,L)     = SO2s(I,J,L) * ( 1d0 - F(I,J,L) )
!            SO2s(I,J,L)     = MAX( SO2s(I,J,L), EPSILON )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!      !-------------------------------
!      ! SO4   (gaseous aerosol) or
!      ! SO4aq (aqueous aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTSO4 .or. N == IDTSO4s .or. N == IDTSO4aq ) THEN
!
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! MSA (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTMSA ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! NH3 (liquid & ice phases)
!      !-------------------------------
!      ELSE IF ( N == IDTNH3 ) THEN
!
!         ! No scavenging at surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute ice to gas ratio for NH3 by co-condensation
!            ! (Eq. 9, Jacob et al, 2000)
!            IF ( C_H2O(I,J,L) > 0d0 ) THEN
!               I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV_NH3
!            ELSE
!               I2G = 0d0
!            ENDIF
!
!            ! Compute liquid to gas ratio for NH3, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 3.3d6,    -4.1d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of NH3 in liquid & ice phases
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G + I2G
!            F_L   = L2G / C_TOT
!            F_I   = I2G / C_TOT
!
!            ! Compute the rate constant K.  The retention factor
!            ! for liquid NH3 is 0.0 for T <= 248 K and 0.05 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * ( F_L + F_I )
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( ( 5d-2 * F_L ) + F_I )
!
!            ELSE
!               K = KC * F_I
!
!            ENDIF
!
!            ! F is the fraction of NH3 scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * BXHEIGHT(I,J,L) / Vud(I,J) )
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! NH4   (gaseous aerosol) or
!      ! NH4aq (aqueous aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTNH4 .or. N == IDTNH4aq ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! NIT / LET / AS / AHS (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTNIT  .or. N == IDTNITs .or.
!     &          N == IDTAS   .or. N == IDTAHS  .or.
!     &          N == IDTLET ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! BC HYDROPHILIC (aerosol) or
!      ! OC HYDROPHILIC (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTBCPI .or. N == IDTOCPI ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! BC HYDROPHOBIC (aerosol) or
!      ! OC HYDROPHOBIC (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTBCPO .or. N == IDTOCPO ) THEN
!
!         ! Force not to be lost in convective updraft for now
!         F    = 0d0
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! DST1/DST2/DST3/DST4 (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTDST1 .or. N == IDTDST2 .or.
!     &          N == IDTDST3 .or. N == IDTDST4 ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! Accum  mode seasalt (aerosol)
!      ! Coarse mode seasalt (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTSALA .or. N == IDTSALC ) THEN
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! ALPH (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTALPH ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for ALPH, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 0.023d0, 0.d0,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of ALPH in liquid phase
!            ! (Eq. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  Assume retention factor
!            ! for liquid ALPH is 0.0 for T <= 248 K and 0.02 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of ALPH scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! LIMO (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTLIMO ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for LIMO, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 0.07d0, 0.d0,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of LIMO in liquid phase
!            ! (Eq. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  Assume retention factor
!            ! for liquid LIMO is 0.0 for T <= 248 K and 0.02 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of LIMO scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! ALCO (liquid phase only)
!      !-------------------------------
!      ELSE IF ( N == IDTALCO ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for ALCO, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 54.d0, 0.d0,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of ALCO in liquid phase
!            ! (Eq. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  Assume retention factor
!            ! for liquid ALCO is 0.0 for T <= 248 K and 0.02 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of ALCO scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-----------------------------------
!      ! SOG[1,2,3,4] (liquid phase only)
!      !-----------------------------------
!      ELSE IF ( N == IDTSOG1 .or. N == IDTSOG2  .or.
!     &          N == IDTSOG3 .or. N == IDTSOG4 ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Start scavenging at level 2
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for GAS1, using
!            ! the appropriate parameters for Henry's law
!            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!            CALL COMPUTE_L2G( 1.0d5, -6.039d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of GAS1 in liquid phase
!            ! (Eq. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  Assume retention factor
!            ! for liquid GAS1 is 0.0 for T <= 248 K and 0.02 for
!            ! 248 K < T < 268 K.  (Eq. 1, Jacob et al, 2000)
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE IF ( T(I,J,L) > 248d0 .and. T(I,J,L) < 268d0 ) THEN
!               K = KC * ( 2d-2 * F_L )
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of GAS1 scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !------------------------------------------
!      ! SOA[1,2,3,4] (aerosol)
!      ! Scavenging efficiency for SOA is 0.8
!      !------------------------------------------
!      ELSE IF ( N == IDTSOA1 .or. N == IDTSOA2  .or.
!     &          N == IDTSOA3 .or. N == IDTSOA4 ) THEN
!         CALL F_AEROSOL( KC, F )
!
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            F(I,J,L) = 0.8d0 * F(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!
!      !------------------------------------------
!      ! SOAG, SOAM (aerosol)
!      ! Scavenging efficiency for SOA is 0.8
!      !------------------------------------------
!      ELSE IF ( N == IDTSOAG .or. N == IDTSOAM ) THEN
!         CALL F_AEROSOL( KC, F )
!
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!            F(I,J,L) = 0.8d0 * F(I,J,L)
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! Hg2 (liquid phase only)
!      !-------------------------------
!      ELSE IF ( IS_Hg2( N ) ) THEN
!
!         ! No scavenging at the surface
!         F(:,:,1) = 0d0
!
!         ! Apply scavenging in levels 2 and higher
!         DO L = 2, LLPAR
!         DO J = 1, JJPAR
!         DO I = 1, IIPAR
!
!            ! Compute liquid to gas ratio for HgCl2, using
!            ! the appropriate parameters for Henry's law
!            ! (Refs: INSERT HERE)
!            !
!            CALL COMPUTE_L2G( 1.0d+14, -8.4d3,
!     &                        T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!            ! Fraction of HgCl2 in liquid phase
!            ! Assume that HgCl2 is not present in ice phase
!            ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!            C_TOT = 1d0 + L2G
!            F_L   = L2G / C_TOT
!
!            ! Compute the rate constant K.  Assume scavenging takes
!            ! place only in warm clouds (retention = 0 where T<268)
!            !
!            IF ( T(I,J,L) >= 268d0 ) THEN
!               K = KC * F_L
!
!            ELSE
!               K = 0d0
!
!            ENDIF
!
!            ! Distance between grid box centers [m]
!            TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!            ! F is the fraction of HgCl2 scavenged out of the updraft
!            ! (Eq. 2, Jacob et al, 2000)
!            F(I,J,L) = 1d0 - EXP( -K * TMP / Vud(I,J) )
!
!         ENDDO
!         ENDDO
!         ENDDO
!
!         ! ND38 index
!         ISOL = GET_ISOL( N )
!
!      !-------------------------------
!      ! HgP (treat like aerosol)
!      !-------------------------------
!      ELSE IF ( IS_HgP( N ) ) THEN
!
!         CALL F_AEROSOL( KC, F )
!         ISOL = GET_ISOL( N )
!
!      !----------------------------
!      ! Insoluble tracer, set F=0
!      !----------------------------
!      ELSE
!         F(:,:,:) = 0d0
!         ISOL     = 0
!
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE COMPUTE_F
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE F_AEROSOL( KC, F )
!!
!!******************************************************************************
!!  Subroutine F_AEROSOL returns the fraction of aerosol scavenged in updrafts
!!  (bmy, 11/7/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) KC (REAL*8) : Conversion rate from cloud condensate to precip [s^-1]
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (2 ) F  (REAL*8) : Fraction of aerosol scavenged in updrafts [unitless]
!!
!!  NOTES:
!!  (1 ) Split off
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD, ONLY : BXHEIGHT
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      REAL*8, INTENT(IN)  :: KC
!      REAL*8, INTENT(OUT) :: F(IIPAR,JJPAR,LLPAR)
!
!      ! Local variables
!      INTEGER             :: I, J, L
!      REAL*8              :: TMP
!
!      !=================================================================
!      ! F_AEROSOL begins here!
!      !
!      ! Aerosol tracers are 100% in the cloud condensate phase, so
!      ! we set K = Kc, and compute F accordingly (cf Jacob et al 2000 )
!      !=================================================================
!
!      ! Turn off scavenging in the first level by setting F = 0
!      F(:,:,1) = 0d0
!
!      ! Apply scavenging in levels 2 and higher
!      DO L = 2, LLPAR
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Distance between grid box centers [m]
!         TMP = 0.5d0 * ( BXHEIGHT(I,J,L-1) + BXHEIGHT(I,J,L) )
!
!         ! (Eq. 2, Jacob et al, 2000, with K = Kc)
!         F(I,J,L) = 1d0 - EXP( -KC * TMP / Vud(I,J) )
!
!      ENDDO
!      ENDDO
!      ENDDO
!
!      ! Return to calling program
!      END SUBROUTINE F_AEROSOL
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_ISOL( N_TEST ) RESULT( VALUE )
!!
!!******************************************************************************
!!  Function GET_ISOL returns the value of ISOL (tracer index for ND38) for
!!  all simulation types.  (bmy, 4/5/04, 7/20/04)
!!
!!
!! NOTES:
!! (1 ) Now initializes a lookup table for faster execution.  Now made into
!!        an EXTERNAL function. (rjp, bmy, 4/5/04)
!! (2 ) Now references N_TRACERS from "tracer_mod.f" (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE TRACER_MOD, ONLY : N_TRACERS
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: N_TEST
!
!      ! Local variables
!      LOGICAL, SAVE       :: FIRST = .TRUE.
!      INTEGER, SAVE       :: NSOL_INDEX(NNPAR)
!      INTEGER             :: I, L, N
!
!      ! Function value
!      INTEGER             :: VALUE
!
!      !=================================================================
!      ! GET_ISOL begins here!
!      !=================================================================
!
!      ! Initialize lookup table on the first call
!      IF ( FIRST ) THEN
!
!         ! Initialize
!         NSOL_INDEX(:) = 0
!
!         ! Loop over tracers
!         DO N = 1, N_TRACERS
!
!            ! Loop over soluble tracers
!            DO L = 1, NSOL
!
!               ! Test if tracer N is among the soluble tracers
!               IF ( IDWETD(L) == N ) THEN
!
!                  ! Save location into the lookup table
!                  NSOL_INDEX(N) = L
!
!                  ! Go to next N
!                  GOTO 100
!               ENDIF
!            ENDDO
!
! 100        CONTINUE
!         ENDDO
!
!         ! Reset first-time flag
!         FIRST = .FALSE.
!      ENDIF
!
!      ! Return value
!      VALUE = NSOL_INDEX(N_TEST)
!
!      ! Return to COMPUTE_F
!      END FUNCTION GET_ISOL
!
!------------------------------------------------------------------------------
!
!      SUBROUTINE RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )
!!
!!******************************************************************************
!!  Subroutine RAINOUT computes RAINFRAC, the fraction of soluble tracer
!!  lost to rainout events in precipitation. (djj, bmy, 2/28/00, 3/5/08)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices
!!  (4  ) N        (INTEGER) : Tracer number
!!  (5  ) K_RAIN   (REAL*8 ) : Rainout rate constant for tracer N [s^-1]
!!  (6  ) DT       (REAL*8 ) : Timestep for rainout event         [s]
!!  (7  ) F        (REAL*8 ) : Fraction of grid box precipitating [unitless]
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (8  ) RAINFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!!
!!  References (see above for full citations):
!!  ============================================================================
!!  (1 ) Jacob et al, 2000
!!  (2 ) Chin et al, 1996
!!
!!  NOTES:
!!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3)
!!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!!        do not carry soluble tracer, so set RAINFRAC = 0. (bmy, 2/28/00)
!!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ,
!!        and CLDICE fields once per dynamic timestep. (bmy, 2/28/00)
!!  (3 ) K_RAIN, the rainout rate constant, and F, the areal fraction of the
!!        grid box undergoing precipitiation, are computed according to
!!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!!  (4 ) Now no longer suppress scavenging of HNO3 and aerosol below 258K.
!!        Updated comments, cosmetic changes.  Now set TK = T(I,J,L) since
!!        T is now sized (IIPAR,JJPAR,LLPAR) in "CMN". (djj, hyl, bmy, 1/24/02)
!!  (5 ) Eliminated obsolete code (bmy, 2/27/02)
!!  (6 ) Now reference T from "dao_mod.f".  Updated comments.  Now bundled
!!        into "wetscav_mod.f". Now refererences "tracerid_mod.f".  Also
!!        removed reference to CMN since we don't need NSRCX. (bmy, 11/8/02)
!!  (7 ) Now updated for carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
!!  (8 ) Now updated for seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!!  (9 ) Now updated for secondary aerosol tracers (rjp, bmy, 7/13/04)
!!  (10) Now treat rainout of mercury aerosol tracers (eck, bmy, 12/9/04)
!!  (11) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!!        statement by grouping blocks together. (cas, bmy, 12/20/04)
!!  (12) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!!  (13) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (14) Change Henry's law constant for Hg2 to 1.0d+14.  Now use functions
!!        IS_Hg2 and IS_HgP to determine if the tracer is a tagged Hg0 or
!!        HgP tracer. (eck, cdh, bmy, 1/6/06)
!!  (15) Updated for SOG4 and SOA4 (dkh, bmy, 5/18/06)
!!  (16) For GEOS-5, suppress rainout when T < 258K (hyl, bmy, 3/5/08)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,      ONLY : T
!      USE ERROR_MOD,    ONLY : ERROR_STOP
!      USE TRACERID_MOD, ONLY : IDTPB,   IDTBE7,   IDTHNO3, IDTH2O2
!      USE TRACERID_MOD, ONLY : IDTCH2O, IDTMP,    IDTSO2,  IDTSO4
!      USE TRACERID_MOD, ONLY : IDTSO4s, IDTSO4aq, IDTMSA,  IDTNH3
!      USE TRACERID_MOD, ONLY : IDTNH4,  IDTNH4aq, IDTNIT,  IDTNITs
!      USE TRACERID_MOD, ONLY : IDTAS,   IDTAHS,   IDTLET,  IDTBCPI
!      USE TRACERID_MOD, ONLY : IDTOCPI, IDTBCPO,  IDTOCPO, IDTDST1
!      USE TRACERID_MOD, ONLY : IDTDST2, IDTDST3,  IDTDST4, IDTSALA
!      USE TRACERID_MOD, ONLY : IDTSALC, IDTALPH,  IDTLIMO, IDTALCO
!      USE TRACERID_MOD, ONLY : IDTSOG1, IDTSOG2,  IDTSOG3, IDTSOG4
!      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2,  IDTSOA3, IDTSOA4
!      USE TRACERID_MOD, ONLY : IS_Hg2,  IS_HgP
!      USE TRACERID_MOD, ONLY : IDTGLYX, IDTMGLY,  IDTGLYC
!      USE TRACERID_MOD, ONLY : IDTSOAG, IDTSOAM
!
!      IMPLICIT NONE
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: I, J, L, N
!      REAL*8,  INTENT(IN)  :: K_RAIN, DT, F
!      REAL*8,  INTENT(OUT) :: RAINFRAC
!
!      ! Local variables
!      REAL*8               :: L2G, I2G, C_TOT, F_L, F_I, K, TK, SO2LOSS
!
!      ! CONV = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
!      REAL*8, PARAMETER    :: CONV = 8.27042925126d-1
!
!      !==================================================================
!      ! RAINOUT begins here!
!      !
!      ! For aerosols, set K = K_RAIN and compute RAINFRAC according
!      ! to Eq. 10 of Jacob et al 2000.  Call function GET_RAINFRAC.
!      !==================================================================
!
!      ! Save the local temperature in TK for convenience
!      TK = T(I,J,L)
!
!#if   defined( GEOS_5 )
!      !------------------------------------------------------------------
!      ! NOTE FROM HONGYU LIU (hyl@nianet.org) -- 3/5/08
!      !
!      ! Lead-210 (210Pb) and Beryllium-7 (7Be) simulations indicate
!      ! that we can improve the GEOS-5 simulation by (1) turning off
!      ! rainout/washout for convective precip (see DO_WETDEP)
!      ! and (2) suppressing rainout for large-scale precip at
!      ! temperatures below 258K.
!      !
!      ! Place an #if block here to set RAINFRAC=0 when T < 258K for
!      ! GEOS-5 met.  This will suppress rainout. (hyl, bmy, 3/5/08)
!      !-------------------------------------------------------------------
!      IF ( TK < 258d0 ) THEN
!         RAINFRAC = 0d0
!         RETURN
!      ENDIF
!#endif
!
!      !------------------------------
!      ! 210Pb and 7Be (aerosol)
!      !------------------------------
!      IF ( N == IDTPb .or. N == IDTBe7 ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! HNO3 (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTHNO3 ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! H2O2 (liquid & ice phases)
!      !------------------------------
!      ELSE IF ( N == IDTH2O2 ) THEN
!
!         ! Compute ice to gas ratio for H2O2 by co-condensation
!         ! (Eq. 9, Jacob et al, 2000)
!         IF ( C_H2O(I,J,L) > 0d0 ) THEN
!            I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
!         ELSE
!            I2G = 0d0
!         ENDIF
!
!         ! Compute liquid to gas ratio for H2O2, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 8.3d4, -7.4d3, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of H2O2 in liquid & ice phases
!         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G + I2G
!         F_L   = L2G / C_TOT
!         F_I   = I2G / C_TOT
!
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid H2O2 is 0.05 for 248 K < T < 268 K, and
!         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * ( F_L + F_I )
!
!         ELSE IF ( TK > 248d0  .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( ( 5d-2 * F_L ) + F_I )
!
!         ELSE
!            K = K_RAIN * F_I
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out H2O2
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! CH2O (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTCH2O ) THEN
!
!         ! Compute liquid to gas ratio for CH2O, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8 and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 3.0d3, -7.2d3, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of CH2O in liquid phase
!         ! NOTE: CH2O does not exist in the ice phase!
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out CH2O
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      ! Update GLYX and MGLY Henry's Law Const calculations (tmf, 9/13/06)
!      !------------------------------
!      ! GLYX (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTGLYX ) THEN
!
!         ! Compute liquid to gas ratio for GLYX, using
!         ! (1) Zhou and Mopper (1990): Kstar298 = 3.6e5 M/atm
!         ! (2) Schweitzer et al. (1998) showed that the temperature dependence
!         ! for CH2O works well for glyoxal, so we use the same H298_R as CH2O
!         CALL COMPUTE_L2G( 3.6d5,   -7.2d3,
!     &                     T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of GLYX in liquid phase
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! assume same retention factor as CH2O
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out GLYX
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! MGLY (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTMGLY ) THEN
!
!         ! Compute liquid to gas ratio for MGLY, using
!         ! the appropriate parameters for Henry's law
!         ! from Betterton and Hoffman 1988): Kstar298 = 3.71d3 M/atm;  H298_R = -7.5d3 K
!         CALL COMPUTE_L2G( 3.7d3,    -7.5d3,
!     &                     T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!
!         ! Fraction of MGLY in liquid phase
!         ! NOTE: CH2O does not exist in the ice phase!
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! assume same retention factor as CH2O
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out MGLY
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! GLYC (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTGLYC ) THEN
!
!         ! Compute liquid to gas ratio for GLYC, using
!         ! the appropriate parameters for Henry's law
!         ! from Betterton and Hoffman 1988): Kstar298 = 4.1d4 M/atm;  H298_R = -4.6d3 K
!         CALL COMPUTE_L2G( 4.1d4,    -4.6d3,
!     &                     T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!
!         ! Fraction of GLYC in liquid phase
!         ! NOTE: CH2O does not exist in the ice phase!
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! assume same retention factor as CH2O
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid CH2O is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out MGLY
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!
!      !------------------------------
!      ! CH3OOH (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTMP ) THEN
!
!         ! Compute liquid to gas ratio for CH3OOH, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 3.1d2, -5.2d3, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of CH3OOH in liquid phase
!         ! NOTE: CH3OOH does not exist in the ice phase!
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid CH3OOH is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out CH3OOH
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! SO2
!      !------------------------------
!      ELSE IF ( N == IDTSO2 ) THEN
!
!         !==============================================================
!         ! NOTE: SO2 and H2O2 are in [v/v] and here RAINFRAC contains
!         ! the amount of SO2 lost due to rainout normalized by the
!         ! total SO2 -- so that in WETDEP routine mulitiplying SO2 in
!         ! [kg] will produce correct amount.  Need to verify this.
!         ! (rjp, 01/16/02)
!         !==============================================================
!
!         ! Treat SO2 as an aerosol
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!         ! We need to save a copy of the original RAINFRAC when
!         ! recalculating these terms for the adjoint .  (dkh, 10/26/05)
!         IF ( LADJ ) THEN
!            RAINFRAC_0(L) = RAINFRAC
!            IF ( L_PRINTFD .and. I == IFD .and. J == JFD .and. L == LFD)
!     &      THEN
!               print*, ' rainfrac_0 = ', rainfrac_0(L)
!            ENDIF
!         ENDIF
!
!         ! Update SO2 and H2O2
!         IF ( SO2s(I,J,L) > EPSILON ) THEN
!
!            ! Limit RAINFRAC
!            SO2LOSS      = MIN( SO2s(I,J,L), H2O2s(I,J,L) )
!            RAINFRAC     = SO2LOSS * RAINFRAC / SO2s(I,J,L)
!            RAINFRAC     = MAX( RAINFRAC, 0d0 )
!
!            ! Update saved H2O2 concentration
!            H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * RAINFRAC )
!            H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )
!
!         ELSE
!            RAINFRAC = 0D0
!
!         ENDIF
!
!         ! Update saved SO2 concentration
!         SO2s(I,J,L) = SO2s(I,J,L) * ( 1.D0 - RAINFRAC )
!         SO2s(I,J,L) = MAX( SO2s(I,J,L), EPSILON )
!
!      !----------------------------
!      ! SO4 and SO4aq (aerosol)
!      !----------------------------
!      ELSE IF ( N == IDTSO4 .or. N == IDTSO4s .or. N == IDTSO4aq ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! MSA (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTMSA ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! NH3 (liquid & ice phases)
!      !------------------------------
!      ELSE IF ( N == IDTNH3 ) THEN
!
!         ! Compute ice to gas ratio for NH3 by co-condensation
!         ! (Eq. 9, Jacob et al, 2000)
!         IF ( C_H2O(I,J,L) > 0d0 ) THEN
!            I2G = ( CLDICE(I,J,L) / C_H2O(I,J,L) ) * CONV
!         ELSE
!            I2G = 0d0
!         ENDIF
!
!         ! Compute liquid to gas ratio for NH3, using
!         ! the appropriate parameters for Henry's law
!         ! (Seinfeld and Pandis, p343 eq. 6.8)
!         ! PH    = 4.5  ! Assumed PH for typical cloud drop
!         ! Hstar = 1.054d11 * (10.**(-PH)) == 3.3d6
!         CALL COMPUTE_L2G( 3.3d6, -4.1d3, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of NH3 in liquid & ice phases
!         ! (Eqs. 4, 5, 6, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G + I2G
!         F_L   = L2G / C_TOT
!         F_I   = I2G / C_TOT
!
!         ! Compute the rate constant K.  The retention factor
!         ! for liquid NH3 is 0.05 for 248 K < T < 268 K, and
!         ! 1.0 for T >= 268 K.  (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * ( F_L + F_I )
!
!         ELSE IF ( TK > 248d0  .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( ( 5d-2 * F_L ) + F_I )
!
!         ELSE
!            K = K_RAIN * F_I
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out NH3
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! NH4 and NH4aq (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTNH4 .or. N == IDTNH4aq ) THEN
!
!         ! NOTE: NH4aq may have a henry's law constant;
!         !       Carine will investigate (cas, bmy, 12/20/04)
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! NIT/AS/AHS/LET (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTNIT .or. N == IDTNITs .or.
!     &          N == IDTAS  .or. N == IDTAHS  .or.
!     &          N == IDTLET ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! BC HYDROPHILIC (aerosol) or
!      ! OC HYDROPHILIC (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTBCPI .or. N == IDTOCPI) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !-------------------------------
!      ! BC HYDROPHOBIC (aerosol) or
!      ! OC HYDROPHOBIC (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTBCPO .or. N == IDTOCPO ) THEN
!
!         ! No rainout
!         RAINFRAC = 0.0D0
!
!      !-------------------------------
!      ! DUST all size bins (aerosol)
!      !-------------------------------
!      ELSE IF ( N == IDTDST1 .or. N == IDTDST2 .or.
!     &          N == IDTDST3 .or. N == IDTDST4 ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! Accum  seasalt (aerosol) or
!      ! Coarse seasalt (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSALA .or. N == IDTSALC ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! ALPH (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTALPH ) THEN
!
!         ! Compute liquid to gas ratio for ALPH, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 0.023d0, 0.d0, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of ALPH in liquid phase
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  Assume that the retention factor
!         ! for liquid ALPH is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out ALPH
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! LIMO (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTLIMO ) THEN
!
!         ! Compute liquid to gas ratio for LIMO, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 0.07d0, 0.d0, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of LIMO in liquid phase
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  Assume that the retention factor
!         ! for liquid LIMO is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out LIMO
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! ALCO (liquid phase only)
!      !------------------------------
!      ELSE IF ( N == IDTALCO ) THEN
!
!         ! Compute liquid to gas ratio for ALCO, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 54.d0, 0.d0, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of ALCO in liquid phase
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  Assume that the retention factor
!         ! for liquid ALCO is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out ALCO
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !----------------------------------
!      ! SOG[1,2,3,4] (liquid phase only)
!      !----------------------------------
!      ELSE IF ( N == IDTSOG1 .or. N == IDTSOG2  .or.
!     &          N == IDTSOG3 .or. N == IDTSOG4 ) THEN
!
!         ! Compute liquid to gas ratio for GAS1, using
!         ! the appropriate parameters for Henry's law
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( 1.0d5, -6.039d3, TK, CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of GAS1 in liquid phase
!         ! (Eqs. 4, 5, Jacob et al, 2000)
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  Assume that the retention factor
!         ! for liquid GAS1 is 0.02 for 248 K < T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            K = K_RAIN * ( 2d-2 * F_L )
!
!         ELSE
!            K = 0d0
!
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out SOG{1,2,3}
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !--------------------------------------
!      ! SOA[1,2,3,4] (aerosol)
!      ! Scavenging efficiency for SOA is 0.8
!      !--------------------------------------
!      ELSE IF ( N == IDTSOA1 .or. N == IDTSOA2  .or.
!     &          N == IDTSOA3 .or. N == IDTSOA4 ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!         RAINFRAC = RAINFRAC * 0.8d0
!
!      !--------------------------------------
!      ! SOAG and SOAM (aerosol)
!      ! Scavenging efficiency for SOA is 0.8
!      !--------------------------------------
!      ELSE IF ( N == IDTSOAG .OR. N == IDTSOAM ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!         RAINFRAC = RAINFRAC * 0.8d0
!
!      !------------------------------
!      ! Hg2 (liquid phase only)
!      !------------------------------
!      ELSE IF ( IS_Hg2( N ) ) THEN
!
!         ! Compute liquid to gas ratio for HgCl2, using
!         ! the appropriate parameters for Henry's law
!         ! (Refs: INSERT HERE)
!         CALL COMPUTE_L2G( 1.0d+14, -8.4d3,
!     &                     T(I,J,L), CLDLIQ(I,J,L), L2G )
!
!         ! Fraction of HgCl2 in liquid phase
!         ! Assume no HgCl2 in the ice phase
!         C_TOT = 1d0 + L2G
!         F_L   = L2G / C_TOT
!
!         ! Compute the rate constant K.  Assume the retention factor
!         ! for liquid HgCl2 is 0 for T < 268 K, and
!         ! 1.0 for T > 268 K. (Eq. 1, Jacob et al, 2000)
!         IF ( TK >= 268d0 ) THEN
!            K = K_RAIN * F_L
!         ELSE
!            K = 0d0
!         ENDIF
!
!         ! Compute RAINFRAC, the fraction of rained-out HgCl2
!         ! (Eq. 10, Jacob et al, 2000)
!         RAINFRAC = GET_RAINFRAC( K, F, DT )
!
!      !------------------------------
!      ! HgP (treat like aerosol)
!      !------------------------------
!      ELSE IF ( IS_HgP( N ) ) THEN
!         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )
!
!      !------------------------------
!      ! ERROR: insoluble tracer!
!      !------------------------------
!      ELSE
!         CALL ERROR_STOP( 'Invalid tracer!', 'RAINOUT (wetscav_mod.f)' )
!
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE RAINOUT
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_RAINFRAC( K, F, DT ) RESULT( RAINFRAC )
!!
!!******************************************************************************
!!  Function GET_RAINFRAC computes the fraction of tracer lost to rainout
!!  according to Jacob et al 2000. (bmy, 11/8/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ===========================================================================
!!  (1 ) K  (REAL*8) : Rainout rate constant              [1/s]
!!  (2 ) DT (REAL*8) : Timestep for rainout event         [s]
!!  (3 ) F  (REAL*8) : Fraction of grid box precipitating [unitless]
!!
!!  NOTES:
!!  (1 ) Now move internal routines GET_RAINFRAC to the module and pass all
!!        arguments explicitly.  This facilitates parallelization on the
!!        Altix platform (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: K, F, DT
!
!      ! Local variables
!      REAL*8             :: RAINFRAC
!
!      !=================================================================
!      ! GET_RAINFRAC begins here!
!      !=================================================================
!
!      ! (Eq. 10, Jacob et al, 2000 )
!      RAINFRAC = F * ( 1 - EXP( -K * DT ) )
!
!      ! Return to RAINOUT
!      END FUNCTION GET_RAINFRAC
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE WASHOUT( I, J, L, N, PP, DT, F, WASHFRAC, AER )
!!
!!******************************************************************************
!!  Subroutine WASHOUT computes WASHFRAC, the fraction of soluble tracer
!!  lost to washout events in precipitation. (djj, bmy, 2/28/00, 5/18/06)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices [unitless]
!!  (4  ) N        (INTEGER) : Tracer number                [unitless]
!!  (5  ) PP       (REAL*8 ) : Precip rate thru at bottom
!!                             of grid box (I,J,L)          [cm3 H2O/cm2 air/s]
!!  (6  ) DT       (REAL*8 ) : Timestep for rainout event   [s]
!!  (7  ) F        (REAL*8 ) : Fraction of grid box
!!                             that is precipitating        [unitless]
!!
!!  Arguments as Output:
!!  ============================================================================
!!  (8  ) WASHFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!!  (9  ) AER      (LOGICAL) : = T if the tracer is an aerosol, =F otherwise
!!
!!  Reference (see above for full citations):
!!  ============================================================================
!!  (1  ) Jacob et al, 2000
!!
!!  NOTES:
!!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3)
!!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!!        do not carry soluble tracers, so set WASHFRAC = 0.
!!  (2 ) K_WASH, the rainout rate constant, and F, the areal fraction of the
!!        grid box undergoing precipitiation, are computed according to
!!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!!  (3 ) Washout is only done for T >= 268 K, when the cloud condensate is
!!        in the liquid phase.
!!  (4 ) T(I+I0,J+J0,L) is now T(I,J,L).  Removed IREF, JREF -- these are
!!        obsolete.  Updated comments. (bmy, 9/27/01)
!!  (5 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!!  (6 ) Now reference BXHEIGHT, T from "dao_mod.f".  Also remove reference
!!        to "CMN_NOX".  Updated comments.  Now bundled into "wetscav_mod.f".
!!        Now also references "tracerid_mod.f".  Added internal routines
!!        WASHFRAC_AEROSOL and WASHFRAC_LIQ_GAS.  Also removed reference to
!!        CMN since we don't need to use NSRCX here. (bmy, 11/6/02)
!!  (7 ) Updated for carbon aerosol and dust tracers (rjp, bmy, 4/5/04)
!!  (8 ) Updated for seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!!  (9 ) Updated for secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!!  (10) Now move internal routines WASHFRAC_AEROSOL and WASHFRAC_LIQ_GAS
!!        to the module and pass all arguments explicitly.  This facilitates
!!        parallelization on the Altix platform (bmy, 7/20/04)
!!  (11) Now handle washout of mercury aerosol tracers (eck, bmy, 12/9/04)
!!  (13) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!!        statement by grouping blocks together (cas, bmy, 12/20/04)
!!  (14) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!!  (15) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (16) Bug fix: Deplete H2O2s the same as SO2s.  Also change Henry's law
!!        constant for Hg2 to 1.0d+14. Now use functions IS_Hg2 and IS_HgP to
!!        determine if a tracer is a tagged Hg0 or HgP tracer.
!!        (dkh, rjp, eck, cdh, bmy, 1/6/06)
!!  (17) Updated for SOG4 and SOA4 (bmy, 5/18/06)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,      ONLY : BXHEIGHT, T
!      USE ERROR_MOD,    ONLY : ERROR_STOP
!      USE TRACERID_MOD, ONLY : IDTPB,    IDTBE7,   IDTHNO3, IDTH2O2
!      USE TRACERID_MOD, ONLY : IDTCH2O,  IDTMP,    IDTSO2,  IDTSO4
!      USE TRACERID_MOD, ONLY : IDTSO4s,  IDTSO4aq, IDTMSA,  IDTNH3
!      USE TRACERID_MOD, ONLY : IDTNH4,   IDTNH4aq, IDTNIT,  IDTNITs
!      USE TRACERID_MOD, ONLY : IDTAS,    IDTAHS,   IDTLET,  IDTBCPI
!      USE TRACERID_MOD, ONLY : IDTOCPI,  IDTBCPO,  IDTOCPO, IDTDST1
!      USE TRACERID_MOD, ONLY : IDTDST2,  IDTDST3,  IDTDST4, IDTSALA
!      USE TRACERID_MOD, ONLY : IDTSALC,  IDTALPH,  IDTLIMO, IDTALCO
!      USE TRACERID_MOD, ONLY : IDTSOG1,  IDTSOG2,  IDTSOG3, IDTSOG4
!      USE TRACERID_MOD, ONLY : IDTSOA1,  IDTSOA2,  IDTSOA3, IDTSOA4
!      USE TRACERID_MOD, ONLY : IS_Hg2,   IS_HgP
!      USE TRACERID_MOD, ONLY : IDTGLYX,  IDTMGLY,  IDTGLYC
!      USE TRACERID_MOD, ONLY : IDTSOAG,  IDTSOAM
!
!#     include "CMN_SIZE"   ! Size parameters
!
!      ! Arguments
!      INTEGER, INTENT(IN)  :: I, J, L, N
!      REAL*8,  INTENT(IN)  :: PP, DT, F
!      REAL*8,  INTENT(OUT) :: WASHFRAC
!      LOGICAL, INTENT(OUT) :: AER
!
!      ! Local variables
!      REAL*8               :: L2G, DZ, TK, SO2LOSS
!
!      ! First order washout rate constant for HNO3, aerosols = 1 cm^-1
!      REAL*8, PARAMETER    :: K_WASH = 1d0
!
!      !=================================================================
!      ! WASHOUT begins here!
!      !
!      ! Call either WASHFRAC_AEROSOL or WASHFRAC_LIQ_GAS to compute the
!      ! fraction of tracer lost to washout according to Jacob et al 2000
!      !=================================================================
!
!      ! TK is Kelvin temperature
!      TK = T(I,J,L)
!
!      ! DZ is the height of the grid box in cm
!      DZ = BXHEIGHT(I,J,L) * 1d2
!
!      !------------------------------
!      ! 210Pb or 7Be (aerosol)
!      !------------------------------
!      IF ( N == IDTPb .or. N == IDTBe7 ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! HNO3 (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTHNO3 ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! H2O2 (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTH2O2 ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 8.3d4, -7.4d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! CH2O (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTCH2O ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 3.0d3, -7.2d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! GLYX (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTGLYX ) THEN
!
!         ! Compute liquid to gas ratio for GLYX, using
!         ! (1) Zhou and Mopper (1990): Kstar298 = 3.6e5 M/atm
!         ! (2) Schweitzer et al. (1998) showed that the temperature dependence for CH2O works well for glyoxal,
!         !      so we use the same H298_R as CH2O
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 3.6d5, -7.2d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! MGLY (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTMGLY ) THEN
!         ! Compute liquid to gas ratio for MGLY, using
!         ! the appropriate parameters for Henry's law
!         ! from Betterton and Hoffman 1988): Kstar298 = 3.71d3 M/atm;  H298_R = -7.5d3 K
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 3.7d3, -7.5d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! GLYC (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTGLYC ) THEN
!         ! Compute liquid to gas ratio for GLYC, using
!         ! the appropriate parameters for Henry's law
!         ! from Betterton and Hoffman 1988): Kstar298 = 4.6d4 M/atm;  H298_R = -4.6d3 K
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 4.1d4, -4.6d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!
!      !------------------------------
!      ! MP (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTMP ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 3.1d2, -5.2d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! SO2 (aerosol treatment)
!      !------------------------------
!      ELSE IF ( N == IDTSO2 ) THEN
!
!         !==============================================================
!         ! NOTE: Even though SO2 is not an aerosol we treat it as SO4 in
!         ! wet scavenging.  When evaporation occurs, it returns to SO4.
!         !==============================================================
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!         ! For adjoint recalculation, need to save the initial value of
!         ! WASHFRAC here.  (dkh, 10/26/05)
!         IF ( LADJ ) THEN
!            WASHFRAC_0(L) = WASHFRAC
!         ENDIF
!
!
!         !==============================================================
!         ! Use the wet-scavenging following [Chin et al, 1996] such
!         ! that a soluble fraction of SO2 is limited by the availability
!         ! of H2O2 in the precipitating grid box.  Then scavenge the
!         ! soluble SO2 at the same rate as sulfate.
!         !==============================================================
!         IF ( TK >= 268d0 .AND. SO2s(I,J,L) > EPSILON ) THEN
!
!            ! Adjust WASHFRAC
!            SO2LOSS  = MIN( SO2s(I,J,L), H2O2s(I,J,L) )
!            WASHFRAC = SO2LOSS * WASHFRAC / SO2s(I,J,L)
!            WASHFRAC = MAX( WASHFRAC, 0d0 )
!
!            ! Deplete H2O2s the same as SO2s (dkh, rjp, bmy, 11/17/05)
!            H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * WASHFRAC )
!            H2O2s(I,J,L) = MAX( H2O2s(I,J,L), EPSILON )
!
!         ELSE
!            WASHFRAC = 0d0
!
!         ENDIF
!
!         ! Update saved SO2 concentration
!         SO2s(I,J,L) = SO2s(I,J,L) * ( 1d0 - WASHFRAC )
!         SO2s(I,J,L) = MAX( SO2s(I,J,L), EPSILON )
!
!      !------------------------------
!      ! SO4 and SO4aq (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSO4 .or. N == IDTSO4s .or. N == IDTSO4aq ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! MSA (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTMSA ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! NH3 (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTNH3 ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 3.3d6, -4.1d3, PP, DT,
!     &                                F,      DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! NH4 and NH4aq (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTNH4 .or. N == IDTNH4aq ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! NIT/AS/AHS/LET (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTNIT  .or. N == IDTNITs .or.
!     &          N == IDTAS   .or. N == IDTAHS  .or.
!     &          N == IDTLET ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! BC HYDROPHILIC (aerosol) or
!      ! OC HYDROPHILIC (aerosol) or
!      ! BC HYDROPHOBIC (aerosol) or
!      ! OC HYDROPHOBIC (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTBCPI .or. N == IDTOCPI  .or.
!     &          N == IDTBCPO .or. N == IDTOCPO ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! DUST all size bins (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTDST1 .or. N == IDTDST2  .or.
!     &          N == IDTDST3 .or. N == IDTDST4 ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! Accum  seasalt (aerosol) or
!      ! Coarse seasalt (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSALA .or. N == IDTSALC ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! ALPH (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTALPH ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 0.023d0, 0.d0, PP, DT,
!     &                                F,       DZ,   TK, K_WASH )
!
!      !------------------------------
!      ! LIMO (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTLIMO ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 0.07d0, 0.d0, PP, DT,
!     &                                F,      DZ,   TK, K_WASH )
!
!      !------------------------------
!      ! ALCO (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( N == IDTALCO ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 54.d0, 0.d0, PP, DT,
!     &                                F,     DZ,   TK, K_WASH )
!
!      !---------------------------------
!      ! SOG[1,2,3,4] (liq & gas phases)
!      !---------------------------------
!      ELSE IF ( N == IDTSOG1 .or. N == IDTSOG2  .or.
!     &          N == IDTSOG3 .or. N == IDTSOG4 ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 1.0d5, -6.039d3, PP, DT,
!     &                                F,      DZ,      TK, K_WASH )
!
!      !------------------------------
!      ! SOA[1,2,3,4] (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSOA1 .or. N == IDTSOA2  .or.
!     &          N == IDTSOA3 .or. N == IDTSOA4 ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! SOAG and SOAM (aerosol)
!      !------------------------------
!      ELSE IF ( N == IDTSOAG .or. N == IDTSOAM ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! Hg2 (liquid & gas phases)
!      !------------------------------
!      ELSE IF ( IS_Hg2( N ) ) THEN
!         AER      = .FALSE.
!         WASHFRAC = WASHFRAC_LIQ_GAS( 1.0d+14, -8.4d3, PP, DT,
!     &                                F,        DZ,    TK, K_WASH )
!
!      !------------------------------
!      ! HgP (treat like aerosol)
!      !------------------------------
!      ELSE IF ( IS_HgP( N ) ) THEN
!         AER      = .TRUE.
!         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!
!      !------------------------------
!      ! ERROR: Insoluble tracer
!      !------------------------------
!      ELSE
!         CALL ERROR_STOP( 'Invalid tracer!', 'WASHOUT (wetscav_mod.f)' )
!
!      ENDIF
!
!      ! Return to calling program
!      END SUBROUTINE WASHOUT
!
!!------------------------------------------------------------------------------
!
!      FUNCTION WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )
!     &         RESULT( WASHFRAC )
!!
!!******************************************************************************
!!  Function WASHFRAC_AEROSOL returns the fraction of soluble aerosol tracer
!!  lost to washout.  (bmy, 11/8/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) TK       (REAL*8 ) : Temperature in grid box            [K]
!!  (2 ) F        (REAL*8 ) : Fraction of grid box
!!                             that is precipitating             [unitless]
!!  (3 ) K_WASH   (REAL*8 ) : 1st order washout rate constant    [1/cm]
!!  (3 ) PP       (REAL*8 ) : Precip rate thru at bottom
!!                             of grid box (I,J,L)           [cm3 H2O/cm2 air/s]
!!
!!  NOTES:
!!  (1 ) WASHFRAC_AEROSOL used to be an internal function to subroutine WASHOUT.
!!        This caused NaN's in the parallel loop on Altix, so we moved it to
!!        the module and now pass Iall arguments explicitly (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: DT, F, K_WASH, PP, TK
!
!      ! Function value
!      REAL*8             :: WASHFRAC
!
!      !=================================================================
!      ! WASHFRAC_AEROSOL begins here!
!      !=================================================================
!
!      ! Washout only happens at or above 268 K
!      IF ( TK >= 268d0 ) THEN
!         WASHFRAC = F * ( 1d0 - EXP( -K_WASH * ( PP / F ) * DT ) )
!      ELSE
!         WASHFRAC = 0d0
!      ENDIF
!
!      ! Return to calling program
!      END FUNCTION WASHFRAC_AEROSOL
!
!!------------------------------------------------------------------------------
!
!      FUNCTION WASHFRAC_LIQ_GAS( Kstar298, H298_R, PP, DT,
!     &                           F,        DZ,     TK, K_WASH )
!     &         RESULT( WASHFRAC )
!!
!!******************************************************************************
!!  Function WASHFRAC_LIQ_GAS returns the fraction of soluble liquid/gas phase
!!  tracer lost to washout. (bmy, 11/8/02, 7/20/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Kstar298 (REAL*8 ) : Eff. Henry's law constant @ 298 K  [moles/atm]
!!  (2 ) H298_R   (REAL*8 ) : Henry's law coefficient            [K]
!!  (3 ) PP       (REAL*8 ) : Precip rate thru at bottom
!!                             of grid box (I,J,L)           [cm3 H2O/cm2 air/s]
!!  (4 ) DT       (REAL*8 ) : Dynamic timestep                   [s]
!!  (5 ) F        (REAL*8 ) : Fraction of grid box
!!                             that is precipitating             [unitless]
!!  (6 ) DZ       (REAL*8 ) : Height of grid box                 [cm]
!!  (7 ) TK       (REAL*8 ) : Temperature in grid box            [K]
!!  (8 ) K_WASH   (REAL*8 ) : 1st order washout rate constant    [1/cm]
!!
!!  NOTES:
!!  (1 ) WASHFRAC_LIQ_GAS used to be an internal function to subroutine WASHOUT.
!!        This caused NaN's in the parallel loop on Altix, so we moved it to
!!        the module and now pass all arguments explicitly (bmy, 7/20/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: Kstar298, H298_R, PP, DT, F, DZ, TK, K_WASH
!
!      ! Local variables
!      REAL*8             :: L2G, LP, WASHFRAC, WASHFRAC_F_14
!
!      !=================================================================
!      ! WASHFRAC_LIQ_GAS begins here!
!      !=================================================================
!
!      ! Suppress washout below 268 K
!      IF ( TK >= 268d0 ) THEN
!
!         ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
!         LP = ( PP * DT ) / ( F * DZ )
!
!         ! Compute liquid to gas ratio for H2O2, using the appropriate
!         ! parameters for Henry's law -- also use rainwater content Lp
!         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
!         CALL COMPUTE_L2G( Kstar298, H298_R, TK, LP, L2G )
!
!         ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
!         WASHFRAC = L2G / ( 1d0 + L2G )
!
!         ! Washout fraction / F from Eq. 14, Jacob et al, 2000
!         WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * DT )
!
!         ! Do not let the Henry's law washout fraction exceed
!         ! ( washout fraction / F ) from Eq. 14 -- this is a cap
!         IF ( WASHFRAC > WASHFRAC_F_14 ) WASHFRAC = WASHFRAC_F_14
!
!      ELSE
!         WASHFRAC = 0d0
!
!      ENDIF
!
!      ! Return to calling program
!      END FUNCTION WASHFRAC_LIQ_GAS
!
!!------------------------------------------------------------------------------

      SUBROUTINE WETDEP_ADJ( LS )
!
!******************************************************************************
!  Subroutine WETDEP_ADJ is the same as WETDEP except
!   - it acts on the adjoint tracers
!   - negative sensitivity values are allowed (i.e., don't call SAFETY)
!   - skip SO2, which is treated separately.
!  (dkh, 10/24/05, 09/28/09)
!
!  Based on WETDEP, which computes the downward mass flux of tracer due to washout
!  and rainout of aerosols and soluble tracers in a column.  The timestep is
!  the dynamic timestep. (hyl, bey, bmy, djj, 4/2/99, 5/24/06)
!
! Notes
! (1 ) All changes have ADJ or adj_group in them.
! (2 ) Completely revised to properly treat adjoint of DSTT (dkh, 03/18/12)
! (3 ) Now support deposition cost function (fp, dkh, 03/04/13)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : IFD, JFD, NFD, LFD
      USE DAO_MOD,           ONLY : BXHEIGHT
      USE DIAG_MOD,          ONLY : AD16, AD17, AD18
      USE DIAG_MOD,          ONLY : CT16, CT17, CT18, AD39
      USE ERROR_MOD,         ONLY : GEOS_CHEM_STOP, IT_IS_NAN
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_WDEP_LS
      USE OCEAN_MERCURY_MOD, ONLY : ADD_Hg2_WD
      USE TIME_MOD,          ONLY : GET_TS_DYN
      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM
      USE TRACERID_MOD,      ONLY : IDTSO2, IDTSO4, IS_Hg2
      USE WETSCAV_MOD,       ONLY : QQ
      USE WETSCAV_MOD,       ONLY : IDWETD
      USE WETSCAV_MOD,       ONLY : LS_K_RAIN
      USE WETSCAV_MOD,       ONLY : LS_F_PRIME
      USE WETSCAV_MOD,       ONLY : PDOWN
      USE WETSCAV_MOD,       ONLY : CONV_F_PRIME
      USE WETSCAV_MOD,       ONLY : RAINOUT
      USE WETSCAV_MOD,       ONLY : WASHOUT
      USE WETSCAV_MOD,       ONLY : NSOL
      ! dkh debug
      USE WETSCAV_MOD,       ONLY : H2O2s
      USE WETSCAV_MOD,       ONLY : SO2s
      USE TRACER_MOD,        ONLY : STT
      USE LOGICAL_ADJ_MOD,   ONLY : LPRINTFD

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! Diagnostic arrays and switches

      ! Arguments
      LOGICAL, INTENT(IN) :: LS

      ! Local Variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: IS_Hg
      LOGICAL             :: AER

      INTEGER             :: I, IDX, J, L, N, NN

      REAL*8              :: Q,     QDOWN,  DT,        DT_OVER_TAU
      REAL*8              :: K,     K_MIN,  K_RAIN,    RAINFRAC
      REAL*8              :: F,     FTOP,   F_PRIME,   WASHFRAC
      REAL*8              :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
      REAL*8              :: ALPHA, ALPHA2, WETLOSS,   TMP
      REAL*8              :: WETLOSS_ADJ

      ! DSTT is the accumulator array of rained-out
      ! soluble tracer for a given (I,J) column
      REAL*8              :: DSTT_ADJ(NSOL,LLPAR,IIPAR,JJPAR)

      REAL*8              :: F_SAVE(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! WETDEP_ADJ begins here!
      !
      ! (1)  I n i t i a l i z e   V a r i a b l e s
      !=================================================================

      ! Is this a mercury simulation with dynamic online ocean?
      IS_Hg = ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN )

      ! Dynamic timestep [s]
      DT    = GET_TS_DYN() * 60d0

      ! Select index for diagnostic arrays -- will archive either
      ! large-scale or convective rainout/washout fractions
      IF ( LS ) THEN
         IDX = 1
      ELSE
         IDX = 2
      ENDIF

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !
      ! just recalculate the values of F and save in F_SAVE
      !=================================================================

      ! initialize
      F_SAVE(:,:,:) = 0d0

#if   !defined( SGI_MIPS )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,           FTOP     )
!$OMP+PRIVATE( F,       F_PRIME,     K_RAIN   )
!$OMP+PRIVATE( L,       Q          )
!$OMP+SCHEDULE( DYNAMIC )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR)
         !==============================================================
         ! Zero variables for this level
         F           = 0d0
         F_PRIME     = 0d0
         K_RAIN      = 0d0
         Q           = 0d0

         ! Start at the top of the atmosphere
         L = LLPAR

         ! If precip forms at (I,J,L), assume it all rains out
         IF ( QQ(L,I,J) > 0d0 ) THEN

            ! Q is the new precip that is forming within grid box (I,J,L)
            Q = QQ(L,I,J)

            ! Compute K_RAIN and F' for either large-scale or convective
            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
            IF ( LS ) THEN
               K_RAIN  = LS_K_RAIN( Q )
               F_PRIME = LS_F_PRIME( Q, K_RAIN )
            ELSE
               K_RAIN  = 1.5d-3
               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
            ENDIF

            ! Set F = F', since there is no FTOP at L = LLPAR
            F = F_PRIME

            F_SAVE(I,J,L) = F

            FTOP      = F

         ENDIF

         !==============================================================
         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
         !==============================================================
         DO L = LLPAR-1, 2, -1

            ! Zero variables for each level
            F           = 0d0
            F_PRIME     = 0d0
            K_RAIN      = 0d0
            Q           = 0d0

            ! Rainout criteria
            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN

               ! Q is the new precip that is forming within grid box (I,J,L)
               Q = QQ(L,I,J)

               ! Compute K_RAIN and F' for either large-scale or convective
               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
               IF ( LS ) THEN
                  K_RAIN  = LS_K_RAIN( Q )
                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
               ELSE
                  K_RAIN  = 1.5d-3
                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
               ENDIF

               ! F is the effective area of precip seen by grid box (I,J,L)
               F = MAX( F_PRIME, FTOP )

               F_SAVE(I,J,L) = F

               ! Save FTOP for next level
               FTOP = F

            !==============================================================
            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
            !==============================================================
            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN

               ! Since no precipitation is forming within grid box (I,J,L),
               ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
               F = FTOP

               F_SAVE(I,J,L) = F

               ! Save FTOP for next level
               FTOP = F

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n
            !===========================================================
            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN

               ! No precipitation at grid box (I,J,L), thus F = 0
               F = 0d0

               F_SAVE(I,J,L) = F

               ! Save FTOP for next level
               FTOP = F
            ENDIF
         ENDDO

         !==============================================================
         ! (7)  W a s h o u t   i n   L e v e l   1
         !==============================================================

         ! We are at the surface, set L = 1
         L = 1

         ! Washout at level 1 criteria
         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN

            ! Since no precipitation is forming within grid box (I,J,L),
            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
            F = FTOP

            F_SAVE(I,J,L) = F

         ENDIF
      ENDDO
      ENDDO
#if   !defined( SGI_MIPS )
!$OMP END PARALLEL DO
#endif

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !
      ! Process rainout / washout by columns.
      !=================================================================

#if   !defined( SGI_MIPS )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,           FTOP,      ALPHA              )
!$OMP+PRIVATE( ALPHA2,  F,           F_PRIME,   GAINED,   K_RAIN   )
!$OMP+PRIVATE( LOST,    MASS_NOWASH, MASS_WASH, RAINFRAC, WASHFRAC )
!$OMP+PRIVATE( WETLOSS, L,           Q,         NN,       N        )
!$OMP+PRIVATE( QDOWN,   AER,         TMP                           )
!$OMP+PRIVATE( WETLOSS_ADJ )
!$OMP+SCHEDULE( DYNAMIC )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero accumulator array
         DO L  = 1, LLPAR
         DO NN = 1, NSOL
            DSTT_ADJ(NN,L,I,J) = 0d0
         ENDDO
         ENDDO


         ! wetdep adj (fp, dkh, 03/04/13)
         IF ( LADJ_WDEP_LS ) THEN
            BOX_DEP  (I,J,:,:) = 0d0
            LOWER_DEP(I,J,:,:) = 1d0
         ENDIF

         !==============================================================
         ! (7)  W a s h o u t   i n   L e v e l   1
         !==============================================================

         ! Zero variables for this level
         F           = 0d0
         Q           = 0d0
         QDOWN       = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0

         L = 1

         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN

            QDOWN = PDOWN(L+1,I,J)

            F = F_SAVE(I,J,L)

            IF ( F > 0d0 ) THEN

               DO NN = 1, NSOL
                  N = IDWETD(NN)

                  CALL WASHOUT( I,     J,  L, N,
     &                          QDOWN, DT, F, WASHFRAC, AER )

                  IF ( AER ) THEN
                     WETLOSS = STT_ADJ(I,J,L,N) * WASHFRAC
                     IF ( LADJ_WDEP_LS ) THEN
                        BOX_DEP(I,J,L,N) = WASHFRAC / DT
                     ENDIF
                  ELSE
                     WETLOSS = STT_ADJ(I,J,L,N) * WASHFRAC * F
                     IF ( LADJ_WDEP_LS ) THEN
                        BOX_DEP(I,J,L,N)=WASHFRAC * F / DT
                     ENDIF
                  ENDIF

                  STT_ADJ(I,J,L,N) = STT_ADJ(I,J,L,N) - WETLOSS

               ENDDO
            ENDIF
         ENDIF

         !==============================================================
         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
         !==============================================================
         DO L = 2, LLPAR-1

            ! Zero variables for each level
            ALPHA       = 0d0
            ALPHA2      = 0d0
            F           = 0d0
            F_PRIME     = 0d0
            GAINED      = 0d0
            K_RAIN      = 0d0
            LOST        = 0d0
            MASS_NOWASH = 0d0
            MASS_WASH   = 0d0
            Q           = 0d0
            QDOWN       = 0d0
            RAINFRAC    = 0d0
            WASHFRAC    = 0d0
            WETLOSS     = 0d0

            WETLOSS_ADJ = 0d0

            ! Rainout criteria
            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN

               Q = QQ(L,I,J)

               IF ( LS ) THEN
                  K_RAIN  = LS_K_RAIN( Q )
                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
               ELSE
                  K_RAIN  = 1.5d-3
                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
               ENDIF

               F = F_SAVE(I,J,L)

               IF ( F > 0d0 ) THEN

                  DO NN = 1, NSOL
                     N = IDWETD(NN)

                     CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                     ! fwd:
                     !DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
                     ! adj:
                     DSTT_ADJ(NN,L+1,I,J) = DSTT_ADJ(NN,L+1,I,J)
     &                                    + DSTT_ADJ(NN,L,I,J)
                     WETLOSS_ADJ          = WETLOSS_ADJ
     &                                    + DSTT_ADJ(NN,L,I,J)
                     DSTT_ADJ(NN,L,I,J)   = 0d0

                     ! fwd:
                     !STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
                     ! adj:
                     WETLOSS_ADJ = WETLOSS_ADJ - STT_ADJ(I,J,L,N)

                     ! fwd:
                     !WETLOSS = STT(I,J,L,N) * RAINFRAC
                     STT_ADJ(I,J,L,N) = STT_ADJ(I,J,L,N)
     &                                + RAINFRAC * WETLOSS_ADJ
                     WETLOSS_ADJ      = 0d0

                     IF ( LADJ_WDEP_LS ) THEN
                        BOX_DEP(I,J,L,N)   = RAINFRAC / DT
                        LOWER_DEP(I,J,L,N) = 1D0
                     ENDIF

                  ENDDO
               ENDIF

            !==============================================================
            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
            !==============================================================
            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN

               QDOWN = PDOWN(L,I,J)
               Q     = QQ(L,I,J)

               F = F_SAVE(I,J,L)

               IF ( F > 0d0 ) THEN

                  DO NN = 1, NSOL
                     N  = IDWETD(NN)

                     CALL WASHOUT( I,     J,  L, N,
     &                             QDOWN, DT, F, WASHFRAC, AER )


                     IF ( AER ) THEN

                        ! fwd:
                        !DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
                        ! adj:
                        DSTT_ADJ(NN,L+1,I,J) = DSTT_ADJ(NN,L+1,I,J)
     &                                       + DSTT_ADJ(NN,L,I,J)
                        WETLOSS_ADJ          = WETLOSS_ADJ
     &                                       + DSTT_ADJ(NN,L,I,J)
                        DSTT_ADJ(NN,L,I,J)   = 0d0

! Treat SO2 separately
!                        IF ( N == IDTSO2 ) THEN
!                            STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
!     &                                        + GAINED * 96D0 / 64D0
!
!                            STT(I,J,L,N)      = STT(I,J,L,N) *
!     &                                          ( 1d0 - WASHFRAC )
!                        ELSE

                            ! fwd:
                            !STT(I,J,L,N)      = STT(I,J,L,N) - WETLOSS
                            ! adj:
                            WETLOSS_ADJ = WETLOSS_ADJ - STT_ADJ(I,J,L,N)

!                        ENDIF

                        ! recalculate ALPHA
                        ALPHA = ( ABS( Q ) * BXHEIGHT(I,J,L) * 100d0 ) /
     &                            PDOWN(L+1,I,J)
                        ALPHA2 = 0.5d0 * ALPHA

                        ! fwd:
                        !GAINED = DSTT(NN,L+1,I,J) * ALPHA2
                        !WETLOSS = STT(I,J,L,N) * WASHFRAC - GAINED
                        STT_ADJ(I,J,L,N)     = STT_ADJ(I,J,L,N)
     &                                       + WASHFRAC
     &                                       * WETLOSS_ADJ
                        DSTT_ADJ(NN,L+1,I,J) = DSTT_ADJ(NN,L+1,I,J)
!----------------------------------------------------------------------
! BUG FIX:
! old:
!     &                                       - ALPHA
! new:
     &                                       - ALPHA2
!----------------------------------------------------------------------
     &                                       * WETLOSS_ADJ
                        WETLOSS_ADJ          = 0d0

                        IF ( LADJ_WDEP_LS ) THEN
                           BOX_DEP(I,J,L,N)  = WASHFRAC / DT
                           LOWER_DEP(I,J,L,N)= 1D0 - ALPHA2
                        ENDIF

                     ELSE


                        ! fwd:
                        !DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
                        ! adj:
                        DSTT_ADJ(NN,L+1,I,J) = DSTT_ADJ(NN,L+1,I,J)
     &                                       + DSTT_ADJ(NN,L,I,J)
                        WETLOSS_ADJ          = WETLOSS_ADJ
     &                                       + DSTT_ADJ(NN,L,I,J)
                        DSTT_ADJ(NN,L,I,J)   = 0d0

                        ! fwd:
                        !STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
                        WETLOSS_ADJ = WETLOSS_ADJ - STT_ADJ(I,J,L,N)

                        ! fwd:
                        !MASS_WASH = ( F*STT(I,J,L,N) ) +DSTT(NN,L+1,I,J)
                        !WETLOSS = MASS_WASH * WASHFRAC -DSTT(NN,L+1,I,J)
                        ! adj:
                        STT_ADJ(I,J,L,N)     = STT_ADJ(I,J,L,N)
     &                                       + F * WASHFRAC
     &                                       * WETLOSS_ADJ
                        DSTT_ADJ(NN,L+1,I,J) = DSTT_ADJ(NN,L+1,I,J)
     &                                       + ( WASHFRAC - 1d0 )
     &                                       * WETLOSS_ADJ
                        WETLOSS_ADJ          = 0d0

                        IF ( LADJ_WDEP_LS ) THEN
                           BOX_DEP(I,J,L,N)   = F * WASHFRAC / DT
                           LOWER_DEP(I,J,L,N) = WASHFRAC
                        ENDIF

                     ENDIF

                  ENDDO
               ENDIF

               FTOP = F

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n
            !===========================================================
            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN

               F = 0d0

               DO NN = 1, NSOL
                  N = IDWETD(NN)

                  ! fwd:
                  !DSTT(NN,L,I,J) = 0d0
                  ! adj:
                  DSTT_ADJ(NN,L,I,J) = 0d0


                  IF ( LADJ_WDEP_LS ) THEN
                     LOWER_DEP(I,J,L,N) = 0D0
                  ENDIF


                  IF ( N == IDTSO2 ) THEN
!                     WETLOSS = -DSTT(NN,L+1,I,J)
!                     STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
!     &                                 - ( WETLOSS * 96d0 / 64d0 )
                  ELSE
                     ! fwd:
                     !WETLOSS = -DSTT(NN,L+1,I,J)
                     !STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                     ! adj:
                     DSTT_ADJ(NN,L+1,I,J) = STT_ADJ(I,J,L,N)

                  ENDIF


               ENDDO

            ENDIF
         ENDDO

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR)
         !==============================================================

         ! Zero variables for this level
         F           = 0d0
         F_PRIME     = 0d0
         K_RAIN      = 0d0
         Q           = 0d0
         RAINFRAC    = 0d0

         WETLOSS_ADJ = 0d0

         ! Start at the top of the atmosphere
         L = LLPAR

         IF ( QQ(L,I,J) > 0d0 ) THEN

            Q = QQ(L,I,J)

            IF ( LS ) THEN
               K_RAIN  = LS_K_RAIN( Q )
               F_PRIME = LS_F_PRIME( Q, K_RAIN )
            ELSE
               K_RAIN  = 1.5d-3
               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
            ENDIF

            F = F_PRIME

            IF ( F > 0d0 ) THEN

               DO NN = 1, NSOL
                  N = IDWETD(NN)

                  CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                  IF ( LADJ_WDEP_LS ) THEN
                     BOX_DEP(I,J,L,N) = RAINFRAC / DT
                  ENDIF

                  WETLOSS_ADJ = 0d0

                  ! fwd:
                  !DSTT(NN,L,I,J) = WETLOSS
                  ! adj:
                  WETLOSS_ADJ        = WETLOSS_ADJ + DSTT_ADJ(NN,L,I,J)
                  DSTT_ADJ(NN,L,I,J) = 0d0

                  ! fwd:
                  !STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
                  ! adj:
                  WETLOSS_ADJ = WETLOSS_ADJ - STT_ADJ(I,J,L,N)

                  ! fwd:
                  !WETLOSS = STT(I,J,L,N) * RAINFRAC
                  ! adj:
                  STT_ADJ(I,J,L,N) = STT_ADJ(I,J,L,N) +
     &                             + WETLOSS_ADJ * RAINFRAC
                  WETLOSS_ADJ      = 0d0

               ENDDO
            ENDIF

         ENDIF

      ENDDO
      ENDDO
#if   !defined( SGI_MIPS )
!$OMP END PARALLEL DO
#endif

      ! Return to calling program
      END SUBROUTINE WETDEP_ADJ

!-----------------------------------------------------------------------------
      SUBROUTINE ADJ_SO2_WETDEP( LS )
!
!******************************************************************************
!  Subroutine ADJ_SO2_WETDEP is the wetdep adjoint for just SO2.  This
!  is treated separately owing to it being nonlinear for this species.
!  (dkh, 10/23/05)
!
!  Notes:
!  (1 ) Updated to GCv8 adjoint (dkh, 09/28/09)
!  (2 ) Remove obsolete STT_MCHK (dkh, 09/30/09)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD, NFD
      USE DAO_MOD,         ONLY : BXHEIGHT
      USE DAO_MOD,         ONLY : T
!      USE DIAG_MOD,       ONLY : AD16, AD17, AD18, CT16, CT17, CT18, AD39
      USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE TIME_MOD,        ONLY : GET_TS_DYN
      USE TRACERID_MOD,    ONLY : IDTSO2, IDTSO4
      USE TRACER_MOD,      ONLY : STT
      !USE WETSCAV_MOD,     ONLY : WASHFRAC_0, RAINFRAC_0
      USE WETSCAV_MOD,     ONLY : SO2s
      USE WETSCAV_MOD,     ONLY : H2O2s
      USE WETSCAV_MOD,     ONLY : QQ
      USE WETSCAV_MOD,     ONLY : LS_K_RAIN
      USE WETSCAV_MOD,     ONLY : LS_F_PRIME
      USE WETSCAV_MOD,     ONLY : CONV_F_PRIME
      USE WETSCAV_MOD,     ONLY : PDOWN
      USE WETSCAV_MOD,     ONLY : RAINOUT
      USE WETSCAV_MOD,     ONLY : WASHOUT
      USE WETSCAV_MOD,     ONLY : SAFETY
      USE WETSCAV_MOD,     ONLY : NSOL
      USE WETSCAV_MOD,     ONLY : NSOLMAX
      USE WETSCAV_MOD,     ONLY : WASHFRAC_FINE_AEROSOL
      USE WETSCAV_MOD,     ONLY : GET_RAINFRAC

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN) :: LS
      !REAL*8  :: IN(10,3)
      !REAL*8  :: OUT(10,3)

      ! Dummy variables used to make adjoint with TAMC
      !REAL*8 :: H2O2s(10,10,10)
      !REAL*8 :: SO2s(10,10,10)
      !REAL*8 :: QQ(10,10,10)
      !REAL*8 :: PDOWN(10,10,10)
      !REAL*8 :: STT(10,10,10,10)
      !INTEGER :: IDTSO2, LLPAR, IDTSO4
      !REAL*8              :: DSTT(10,10)

      ! Local Variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: AER

      INTEGER             :: I, IDX, J, L, N, NN

      REAL*8              :: Q,     QDOWN,  DT,        DT_OVER_TAU
      REAL*8              :: K,     K_MIN,  K_RAIN,    RAINFRAC
      REAL*8              :: F,     FTOP,   F_PRIME,   WASHFRAC
      REAL*8              :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
      REAL*8              :: ALPHA, ALPHA2, WETLOSS,   L_PLUS_W
      REAL*8              :: XDSTT, TMP

      ! DSTT is the accumulator array of rained-out
      ! soluble tracer for a given (I,J) column
      REAL*8              :: DSTT(NSOLMAX,LLPAR)

      ! Checkpointing vectors. Label them MCHK because we store
      ! them in Memory.
      REAL*8 :: F_MCHK(LLPAR)
      REAL*8 :: SO2_MCHK(LLPAR)
      REAL*8 :: SO4_MCHK(LLPAR)
      REAL*8 :: H2O2s_MCHK(LLPAR)
      REAL*8 :: SO2s_MCHK(LLPAR)
      REAL*8 :: ALPHA_MCHK(LLPAR)
      REAL*8 :: RAINFRAC_MCHK(LLPAR)
      REAL*8 :: WASHFRAC_MCHK(LLPAR)
      REAL*8 :: WASHFRAC_0(LLPAR)
      REAL*8 :: RAINFRAC_0(LLPAR)


C==============================================
C define local variables
C==============================================
      ! Replace adh2o2s and adso2s with ADJ_H2O2s and ADJ_SO2s, which
      !  are module variables.
      ! Replace adstt with ADJ_STT.
      !real*8 addstt(10,10)
      REAL*8 ADDSTT(1,LLPAR)
      real*8 adgained
      !real*8 ADJ_H2O2s(10,10,10)
      real*8 adrainfrac
      !real*8 ADJ_SO2s(10,10,10)
      !real*8 adstt(10,10,10,10)
      real*8 adwashfrac
      real*8 adwetloss

      REAL*8, PARAMETER   :: NEG_SMALL = -1.0D-10

      !=================================================================
      ! ADJ_SO2_WETDEP begins here!
      !
      ! (1)  I n i t i a l i z e   V a r i a b l e s
      !=================================================================


      !STT(10,10,:,IDTSO2) = IN(:,1)
      !H2O2s(10,10,:)      = IN(:,2)
      !SO2s(10,10,:)       = IN(:,3)

      ! Reset checkpointing arrays
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( L )
      DO L = 1, LLPAR

         WASHFRAC_0(L)    = 0d0
         RAINFRAC_0(L)    = 0d0
         F_MCHK(L)        = 0d0
         SO2_MCHK(L)      = 0d0
         SO4_MCHK(L)      = 0d0
         H2O2s_MCHK(L)    = 0d0
         SO2s_MCHK(L)     = 0d0
         ALPHA_MCHK(L)    = 0d0
         RAINFRAC_MCHK(L) = 0d0
         WASHFRAC_MCHK(L) = 0d0

      ENDDO
!$OMP END PARALLEL DO

      ! Dynamic timestep [s]
      DT = GET_TS_DYN() * 60d0

!      ! Select index for diagnostic arrays -- will archive either
!      ! large-scale or convective rainout/washout fractions
!      IF ( LS ) THEN
!         IDX = 1
!      ELSE
!         IDX = 2
!      ENDIF

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !
      ! Process rainout / washout by columns.
      !=================================================================

#if   !defined( SGI_MIPS )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,           FTOP,      DSTT,     ALPHA    )
!$OMP+PRIVATE( ALPHA2,  F,           F_PRIME,   GAINED,   K_RAIN   )
!$OMP+PRIVATE( LOST,    MASS_NOWASH, MASS_WASH, RAINFRAC, WASHFRAC )
!$OMP+PRIVATE( WETLOSS, L,           Q,         NN,       N        )
!$OMP+PRIVATE( XDSTT,   QDOWN,       AER,       TMP                )
!$OMP+PRIVATE( ADDSTT,  ADGAINED,    ADRAINFRAC, ADWASHFRAC        )
!$OMP+PRIVATE( ADWETLOSS )
!$OMP+PRIVATE( F_MCHK,  H2O2s_MCHK, SO2s_MCHK                      )
!$OMP+PRIVATE( SO2_MCHK                                            )
!$OMP+PRIVATE( ALPHA_MCHK, RAINFRAC_MCHK,  WASHFRAC_MCHK           )
!$OMP+PRIVATE( WASHFRAC_0, RAINFRAC_0                              )
!$OMP+SCHEDULE( DYNAMIC )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero arrays for this column
         FTOP      = 0d0
         DSTT(:,:) = 0d0

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR)
         !==============================================================

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0


         ! Set a few things for the entire routine that are just for
         ! SO2
         NN  = 1  ! The dimension of DSTT is (1,LLPAR)
         AER = .TRUE.

         ! Start at the top of the atmosphere
         L = LLPAR

         ! Checkpt H2O2s, SO2s, STT
         H2O2s_MCHK(L) = H2O2s(I,J,L)
         SO2s_MCHK(L)  = SO2s(I,J,L)
         SO2_MCHK(L)   = STT(I,J,L,IDTSO2)
         SO4_MCHK(L)   = STT(I,J,L,IDTSO4)

         ! If precip forms at (I,J,L), assume it all rains out
         IF ( QQ(L,I,J) > 0d0 ) THEN

            ! Q is the new precip that is forming within grid box (I,J,L)
            Q = QQ(L,I,J)

            ! Compute K_RAIN and F' for either large-scale or convective
            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
            IF ( LS ) THEN
               K_RAIN  = LS_K_RAIN( Q )
               F_PRIME = LS_F_PRIME( Q, K_RAIN )
            ELSE
               K_RAIN  = 1.5d-3
               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
            ENDIF

            ! Set F = F', since there is no FTOP at L = LLPAR
            F = F_PRIME

            ! Only compute rainout if F > 0.
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN

!               ! ND16 diagnostic...save LS and Conv fractions
!               IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!               ENDIF

!               ! ND17 diagnostic...increment counter
!               IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                  CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
!               ENDIF

               ! Loop over soluble tracers and/or aerosol tracers
!               DO NN = 1, NSOL
!                  N = IDWETD(NN)
                  N = IDTSO2

                  ! save a copy of the initial value of RAINFRAC (dkh, 10/08/09)
                  RAINFRAC_0(L) = GET_RAINFRAC( K_RAIN, F, DT )
#if   defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA2 )
                  IF ( T(I,J,L) < 258d0 ) THEN
                     RAINFRAC_0(L) = 0d0
                  ENDIF
#endif

                  ! Call subroutine RAINOUT to compute the fraction
                  ! of tracer lost to rainout in grid box (I,J,L=LLPAR)
                  CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                   ! use this for generating adjoint (dkh, 10/24/05)
!                  CALL RAINOUT_SO2(I,J,L,RAINFRAC, H2O2s, SO2s )

                  ! WETLOSS is the amount of soluble tracer
                  ! lost to rainout in grid box (I,J,L=LLPAR)
                  WETLOSS = STT(I,J,L,N) * RAINFRAC

                  ! Remove rainout losses in grid box (I,J,L=LLPAR) from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                  ! DSTT is an accumulator array for rained-out tracers.
                  ! The tracers in DSTT are in the liquid phase and will
                  ! precipitate to the levels below until a washout occurs.
                  ! Initialize DSTT at (I,J,L=LLPAR) with WETLOSS.
                  DSTT(NN,L) = WETLOSS

!                  ! ND17 diagnostic...LS and conv rainout fractions [unitless]
!                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                     AD17(I,J,L,NN,IDX) =
!     &                    AD17(I,J,L,NN,IDX) + RAINFRAC / F
!                  ENDIF

!                  ! ND39 diag - save rainout losses in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     XDSTT          = WETLOSS / DT
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
!                  ENDIF

                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN

                     IF ( STT(I,J,L,N) < NEG_SMALL ) THEN

                     CALL SAFETY( I, J, L, N, 4,
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),
     &                            ALPHA,      ALPHA2,       RAINFRAC,
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH,
     &                            WETLOSS,    GAINED,       LOST,
     &                            DSTT(NN,:), STT(I,J,:,N) )

                     END IF
                  ENDIF
               !ENDDO
           ENDIF

            ! Save FTOP for the next lower level
            FTOP = F
         ENDIF

         ! Checkpt
         F_MCHK(L)    = F
         RAINFRAC_MCHK(L) = RAINFRAC

         !==============================================================
         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
         !
         !==============================================================
         DO L = LLPAR-1, 2, -1


            ! Checkpt H2O2s, SO2s, STT
            H2O2s_MCHK(L) = H2O2s(I,J,L)
            SO2s_MCHK(L)  = SO2s(I,J,L)
            SO2_MCHK(L)   = STT(I,J,L,IDTSO2)
            SO4_MCHK(L)   = STT(I,J,L,IDTSO4)

            ! Zero variables for each level
            ALPHA       = 0d0
            ALPHA2      = 0d0
            F           = 0d0
            F_PRIME     = 0d0
            GAINED      = 0d0
            K_RAIN      = 0d0
            LOST        = 0d0
            MASS_NOWASH = 0d0
            MASS_WASH   = 0d0
            RAINFRAC    = 0d0
            WASHFRAC    = 0d0
            WETLOSS     = 0d0

            ! Rainout criteria
            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN

               ! Q is the new precip that is forming within grid box (I,J,L)
               Q = QQ(L,I,J)

               ! Compute K_RAIN and F' for either large-scale or convective
               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
               IF ( LS ) THEN
                  K_RAIN  = LS_K_RAIN( Q )
                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
               ELSE
                  K_RAIN  = 1.5d-3
                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
               ENDIF

               ! F is the effective area of precip seen by grid box (I,J,L)
               F = MAX( F_PRIME, FTOP )

               ! Only compute rainout if F > 0.
               ! This helps to eliminate unnecessary CPU cycles.
               IF ( F > 0d0 ) THEN

!                  ! ND16 diagnostic...save F
!                  IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!                  ENDIF

!                  ! ND17 diagnostic...increment counter
!                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                     CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
!                  ENDIF

                  ! Loop over soluble tracers and/or aerosol tracers
!                  DO NN = 1, NSOL
!                     N = IDWETD(NN)
                     N = IDTSO2

                     ! save a copy of the initial value of RAINFRAC (dkh, 10/08/09)
                     RAINFRAC_0(L) = GET_RAINFRAC( K_RAIN, F, DT )
#if   defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA2 )
                     IF ( T(I,J,L) < 258d0 ) THEN
                        RAINFRAC_0(L) = 0d0
                     ENDIF
#endif

                     ! Call subroutine RAINOUT to comptue the fraction
                     ! of tracer lost to rainout in grid box (I,J,L)
                     CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )

                     ! use this for generating adjoint (dkh, 10/24/05)
!                     CALL RAINOUT_SO2( I, J, L,RAINFRAC, H2O2s, SO2s )

                     ! WETLOSS is the amount of tracer in grid box
                     ! (I,J,L) that is lost to rainout.
                     WETLOSS = STT(I,J,L,N) * RAINFRAC

                     ! dkh debug
                     IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
                        print*, ' RAINFRAC at adj so2 = ', RAINFRAC
                        print*, ' K_RAIN   at adj so2 = ', K_RAIN
                        print*, ' F        at adj so2 = ', F
                        print*, ' WETLOSS = ', WETLOSS
                        print*, ' STT(FD) before = ', STT(I,J,L,NFD)
                        print*, ' DSTT(LFD)=', DSTT(NN,L)
                        print*, ' DSTT(LFD+1)=', DSTT(NN,L+1)
                     ENDIF


                     ! Subtract the rainout loss in grid box (I,J,L) from STT
                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                     ! Add to DSTT the tracer lost to rainout in grid box
                     ! (I,J,L) plus the tracer lost to rainout from grid box
                     ! (I,J,L+1), which has by now precipitated down into
                     ! grid box (I,J,L).  DSTT will continue to accumulate
                     ! rained out tracer in this manner until a washout
                     ! event occurs.
                     DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS

!                     ! ND17 diagnostic...rainout fractions [unitless]
!                     IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                        AD17(I,J,L,NN,IDX) =
!     &                       AD17(I,J,L,NN,IDX) + RAINFRAC / F
!                     ENDIF
!
!                     ! ND39 diag -- save rainout losses in [kg/s]
!                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                        XDSTT          = WETLOSS / DT
!                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
!                     ENDIF

                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 ) THEN

                        IF ( STT(I,J,L,N) < NEG_SMALL ) THEN

                        CALL SAFETY( I, J, L, N, 5,
     &                               LS,           PDOWN(L,I,J),
     &                               QQ(L,I,J),    ALPHA,
     &                               ALPHA2,       RAINFRAC,
     &                               WASHFRAC,     MASS_WASH,
     &                               MASS_NOWASH,  WETLOSS,
     &                               GAINED,       LOST,
     &                               DSTT(NN,:),   STT(I,J,:,N) )

                        ENDIF
                     ENDIF
!                  ENDDO
               ENDIF

               ! Save FTOP for next level
               FTOP = F

            !==============================================================
            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
            !
            !==============================================================
            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN

               ! QDOWN is the precip leaving thru the bottom of box (I,J,L)
               ! Q     is the new precip that is forming within box (I,J,L)
               QDOWN = PDOWN(L,I,J)
               Q     = QQ(L,I,J)

               ! Since no precipitation is forming within grid box (I,J,L),
               ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
               F = FTOP

               ! Only compute washout if F > 0.
               ! This helps to eliminate needless CPU cycles.
               IF ( F > 0d0 ) THEN

!                  ! ND16 diagnostic...save F (fraction of grid box raining)
!                  IF ( ND16 > 0d0 .and. L <= LD16 ) THEN
!                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!                  ENDIF

                  ! ND18 diagnostic...increment counter
!                  IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                     CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
!                  ENDIF
!
                  ! Loop over soluble tracers and/or aerosol tracers
!                  DO NN = 1, NSOL
!                     N  = IDWETD(NN)

                     N = IDTSO2

                     ! save a copy of the initial value of WASHFRAC (dkh, 10/08/09)
                     WASHFRAC_0(L) =
     &                WASHFRAC_FINE_AEROSOL( DT, F, QDOWN, T(I,J,L) )


                     ! Call WASHOUT to compute the fraction of
                     ! tracer lost to washout in grid box (I,J,L)
                     CALL WASHOUT( I,     J,  L, N,
     &                             QDOWN, DT, F, WASHFRAC, AER )
                     ! Use this for generating adjoint (dkh, 10/24/05)
!                     CALL WASHOUT_SO2(I, J, L, WASHFRAC, H2O2s, SO2s )


                     !=====================================================
                     ! Washout of aerosol tracers --
                     ! this is modeled as a kinetic process
                     !=====================================================
!                     IF ( AER ) THEN

                        ! ALPHA is the fraction of the raindrops that
                        ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                        ALPHA = ( ABS( Q ) * BXHEIGHT(I,J,L) * 100d0 ) /
     &                            PDOWN(L+1,I,J)

                        ! ALPHA2 is the fraction of the rained-out aerosols
                        ! that gets resuspended in grid box (I,J,L)
                        ALPHA2 = 0.5d0 * ALPHA

                        ! GAINED is the rained out aerosol coming down from
                        ! grid box (I,J,L+1) that will evaporate and re-enter
                        ! the atmosphere in the gas phase in grid box (I,J,L).
                        GAINED = DSTT(NN,L+1) * ALPHA2

                        ! Amount of aerosol lost to washout in grid box
                        ! (qli, bmy, 10/29/02)
                        WETLOSS = STT(I,J,L,N) * WASHFRAC - GAINED

                        ! Remove washout losses in grid box (I,J,L) from STT.
                        ! Add the aerosol that was reevaporated in (I,J,L).
                        ! SO2 in sulfate chemistry is wet-scavenged on the
                        ! raindrop and converted to SO4 by aqeuous chem.
                        ! If evaporation occurs then SO2 comes back as SO4
                        ! (rjp, bmy, 3/23/03)
!                        IF ( N == IDTSO2 ) THEN
                            STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
     &                                        + GAINED * 96D0 / 64D0

                            STT(I,J,L,N)      = STT(I,J,L,N)
     &                                        - ( WETLOSS + GAINED )
!                        ELSE
!                            STT(I,J,L,N)      = STT(I,J,L,N) - WETLOSS
!                        ENDIF

                        ! LOST is the rained out aerosol coming down from
                        ! grid box (I,J,L+1) that will remain in the liquid
                        ! phase in grid box (I,J,L) and will NOT re-evaporate.
                        LOST = DSTT(NN,L+1) - GAINED

                        ! Add the washed out tracer from grid box (I,J,L) to
                        ! DSTT.  Also add the amount of tracer coming down
                        ! from grid box (I,J,L+1) that does NOT re-evaporate.
                        DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS

!                        ! ND18 diagnostic...divide washout fraction by F
!                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                           AD18(I,J,L,NN,IDX) =
!     &                          AD18(I,J,L,NN,IDX) + WASHFRAC / F
!                        ENDIF

                     !=====================================================
                     ! Washout of non-aerosol tracers
                     ! This is modeled as an equilibrium process
                     !=====================================================
!                     ELSE
!
!                        ! MASS_NOWASH is the amount of non-aerosol tracer in
!                        ! grid box (I,J,L) that is NOT available for washout.
!                        MASS_NOWASH = ( 1d0 - F ) * STT(I,J,L,N)
!
!                        ! MASS_WASH is the total amount of non-aerosol tracer
!                        ! that is available for washout in grid box (I,J,L).
!                        ! It consists of the mass in the precipitating
!                        ! part of box (I,J,L), plus the previously rained-out
!                        ! tracer coming down from grid box (I,J,L+1).
!                        ! (Eq. 15, Jacob et al, 2000).
!                        MASS_WASH = ( F * STT(I,J,L,N) ) + DSTT(NN,L+1)
!
!                        ! WETLOSS is the amount of tracer mass in
!                        ! grid box (I,J,L) that is lost to washout.
!                        ! (Eq. 16, Jacob et al, 2000)
!                        WETLOSS = MASS_WASH * WASHFRAC - DSTT(NN,L+1)
!
!                        ! The tracer left in grid box (I,J,L) is what was
!                        ! in originally in the non-precipitating fraction
!                        ! of the box, plus MASS_WASH, less WETLOSS.
!                        STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!
!                        ! Add washout losses in grid box (I,J,L) to DSTT
!                        DSTT(NN,L) = DSTT(NN,L+1) + WETLOSS
!
!                        ! ND18 diagnostic...we don't have to divide the
!                        ! washout fraction by F since this is accounted for.
!                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                           AD18(I,J,L,NN,IDX) =
!     &                          AD18(I,J,L,NN,IDX) + WASHFRAC
!                        ENDIF
!                     ENDIF
!
!                     ! ND39 diag -- save rainout losses in [kg/s]
!                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                        XDSTT          = WETLOSS / DT
!                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
!                     ENDIF
!
!                     ! Negative tracer...call subroutine SAFETY
!                     IF ( STT(I,J,L,N) < 0d0 ) THEN
!                        CALL SAFETY( I, J, L, N, 6,
!     &                               LS,           PDOWN(L,I,J),
!     &                               QQ(L,I,J),    ALPHA,
!     &                               ALPHA2,       RAINFRAC,
!     &                               WASHFRAC,     MASS_WASH,
!     &                               MASS_NOWASH,  WETLOSS,
!     &                               GAINED,       LOST,
!     &                               DSTT(NN,:),   STT(I,J,:,N) )
!                     ENDIF
!                  ENDDO
               ENDIF

               ! Save FTOP for next level
               FTOP = F

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n
            !
            ! If there is no precipitation leaving grid box (I,J,L),
            ! then  set F, the effective area of precipitation in grid
            ! box (I,J,L), to zero.
            !
            ! Also, all of the previously rained-out tracer that is now
            ! coming down from grid box (I,J,L+1) will evaporate and
            ! re-enter the atmosphere in the gas phase in grid box
            ! (I,J,L).  This is called "resuspension".
            !===========================================================
            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN

               ! No precipitation at grid box (I,J,L), thus F = 0
               F = 0d0

               ! Loop over soluble tracers and/or aerosol tracers
!               DO NN = 1, NSOL
!                  N = IDWETD(NN)
                  N = IDTSO2

                  ! WETLOSS is the amount of tracer in grid box (I,J,L)
                  ! that is lost to rainout. (qli, bmy, 10/29/02)
                  WETLOSS = -DSTT(NN,L+1)

                  ! All of the rained-out tracer coming from grid box
                  ! (I,J,L+1) goes back into the gas phase at (I,J,L)
                  ! In evap, SO2 comes back as SO4 (rjp, bmy, 3/23/03)
!                  IF ( N == IDTSO2 ) THEN
                     STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
     &                                 - ( WETLOSS * 96d0 / 64d0 )
!                  ELSE
!                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!                  ENDIF

                  ! There is nothing rained out/washed out in grid box
                  ! (I,J,L), so set DSTT at grid box (I,J,L) to zero.
		  DSTT(NN,L) = 0d0

!                  ! ND39 diag -- save rainout losses in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     XDSTT          = WETLOSS / DT
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
!                  ENDIF

                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN

                     IF ( STT(I,J,L,N) < NEG_SMALL ) THEN
                     CALL SAFETY( I, J, L, N, 7,
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),
     &                            ALPHA,      ALPHA2,       RAINFRAC,
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH,
     &                            WETLOSS,    GAINED,       LOST,
     &                            DSTT(NN,:), STT(I,J,:,N) )
                     ENDIF
                  ENDIF
!               ENDDO

               ! Save FTOP for next level
               FTOP = F
            ENDIF

            ! Checkpt
            F_MCHK(L) = F
            WASHFRAC_MCHK(L) = WASHFRAC
            RAINFRAC_MCHK(L) = RAINFRAC
            ALPHA_MCHK(L)    = ALPHA

         ENDDO                  ! L

         !==============================================================
         ! (7)  W a s h o u t   i n   L e v e l   1
         !
         ! Assume all of the tracer precipitating down from grid box
         ! (I,J,L=2) to grid box (I,J,L=1) gets washed out in grid box
         ! (I,J,L=1).
         !==============================================================

         ! BUG FIX? This should be at L = 1  (dkh, 10/08/09)
         !! Checkpt H2O2s, SO2s, STT
         !H2O2s_MCHK(L) = H2O2s(I,J,L)
         !SO2s_MCHK(L)  = SO2s(I,J,L)
         !SO2_MCHK(L)   = STT(I,J,L,IDTSO2)
         !SO4_MCHK(L)   = STT(I,J,L,IDTSO4)

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0


         ! We are at the surface, set L = 1
         L = 1

         ! BUG FIX This should be here at L = 1  (dkh, 10/08/09)
         ! Checkpt H2O2s, SO2s, STT
         H2O2s_MCHK(L) = H2O2s(I,J,L)
         SO2s_MCHK(L)  = SO2s(I,J,L)
         SO2_MCHK(L)   = STT(I,J,L,IDTSO2)
         SO4_MCHK(L)   = STT(I,J,L,IDTSO4)

         ! Washout at level 1 criteria
         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN

            ! QDOWN is the precip leaving thru the bottom of box (I,J,L+1)
            QDOWN = PDOWN(L+1,I,J)

            ! Since no precipitation is forming within grid box (I,J,L),
            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
            F = FTOP

            ! Only compute washout if F > 0.
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN

!               ! ND16 diagnostic...save F
!               IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!               ENDIF
!
!               ! ND18 diagnostic...increment counter
!               IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                  CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
!               ENDIF

!               ! Loop over soluble tracers and/or aerosol tracers
!               DO NN = 1, NSOL
!                  N = IDWETD(NN)
                  N = IDTSO2

                  ! save a copy of the initial value of WASHFRAC (dkh, 10/08/09)
                  WASHFRAC_0(L) =
     &                WASHFRAC_FINE_AEROSOL( DT, F, QDOWN, T(I,J,L) )


                  ! Call WASHOUT to compute the fraction of tracer
                  ! in grid box (I,J,L) that is lost to washout.
                  CALL WASHOUT( I,     J,  L, N,
     &                          QDOWN, DT, F, WASHFRAC, AER )

                  ! use this for generating adjoint (dkh, 10/24/05)
!                  CALL WASHOUT_SO2( I,     J,  L, WASHFRAC, H2O2s, SO2s)

                  ! NOTE: for HNO3 and aerosols, there is an F factor
                  ! already present in WASHFRAC.  For other soluble
                  ! gases, we need to multiply by the F (hyl, bmy, 10/27/00)
!                  IF ( AER ) THEN
                     WETLOSS = STT(I,J,L,N) * WASHFRAC
!                  ELSE
!                     WETLOSS = STT(I,J,L,N) * WASHFRAC * F
!                  ENDIF

                  ! Subtract WETLOSS from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

!                  ! ND18 diagnostic...LS and conv washout fractions [unitless]
!                  IF ( ND18 > 0 .and. L <= LD18 ) THEN
!
!                     ! Only divide WASHFRAC by F for aerosols, since
!                     ! for non-aerosols this is already accounted for
!                     IF ( AER ) THEN
!                        TMP = WASHFRAC / F
!                     ELSE
!                        TMP = WASHFRAC
!                     ENDIF
!
!                     AD18(I,J,L,NN,IDX) = AD18(I,J,L,NN,IDX) + TMP
!                  ENDIF

!                  ! ND39 diag -- save washout loss in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     XDSTT          = WETLOSS / DT
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + XDSTT
!                  ENDIF

!                  !-----------------------------------------------------
!                  ! Dirty kludge to prevent wet deposition from removing
!                  ! stuff from stratospheric boxes -- this can cause
!                  ! negative tracer (rvm, bmy, 6/21/00)
!                  !
!                  IF ( STT(I,J,L,N) < 0d0 .and. L > 23 ) THEN
!                      WRITE ( 6, 101 ) I, J, L, N, 7
! 101                  FORMAT( 'WETDEP - STT < 0 at ', 3i4,
!     &                        ' for tracer ', i4, 'in area ', i4 )
!                      PRINT*, 'STT:', STT(I,J,:,N)
!                      STT(I,J,L,N) = 0d0
!                  ENDIF
!                  !-----------------------------------------------------
!
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     IF ( STT(I,J,L,N) < NEG_SMALL ) THEN
                     CALL SAFETY( I, J, L, N, 8,
     &                            LS,         PDOWN(L,I,J), QQ(L,I,J),
     &                            ALPHA,      ALPHA2,       RAINFRAC,
     &                            WASHFRAC,   MASS_WASH,    MASS_NOWASH,
     &                            WETLOSS,    GAINED,       LOST,
     &                            DSTT(NN,:), STT(I,J,:,N) )
                     ENDIF
                  ENDIF
!               ENDDO
            ENDIF
         ENDIF

         ! Checkpt
         F_MCHK(L)        = F
         WASHFRAC_MCHK(L) = WASHFRAC

      ! dkh debug
      IF ( LPRINTFD .and. I == IFD .and. J == JFD ) THEN
         WRITE(6,*) ' WETD CHK variables after WETDEP(T) SO2'
         print*, ' H2O2s(FD) aft so2 = ', H2O2s(IFD,JFD,LFD)
         print*, ' SO2s(FD)  aft so2 = ', SO2s(IFD,JFD,LFD)
         print*, ' SO4(FD)   aft so2 = ', STT(IFD,JFD,LFD,IDTSO4)
         print*, ' STT(FD)   aft so2 = ', STT(IFD,JFD,LFD,NFD)
      ENDIF


         !============================================================
         ! Adjoint begins here
         !============================================================

         ! Set N for adjoint of SO2
         N = IDTSO2
C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
         ! Don't reset arrays that are involved in other adjoint routines
         ! such as ADH2O2s, ADSO2s and ADSTT.
         !do ip2 = 1, 10
         !  do ip1 = 1, 10
         !    addstt(ip1,ip2) = 0.
         !  end do
         !end do
         ADDSTT(:,:) = 0d0
         adgained = 0.d0
         !do ip3 = 1, 10
         !  do ip2 = 1, 10
         !    do ip1 = 1, 10
         !      ADJ_H2O2s(ip1,ip2,ip3) = 0.
         !    end do
         !  end do
         !end do
         adrainfrac = 0.
         !do ip3 = 1, 10
         !  do ip2 = 1, 10
         !    do ip1 = 1, 10
         !      ADJ_SO2s(ip1,ip2,ip3) = 0.
         !    end do
         !  end do
         !end do
         !do ip4 = 1, 10
         !  do ip3 = 1, 10
         !    do ip2 = 1, 10
         !      do ip1 = 1, 10
         !        adstt(ip1,ip2,ip3,ip4) = 0.
         !      end do
         !    end do
         !  end do
         !end do
         adwashfrac = 0.
         adwetloss = 0.

         !============================================================
         ! Adjoint of 7
         !============================================================
         L = 1

         F                  = F_MCHK(L)
         H2O2s(I,J,L)       = H2O2s_MCHK(L)
         SO2s(I,J,L)        = SO2s_MCHK(L)
         !STT(I,J,L,IDTSO4)  = SO2_MCHK(L)

         if (pdown(l+1,i,j) .gt. 0.d0) then

            !f = ftop
            if (f .gt. 0.d0) then
               !n = idtso2
               !call washout_so2( i,j,l,washfrac,h2o2s,so2s )
               WASHFRAC = WASHFRAC_MCHK(L)
               ! Reference the working adjoint array
               !adwetloss = adwetloss-adstt(i,j,l,n)
               adwetloss = adwetloss-STT_ADJ(i,j,l,n)
               ! Reference the working adjoint array
               !adstt(i,j,l,n) = adstt(i,j,l,n)+adwetloss*washfrac
               STT_ADJ(i,j,l,n) = STT_ADJ(i,j,l,n)+adwetloss*washfrac
               ! Use original SO2 stored in memory
               !adwashfrac = adwashfrac+adwetloss*stt(i,j,l,STT2ADJ(n))
               adwashfrac = adwashfrac+adwetloss*SO2_MCHK(L)
               adwetloss = 0.d0
               !call adresto( 'memory_2_wetdep_so2_h2o2s',25,h2o2s,8,1000,1 )
               !call adresto( 'memory_2_wetdep_so2_so2s',24,so2s,8,1000,1 )
               ! Do not reset. Instead, we have checkpted the final value
               ! and pass this instead. Also, don't have to arguments that
               ! are actually module variables.
               !washfrac = 0.d0
!               call adwashout_so2( i,j,l,washfrac,h2o2s,so2s,adwashfrac,
!!     $adh2o2s,adso2s )
               ! Now make WASHFRAC_0 a local variable
               !CALL ADWASHOUT_SO2( I, J, L, WASHFRAC, ADWASHFRAC )
               CALL ADWASHOUT_SO2( I, J, L, WASHFRAC, ADWASHFRAC,
     &                             WASHFRAC_0                       )
            endif
         endif
         adwashfrac = 0.d0

         do l = 2, llpar-1
!        call adresto( 'memory_1_wetdep_so2_h2o2s',25,h2o2s,8,1000,1+
!     $llpar-1-l )
!        call adresto( 'memory_1_wetdep_so2_so2s',24,so2s,8,1000,1+llpar-
!     $1-l )
!        call adresto( 'memory_1_wetdep_so2_stt',23,stt,8,10000,1+llpar-
!     $1-l )
!        call adresto( 'memory_1_wetdep_so2_ftop',24,ftop,8,1,1+llpar-1-
!     $l )
            !alpha = 0.d0
            !f = 0.d0
            !rainfrac = 0.d0
            !washfrac = 0.d0

            F                  = F_MCHK(L)
            H2O2s(I,J,L)       = H2O2s_MCHK(L)
            SO2s(I,J,L)        = SO2s_MCHK(L)
            !STT(I,J,L,IDTSO4)  = STT_MCHK(L)

            ! just to be safe
            ADRAINFRAC         = 0d0
            ADWETLOSS          = 0d0
            ADGAINED           = 0d0

            !============================================================
            ! Adjoint of 4
            !============================================================
            if (pdown(l,i,j) .gt. 0.d0 .and. qq(l,i,j) .gt. 0.d0) then
              if (f .gt. 0.d0) then
                !n = idtso2
                !call rainout_so2( i,j,l,rainfrac,h2o2s,so2s )
                RAINFRAC = RAINFRAC_MCHK(L)
                adwetloss = adwetloss+addstt(nn,l)
                addstt(nn,l+1) = addstt(nn,l+1)+addstt(nn,l)
                addstt(nn,l) = 0.
                ! Reference the working adjoint array
                !adwetloss = adwetloss-adstt(i,j,l,n)
                adwetloss = adwetloss-STT_ADJ(i,j,l,n)
                ! Use original SO2 stored in memory
                !adrainfrac = adrainfrac+adwetloss*stt(i,j,l,STT2ADJ(n))
                adrainfrac = adrainfrac+adwetloss*SO2_MCHK(L)
               ! Reference the working adjoint array
                !adstt(i,j,l,n) = adstt(i,j,l,n)+adwetloss*rainfrac
                ! dkh debug
                IF ( I == IFD .and. J == JFD .and. L == LFD
     &                    .and. N == NFD .and. LPRINTFD ) THEN
                   print*, ' STT_ADJ(FD) before = ', STT_ADJ(I,J,L,N)
                   print*, ' ADWETLOSS = ', ADWETLOSS
                   print*, ' RAINFRAC = ', RAINFRAC
                   print*, ' H2O2s = ', H2O2s(I,J,L)
                   print*, ' SO2s = ', SO2s(I,J,L)
                   print*, ' STT(FD) =', STT(IFD,JFD,LFD,NFD)
                   print*, ' adrainfrac = ', adrainfrac
                   print*, ' SO2_MCHK   = ', SO2_MCHK(L)
                ENDIF
                STT_ADJ(i,j,l,n) = STT_ADJ(i,j,l,n)+adwetloss*rainfrac
                adwetloss = 0.
                IF ( I == IFD .and. J == JFD .and. L == LFD
     &               .and. N == NFD .and. LPRINTFD ) THEN
                   print*, ' STT_ADJ(FD) after = ', STT_ADJ(I,J,L,N)
                ENDIF
!                call adresto( 'memory_1_wetdep_so2_h2o2s',25,h2o2s,8,1000,1+
!     $llpar-1-l )
!                call adresto( 'memory_1_wetdep_so2_so2s',24,so2s,8,1000,1+
!     $llpar-1-l )
                ! Do not reset. Instead, we have checkpted the final value
                ! of RAINFRAC and pass this instead. Also, don't have to arguments that
                ! are actually module variables.
                !rainfrac = 0.d0
!                call adrainout_so2( i,j,l,rainfrac,h2o2s,so2s,adrainfrac,
!     $adh2o2ss,adso2s )
                ! Now make RAINFRAC_0 a local variable
                !CALL ADRAINOUT_SO2( I, J, L, RAINFRAC, ADRAINFRAC )
                CALL ADRAINOUT_SO2( I, J, L, RAINFRAC, ADRAINFRAC,
     &             RAINFRAC_0 )
                IF ( I == IFD .and. J == JFD .and. L == LFD
     &                    .and. N == NFD .and. LPRINTFD ) THEN
                 print*, ' SO2s_ADJ(FD) = ', SO2s_ADJ(I,J,L)
                ENDIF
              endif
            !============================================================
            ! Adjoint of 5
            !============================================================
            else if (pdown(l,i,j) .gt. 0.d0 .and. qq(l,i,j) .le. 0.d0)
     &        then
              !f = ftop
              if (f .gt. 0.d0) then
                !n = idtso2
                !call washout_so2( i,j,l,washfrac,h2o2s,so2s )
                WASHFRAC = WASHFRAC_MCHK(L)
                ALPHA    = ALPHA_MCHK(L)
                alpha2 = 0.5d0*alpha
                adwetloss = adwetloss+addstt(nn,l)
                addstt(nn,l+1) = addstt(nn,l+1)+addstt(nn,l)
                addstt(nn,l) = 0.
                ! Reference the working adjoint array
                !adgained = adgained-adstt(i,j,l,n)
                !adwetloss = adwetloss-adstt(i,j,l,n)
                !adgained = adgained+1.5d0*adstt(i,j,l,idtso4)
                adgained = adgained-STT_ADJ(i,j,l,n)
                adwetloss = adwetloss-STT_ADJ(i,j,l,n)
                adgained = adgained+1.5d0*STT_ADJ(i,j,l,IDTSO4)
                adgained = adgained-adwetloss
                ! Reference the working adjoint array
                !adstt(i,j,l,n) = adstt(i,j,l,n)+adwetloss*washfrac
                STT_ADJ(i,j,l,n) = STT_ADJ(i,j,l,n)+adwetloss*washfrac
                ! Use the original SO2 stored in memory
                !adwashfrac = adwashfrac+adwetloss*stt(i,j,l,STT2ADJ(n))
                adwashfrac = adwashfrac+adwetloss*SO2_MCHK(L)
                adwetloss = 0.
                addstt(nn,l+1) = addstt(nn,l+1)+adgained*alpha2
                adgained = 0.
!            call adresto( 'memory_1_wetdep_so2_h2o2s',25,h2o2s,8,1000,1+
!     $llpar-1-l )
!            call adresto( 'memory_1_wetdep_so2_so2s',24,so2s,8,1000,1+
!     $llpar-1-l )
                ! Do not reset. Instead, we have checkpted the final value
                ! and pass this instead. Also, don't have to arguments that
                ! are actually module variables.
                !washfrac = 0.d0
!                call adwashout_so2( i,j,l,washfrac,h2o2s,so2s,adwashfrac,
!     $adh2o2ss,adso2s )
                ! Now make WASHFRAC_0 a local variable
                !CALL ADWASHOUT_SO2( I, J, L, WASHFRAC, ADWASHFRAC )
                CALL ADWASHOUT_SO2( I, J, L, WASHFRAC, ADWASHFRAC,
     &                              WASHFRAC_0                     )
              endif
            !============================================================
            ! Adjoint of 6
            !============================================================
            else if (abs(pdown(l,i,j)) .lt. 1.d-30) then
              addstt(nn,l) = 0.
              ! Reference the working adjoint array
              !adwetloss = adwetloss-1.5d0*adstt(i,j,l,idtso4)
              adwetloss = adwetloss-1.5d0*STT_ADJ(i,j,l,IDTSO4)
              addstt(nn,l+1) = addstt(nn,l+1)-adwetloss
              adwetloss = 0.
            endif
            adwashfrac = 0.
            adrainfrac = 0.
         end do  ! L

         !============================================================
         ! Adjoint of 3
         !============================================================
         l = llpar

         F                  = F_MCHK(L)
         H2O2s(I,J,L)       = H2O2s_MCHK(L)
         SO2s(I,J,L)        = SO2s_MCHK(L)
         !STT(I,J,L,IDTSO4)  = STT_MCHK(L)

         if (qq(l,i,j) .gt. 0.d0) then
           if (f .gt. 0.d0) then
             adwetloss = adwetloss+addstt(nn,l)
             addstt(nn,l) = 0.
             ! Reference the working adjoint array
             !adwetloss = adwetloss-adstt(i,j,l,n)
             adwetloss = adwetloss-STT_ADJ(i,j,l,n)
!             stt(10,10,:,idtso2) = in(:,1)
!             h2o2s(10,10,:) = in(:,2)
!             so2s(10,10,:) = in(:,3)
             ! Use original SO2 stored in memory
             !adrainfrac = adrainfrac+adwetloss*stt(i,j,l,STT2ADJ(n))
             adrainfrac = adrainfrac+adwetloss*SO2_MCHK(L)
             ! Why didn't TAMC have a call to rainout_so2 here?
             RAINFRAC = RAINFRAC_MCHK(L)
             ! Reference the working adjoint array
             !adstt(i,j,l,n) = adstt(i,j,l,n)+adwetloss*rainfrac
             STT_ADJ(i,j,l,n) = STT_ADJ(i,j,l,n)+adwetloss*rainfrac
             adwetloss = 0.
!             h2o2s(10,10,:) = in(:,2)
!             so2s(10,10,:) = in(:,3)
             ! Do not reset. Instead, we have checkpted the final value
             ! of RAINFRAC and pass this instead. Also, don't have to arguments that
             ! are actually module variables.
             !rainfrac = 0.d0
!             call adrainout_so2( i,j,l,rainfrac,h2o2s,so2s,adrainfrac,
!     $ADJ_H2O2s,ADJ_SO2s )
             ! Now make RAINFRAC_0 a local variable (dkh, 10/08/09)
             !CALL ADRAINOUT_SO2(I, J, L, RAINFRAC, ADRAINFRAC )
             CALL ADRAINOUT_SO2(I, J, L, RAINFRAC, ADRAINFRAC,
     &                          RAINFRAC_0                     )
           endif
         endif
         !dkh debug
         ADRAINFRAC = 0d0

      ENDDO
      ENDDO
#if   !defined( SGI_MIPS )
!$OMP END PARALLEL DO
#endif


!       OUT(:,1) = STT(10,10,:,IDTSO2)
!       OUT(:,2) = H2O2s(10,10,:)
!       OUT(:,3) = SO2s(10,10,:)

      ! Return to calling program
      END SUBROUTINE ADJ_SO2_WETDEP

!------------------------------------------------------------------------------
C***************************************************************
C                           DISCLAIMER
C
C   This file was generated by TAMC version 5.3.2
C
C   THE AUTHOR DOES NOT MAKE  ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR  RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS,
C   OR USEFULNESS  OF ANY INFORMATION OR PROCESS  DISCLOSED, OR REPRESENTS
C   THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C   THIS CODE IS FOR NON-PROFIT-ORIENTED ACADEMIC RESEARCH AND EDUCATION
C   ONLY.  ANY COMMERCIAL OR OTHER PROFIT-ORIENTED USE OR  EVALUATION IS
C   STRICTLY  FORBIDDEN.  PASSING  THIS CODE  TO  ANY THIRD  PARTY IS  NOT
C   ALLOWED.
C
C   FOR COMMERCIAL OR  OTHER PROFIT-ORIENTED APPLICATIONS PLEASE CONTACT
C   info@FastOpt.com
C
      ! We don't actually need to pass H2O2s and SO2s  or their
      ! adjoints as arguments because they're module variables. .
!      subroutine adrainout_so2( i, j, l, rainfrac, h2o2s, so2s,
!     $adrainfrac, adjh2o2s, adso2s )
      SUBROUTINE ADRAINOUT_SO2( I,          J,       L,
     &                          RAINFRAC,   ADRAINFRAC,
     &                          RAINFRAC_0               )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
! Modifications in CAPS by dkh.
!
! Notes:
! (1 ) Updated to GCv8 (dkh, 09/29/09)
!
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE WETSCAV_MOD,     ONLY : SO2s
      USE WETSCAV_MOD,     ONLY : H2O2s
      ! Now make this a local variable
      !USE WETSCAV_MOD,     ONLY : RAINFRAC_0
      USE WETSCAV_MOD,     ONLY : EPSILON


      implicit none
#     include "CMN_SIZE"   ! LLPAR
C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      ! make this a module var
      !real*8 adh2o2s(10,10,10)
      real*8 adrainfrac
      ! make this a module var
      !real*8 adso2s(10,10,10)
      ! This is a module variable
      !real*8 h2o2s(10,10,10)
      integer i
      integer j
      integer l
      real*8 rainfrac
      ! This is a module variable
      !real*8 so2s(10,10,10)

      REAL*8 RAINFRAC_0(LLPAR)

C==============================================
C define local variables
C==============================================
      real*8 adso2loss
      ! EPSILON is actually a module variable
      !real*8 epsilon
      !real*8 h2o2sh(10,10,10)
      REAL*8 H2O2sh
      !integer ip1
      !integer ip2
      !integer ip3
      real*8 rainfrach
      real*8 so2loss
      !real*8 so2sh(10,10,10)
      REAL*8 SO2sh

C----------------------------------------------
C SAVE ARGUMENTS
C----------------------------------------------
      ! Only save local scalar values
      !do ip3 = 1, 10
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      so2sh(ip1,ip2,ip3) = so2s(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      ! dkh debug
      !SO2s(I,J,L) = 1.5995817D-9
      !H2O2s(I,J,L) = 6.163918D-10

      SO2sh = SO2s(I,J,L)
      rainfrach = rainfrac
      !do ip3 = 1, 10
      ! Only save local scalar values
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      h2o2sh(ip1,ip2,ip3) = h2o2s(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      H2O2sh = H2O2s(I,J,L)

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adso2loss = 0.
C ROUTINE BODY
C----------------------------------------------
      if (so2s(i,j,l) .gt. epsilon) then
        so2loss = min(so2s(i,j,l),h2o2s(i,j,l))
        ! Rather than recompute RAINFRAC, which would mean recomputing K_RAIN,
        ! we'll just use the checkpointed value.  So RAINFRAC is actually
        ! an input of adrainout_so2.
        !rainfrac = so2loss*rainfrac/so2s(i,j,l)
        !rainfrac = max(rainfrac,0.d0)
      else
        !rainfrac = 0.d0
      endif
      so2s(i,j,l) = so2s(i,j,l)*(1.d0-rainfrac)
      SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)
     &                  *(0.5d0+sign(0.5,so2s(i,j,l)-epsilon))

      ! dkh debug
      IF ( LPRINTFD .and. I == IFD .and. J == JFD .and. L == LFD )
     &    print*, ' SO2s_ADJ 0: ', SO2s_ADJ(I,J,L),
     &            ' so2s (final) = ', so2s(I,J,L),
     &            ' adrainfrac (initial) = ', adrainfrac

      ! Only hold the local scalar value
      !do ip3 = 1, 10
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      so2s(ip1,ip2,ip3) = so2sh(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      SO2s(I,J,L) = SO2sh
      adrainfrac = adrainfrac-SO2s_ADJ(i,j,l)*so2s(i,j,l)
      SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)*(1.d0-rainfrac)
      IF ( LPRINTFD .and. I == IFD .and. J == JFD .and. L == LFD )
     &    print*, ' SO2s_ADJ 1: ', SO2s_ADJ(I,J,L),
     &            ' adrainfrac (1) = ', adrainfrac
      if (so2s(i,j,l) .gt. epsilon) then
        h2o2s(i,j,l) = h2o2s(i,j,l)-so2s(i,j,l)*rainfrac
        H2O2s_ADJ(i,j,l) = H2O2s_ADJ(i,j,l)*(0.5d0+sign(0.5,h2o2s(i,j,l)
     $-epsilon))
        adrainfrac = adrainfrac-H2O2s_ADJ(i,j,l)*so2s(i,j,l)
        SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)-H2O2s_ADJ(i,j,l)*rainfrac
      IF ( LPRINTFD .and.  I == IFD .and. J == JFD .and. L == LFD )
     &    print*, ' SO2s_ADJ 2: ', SO2s_ADJ(I,J,L),
     &            ' adrainfrac (2) = ', adrainfrac
        ! At this point rainfrac is already at its final value
        ! so don't need to recalculate it.  Also, rainfrach IS
        ! the final value.
        !rainfrac = rainfrach
        !rainfrac = so2loss*rainfrac/so2s(i,j,l)
        adrainfrac = adrainfrac*(0.5d0+sign(0.5d0,rainfrac-0.d0))
        ! Here we need the original RAINFRAC, so restore using RAINFRAC_0
        ! before rescaling because we didn't checkpt that.
        !rainfrac = rainfrach
        RAINFRAC = RAINFRAC_0(L)
        adso2loss = adso2loss+adrainfrac*(rainfrac/so2s(i,j,l))
        SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)-adrainfrac*(so2loss*rainfrac/
     $(so2s(i,j,l)*so2s(i,j,l)))

        ! dkh debug
      IF ( LPRINTFD .and.  I == IFD .and. J == JFD .and. L == LFD )
     &    print*, ' SO2s_ADJ 3: ', SO2s_ADJ(I,J,L),
     &            ' adrainfrac = ', adrainfrac,
     &            ' so2loss = ', so2loss,
     &            ' rainfrac = ', rainfrac,
     &            ' so2s^2   = ', so2s(I,J,L)*so2s(I,J,L),
     &            ' so2s   = ', so2s(I,J,L)

        adrainfrac = adrainfrac*(so2loss/so2s(i,j,l))
        ! Only hold the local scalar value
        !do ip3 = 1, 10
        !  do ip2 = 1, 10
        !    do ip1 = 1, 10
        !      h2o2s(ip1,ip2,ip3) = h2o2sh(ip1,ip2,ip3)
        !    end do
        !  end do
        !end do
        H2O2s(I,J,L) = H2O2sh
        H2O2s_ADJ(i,j,l) =
     &H2O2s_ADJ(i,j,l)+adso2loss*(0.5d0-sign(0.5d0,h2o2s(i,
     $j,l)-so2s(i,j,l)))
        SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)
     &+adso2loss*(0.5d0+sign(0.5d0,h2o2s(i,j,
     $l)-so2s(i,j,l)))

      ! dkh debug
      IF ( LPRINTFD .and.  I == IFD .and. J == JFD .and. L == LFD )
     &    print*, ' SO2s_ADJ 4: ', SO2s_ADJ(I,J,L),
     &            ' adso2loss = ', adso2loss,
     &            ' so2s ', so2s(I,J,L),
     &            ' h2o2s ', h2o2s(I,J,L),
     &            ' adrainfrac ', adrainfrac

        adso2loss = 0.d0
      else
        adrainfrac = 0.d0
      endif

      !end
      END SUBROUTINE ADRAINOUT_SO2

C***************************************************************

      ! We don't actually need to pass H2O2s and SO2s  or their
      ! adjoints as arguments because they're module variables. .
!      subroutine adwashout_so2( i, j, l, washfrac, h2o2s, so2s,
!     $adwashfrac, adh2o2s, adso2s )
      SUBROUTINE ADWASHOUT_SO2( I,        J,    L,
     &                          WASHFRAC, ADWASHFRAC,
     &                          WASHFRAC_0             )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
! Modifications in CAPS by dkh.
!
! Notes:
! (1 ) Updated to GCv8 (dkh, 09/29/09)
!
C***************************************************************
C***************************************************************
      ! References to f90 mofules
      USE DAO_MOD,     ONLY : T
      USE WETSCAV_MOD, ONLY : H2O2s
      USE WETSCAV_MOD, ONLY : SO2s
      ! Now make this a local variable
      !USE WETSCAV_MOD, ONLY : WASHFRAC_0
      USE WETSCAV_MOD, ONLY : EPSILON


C==============================================
C all entries are defined explicitly
C==============================================
      implicit none

#     include "CMN_SIZE"    ! LLPAR

C==============================================
C define common blocks
C==============================================
C==============================================
C define arguments
C==============================================
      ! Comment out those that are module variables.
      !real*8 adh2o2s(10,10,10)
      !real*8 adso2s(10,10,10)
      real*8 adwashfrac
      !real*8 h2o2s(10,10,10)
      integer i
      integer j
      integer l
      !real*8 so2s(10,10,10)
      real*8 washfrac

      REAL*8 WASHFRAC_0(LLPAR)

C==============================================
C define local variables
C==============================================
      real*8 adso2loss
      ! Epsilon is actually a module variable
      !real*8 epsilon
      ! Only save scalar of local value
      !real*8 h2o2sh(10,10,10)
      REAL*8 :: H2O2sh
      ! Don't need TAMC's dummy indicies
      !integer ip1
      !integer ip2
      !integer ip3
      real*8 so2loss
      ! Only save scalar of local value
      !real*8 so2sh(10,10,10)
      REAL*8 :: SO2sh
      real*8 tk
      real*8 washfrach

C----------------------------------------------
C SAVE ARGUMENTS
C----------------------------------------------
      washfrach = washfrac
      !do ip3 = 1, 10
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      so2sh(ip1,ip2,ip3) = so2s(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      SO2sh = SO2s(I,J,L)
      !do ip3 = 1, 10
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      h2o2sh(ip1,ip2,ip3) = h2o2s(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      H2O2sh = H2O2s(I,J,L)

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adso2loss = 0.d0

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
      ! TK is Kelvin temperature
      TK = T(I,J,L)

      if (tk .ge. 268.d0 .and. so2s(i,j,l) .gt. epsilon) then
        so2loss = min(so2s(i,j,l),h2o2s(i,j,l))
        ! Rather than recalculate WASHFRAC, which would involve
        ! first recalculating WASHFRAC_AEROSOL,  we pass the final checkpted
        ! value of WASHFRAC to this routine.
        !washfrac = so2loss*washfrac/so2s(i,j,l)
        !washfrac = max(washfrac,0.d0)
      else
        !washfrac = 0.d0
      endif
      so2s(i,j,l) = so2s(i,j,l)*(1.d0-washfrac)
      SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)
     &                  *(0.5d0+sign(0.5d0,so2s(i,j,l)-epsilon))
      !do ip3 = 1, 10
      !  do ip2 = 1, 10
      !    do ip1 = 1, 10
      !      so2s(ip1,ip2,ip3) = so2sh(ip1,ip2,ip3)
      !    end do
      !  end do
      !end do
      SO2s(I,J,L) = SO2sh
      adwashfrac = adwashfrac-SO2s_ADJ(i,j,l)*so2s(i,j,l)
      SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)*(1.d0-washfrac)
      if (tk .ge. 268.d0 .and. so2s(i,j,l) .gt. epsilon) then
        !!! Upgrade to GCv8, no longer use this condition (dkh, 10/07/09)
        !!!if (so2s(i,j,l) .lt. h2o2s(i,j,l)) then
          !!! Upgrade to GCv8, now multiply so2s by final washfrac
          !!!h2o2s(i,j,l) = h2o2s(i,j,l)-so2s(i,j,l)
          h2o2s(i,j,l) = h2o2s(i,j,l)-so2s(i,j,l) * WASHFRAC
          H2O2s_ADJ(i,j,l) =
     &H2O2s_ADJ(i,j,l)*(0.5d0+sign(0.5d0,h2o2s(i,j,l)-
     $epsilon))
          !!! Upgrade to GCv8, now multiply H2O2s_ADJ by WASHFRAC
          !!!SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)-H2O2s_ADJ(i,j,l)
          SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)-H2O2s_ADJ(i,j,l) * WASHFRAC

          !!! and add this line:
          ADWASHFRAC = ADWASHFRAC-H2O2s_ADJ(i,j,l)*SO2s(i,j,l)

        !!!else
        !!!  H2O2s_ADJ(i,j,l) = 0.d0
        !!!endif
        ! We don't need to recalculate the final WASHFRAC because
        ! that value has been chekpted.
        !washfrac = washfrach
        !washfrac = so2loss*washfrac/so2s(i,j,l)
        adwashfrac = adwashfrac*(0.5d0+sign(0.5d0,washfrac-0.d0))
        ! Here we need the original WASHFRAC, so restore using WASHFRAC_0
        !washfrac = washfrach
        WASHFRAC = WASHFRAC_0(L)
        adso2loss = adso2loss+adwashfrac*(washfrac/so2s(i,j,l))
        SO2s_ADJ(i,j,l) = SO2s_ADJ(i,j,l)-adwashfrac*(so2loss*washfrac/
     $(so2s(i,j,l)*so2s(i,j,l)))
        adwashfrac = adwashfrac*(so2loss/so2s(i,j,l))
        !do ip3 = 1, 10
        !  do ip2 = 1, 10
        !    do ip1 = 1, 10
        !      h2o2s(ip1,ip2,ip3) = h2o2sh(ip1,ip2,ip3)
        !    end do
        !  end do
        !end do
        H2O2s(I,J,L) = H2O2sh
        H2O2s_ADJ(i,j,l) = H2O2s_ADJ(i,j,l)
     &+adso2loss*(0.5d0-sign(0.5d0,h2o2s(i,
     $j,l)-so2s(i,j,l)))
        SO2s_ADJ(i,j,l) =
     &SO2s_ADJ(i,j,l)+adso2loss*(0.5d0+sign(0.5d0,h2o2s(i,j,
     $l)-so2s(i,j,l)))
        adso2loss = 0.d0
      else
        adwashfrac = 0.d0
      endif

      !end
      END SUBROUTINE ADWASHOUT_SO2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!
!      SUBROUTINE RECALC_SOX_WETDEP( LS )
!!
!!******************************************************************************
!! !  Subroutine RECALC_SOX_WETDEP is a copy of WETDEP that only acts on
!!  SO2 and SO4. See WETDEP for notes. We only need to recompute values
!!  of these two species in order to do adjoint of wetdep for convective precip
!!  of SO2.  (dkh, 10/23/05)
!!
!! NOTEs
!! (1 ) Call GET_NN which returns the NN values for SO2 and SO4. (dkh, 10/23/05)
!! (2 ) Updated to v8-02-01, adj_group (dkh, ks, mak, cs  06/08/09)
!!
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE DAO_MOD,           ONLY : BXHEIGHT
!      USE DIAG_MOD,          ONLY : AD16, AD17, AD18
!      USE DIAG_MOD,          ONLY : CT16, CT17, CT18, AD39
!      USE ERROR_MOD,         ONLY : GEOS_CHEM_STOP, IT_IS_NAN
!      USE LOGICAL_MOD,       ONLY : LDYNOCEAN
!      USE OCEAN_MERCURY_MOD, ONLY : ADD_Hg2_WD
!      USE TIME_MOD,          ONLY : GET_TS_DYN
!      USE TRACER_MOD,        ONLY : ITS_A_MERCURY_SIM, STT
!      USE TRACERID_MOD,      ONLY : IDTSO2, IDTSO4, IS_Hg2
!
!      IMPLICIT NONE
!
!#     include "CMN_SIZE"  ! Size parameters
!#     include "CMN_DIAG"  ! Diagnostic arrays and switches
!
!      ! Arguments
!      LOGICAL, INTENT(IN) :: LS
!
!      ! Local Variables
!      LOGICAL, SAVE       :: FIRST = .TRUE.
!      LOGICAL             :: IS_Hg
!      LOGICAL             :: AER
!
!      INTEGER             :: I, IDX, J, L, N, NN
!
!      REAL*8              :: Q,     QDOWN,  DT,        DT_OVER_TAU
!      REAL*8              :: K,     K_MIN,  K_RAIN,    RAINFRAC
!      REAL*8              :: F,     FTOP,   F_PRIME,   WASHFRAC
!      REAL*8              :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
!      REAL*8              :: ALPHA, ALPHA2, WETLOSS,   TMP
!
!      ! DSTT is the accumulator array of rained-out
!      ! soluble tracer for a given (I,J) column
!      REAL*8              :: DSTT(NSOL,LLPAR,IIPAR,JJPAR)
!
!      !=================================================================
!      ! WETDEP begins here!
!      !
!      ! (1)  I n i t i a l i z e   V a r i a b l e s
!      !=================================================================
!
!      ! Is this a mercury simulation with dynamic online ocean?
!      IS_Hg = ( ITS_A_MERCURY_SIM() .and. LDYNOCEAN )
!
!      ! Dynamic timestep [s]
!      DT    = GET_TS_DYN() * 60d0
!
!      ! Select index for diagnostic arrays -- will archive either
!      ! large-scale or convective rainout/washout fractions
!      IF ( LS ) THEN
!         IDX = 1
!      ELSE
!         IDX = 2
!      ENDIF
!
!      !=================================================================
!      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
!      !
!      ! Process rainout / washout by columns.
!      !=================================================================
!
!#if   !defined( SGI_MIPS )
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I,       J,           FTOP,      ALPHA              )
!!$OMP+PRIVATE( ALPHA2,  F,           F_PRIME,   GAINED,   K_RAIN   )
!!$OMP+PRIVATE( LOST,    MASS_NOWASH, MASS_WASH, RAINFRAC, WASHFRAC )
!!$OMP+PRIVATE( WETLOSS, L,           Q,         NN,       N        )
!!$OMP+PRIVATE( QDOWN,   AER,         TMP                           )
!!$OMP+SCHEDULE( DYNAMIC )
!#endif
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Zero FTOP
!         FTOP = 0d0
!
!         ! Zero accumulator array
!         DO L  = 1, LLPAR
!         DO NN = 1, NSOL
!            DSTT(NN,L,I,J) = 0d0
!         ENDDO
!         ENDDO
!
!         !==============================================================
!         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR)
!         !
!         ! Assume that rainout is happening in the top layer if
!         ! QQ(LLPAR,I,J) > 0.  In other words, if any precipitation
!         ! forms in grid box (I,J,LLPAR), assume that all of it falls
!         ! down to lower levels.
!         !
!         ! Soluble gases/aerosols are incorporated into the raindrops
!         ! and are completely removed from grid box (I,J,LLPAR).  There
!         ! is no evaporation and "resuspension" of aerosols during a
!         ! rainout event.
!         !
!         ! For large-scale (a.k.a. stratiform) precipitation, the first
!         ! order rate constant for rainout in the grid box (I,J,L=LLPAR)
!         ! (cf. Eq. 12, Jacob et al, 2000) is given by:
!         !
!         !                        Q
!         !    K_RAIN = K_MIN + -------    [units: s^-1]
!         !                      L + W
!         !
!         ! and the areal fraction of grid box (I,J,L=LLPAR) that
!         ! is actually experiencing large-scale precipitation
!         ! (cf. Eq. 11, Jacob et al, 2000) is given by:
!         !
!         !                  Q
!         !    F' =  -------------------   [unitless]
!         !           K_RAIN * ( L + W )
!         !
!         ! Where:
!         !
!         !    K_MIN  = minimum value for K_RAIN
!         !           = 1.0e-4 [s^-1]
!         !
!         !    L + W  = condensed water content in cloud
!         !           = 1.5e-6 [cm3 H2O/cm3 air]
!         !
!         !    Q = QQ = rate of precipitation formation
!         !             [ cm3 H2O / cm3 air / s ]
!         !
!         ! For convective precipitation, K_RAIN = 5.0e-3 [s^-1], and the
!         ! expression for F' (cf. Eq. 13, Jacob et al, 2000) becomes:
!         !
!         !                                  { DT        }
!         !                    FMAX * Q * MIN{ --- , 1.0 }
!         !                                  { TAU       }
!         !  F' = ------------------------------------------------------
!         !               { DT        }
!         !        Q * MIN{ --- , 1.0 }  +  FMAX * K_RAIN * ( L + W )
!         !               { TAU       }
!         !
!         ! Where:
!         !
!         !    Q = QQ = rate of precipitation formation
!         !             [cm3 H2O/cm3 air/s]
!         !
!         !    FMAX   = maximum value for F'
!         !           = 0.3
!         !
!         !    DT     = dynamic time step from the CTM [s]
!         !
!         !    TAU    = duration of rainout event
!         !           = 1800 s (30 min)
!         !
!         !    L + W  = condensed water content in cloud
!         !           = 2.0e-6 [cm3 H2O/cm3 air]
!         !
!         ! K_RAIN and F' are needed to compute the fraction of tracer
!         ! in grid box (I,J,L=LLPAR) lost to rainout.  This is done in
!         ! module routine RAINOUT.
!         !==============================================================
!
!         ! Zero variables for this level
!         ALPHA       = 0d0
!         ALPHA2      = 0d0
!         F           = 0d0
!         F_PRIME     = 0d0
!         GAINED      = 0d0
!         K_RAIN      = 0d0
!         LOST        = 0d0
!         Q           = 0d0
!         QDOWN       = 0d0
!         MASS_NOWASH = 0d0
!         MASS_WASH   = 0d0
!         RAINFRAC    = 0d0
!         WASHFRAC    = 0d0
!         WETLOSS     = 0d0
!
!         ! Start at the top of the atmosphere
!         L = LLPAR
!
!         ! If precip forms at (I,J,L), assume it all rains out
!         IF ( QQ(L,I,J) > 0d0 ) THEN
!
!            ! Q is the new precip that is forming within grid box (I,J,L)
!            Q = QQ(L,I,J)
!
!            ! Compute K_RAIN and F' for either large-scale or convective
!            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
!            IF ( LS ) THEN
!               K_RAIN  = LS_K_RAIN( Q )
!               F_PRIME = LS_F_PRIME( Q, K_RAIN )
!            ELSE
!               K_RAIN  = 1.5d-3
!               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
!            ENDIF
!
!            ! Set F = F', since there is no FTOP at L = LLPAR
!            F = F_PRIME
!
!            ! Only compute rainout if F > 0.
!            ! This helps to eliminate unnecessary CPU cycles.
!            IF ( F > 0d0 ) THEN
!
!               ! ND16 diagnostic...save LS and Conv fractions
!               IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!               ENDIF
!
!               ! ND17 diagnostic...increment counter
!               IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                  CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
!               ENDIF
!
!               ! Loop over soluble tracers and/or aerosol tracers
!               !DO NN = 1, NSOL
!               ! Loop over just SO2 then SO4. (dkh, 10/23/05)
!               DO NNN = 1, 2
!                  NN = GET_NN(NNN)
!                  N = IDWETD(NN)
!
!                  ! Call subroutine RAINOUT to compute the fraction
!                  ! of tracer lost to rainout in grid box (I,J,L=LLPAR)
!                  CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )
!
!                  ! WETLOSS is the amount of soluble tracer
!                  ! lost to rainout in grid box (I,J,L=LLPAR)
!                  WETLOSS = STT(I,J,L,N) * RAINFRAC
!
!                  ! Remove rainout losses in grid box (I,J,L=LLPAR) from STT
!                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!
!                  ! DSTT is an accumulator array for rained-out tracers.
!                  ! The tracers in DSTT are in the liquid phase and will
!                  ! precipitate to the levels below until a washout occurs.
!                  ! Initialize DSTT at (I,J,L=LLPAR) with WETLOSS.
!                  DSTT(NN,L,I,J) = WETLOSS
!
!                  ! ND17 diagnostic...LS and conv rainout fractions [unitless]
!                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                     AD17(I,J,L,NN,IDX) =
!     &                    AD17(I,J,L,NN,IDX) + RAINFRAC / F
!                  ENDIF
!
!                  ! ND39 diag - save rainout losses in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                  ENDIF
!
!                  ! Negative tracer...call subroutine SAFETY
!                  IF ( STT(I,J,L,N) < 0d0 ) THEN
!                     CALL SAFETY( I, J, L, N, 3,
!     &                            LS,             PDOWN(L,I,J),
!     &                            QQ(L,I,J),      ALPHA,
!     &                            ALPHA2,         RAINFRAC,
!     &                            WASHFRAC,       MASS_WASH,
!     &                            MASS_NOWASH,    WETLOSS,
!     &                            GAINED,         LOST,
!     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
!                  ENDIF
!               ENDDO
!            ENDIF
!
!            ! Save FTOP for the next lower level
!            FTOP = F
!         ENDIF
!
!         !==============================================================
!         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
!         !
!         ! Rainout occurs when there is more precipitation in grid box
!         ! (I,J,L) than in grid box (I,J,L+1).  In other words, rainout
!         ! occurs when the amount of rain falling through the bottom of
!         ! grid box (I,J,L) is more than the amount of rain coming in
!         ! through the top of grid box (I,J,L).
!         !
!         ! Thus ( PDOWN(L,I,J) > 0 and QQ(L,I,J) > 0 ) is the
!         ! criterion for Rainout.
!         !
!         ! Soluble gases/aerosols are incorporated into the raindrops
!         ! and are completely removed from grid box (I,J,L).  There is
!         ! no evaporation and "resuspension" of aerosols during a
!         ! rainout event.
!         !
!         ! Compute K_RAIN and F' for grid box (I,J,L) exactly as
!         ! described above in Section (4).  K_RAIN and F' depend on
!         ! whether we have large-scale or convective precipitation.
!         !
!         ! F' is the areal fraction of grid box (I,J,L) that is
!         ! precipitating.  However, the effective area of precipitation
!         ! that layer L sees (cf. Eqs. 11-13, Jacob et al, 2000) is
!         ! given by:
!         !
!         !                   F = MAX( F', FTOP )
!         !
!         ! where FTOP = F' at grid box (I,J,L+1), that is, for the grid
!         ! box immediately above the current grid box.
!         !
!         ! Therefore, the effective area of precipitation in grid box
!         ! (I,J,L) depends on the area of precipitation in the grid
!         ! boxes above it.
!         !
!         ! Having computed K_RAIN and F for grid box (I,J,L), call
!         ! routine RAINOUT to compute the fraction of tracer lost to
!         ! rainout conditions.
!         !==============================================================
!         DO L = LLPAR-1, 2, -1
!
!            ! Zero variables for each level
!            ALPHA       = 0d0
!            ALPHA2      = 0d0
!            F           = 0d0
!            F_PRIME     = 0d0
!            GAINED      = 0d0
!            K_RAIN      = 0d0
!            LOST        = 0d0
!            MASS_NOWASH = 0d0
!            MASS_WASH   = 0d0
!            Q           = 0d0
!            QDOWN       = 0d0
!            RAINFRAC    = 0d0
!            WASHFRAC    = 0d0
!            WETLOSS     = 0d0
!
!            ! Rainout criteria
!            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN
!
!               ! Q is the new precip that is forming within grid box (I,J,L)
!               Q = QQ(L,I,J)
!
!               ! Compute K_RAIN and F' for either large-scale or convective
!               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000)
!               IF ( LS ) THEN
!                  K_RAIN  = LS_K_RAIN( Q )
!                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
!               ELSE
!                  K_RAIN  = 1.5d-3
!                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
!               ENDIF
!
!               ! F is the effective area of precip seen by grid box (I,J,L)
!               F = MAX( F_PRIME, FTOP )
!
!               ! Only compute rainout if F > 0.
!               ! This helps to eliminate unnecessary CPU cycles.
!               IF ( F > 0d0 ) THEN
!
!                  ! ND16 diagnostic...save F
!                  IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!                  ENDIF
!
!                  ! ND17 diagnostic...increment counter
!                  IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                     CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
!                  ENDIF
!
!                  ! Loop over soluble tracers and/or aerosol tracers
!                  !DO NN = 1, NSOL
!                  ! Loop over just SO2 then SO4. (dkh, 10/23/05)
!                  DO NNN = 1, 2
!                     NN = GET_NN(NNN)
!                     N = IDWETD(NN)
!
!                     ! Call subroutine RAINOUT to comptue the fraction
!                     ! of tracer lost to rainout in grid box (I,J,L)
!                     CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC )
!
!                     ! WETLOSS is the amount of tracer in grid box
!                     ! (I,J,L) that is lost to rainout.
!                     WETLOSS = STT(I,J,L,N) * RAINFRAC
!
!                     ! For the mercury simulation, we need to archive the
!                     ! amt of Hg2 [kg] that is scavenged out of the column.
!                     ! Also for tagged Hg2. (sas, cdh, bmy, 1/6/06)
!                     IF ( IS_Hg .and. IS_Hg2( N ) ) THEN
!                        CALL ADD_Hg2_WD( I, J, N, WETLOSS )
!                     ENDIF
!
!                     ! Subtract the rainout loss in grid box (I,J,L) from STT
!                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!
!                     ! Add to DSTT the tracer lost to rainout in grid box
!                     ! (I,J,L) plus the tracer lost to rainout from grid box
!                     ! (I,J,L+1), which has by now precipitated down into
!                     ! grid box (I,J,L).  DSTT will continue to accumulate
!                     ! rained out tracer in this manner until a washout
!                     ! event occurs.
!                     DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
!
!                     ! ND17 diagnostic...rainout fractions [unitless]
!                     IF ( ND17 > 0 .and. L <= LD17 ) THEN
!                        AD17(I,J,L,NN,IDX) =
!     &                       AD17(I,J,L,NN,IDX) + RAINFRAC / F
!                     ENDIF
!
!                     ! ND39 diag -- save rainout losses in [kg/s]
!                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                     ENDIF
!
!                     ! Negative tracer...call subroutine SAFETY
!                     IF ( STT(I,J,L,N) < 0d0 .or.
!     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!                        CALL SAFETY( I, J, L, N, 4,
!     &                               LS,             PDOWN(L,I,J),
!     &                               QQ(L,I,J),      ALPHA,
!     &                               ALPHA2,         RAINFRAC,
!     &                               WASHFRAC,       MASS_WASH,
!     &                               MASS_NOWASH,    WETLOSS,
!     &                               GAINED,         LOST,
!     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
!                     ENDIF
!                  ENDDO
!               ENDIF
!
!               ! Save FTOP for next level
!               FTOP = F
!
!            !==============================================================
!            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
!            !
!            ! Washout occurs when we have evaporation (or no precipitation
!            ! at all) at grid box (I,J,L), but have rain coming down from
!            ! grid box (I,J,L+1).
!            !
!            ! Thus PDOWN(L,I,J) > 0 and QQ(L,I,J) <= 0 is the criterion
!            ! for Washout.  Also recall that QQ(L,I,J) < 0 denotes
!            ! evaporation and not precipitation.
!            !
!            ! A fraction ALPHA of the raindrops falling down from grid
!            ! box (I,J,L+1) to grid box (I,J,L) will evaporate along the
!            ! way.  ALPHA is given by:
!            !
!            !            precip leaving (I,J,L+1) - precip leaving (I,J,L)
!            !  ALPHA = ---------------------------------------------------
!            !                     precip leaving (I,J,L+1)
!            !
!            !
!            !                    -QQ(L,I,J) * DZ(I,J,L)
!            !        =         --------------------------
!            !                        PDOWN(L+1,I,J)
!            !
!            ! We assume that a fraction ALPHA2 = 0.5 * ALPHA of the
!            ! previously rained-out aerosols and HNO3 coming down from
!            ! level (I,J,L+1) will evaporate and re-enter the atmosphere
!            ! in the gas phase in grid box (I,J,L).  This process is
!            ! called "resuspension".
!            !
!            ! For non-aerosol species, the amount of previously rained
!            ! out mass coming down from grid box (I,J,L+1) to grid box
!            ! (I,J,L) is figured into the total mass available for
!            ! washout in grid box (I,J,L).  We therefore do not have to
!            ! use the fraction ALPHA2 to compute the resuspension.
!            !
!            ! NOTE from Hongyu Liu about ALPHA (hyl, 2/29/00)
!            ! =============================================================
!            ! If our QQ field was perfect, the evaporated amount in grid
!            ! box (I,J,L) would be at most the total rain amount coming
!            ! from above (i.e. PDOWN(I,J,L+1) ). But this is not true for
!            ! the MOISTQ field we are using.  Sometimes the evaporation in
!            ! grid box (I,J,L) can be more than the rain amount from above.
!            ! The reason is our "evaporation" also includes the effect of
!            ! cloud detrainment.  For now we cannot find a way to
!            ! distinguish betweeen the two. We then decided to release
!            ! aerosols in both the detrained air and the evaporated air.
!            !
!            ! Therefore, we should use this term in the numerator:
!            !
!            !                -QQ(I,J,L) * BXHEIGHT(I,J,L)
!            !
!            ! instead of the term:
!            !
!            !                PDOWN(L+1)-PDOWN(L)
!            !
!            ! Recall that in make_qq.f we have restricted PDOWN to
!            ! positive values, otherwise, QQ would be equal to
!            ! PDOWN(L+1)-PDOWN(L).
!            !==============================================================
!            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN
!
!               ! QDOWN is the precip leaving thru the bottom of box (I,J,L)
!               ! Q     is the new precip that is forming within box (I,J,L)
!               QDOWN = PDOWN(L,I,J)
!               Q     = QQ(L,I,J)
!
!               ! Since no precipitation is forming within grid box (I,J,L),
!               ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
!               F = FTOP
!
!               ! Only compute washout if F > 0.
!               ! This helps to eliminate needless CPU cycles.
!               IF ( F > 0d0 ) THEN
!
!                  ! ND16 diagnostic...save F (fraction of grid box raining)
!                  IF ( ND16 > 0d0 .and. L <= LD16 ) THEN
!                     AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                     CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!                  ENDIF
!
!                  ! ND18 diagnostic...increment counter
!                  IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                     CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
!                  ENDIF
!
!                  ! Loop over soluble tracers and/or aerosol tracers
!                  !DO NN = 1, NSOL
!                  ! Loop over just SO2 then SO4. (dkh, 10/23/05)
!                  DO NNN = 1, 2
!                     NN = GET_NN(NNN)
!                     N  = IDWETD(NN)
!
!                     ! Call WASHOUT to compute the fraction of
!                     ! tracer lost to washout in grid box (I,J,L)
!                     CALL WASHOUT( I,     J,  L, N,
!     &                             QDOWN, DT, F, WASHFRAC, AER )
!
!                     !=====================================================
!                     ! Washout of aerosol tracers --
!                     ! this is modeled as a kinetic process
!                     !=====================================================
!                     IF ( AER ) THEN
!
!                        ! ALPHA is the fraction of the raindrops that
!                        ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
!                        ALPHA = ( ABS( Q ) * BXHEIGHT(I,J,L) * 100d0 ) /
!     &                            PDOWN(L+1,I,J)
!
!                        ! ALPHA2 is the fraction of the rained-out aerosols
!                        ! that gets resuspended in grid box (I,J,L)
!                        ALPHA2 = 0.5d0 * ALPHA
!
!                        ! GAINED is the rained out aerosol coming down from
!                        ! grid box (I,J,L+1) that will evaporate and re-enter
!                        ! the atmosphere in the gas phase in grid box (I,J,L).
!                        GAINED = DSTT(NN,L+1,I,J) * ALPHA2
!
!                        ! Amount of aerosol lost to washout in grid box
!                        ! (qli, bmy, 10/29/02)
!                        WETLOSS = STT(I,J,L,N) * WASHFRAC - GAINED
!
!                        ! Remove washout losses in grid box (I,J,L) from STT.
!                        ! Add the aerosol that was reevaporated in (I,J,L).
!                        ! SO2 in sulfate chemistry is wet-scavenged on the
!                        ! raindrop and converted to SO4 by aqeuous chem.
!                        ! If evaporation occurs then SO2 comes back as SO4
!                        ! (rjp, bmy, 3/23/03)
!                        IF ( N == IDTSO2 ) THEN
!                            STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
!     &                                        + GAINED * 96D0 / 64D0
!
!                            STT(I,J,L,N)      = STT(I,J,L,N) *
!     &                                          ( 1d0 - WASHFRAC )
!                        ELSE
!                            STT(I,J,L,N)      = STT(I,J,L,N) - WETLOSS
!                        ENDIF
!
!                        ! LOST is the rained out aerosol coming down from
!                        ! grid box (I,J,L+1) that will remain in the liquid
!                        ! phase in grid box (I,J,L) and will NOT re-evaporate.
!                        LOST = DSTT(NN,L+1,I,J) - GAINED
!
!                        ! Add the washed out tracer from grid box (I,J,L) to
!                        ! DSTT.  Also add the amount of tracer coming down
!                        ! from grid box (I,J,L+1) that does NOT re-evaporate.
!                        DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
!                        ! Maybe it should be this ????
!                        !DSTT(NN,L,I,J) = LOST + WETLOSS
!
!                        ! ND18 diagnostic...divide washout fraction by F
!                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                           AD18(I,J,L,NN,IDX) =
!     &                          AD18(I,J,L,NN,IDX) + WASHFRAC / F
!                        ENDIF
!
!                     !=====================================================
!                     ! Washout of non-aerosol tracers
!                     ! This is modeled as an equilibrium process
!                     !=====================================================
!                     ELSE
!
!                        ! MASS_NOWASH is the amount of non-aerosol tracer in
!                        ! grid box (I,J,L) that is NOT available for washout.
!                        MASS_NOWASH = ( 1d0 - F ) * STT(I,J,L,N)
!
!                        ! MASS_WASH is the total amount of non-aerosol tracer
!                        ! that is available for washout in grid box (I,J,L).
!                        ! It consists of the mass in the precipitating
!                        ! part of box (I,J,L), plus the previously rained-out
!                        ! tracer coming down from grid box (I,J,L+1).
!                        ! (Eq. 15, Jacob et al, 2000).
!                        MASS_WASH = ( F*STT(I,J,L,N) ) +DSTT(NN,L+1,I,J)
!
!                        ! WETLOSS is the amount of tracer mass in
!                        ! grid box (I,J,L) that is lost to washout.
!                        ! (Eq. 16, Jacob et al, 2000)
!                        WETLOSS = MASS_WASH * WASHFRAC -DSTT(NN,L+1,I,J)
!
!                        ! The tracer left in grid box (I,J,L) is what was
!                        ! in originally in the non-precipitating fraction
!                        ! of the box, plus MASS_WASH, less WETLOSS.
!                        STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!
!                        ! Add washout losses in grid box (I,J,L) to DSTT
!                        DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS
!
!                        ! For the mercury simulation, we need to archive the
!                        ! amt of Hg2 [kg] that is scavenged out of the column.
!                        ! Also for tagged Hg2. (sas, cdh, bmy, 1/6/06)
!                        IF ( IS_Hg .and. IS_Hg2( N ) ) THEN
!                           CALL ADD_Hg2_WD( I, J, N, WETLOSS )
!                        ENDIF
!
!                        ! ND18 diagnostic...we don't have to divide the
!                        ! washout fraction by F since this is accounted for.
!                        IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                           AD18(I,J,L,NN,IDX) =
!     &                          AD18(I,J,L,NN,IDX) + WASHFRAC
!                        ENDIF
!                     ENDIF
!
!                     ! ND39 diag -- save rainout losses in [kg/s]
!                     IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                        AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                     ENDIF
!
!                     ! Negative tracer...call subroutine SAFETY
!                     IF ( STT(I,J,L,N) < 0d0 .or.
!     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!                        CALL SAFETY( I, J, L, N, 5,
!     &                               LS,             PDOWN(L,I,J),
!     &                               QQ(L,I,J),      ALPHA,
!     &                               ALPHA2,         RAINFRAC,
!     &                               WASHFRAC,       MASS_WASH,
!     &                               MASS_NOWASH,    WETLOSS,
!     &                               GAINED,         LOST,
!     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
!                     ENDIF
!                  ENDDO
!               ENDIF
!
!               ! Save FTOP for next level
!               FTOP = F
!
!            !===========================================================
!            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n
!            !
!            ! If there is no precipitation leaving grid box (I,J,L),
!            ! then  set F, the effective area of precipitation in grid
!            ! box (I,J,L), to zero.
!            !
!            ! Also, all of the previously rained-out tracer that is now
!            ! coming down from grid box (I,J,L+1) will evaporate and
!            ! re-enter the atmosphere in the gas phase in grid box
!            ! (I,J,L).  This is called "resuspension".
!            !===========================================================
!            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN
!
!               ! No precipitation at grid box (I,J,L), thus F = 0
!               F = 0d0
!
!               ! Loop over soluble tracers and/or aerosol tracers
!               !DO NN = 1, NSOL
!               ! Loop over just SO2 then SO4. (dkh, 10/23/05)
!               DO NNN = 1, 2
!                  NN = GET_NN(NNN)
!                  N = IDWETD(NN)
!
!                  ! WETLOSS is the amount of tracer in grid box (I,J,L)
!                  ! that is lost to rainout. (qli, bmy, 10/29/02)
!                  WETLOSS = -DSTT(NN,L+1,I,J)
!
!                  ! For the mercury simulation, we need to archive the
!                  ! amt of Hg2 [kg] that is scavenged out of the column.
!                  ! Also for tagged Hg2. (sas, cdh, bmy, 1/6/06)
!                  IF ( IS_Hg .and. IS_Hg2( N ) ) THEN
!                     CALL ADD_Hg2_WD( I, J, N, WETLOSS )
!                  ENDIF
!
!                  ! All of the rained-out tracer coming from grid box
!                  ! (I,J,L+1) goes back into the gas phase at (I,J,L)
!                  ! In evap, SO2 comes back as SO4 (rjp, bmy, 3/23/03)
!                  IF ( N == IDTSO2 ) THEN
!                     STT(I,J,L,IDTSO4) = STT(I,J,L,IDTSO4)
!     &                                 - ( WETLOSS * 96d0 / 64d0 )
!                  ELSE
!                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!                  ENDIF
!
!                  ! There is nothing rained out/washed out in grid box
!                  ! (I,J,L), so set DSTT at grid box (I,J,L) to zero.
!		  DSTT(NN,L,I,J) = 0d0
!
!                  ! ND39 diag -- save rainout losses in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                  ENDIF
!
!                  ! Negative tracer...call subroutine SAFETY
!                  IF ( STT(I,J,L,N) < 0d0 ) THEN
!                     CALL SAFETY( I, J, L, N, 6,
!     &                            LS,             PDOWN(L,I,J),
!     &                            QQ(L,I,J),      ALPHA,
!     &                            ALPHA2,         RAINFRAC,
!     &                            WASHFRAC,       MASS_WASH,
!     &                            MASS_NOWASH,    WETLOSS,
!     &                            GAINED,         LOST,
!     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
!                  ENDIF
!               ENDDO
!
!               ! Save FTOP for next level
!               FTOP = F
!            ENDIF
!         ENDDO
!
!         !==============================================================
!         ! (7)  W a s h o u t   i n   L e v e l   1
!         !
!         ! Assume all of the tracer precipitating down from grid box
!         ! (I,J,L=2) to grid box (I,J,L=1) gets washed out in grid box
!         ! (I,J,L=1).
!         !==============================================================
!
!         ! Zero variables for this level
!         ALPHA       = 0d0
!         ALPHA2      = 0d0
!         F           = 0d0
!         F_PRIME     = 0d0
!         GAINED      = 0d0
!         K_RAIN      = 0d0
!         LOST        = 0d0
!         MASS_NOWASH = 0d0
!         MASS_WASH   = 0d0
!         Q           = 0d0
!         QDOWN       = 0d0
!         RAINFRAC    = 0d0
!         WASHFRAC    = 0d0
!         WETLOSS     = 0d0
!
!         ! We are at the surface, set L = 1
!         L = 1
!
!         ! Washout at level 1 criteria
!         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN
!
!            ! QDOWN is the precip leaving thru the bottom of box (I,J,L+1)
!            QDOWN = PDOWN(L+1,I,J)
!
!            ! Since no precipitation is forming within grid box (I,J,L),
!            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
!            F = FTOP
!
!            ! Only compute washout if F > 0.
!            ! This helps to eliminate unnecessary CPU cycles.
!            IF ( F > 0d0 ) THEN
!
!               ! ND16 diagnostic...save F
!               IF ( ND16 > 0 .and. L <= LD16 ) THEN
!                  AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
!                  CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
!               ENDIF
!
!               ! ND18 diagnostic...increment counter
!               IF ( ND18 > 0 .and. L <= LD18 ) THEN
!                  CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
!               ENDIF
!
!               ! Loop over soluble tracers and/or aerosol tracers
!               !DO NN = 1, NSOL
!               ! Loop over just SO2 then SO4. (dkh, 10/23/05)
!               DO NNN = 1, 2
!                  NN = GET_NN(NNN)
!                  N = IDWETD(NN)
!
!                  ! Call WASHOUT to compute the fraction of tracer
!                  ! in grid box (I,J,L) that is lost to washout.
!                  CALL WASHOUT( I,     J,  L, N,
!     &                          QDOWN, DT, F, WASHFRAC, AER )
!
!                  ! NOTE: for HNO3 and aerosols, there is an F factor
!                  ! already present in WASHFRAC.  For other soluble
!                  ! gases, we need to multiply by the F (hyl, bmy, 10/27/00)
!                  IF ( AER ) THEN
!                     WETLOSS = STT(I,J,L,N) * WASHFRAC
!                  ELSE
!                     WETLOSS = STT(I,J,L,N) * WASHFRAC * F
!                  ENDIF
!
!                  ! Subtract WETLOSS from STT
!                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS
!
!                  ! For the mercury simulation, we need to archive the
!                  ! amt of Hg2 [kg] that is scavenged out of the column.
!                  ! Also for tagged Hg2. (sas, cdh, bmy, 1/6/06)
!                  IF ( IS_Hg .and. IS_Hg2( N ) ) THEN
!                     CALL ADD_Hg2_WD( I, J, N, WETLOSS )
!                  ENDIF
!
!                  ! ND18 diagnostic...LS and conv washout fractions [unitless]
!                  IF ( ND18 > 0 .and. L <= LD18 ) THEN
!
!                     ! Only divide WASHFRAC by F for aerosols, since
!                     ! for non-aerosols this is already accounted for
!                     IF ( AER ) THEN
!                        TMP = WASHFRAC / F
!                     ELSE
!                        TMP = WASHFRAC
!                     ENDIF
!
!                     AD18(I,J,L,NN,IDX) = AD18(I,J,L,NN,IDX) + TMP
!                  ENDIF
!
!                  ! ND39 diag -- save washout loss in [kg/s]
!                  IF ( ND39 > 0 .and. L <= LD39 ) THEN
!                     AD39(I,J,L,NN) = AD39(I,J,L,NN) + WETLOSS / DT
!                  ENDIF
!
!                  !-----------------------------------------------------
!                  ! Dirty kludge to prevent wet deposition from removing
!                  ! stuff from stratospheric boxes -- this can cause
!                  ! negative tracer (rvm, bmy, 6/21/00)
!                  !
!                  IF ( STT(I,J,L,N) < 0d0 .and. L > 23 ) THEN
!                      WRITE ( 6, 101 ) I, J, L, N, 7
! 101                  FORMAT( 'WETDEP - STT < 0 at ', 3i4,
!     &                        ' for tracer ', i4, 'in area ', i4 )
!                      PRINT*, 'STT:', STT(I,J,:,N)
!                      STT(I,J,L,N) = 0d0
!                  ENDIF
!                  !-----------------------------------------------------
!
!                  ! Negative tracer...call subroutine SAFETY
!                  IF ( STT(I,J,L,N) < 0d0 ) THEN
!                     CALL SAFETY( I, J, L, N, 7,
!     &                            LS,             PDOWN(L,I,J),
!     &                            QQ(L,I,J),      ALPHA,
!     &                            ALPHA2,         RAINFRAC,
!     &                            WASHFRAC,       MASS_WASH,
!     &                            MASS_NOWASH,    WETLOSS,
!     &                            GAINED,         LOST,
!     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
!                  ENDIF
!               ENDDO
!            ENDIF
!         ENDIF
!      ENDDO
!      ENDDO
!#if   !defined( SGI_MIPS )
!!$OMP END PARALLEL DO
!#endif
!
!      ! Return to calling program
!      END SUBROUTINE WETDEP
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_NN( NNN ) RESULT( NN )
!!
!!******************************************************************************
!!  Function GET_NN returns the value of NN for SO2 or SO4. (dkh, 10/23/05)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) NNN      (INTEGER) : Index of loop over SOX species
!!
!!  Function value:
!!  ============================================================================
!!  (2 ) NN       (INTEGER) : Index of SOX species in IDWETD
!!
!!  NOTES:
!! (1 ) Added by adj_group                 (dkh, ks, mak, cs  06/08/09)
!!
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE TRACERID_MOD,  ONLY : IDTSO2, IDTSO4
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: NNN
!
!      ! Function value
!      INTEGER             :: NN
!
!      ! Local variables
!      INTEGER             :: N
!
!      !==================================================================
!      ! GET_NN begins here!
!      !==================================================================
!
!      SELECT CASE( NNN )
!
!        ! Look for NN coresponding to SO2
!        CASE ( 1 )
!
!           DO N = 1, NSOLMAX
!             IF ( IDWETD( N ) == IDTSO2 ) THEN
!                NN = N
!                RETURN
!             ENDIF
!           ENDDO
!
!        ! Look for NN coresponding to SO4
!        CASE ( 2 )
!
!           DO N = 1, NSOLMAX
!             IF ( IDWETD( N ) == IDTSO4 ) THEN
!                NN = N
!                RETURN
!             ENDIF
!           ENDDO
!
!       END SELECT
!
!      ! Return
!      END FUNCTION GET_NN
!
!!------------------------------------------------------------------------------
!
!
!      FUNCTION LS_K_RAIN( Q ) RESULT( K_RAIN )
!!
!!******************************************************************************
!!  Function LS_K_RAIN computes K_RAIN, the first order rainout rate constant
!!  for large-scale (a.k.a. stratiform) precipitation (bmy, 3/18/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Q      (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!!
!!  Function value:
!!  ============================================================================
!!  (2 ) K_RAIN (REAL*8) : 1st order rainout rate constant [s-1]
!!
!!  NOTES:
!!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: Q
!
!      ! Function value
!      REAL*8             :: K_RAIN
!
!      !==================================================================
!      ! LS_K_RAIN begins here!
!      !==================================================================
!
!      ! Compute rainout rate constant K in s^-1 (Eq. 12, Jacob et al, 2000).
!      ! 1.0d-4 = K_MIN, a minimum value for K_RAIN
!      ! 1.5d-6 = L + W, the condensed water content (liq + ice) in the cloud
!      K_RAIN = 1.0d-4 + ( Q / 1.5d-6 )
!
!      ! Return to WETDEP
!      END FUNCTION LS_K_RAIN
!
!!------------------------------------------------------------------------------
!
!      FUNCTION LS_F_PRIME( Q, K_RAIN ) RESULT( F_PRIME )
!!
!!******************************************************************************
!!  Function LS_F_PRIME computes F', the fraction of the grid box that is
!!  precipitating during large scale (a.k.a. stratiform) precipitation.
!!  (bmy, 3/18/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Q       (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!!
!!  Function value:
!!  ============================================================================
!!  (3 ) F_PRIME (REAL*8) : Fraction of grid box undergoing LS precip [unitless]
!!
!!  NOTES:
!!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: Q, K_RAIN
!
!      ! Function value
!      REAL*8             :: F_PRIME
!
!      !=================================================================
!      ! LS_F_PRIME begins here!
!      !=================================================================
!
!      ! Compute F', the area of the grid box undergoing precipitation
!      ! 1.5d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
!      F_PRIME = Q / ( K_RAIN * 1.5d-6 )
!
!      ! Return to WETDEP
!      END FUNCTION LS_F_PRIME
!
!!------------------------------------------------------------------------------
!
!      FUNCTION CONV_F_PRIME( Q, K_RAIN, DT ) RESULT( F_PRIME )
!!
!!******************************************************************************
!!  Function CONV_F_PRIME computes F', the fraction of the grid box that is
!!  precipitating during convective precipitation. (bmy, 3/18/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!!
!!  Function value:
!!  ============================================================================
!!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!!
!!  NOTES:
!!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!!******************************************************************************
!!
!      ! Arguments
!      REAL*8, INTENT(IN) :: Q, K_RAIN, DT
!
!      ! Local variables
!      REAL*8             :: TIME
!
!      ! Function value
!      REAL*8             :: F_PRIME
!
!      !=================================================================
!      ! CONV_F_PRIME begins here!
!      !=================================================================
!
!      ! Assume the rainout event happens in 30 minutes (1800 s)
!      ! Compute the minimum of DT / 1800s and 1.0
!      TIME = MIN( DT / 1800d0, 1d0 )
!
!      ! Compute F' for convective precipitation (Eq. 13, Jacob et al, 2000)
!      ! 0.3  = FMAX, the maximum value of F' for convective precip
!      ! 2d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
!      F_PRIME = ( 0.3d0 * Q * TIME ) /
!     &          ( ( Q * TIME ) + ( 0.3d0 * K_RAIN * 2d-6 ) )
!
!      ! Return to WETDEP
!      END FUNCTION CONV_F_PRIME
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE SAFETY( I,         J,           L,        N,
!     &                   A,         LS,          PDOWN,    QQ,
!     &                   ALPHA,     ALPHA2,      RAINFRAC, WASHFRAC,
!     &                   MASS_WASH, MASS_NOWASH, WETLOSS,  GAINED,
!     &                   LOST,      DSTT,        STT )
!!
!!******************************************************************************
!!  Subroutine SAFETY stops the run with debug output and an error message
!!  if negative tracers are found. (bmy, 3/18/04)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!!
!!  Function value:
!!  ============================================================================
!!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!!
!!  NOTES:
!!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
!#     include "CMN_SIZE"
!
!      ! Arguments
!      LOGICAL, INTENT(IN) :: LS
!      INTEGER, INTENT(IN) :: I, J, L, N, A
!      REAL*8,  INTENT(IN) :: PDOWN,    QQ,       ALPHA,     ALPHA2
!      REAL*8,  INTENT(IN) :: RAINFRAC, WASHFRAC, MASS_WASH, MASS_NOWASH
!      REAL*8,  INTENT(IN) :: WETLOSS,  GAINED,   LOST,      DSTT(LLPAR)
!      REAL*8,  INTENT(IN) :: STT(LLPAR)
!
!      !=================================================================
!      ! SAFETY begins here!
!      !=================================================================
!
!      ! Print line
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!      ! Write error message and stop the run
!      WRITE ( 6, 100 ) I, J, L, N, A
! 100  FORMAT( 'WETDEP - STT < 0 at ', 3i4, ' for tracer ', i4,
!     &        ' in area ', i4 )
!
!      PRINT*, 'LS          : ', LS
!      PRINT*, 'PDOWN       : ', PDOWN
!      PRINT*, 'QQ          : ', QQ
!      PRINT*, 'ALPHA       : ', ALPHA
!      PRINT*, 'ALPHA2      : ', ALPHA2
!      PRINT*, 'RAINFRAC    : ', RAINFRAC
!      PRINT*, 'WASHFRAC    : ', WASHFRAC
!      PRINT*, 'MASS_WASH   : ', MASS_WASH
!      PRINT*, 'MASS_NOWASH : ', MASS_NOWASH
!      PRINT*, 'WETLOSS     : ', WETLOSS
!      PRINT*, 'GAINED      : ', GAINED
!      PRINT*, 'LOST        : ', LOST
!      PRINT*, 'DSTT(NN,:)  : ', DSTT(:)
!      PRINT*, 'STT(I,J,:N) : ', STT(:)
!
!      ! Print line
!      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!
!      ! Deallocate memory and stop
!      CALL GEOS_CHEM_STOP
!
!      ! Return to WETDEP
!      END SUBROUTINE SAFETY
!
!!------------------------------------------------------------------------------
!
!      SUBROUTINE WETDEPID
!!
!!******************************************************************************
!!  Subroutine WETDEPID sets up the index array of soluble tracers used in
!!  the WETDEP routine above (bmy, 11/8/02, 5/18/06)
!!
!!  NOTES:
!!  (1 ) Now references "tracerid_mod.f".  Also references "CMN" in order to
!!        pass variables NSRCX and NTRACE. (bmy, 11/8/02)
!!  (2 ) Updated for carbon aerosol & dust tracers (rjp, bmy, 4/5/04)
!!  (3 ) Updated for seasalt aerosol tracers.  Also added fancy output.
!!        (rjp, bec, bmy, 4/20/04)
!!  (4 ) Updated for secondary organic aerosol tracers (bmy, 7/13/04)
!!  (5 ) Now references N_TRACERS, TRACER_NAME, TRACER_MW_KG from
!!        "tracer_mod.f".  Removed reference to NSRCX.  (bmy, 7/20/04)
!!  (6 ) Updated for mercury aerosol tracers (eck, bmy, 12/9/04)
!!  (7 ) Updated for AS, AHS, LET, NH4aq, SO4aq (cas, bmy, 12/20/04)
!!  (8 ) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!!  (10) Now use IS_Hg2 and IS_HgP to determine if a tracer is a tagged Hg2
!!        or HgP tracer (bmy, 1/6/06)
!!  (11) Now added SOG4 and SOA4 (dkh, bmy, 5/18/06)
!!******************************************************************************
!!
!      ! References To F90 modules
!      USE ERROR_MOD,    ONLY : ERROR_STOP
!      USE TRACER_MOD,   ONLY : N_TRACERS, TRACER_NAME, TRACER_MW_G
!      USE TRACERID_MOD, ONLY : IDTPB,     IDTBE7,    IDTHNO3, IDTH2O2
!      USE TRACERID_MOD, ONLY : IDTCH2O,   IDTMP,     IDTSO2,  IDTSO4
!      USE TRACERID_MOD, ONLY : IDTSO4s,   IDTSO4aq,  IDTMSA,  IDTNH3
!      USE TRACERID_MOD, ONLY : IDTNH4,    IDTNH4aq,  IDTNIT,  IDTNITs
!      USE TRACERID_MOD, ONLY : IDTAS,     IDTAHS,    IDTLET,  IDTBCPI
!      USE TRACERID_MOD, ONLY : IDTOCPI,   IDTBCPO,   IDTOCPO, IDTDST1
!      USE TRACERID_MOD, ONLY : IDTDST2,   IDTDST3,   IDTDST4, IDTSALA
!      USE TRACERID_MOD, ONLY : IDTSALC,   IDTALPH,   IDTLIMO, IDTALCO
!      USE TRACERID_MOD, ONLY : IDTSOG1,   IDTSOG2,   IDTSOG3, IDTSOG4
!      USE TRACERID_MOD, ONLY : IDTSOA1,   IDTSOA2,   IDTSOA3, IDTSOA4
!      USE TRACERID_MOD, ONLY : IS_Hg2,    IS_HgP
!      USE TRACERID_MOD, ONLY : IDTGLYX,   IDTMGLY,   IDTGLYC
!      USE TRACERID_MOD, ONLY : IDTSOAG,   IDTSOAM
!
!#     include "CMN_SIZE"  ! Size parameters
!
!      ! Local variables
!      INTEGER :: N, NN
!
!      !=================================================================
!      ! WETDEPID begins here!
!      !=================================================================
!
!      ! Zero NSOL
!      NSOL = 0
!
!      ! Sort soluble tracers into IDWETD
!      DO N = 1, N_TRACERS
!
!         !-----------------------------
!         ! Rn-Pb-Be tracers
!         !-----------------------------
!         IF ( N == IDTPB ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTPB
!
!         ELSE IF ( N == IDTBE7 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTBE7
!
!         !-----------------------------
!         ! Full chemistry tracers
!         !-----------------------------
!         ELSE IF ( N == IDTHNO3 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTHNO3
!
!         ELSE IF ( N == IDTH2O2 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTH2O2
!
!         ELSE IF ( N == IDTCH2O ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTCH2O
!
!         ELSE IF ( N == IDTMP ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTMP
!
!         ELSE IF ( N == IDTGLYX ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTGLYX
!
!         ELSE IF ( N == IDTMGLY ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTMGLY
!
!         ELSE IF ( N == IDTGLYC ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTGLYC
!
!         !-----------------------------
!         ! Sulfate aerosol tracers
!         !-----------------------------
!         ELSE IF ( N == IDTSO2 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSO2
!
!         ELSE IF ( N == IDTSO4 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSO4
!
!         ELSE IF ( N == IDTSO4s ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSO4s
!
!         ELSE IF ( N == IDTMSA ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTMSA
!
!         ELSE IF ( N == IDTNH3 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTNH3
!
!         ELSE IF ( N == IDTNH4 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTNH4
!
!         ELSE IF ( N == IDTNIT ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTNIT
!
!         ELSE IF ( N == IDTNITs ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTNITs
!
!         !-----------------------------
!         ! Crystal & Aqueous aerosols
!         !-----------------------------
!         ELSE IF ( N == IDTAS ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTAS
!
!         ELSE IF ( N == IDTAHS ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTAHS
!
!         ELSE IF ( N == IDTLET ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTLET
!
!         ELSE IF ( N == IDTNH4aq ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTNH4aq
!
!         ELSE IF ( N == IDTSO4aq ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSO4aq
!
!         !-----------------------------
!         ! Carbon & SOA aerosol tracers
!         !-----------------------------
!         ELSE IF ( N == IDTBCPI ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTBCPI
!
!         ELSE IF ( N == IDTOCPI ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTOCPI
!
!         ELSE IF ( N == IDTBCPO ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTBCPO
!
!         ELSE IF ( N == IDTOCPO ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTOCPO
!
!         ELSE IF ( N == IDTALPH ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTALPH
!
!         ELSE IF ( N == IDTLIMO ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTLIMO
!
!         ELSE IF ( N == IDTALCO ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTALCO
!
!         ELSE IF ( N == IDTSOG1 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOG1
!
!         ELSE IF ( N == IDTSOG2 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOG2
!
!         ELSE IF ( N == IDTSOG3 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOG3
!
!         ELSE IF ( N == IDTSOG4 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOG4
!
!         ELSE IF ( N == IDTSOA1 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOA1
!
!         ELSE IF ( N == IDTSOA2 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOA2
!
!         ELSE IF ( N == IDTSOA3 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOA3
!
!         ELSE IF ( N == IDTSOA4 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOA4
!
!         ELSE IF ( N == IDTSOAG ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOAG
!
!         ELSE IF ( N == IDTSOAM ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSOAM
!
!         !-----------------------------
!         ! Dust aerosol tracers
!         !-----------------------------
!         ELSE IF ( N == IDTDST1 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTDST1
!
!         ELSE IF ( N == IDTDST2 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTDST2
!
!         ELSE IF ( N == IDTDST3 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTDST3
!
!         ELSE IF ( N == IDTDST4 ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTDST4
!
!         !-----------------------------
!         ! Seasalt aerosol tracers
!         !-----------------------------
!         ELSE IF ( N == IDTSALA ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSALA
!
!         ELSE IF ( N == IDTSALC ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = IDTSALC
!
!         !-----------------------------
!         ! Total and tagged Hg tracers
!         !-----------------------------
!         ELSE IF ( IS_Hg2( N ) ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = N
!
!         ELSE IF ( IS_HgP( N ) ) THEN
!            NSOL         = NSOL + 1
!            IDWETD(NSOL) = N
!
!         ENDIF
!      ENDDO
!
!      ! Error check: Make sure that NSOL is less than NSOLMAX
!      IF ( NSOL > NSOLMAX ) THEN
!         CALL ERROR_STOP( 'NSOL > NSOLMAX!', 'WETDEPID (wetscav_mod.f)')
!      ENDIF
!
!      ! Also check to see if NSOL is larger than the maximum
!      ! number of soluble tracers for a particular simulation
!      IF ( NSOL > GET_WETDEP_NMAX() ) THEN
!         CALL ERROR_STOP( 'NSOL > NMAX', 'WETDEPID (wetscav_mod.f)')
!      ENDIF
!
!      !=================================================================
!      ! Echo list of soluble tracers to the screen
!      !=================================================================
!      WRITE( 6, '(/,a,/)' ) 'WETDEPID: List of soluble tracers: '
!      WRITE( 6, '(a)  '   ) '  #             Name  Tracer Mol Wt'
!      WRITE( 6, '(a)'     ) '                      Number g/mole'
!      WRITE( 6, '(a)'     ) REPEAT( '-', 36 )
!
!      DO NN = 1, NSOL
!         N = IDWETD(NN)
!         WRITE( 6, '(i3,3x,a14,3x,i3,3x,f6.1)' )
!     &        NN, TRIM( TRACER_NAME(N) ), N, TRACER_MW_G(N)
!      ENDDO
!
!      ! Return to calling program
!      END SUBROUTINE WETDEPID
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_WETDEP_NMAX() RESULT ( NMAX )
!!
!!******************************************************************************
!!  Function GET_WETDEP_NMAX returns the maximum number of soluble tracers
!!  for a given type of simulation.  Primarily used for allocation of
!!  diagnostic arrays. (bmy, 12/2/02, 5/18/06)
!!
!!  NOTES:
!!  (1 ) Modified to include carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
!!  (2 ) Modified to include seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!!  (3 ) Modified to include 2ndary organic aerosol tracers (rjp, bmy, 7/13/04)
!!  (4 ) Now references ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM, and
!!        ITS_A_RnPbBe_SIM from "tracer_mod.f".  Now references LCARB, LDUST,
!!        LSOA, LSSALT, LSULF from "logical_mod.f". (bmy, 7/20/04)
!!  (5 ) Modified to include mercury aerosol tracers (eck, bmy, 12/14/04)
!!  (6 ) Modified for AS, AHS, LET, NH4aq, SO4aq (cas, bmy, 12/20/04)
!!  (7 ) Modified for SO4s, NITs (bec, bmy, 4/25/05)
!!  (8 ) Modified for SOG4, SOA4 (dkh, bmy, 5/18/06)
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE LOGICAL_MOD, ONLY : LCARB,  LDUST, LSOA
!      USE LOGICAL_MOD, ONLY : LSSALT, LSULF, LSPLIT, LCRYST
!      USE TRACER_MOD,  ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM
!      USE TRACER_MOD,  ONLY : ITS_A_RnPbBe_SIM,   ITS_A_MERCURY_SIM
!      USE TRACERID_MOD, ONLY : IDTSOAG,  IDTSOAM
!
!#     include "CMN_SIZE"   ! Size Parameters
!
!      ! Function value
!      INTEGER :: NMAX
!
!      !=================================================================
!      ! GET_WETDEP_NMAX begins here!
!      !
!      ! NOTE: If you add tracers to a simulation, update as necessary
!      !=================================================================
!      IF ( ITS_A_FULLCHEM_SIM() ) THEN
!
!         !-----------------------
!         ! Fullchem simulation
!         !-----------------------
!         NMAX = 7                               ! HNO3, H2O2, CH2O, MP,
!                                                ! GLYX, MGLY, GLYC
!         IF ( LSULF )    NMAX = NMAX + 8        ! SO2, SO4, MSA, NH3, NH4, NIT
!         IF ( LDUST  )   NMAX = NMAX + NDSTBIN  ! plus # of dust bins
!         IF ( LSSALT )   NMAX = NMAX + 2        ! plus 2 seasalts
!
!         IF ( LSOA ) THEN
!            IF ( LCARB ) NMAX = NMAX + 15       ! carbon + SOA aerosols
!            IF ( IDTSOAG /= 0 ) NMAX = NMAX + 1 ! SOAG deposition
!            IF ( IDTSOAM /= 0 ) NMAX = NMAX + 1 ! SOAM deposition
!         ELSE
!            IF ( LCARB ) NMAX = NMAX + 4        ! just carbon aerosols
!         ENDIF
!
!      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
!
!         !-----------------------
!         ! Offline simulation
!         !-----------------------
!         NMAX = 0
!         IF ( LSULF  )   NMAX = NMAX + 9        ! add 9 sulfur species
!         IF ( LCRYST )   NMAX = NMAX + 5        ! add 5 cryst & aq species
!         IF ( LDUST  )   NMAX = NMAX + NDSTBIN  ! Add number of dust bins
!         IF ( LSSALT )   NMAX = NMAX + 2        ! plus 2 seasalts
!
!         IF ( LSOA ) THEN
!            IF ( LCARB ) NMAX = NMAX + 15       ! carbon + SOA aerosols
!         ELSE
!            IF ( LCARB ) NMAX = NMAX + 4        ! just carbon aerosols
!         ENDIF
!
!      ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
!
!         !-----------------------
!         ! Rn-Pb-Be simulation
!         !-----------------------
!         NMAX = 2                               ! 210Pb, 7Be
!
!      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN
!
!         !-----------------------
!         ! Mercury simulation
!         !-----------------------
!         NMAX = 2                               ! Hg2, HgP
!         IF ( LSPLIT ) NMAX = NMAX + 14         ! Tagged tracers
!
!      ELSE
!
!         !-----------------------
!         ! Everything else
!         !-----------------------
!         NMAX = 0
!
!      ENDIF
!
!      ! Return to calling program
!      END FUNCTION GET_WETDEP_NMAX
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_WETDEP_NSOL() RESULT( N_SOLUBLE )
!!
!!******************************************************************************
!!  Function GET_WETDEP_NSOL returns NSOL (# of soluble tracers) to a calling
!!  program outside WETSCAV_MOD.  This is so that we can keep NSOL declared
!!  as a PRIVATE variable. (bmy, 1/10/03)
!!
!!  NOTES:
!!******************************************************************************
!!
!      ! Function value
!      INTEGER :: N_SOLUBLE
!
!      !=================================================================
!      ! GET_WETDEP_NSOL begins here!
!      !=================================================================
!
!      ! Get the # of soluble tracers
!      N_SOLUBLE = NSOL
!
!      ! Return to calling program
!      END FUNCTION GET_WETDEP_NSOL
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GET_WETDEP_IDWETD( NWET ) RESULT( N )
!!
!!******************************************************************************
!!  Function GET_WETDEP_IDWETD returns the tracer number of wet deposition
!!  species  NWET.  This is meant to be called outside of WETSCAV_MOD so that
!!  IDWETD can be kept as a PRIVATE variable. (bmy, 1/10/03)
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) NWET (INTEGER) : Wet deposition species N
!!
!!  NOTES:
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE ERROR_MOD, ONLY : ERROR_STOP
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: NWET
!
!      ! Function value
!      INTEGER             :: N
!
!      !=================================================================
!      ! GET_WETDEP_IDWETD begins here!
!      !=================================================================
!
!      ! Make sure NWET is valid
!      IF ( NWET < 1 .or. NWET > NSOLMAX ) THEN
!         CALL ERROR_STOP( 'Invalid value of NWET!',
!     &                    'GET_N_WETDEP (wetscav_mod.f)' )
!      ENDIF
!
!      ! Get the tracer # for wet deposition species N
!      N = IDWETD(NWET)
!
!      ! Return to calling program
!      END FUNCTION GET_WETDEP_IDWETD
!
!!------------------------------------------------------------------------------

      SUBROUTINE WETSCAV_ADJ_FORCE()
!
!******************************************************************************
!  Subroutine WETSCAV_ADJ_FORCE calculates adjoint forcing for sensitivity of
!  wetdeposition.  Also works for wetdep 4DVAR.  (fp, dkh, 03/04/13)
!
!  NOTES
!
!******************************************************************************
!

      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : TRACER_IND, NOBS, DEP_UNIT
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD, ONLY : NHX_ADJ_FORCE
      USE ADJ_ARRAYS_MOD, ONLY : GET_CF_REGION, ADJOINT_AREA_M2
      USE ADJ_ARRAYS_MOD, ONLY : NSPAN
      USE ADJ_ARRAYS_MOD, ONLY : TRACER_IND
      USE ADJ_ARRAYS_MOD, ONLY : TR_WDEP_CONV
      USE TIME_MOD,       ONLY : GET_TS_DYN
      USE TRACER_MOD,     ONLY : TRACER_MW_G, TRACER_NAME
      USE TRACER_MOD,     ONLY : N_TRACERS
      USE GRID_MOD,       ONLY : GET_AREA_M2


#     include "CMN_SIZE"       ! Size parameters
#     include "define_adj.h"   ! Obs operators


      INTEGER                 :: N, I, J, L, LL, NFORCE
      REAL*8                  :: ADD_STT_ADJ
      REAL*8                  :: ADJ_WET
      REAL*8, SAVE            :: NTSDYN
      REAL*8, SAVE            :: CONV_TIME
      LOGICAL, SAVE           :: FIRST = .TRUE.

      !=================================================================
      ! WETSCAV_ADJ_FORCE begins here!
      !=================================================================

      IF ( FIRST ) THEN

         NTSDYN = NSPAN / ( GET_TS_DYN() / 60D0 )

         !default is kg/s
         !CONV_AREA = 1d0
         CONV_TIME = 1D0 / NTSDYN

         ! print tracer names (all steps fp 06/03/2013)
         DO N = 1, NOBS
            WRITE(*,*) 'Forcing ', TRACER_NAME(TRACER_IND(N)),
     &           ' in ls wetdep (' ,TRIM( DEP_UNIT ),')'
         ENDDO

         FIRST = .FALSE.

      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ADJ_WET, LL, NFORCE, ADD_STT_ADJ )
      DO N = 1, NOBS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! weight is only used for sensitivity studies
         IF ( GET_CF_REGION(I,J,L)   .GT. 0 ) THEN

            NFORCE  = TRACER_IND(N)
            ADJ_WET = BOX_DEP(I,J,L,NFORCE)

            DO LL = 2, L - 1
               ADJ_WET = ADJ_WET * LOWER_DEP(I,J,LL,NFORCE)
            ENDDO

            ADD_STT_ADJ = ADJ_WET
     &                  * GET_CF_REGION(I,J,L)
     &                  * CONV_TIME
     &                  * TR_WDEP_CONV(J,NFORCE)

            ! default unit is kg/s
            STT_ADJ(I,J,L,NFORCE) = STT_ADJ(I,J,L,NFORCE)
     &                            + ADD_STT_ADJ

          ENDIF

      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! return to calling routine
      END SUBROUTINE WETSCAV_ADJ_FORCE
!----------------------------------------------------------------------

      SUBROUTINE INIT_WETSCAV_ADJ
!
!******************************************************************************
!  Subroutine INIT_WETSCAV_ADJ initializes updraft velocity, cloud liquid water
!  content, cloud ice content, and mixing ratio of water fields, which
!  are used in the wet scavenging routines. (bmy, 2/23/00, 3/7/05)
!
!  Need the forward arrays for recalculation.  Also init adjoint arrays. (dkh, 09/30/09)
!
!  NOTES:
!  (1 ) References "e_ice.f" -- routine to compute Eice(T).
!  (2 ) Vud, CLDLIQ, CLDICE, C_H2O are all independent of tracer, so we
!        can compute them once per timestep, before calling the cloud
!        convection and wet deposition routines.
!  (3 ) Set C_H2O = 0 below -120 Celsius.  E_ICE(T) has a lower limit of
!        -120 Celsius, so temperatures lower than this will cause a stop
!        with an error message. (bmy, 6/15/00)
!  (4 ) Replace {IJL}GLOB with IIPAR,JJPAR,LLPAR.  Also rename PW to P.
!        Remove IREF, JREF, these are obsolete.  Now reference IS_WATER
!        from "dao_mod.f" to determine water boxes.
!  (5 ) Removed obsolete code from 9/01.  Updated comments and made
!        cosmetic changes. (bmy, 10/24/01)
!  (6 ) Now use routine GET_PCENTER from "pressure_mod.f" to compute the
!        pressure at the midpoint of grid box (I,J,L).  Also removed P and
!        SIG from the argument list (dsa, bdf, bmy, 8/20/02)
!  (7 ) Now reference T from "dao_mod.f".  Updated comments.  Now allocate
!        Vud, C_H2O, CLDLIQ and CLDICE here on the first call.  Now references
!        ALLOC_ERR from "error_mod.f".  Now set H2O2s and SO2s to the initial
!        values from for the first call to COMPUTE_F .  Now call WETDEPID
!        on the first call to initialize the wetdep index array. (bmy, 1/27/03)
!  (8 ) Now references STT from "tracer_mod.f".  Also now we call WETDEPID
!        from "input_mod.f" (bmy, 7/20/04)
!  (9 ) Now references new function E_ICE, which is an analytic function of
!        Kelvin temperature instead of Celsius. (bmy, 3/7/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,         ONLY : T, IS_WATER
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_WDEP_LS
      USE PRESSURE_MOD,    ONLY : GET_PCENTER
      USE TRACER_MOD,      ONLY : STT
      USE TRACERID_MOD,    ONLY : IDTH2O2, IDTSO2
      USE TRACER_MOD,      ONLY : N_TRACERS !fp wetdep ls


#     include "CMN_SIZE"        ! Size parameters

      ! Local variables
      INTEGER                  :: I, J, L, AS
      REAL*8                   :: PL, TK
      LOGICAL, SAVE            :: FIRST = .TRUE.

      !=================================================================
      ! INIT_WETSCAV_ADJ begins here!
      !=================================================================
      IF ( FIRST ) THEN

!         ! Allocate Vud on first call
!         ALLOCATE( Vud( IIPAR, JJPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'Vud' )
!         Vud = 0d0
!
!         ! Allocate C_H2O on first call
!         ALLOCATE( C_H2O( IIPAR, JJPAR, LLPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'C_H2O' )
!         C_H2O = 0d0
!
!         ! Allocate CLDLIQ on first call
!         ALLOCATE( CLDLIQ( IIPAR, JJPAR, LLPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDLIQ' )
!         CLDLIQ = 0d0
!
!         ! Allocate CLDICE on first call
!         ALLOCATE( CLDICE( IIPAR, JJPAR, LLPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDICE' )
!         CLDICE = 0d0
!
!         ! Allocate H2O2s for wet deposition
!         ALLOCATE( H2O2s( IIPAR, JJPAR, LLPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2s' )
!
!         ! Set H2O2s to the initial H2O2 from STT, so that we will have
!         ! nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
!         IF ( IDTH2O2 > 0 ) THEN
!            H2O2s = STT(:,:,:,IDTH2O2)
!         ELSE
!            H2O2s = 0d0
!         ENDIF
!
!         ! Allocate SO2s for wet deposition
!         ALLOCATE( SO2s( IIPAR, JJPAR, LLPAR ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2s' )
!
!         ! Set SO2s to the initial SO2 from STT, so that we will have
!         ! nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
!         IF ( IDTSO2 > 0 ) THEN
!            SO2s = STT(:,:,:,IDTSO2)
!         ELSE
!            SO2s = 0d0
!         ENDIF

         IF ( IDTSO2 > 0 ) THEN
            ! Allocate SO2s_ADJ for wet deposition
            ALLOCATE( SO2s_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2s_ADJ' )
            SO2s_ADJ = 0d0
         ENDIF

         IF ( IDTH2O2 > 0 ) THEN
            ! Allocate H2O2s_ADJ for wet deposition
            ALLOCATE( H2O2s_ADJ( IIPAR, JJPAR, LLPAR ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'H2O2s_ADJ' )
            H2O2s_ADJ = 0d0
         ENDIF

         IF ( LADJ_WDEP_LS ) THEN

            ALLOCATE( BOX_DEP( IIPAR, JJPAR, LLPAR, N_TRACERS), STAT=AS)
            IF( AS/=0) CALL ALLOC_ERR('BOX_DEP')
            BOX_DEP = 0d0

            ALLOCATE( LOWER_DEP( IIPAR, JJPAR, LLPAR, N_TRACERS),
     &         STAT=AS)
            IF( AS/=0) CALL ALLOC_ERR('LOWER_DEP')
            LOWER_DEP = 1d0

         ENDIF

         ! Reset flag
         FIRST = .FALSE.
      ENDIF

!      !=================================================================
!      ! Compute Vud, CLDLIQ, CLDICE, C_H2O, following Jacob et al, 2000.
!      !=================================================================
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, L, TK, PL )
!!$OMP+SCHEDULE( DYNAMIC )
!      DO L = 1, LLPAR
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Compute Temp [K] and Pressure [hPa]
!         TK = T(I,J,L)
!         PL = GET_PCENTER(I,J,L)
!
!         !==============================================================
!         ! Compute Vud -- 5 m/s over oceans, 10 m/s over land (or ice?)
!         ! Assume Vud is the same at all altitudes; the array can be 2-D
!         !==============================================================
!         IF ( L == 1 ) THEN
!            IF ( IS_WATER( I, J ) ) THEN
!               Vud(I,J) = 5d0
!            ELSE
!               Vud(I,J) = 10d0
!            ENDIF
!         ENDIF
!
!         !==============================================================
!         ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air],
!         ! is a function of the local Kelvin temperature:
!         !
!         !    CLDLIQ = 2e-6                    [     T >= 268 K    ]
!         !    CLDLIQ = 2e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
!         !    CLDLIQ = 0                       [     T <= 248 K    ]
!         !==============================================================
!         IF ( TK >= 268d0 ) THEN
!            CLDLIQ(I,J,L) = 2d-6
!
!         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
!            CLDLIQ(I,J,L) = 2d-6 * ( ( TK - 248d0 ) / 20d0 )
!
!         ELSE
!            CLDLIQ(I,J,L) = 0d0
!
!         ENDIF
!
!         !=============================================================
!         ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
!         !
!         !    CLDICE = 2e-6 - CLDLIQ
!         !=============================================================
!         CLDICE(I,J,L) = 2d-6 - CLDLIQ(I,J,L)
!
!         !=============================================================
!         ! C_H2O is given by Dalton's Law as:
!         !
!         !       C_H2O = Eice( Tk(I,J,L) ) / P(I,J,L)
!         !
!         ! where P(L) = pressure in grid box (I,J,L)
!         !
!         ! and   Tk(I,J,L) is the Kelvin temp. of grid box (I,J,L).
!         !
!         ! and   Eice( Tk(I,J,L) ) is the saturation vapor pressure
!         !       of ice [hPa] at temperature Tk(I,J,L) -- computed in
!         !       routine E_ICE above.
!         !==============================================================
!         C_H2O(I,J,L) = E_ICE( TK ) / PL
!
!      ENDDO
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE INIT_WETSCAV_ADJ

!------------------------------------------------------------------------------


      SUBROUTINE ADJ_INIT_WETSCAV(  )
!
!******************************************************************************
!  Subroutine ADJ_INIT_WETSCAV is the adjoint of INIT_WETSCAV, passing adjoint
!  of H2O2s and SO2s (H2O2s_ADJ, SO2s_ADJ) to the total adjoint tracer array,
!  STT_ADJ (dkh, 10/??/05)
!
!  Module variables as Input:
!  ============================================================================
!  (1 ) H2O2s_ADJ
!  (2 ) SO2s_ADJ
!  (3 ) STT_ADJ
!
!  Module variables as Output:
!  ============================================================================
!  (1 ) STT_ADJ
!
!  NOTES:
!  (1 ) Updated to GCv8 (dkh, 10/04/09)
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE TIME_MOD,       ONLY : GET_NHMS, GET_NHMSb
      USE TIME_MOD,       ONLY : GET_NYMD, GET_NYMDb
      USE TRACERID_MOD,   ONLY : IDTSO2,   IDTH2O2

#     include "CMN_SIZE"       ! Size params

      ! Local variables
      INTEGER                 :: I, J, L

      !=================================================================
      ! ADJ_INIT_WETSCAV begins here!
      !=================================================================

      ! We only need to do this at the initial time step corresponding to
      ! when H2O2s and SO2s were allocated and initialized.
      IF ( GET_NHMSb() == GET_NHMS() .and.
     &     GET_NYMDb() == GET_NYMD() .and.
     &     IDTSO2      >  0          .and.
     &     IDTH2O2     >  0                ) THEN

         ! dkh debug
         print*, ' do ADJ_INIT_WETSCAV '

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

             STT_ADJ(I,J,L,IDTH2O2) = STT_ADJ(I,J,L,IDTH2O2)
     &                                + H2O2s_ADJ(I,J,L)
             STT_ADJ(I,J,L,IDTSO2)  = STT_ADJ(I,J,L,IDTSO2)
     &                                + SO2s_ADJ(I,J,L)

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      ! Return to calling program
      END SUBROUTINE ADJ_INIT_WETSCAV

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_WETSCAV_ADJ

      !=================================================================
      ! Subroutine CLEANUP_WETSCAV deallocates arrays for
      ! wet scavenging / wet deposition
      !=================================================================
!      IF ( ALLOCATED( Vud    ) ) DEALLOCATE( Vud    )
!      IF ( ALLOCATED( C_H2O  ) ) DEALLOCATE( C_H2O  )
!      IF ( ALLOCATED( CLDLIQ ) ) DEALLOCATE( CLDLIQ )
!      IF ( ALLOCATED( CLDICE ) ) DEALLOCATE( CLDICE )
!      IF ( ALLOCATED( PDOWN  ) ) DEALLOCATE( PDOWN  )
!      IF ( ALLOCATED( QQ     ) ) DEALLOCATE( QQ     )
!      IF ( ALLOCATED( H2O2s  ) ) DEALLOCATE( H2O2s  )
!      IF ( ALLOCATED( SO2s   ) ) DEALLOCATE( SO2s   )
      IF ( ALLOCATED( SO2s_ADJ   ) ) DEALLOCATE( SO2s_ADJ   )
      IF ( ALLOCATED( H2O2s_ADJ  ) ) DEALLOCATE( H2O2s_ADJ  )
      IF ( ALLOCATED( BOX_DEP    ) ) DEALLOCATE( BOX_DEP    )
      IF ( ALLOCATED( LOWER_DEP  ) ) DEALLOCATE( LOWER_DEP  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_WETSCAV_ADJ

!-----------------------------------------------------------------------------

      ! End of module
      END MODULE WETSCAV_ADJ_MOD
