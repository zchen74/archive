! $Id: schem_adj.f,v 1.1 2010/05/07 20:39:47 daven Exp $
      SUBROUTINE SCHEM_ADJ
!
!******************************************************************************
!  Subroutine SCHEM_ADJ performs adjoint of strat chem.  (dkh, 05/02/10)
!
!  Based on forward model routine SCHEM (qli, bmy, 11/20/1999, 10/25/05).
!
!  NOTES:
!  (1 ) Use ITS_A_NEW_MONTH instead of older method (dkh, 05/02/10)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT,     GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,         READ_BPCH2
      USE DAO_MOD,        ONLY : AD, T
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE ERROR_MOD,      ONLY : ALLOC_ERR
      USE TIME_MOD,       ONLY : GET_MONTH,        GET_TAU
      USE TIME_MOD,       ONLY : GET_TS_CHEM,      TIMESTAMP_STRING
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : N_TRACERS
      USE TRACER_MOD,     ONLY : TRACER_MW_KG,     XNUMOLAIR
      USE TRACERID_MOD,   ONLY : IDTACET, IDTALD2, IDTALK4, IDTC2H6
      USE TRACERID_MOD,   ONLY : IDTC3H8, IDTCH2O, IDTH2O2, IDTHNO4
      USE TRACERID_MOD,   ONLY : IDTISOP, IDTMACR, IDTMEK,  IDTMP
      USE TRACERID_MOD,   ONLY : IDTMVK,  IDTPMN,  IDTPRPE, IDTR4N2
      USE TRACERID_MOD,   ONLY : IDTRCHO
      USE TRANSFER_MOD,   ONLY : TRANSFER_ZONAL
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, ITS_IN_THE_STRAT

      IMPLICIT NONE

#     include "CMN_SIZE"        ! Size parameters

      ! Local variables
      LOGICAL, SAVE             :: FIRST = .TRUE.

      INTEGER                   :: I, IOS, J, L, N, NN, LMIN
      INTEGER, SAVE             :: MONTHSAVE = 0

      ! Number of photolysis species (currently is 13)
      INTEGER, PARAMETER        :: NSPHOTO = 13

      ! Tracers that undergo photolysis loss in the stratosphere
      INTEGER                   :: SPHOTOID(NSPHOTO) = (/
     &                               3,  8,  9, 10, 11, 12, 13,
     &                              14, 17, 20, 22, 23, 24 /)

      ! Character variables
      CHARACTER(LEN=16 )        :: STAMP
      CHARACTER(LEN=255)        :: FILENAME

      ! REAL*4 arrays -- for reading from binary data files
      REAL*4                    :: ARRAY(1,JGLOB,LGLOB)
      REAL*4, ALLOCATABLE, SAVE :: STRATOH(:,:)
      REAL*4, ALLOCATABLE, SAVE :: SJVALUE(:,:,:)
      REAL*4, ALLOCATABLE, SAVE :: COPROD(:,:)
      REAL*4, ALLOCATABLE, SAVE :: COLOSS(:,:)

      ! REAL*8 variables
      REAL*8                    :: k0,     k1,     k2,  k3, XTAU
      REAL*8                    :: DTCHEM, RDLOSS, T1L, M,  TK, RC

      ! External functions
      REAL*8, EXTERNAL          :: BOXVL

      !=================================================================
      ! SCHEM_ADJ begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - SCHEM_ADJ: Strat chemistry adjoint at ', a )

      !=================================================================
      ! If it is the first call to SCHEM, allocate arrays for reading
      ! data. These arrays are declared SAVE so they will be preserved
      ! between calls.
      !=================================================================
      IF ( FIRST ) THEN
         ALLOCATE( STRATOH( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'STRATOH' )
         STRATOH = 0e0

         ALLOCATE( SJVALUE( JJPAR, LLPAR, NSPHOTO ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'SJVALUE' )
         SJVALUE = 0e0

         ALLOCATE( COPROD( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
         COPROD = 0e0

         ALLOCATE( COLOSS( JJPAR, LLPAR ), STAT=IOS )
         IF ( IOS /= 0 ) CALL ALLOC_ERR( 'COLOSS' )
         COLOSS = 0e0
      ENDIF

      !=================================================================
      ! If it is a new month (or the first call to SCHEM),
      ! do the following:
      !
      ! (1) Read archived J-values and store in SJVALUE
      ! (2) Read archived CO production rates and store in COPROD
      ! (3) Read archived CO loss rates and store in COLOSS
      !
      ! NOTES
      ! (a) All of the above-mentioned data are stored in binary punch
      !     files, for ease of use.
      !
      ! (b) STRATOH, SJVALUE, CO_PROD, and CO_LOSS are now declared
      !     as both ALLOCATABLE and SAVE.  If SCHEM is called, then
      !     data will be declared for these arrays, and the values in
      !     these arrays will be preserved between calls.
      !
      ! (c) If SCHEM is never called (i.e. if you are running another
      !     type of chemistry simulation), then memory never gets
      !     allocated to STRATOH, SJVALUE, CO_PROD, and CO_LOSS.
      !     This saves on computational resources.
      !=================================================================
      ! adj_group:  now use ITS_A_NEW_MONTH
      !IF ( GET_MONTH() /= MONTHSAVE .or. FIRST ) THEN
      !   MONTHSAVE = GET_MONTH()
      IF ( ITS_A_NEW_MONTH() ) THEN

         ! TAU value at the beginning of this month
         XTAU = GET_TAU0( GET_MONTH(), 1, 1985 )

         !==============================================================
         ! Read this month's OH
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'stratOH_200203/stratOH.' //
     &              GET_NAME_EXT()   // '.'                       //
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,
     &                    XTAU,      1,         JGLOB,
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATOH )

         !==============================================================
         ! Read in monthly mean archived J-values
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'stratjv_200203/stratjv.' //
     &              GET_NAME_EXT()   // '.'                       //
     &              GET_RES_EXT()

         DO NN = 1, NSPHOTO
            N = SPHOTOID(NN)

            ! Read data
            CALL READ_BPCH2( FILENAME, 'JV-MAP-$', N,
     &                       XTAU,      1,         JGLOB,
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), SJVALUE(:,:,NN) )
         ENDDO

         !==============================================================
         ! Read in CO production rates
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COprod.' //
     &              GET_NAME_EXT()   // '.'                      //
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 9,
     &                    XTAU,      1,         JGLOB,
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COPROD )

         !==============================================================
         ! Read in CO loss rates
         !==============================================================
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/COloss.' //
     &              GET_NAME_EXT()   // '.'                      //
     &              GET_RES_EXT()

         ! Read data
         CALL READ_BPCH2( FILENAME, 'PORL-L=$', 10,
     &                    XTAU,      1,         JGLOB,
     &                    LGLOB,     ARRAY,     QUIET=.TRUE. )

         ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
         CALL TRANSFER_ZONAL( ARRAY(1,:,:), COLOSS )

      ENDIF

      !=================================================================
      ! Do photolysis for selected tracers with this
      ! month's archived J-values
      !=================================================================

      ! Get the minimum level extent of the ann mean tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NN )
!$OMP+SCHEDULE( DYNAMIC )
      DO NN = 1, NSPHOTO
         N = SPHOTOID(NN)

         DO L = LMIN, LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Only proceed for stratospheric boxes
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

               ! fwd code:
               !STT(I,J,L,N) = STT(I,J,L,N) *
               !               EXP( -SJVALUE(J,L,NN) * DTCHEM )
               ! adj code:
               STT_ADJ(I,J,L,N) = STT_ADJ(I,J,L,N) *
     &                            EXP( -SJVALUE(J,L,NN) * DTCHEM )
            ENDIF

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !print*, 'In schem, done with photolysis'

      !=================================================================
      ! CO is special --
      ! use archived P, L rates for CO chemistry in stratosphere
      !=================================================================
      CALL CO_STRAT_PL_ADJ( COPROD, COLOSS )

      !=================================================================
      ! Reaction with OH -- compute rate constants for each tracer
      !=================================================================
      !print*, 'In schem, before reaction with OH'

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, M, TK, RC, k0, k1, RDLOSS, T1L )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1,    N_TRACERS
      DO L = LMIN, LLPAR
      DO J = 1,    JJPAR
      DO I = 1,    IIPAR

         ! Only proceed for stratospheric boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN

            ! Density of air at grid box (I,J,L) in molec/cm3
            M = AD(I,J,L) / BOXVL(I,J,L) * XNUMOLAIR

            ! Temperature at grid box (I,J,L) in K
            TK = T(I,J,L)

            ! Select proper reaction rate w/ OH for the given tracer
            ! Some rates are temperature or density dependent
            IF ( N == IDTALK4 ) THEN
               RC = 8.20D-12 * EXP(  -300.D0 / TK )

            ELSE IF ( N == IDTISOP ) THEN
               RC = 2.55D-11 * EXP(   410.D0 / TK )

            ELSE IF ( N == IDTH2O2 ) THEN
               RC = 2.90D-12 * EXP(  -160.D0 / TK )

            ELSE IF ( N == IDTACET ) THEN
               RC = 1.70D-12 * EXP(  -600.D0 / TK )

            ELSE IF ( N == IDTMEK  ) THEN
               RC = 2.92D-13 * EXP(   414.D0 / TK )

            ELSE IF ( N == IDTALD2 ) THEN
               RC = 1.40D-12 * EXP( -1860.D0 / TK )

            ELSE IF ( N == IDTRCHO ) THEN
               RC = 2.00D-11

            ELSE IF ( N == IDTMVK  ) THEN
               RC = 4.13D-12 * EXP(   452.D0 / TK )

            ELSE IF ( N == IDTMACR ) THEN
               RC = 1.86D-11 * EXP(  -175.D0 / TK )

            ELSE IF ( N == IDTPMN  ) THEN
               RC = 3.60D-12

            ELSE IF ( N == IDTR4N2 ) THEN
               RC = 1.30D-12

            ELSE IF ( N == IDTPRPE ) THEN
               k0 = 8.0D-27 * ( 300.D0 / TK )**3.5
               k1 = 3.0D-11

               RC = k1 * k0 * M / ( k1 + k0*M )
               RC = RC * 0.5 ** (1 / ( 1 + LOG10( k0*M/k1 )**2 ) )

            ELSE IF ( N == IDTC3H8 ) THEN
               RC = 8.00D-12 * EXP(  -590.D0 / TK )

            ELSE IF ( N == IDTCH2O ) THEN
               RC = 1.00D-12

            ELSE IF ( N == IDTC2H6 ) THEN
               RC =  7.9D-12 * EXP( -1030.D0 / TK )

            ELSE IF ( N == IDTHNO4 ) THEN
               RC = 1.30D-12 * EXP(   380.D0 / TK )

            ELSE IF ( N == IDTMP ) THEN
               RC = 1.14D-12 * EXP(   200.D0 / TK )

            ELSE
               RC = 0d0

            ENDIF

            ! Compute loss with OH based on the rate constants from above
            ! Cap RDLOSS so that it does not exceed 1.0 (bmy, 5/4/00)
            RDLOSS       = RC * STRATOH(J,L) * DTCHEM
            RDLOSS       = MIN( RDLOSS, 1d0 )

            ! T1L is the absolute amount of STT lost to rxn with OH
            ! Subtract T1L from STT
            ! fwd code:
            !T1L          = STT(I,J,L,N) * RDLOSS
            !STT(I,J,L,N) = STT(I,J,L,N) - T1L
            ! adj code:
            STT_ADJ(I,J,L,N) = STT_ADJ(I,J,L,N) * ( 1D0 - RDLOSS )


            ! Oxidation of PRPE as source of ACET with 80% yield
            IF ( N == IDTPRPE ) THEN
               ! fwd code:
               !STT(I,J,L,IDTACET) = STT(I,J,L,IDTACET) +
               !     0.8d0 * T1L *
               !     TRACER_MW_KG(IDTACET) / TRACER_MW_KG(IDTPRPE)
               ! adj code:
               STT_ADJ(I,J,L,IDTACET) = STT_ADJ(I,J,L,IDTACET) *
     &              0.8d0 * RDLOSS *
     &              TRACER_MW_KG(IDTACET) / TRACER_MW_KG(IDTPRPE)
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Set FIRST = .FALSE. -- we have been thru SCHEM at least once now
      FIRST = .FALSE.

      ! Return to calling program
      END SUBROUTINE SCHEM_ADJ
