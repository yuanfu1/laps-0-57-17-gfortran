MODULE STMAS

!doc=============================================================================================
!doc
!doc This module defines constants, variables and processes of STMAS surface analysis.
!doc It offers uni-variate and multi-variate analysis and balance options.
!doc
!doc History: April. 2010 by Yuanfu Xie.
!doc
!doc This module is modified from STMASFC_module.f90 with FC dropped. 
!doc
!doc History: May, 2011 by Hongli Jiang
!doc
!doc=============================================================================================

  IMPLICIT NONE

  !--------------
  ! Constants:
  !--------------
  
  INTEGER, PARAMETER :: STMAS_maxvars = 6         ! Max analysis vars allowed
  INTEGER, PARAMETER :: STMAS_maxdim = 4          ! Max dimension HJ mod. 6/7/2011
  
  REAL, PARAMETER :: STMAS_invalid = -1.0E10      ! Invalid value
  REAL, PARAMETER :: T00 = 273.16                 ! Absolute 0,was named temp_0. HJ 8/31/2011
  REAL, PARAMETER :: knt2ms = 0.51444444          ! Knot to m/s
  REAL, PARAMETER :: mile2m = 1609.0              ! Mile to meter
  REAL, PARAMETER :: mb2pas = 100.0               ! Mb 2 pascal
  REAL, PARAMETER :: inch2m = 0.0254              ! Inch to meter
  REAL, PARAMETER :: gascnt = 287.0               ! Gas constant 
                                                  ! for dry air
  REAL, PARAMETER :: CP = 1004.0                  ! Specific heat at constant pres, was named spheat

  ! Lapse rates: see LAPS mdatlaps.f under sfc: lapses(1): temp; lapses(2): dewp
  REAL, PARAMETER :: lapses(2) = (/-0.01167, -0.007/)         
  
! HJ add 5/12/2011
  real, parameter :: P00 = 1.e5      ! standard pressure 
  real, parameter :: G = 9.8         ! gravity acceleration
  real, parameter :: kappa = 0.287   ! gascnt/CP=287/1004
  real, parameter :: kapp1 = -0.713  ! 1-kappa
  real, parameter :: kapp2 = -1.713  ! 2-kappa
  real, parameter :: rpk=10.541      ! gascnt/p00^kappa
!
! HJ add: ZTOP, 6/29/2011
  real, parameter :: ZTOP=20000.     ! domain top, dz=ztop/(nz-1)=500. 

  !----------------
  ! Variables:
  !----------------

  CHARACTER*4 :: STMAS_varnames(STMAS_maxvars) ! Analysis variable names
  CHARACTER*4 :: STMASFC_varnames(1)           ! Sfc pressure only HJ add: 6/22/2011

  INTEGER :: STMAS_maxobs                      ! Max numbers of observations read in stmasfc_vars.nl
  INTEGER :: STMAS_numvars                     ! Number of analysis variables
  INTEGER :: STMAS_time_window(2)              ! Seconds before (-) and after of current time
  INTEGER :: STMAS_debugging                   ! Option for printing debugging information
  INTEGER :: STMAS_start_grdpts(STMAS_maxdim)  ! Multigrid start levels in x, y, Z and t
  INTEGER :: STMAS_numlevels                   ! Total number of multigrid levels
  INTEGER :: STMAS_maxlevels(STMAS_maxdim)     ! Maximum number levels in each dimension 
  INTEGER :: STMAS_iters                       ! Total number of multigrid levels
  INTEGER :: STMAS_success                     ! Flag for a successful run
  INTEGER :: LAPS_i4time                       ! The current time at i4 format of LAPS

  REAL :: STMAS_domain(2,STMAS_maxdim)         ! Analysis domain, WAS domain(2,3), x,y,Z,t, HJ mod 6/7/2011
  REAL :: STMAS_radius(6,STMAS_maxvars)        ! Obs influence radius in x,y,t,z,flow, and land
  REAL :: STMAS_inc(STMAS_maxvars)             ! Increment sizes used to control obs mapping following bkgd
  REAL :: STMAS_thresholds(STMAS_maxvars)      ! QC threshold check values
  REAL :: STMAS_stddev(STMAS_maxvars)          ! QC factors of std dev check
  REAL :: STMAS_penal(STMAS_maxvars)           ! penality coeffient
  REAL :: STMAS_smooth(STMAS_maxvars)          ! smoothing parameters
  REAL :: LAPS_rdplvl                          ! LAPS reduced pressure level


  !----------------
  ! Dynamic arrays:
  !----------------

  !REAL,ALLOCATABLE :: obserrors(:,:)          ! Each error of obs uses a column.
  !REAL,ALLOCATABLE :: obsxyzt(:,:)            ! Observation locations
!
  
  !----------------
  ! Data type:
  !----------------

  TYPE STMAS_obs                              ! Raw observations
    INTEGER :: numobs,numvars                 ! Number of obs and vars
    CHARACTER*10,ALLOCATABLE :: stnames(:)    ! Station name
    CHARACTER*6, ALLOCATABLE :: types(:)      ! Observation type
    REAL :: weights(20)                       ! Observation weights in cost
    REAL,ALLOCATABLE :: lat(:),lon(:),time(:) ! Observation location and time
    REAL,ALLOCATABLE :: xyzt(:,:)             ! Grid locations: 4*numobs, HJ changed from xytz to xyzt
    REAL,ALLOCATABLE :: value(:,:),error(:,:) ! Obs value and error: numobs*numvars
    REAL,ALLOCATABLE :: bkgd(:,:)             ! Background obs-values at this site
                                              ! i.e. obs-value from backgrounds after 
                                              ! forward operation: numobs*numvars
    REAL,ALLOCATABLE :: topo(:),land(:)       ! Topography and land at obs site
  END TYPE STMAS_obs

  TYPE STMAS_gridded_obs                      ! Observations mapped on to current multigrid
    INTEGER :: numobs,numvars                 ! Number of obs and vars
    INTEGER,ALLOCATABLE :: ixyzt(:,:)         ! Observation location and time. 
                                              ! HJ mod: changed from ixyt(:,:) 6/10/2011
    REAL,   ALLOCATABLE :: value(:,:),error(:,:) ! Obs value and error: numobs*numvars
  END TYPE STMAS_gridded_obs

  TYPE STMAS_bkgd                             ! STMAS background grid
    INTEGER :: numgrid(4)                     ! Number of multigrid points  HJ mod: changed from 3 to 4. 6/8/2011
    INTEGER :: nfinest(4)                     ! Number of the finest multigrid points 
    INTEGER :: incr(4)                        ! Grid number increments (dims)
    REAL    :: gridspc(4)                     ! Grid spacing
    REAL,ALLOCATABLE :: anal(:,:,:,:,:)       ! Saving the analysis HJ mod: from 4 to 5. 6/8/2011
    REAL,ALLOCATABLE :: bkgd(:,:,:,:,:)       ! Background
    REAL,ALLOCATABLE :: lat(:,:),lon(:,:)     ! lat, long
    REAL,ALLOCATABLE :: topo(:,:),land(:,:)   ! Topography and land-water factor
    REAL,ALLOCATABLE :: mapf(:,:)             ! Mapping factor
    REAL,ALLOCATABLE :: grdt(:,:,:,:,:)       ! Saving the gradint HJ mod: try to define grdt
                                              ! in the same structure as in anal. Using grdt
                                              ! instead of grad to avoid confusion since grad
                                              ! is defined in STMAS_minimizer.f90 6/15/2011
!------------
! Define vertical grid spacing ZZ here for now. Horizontal and time
! grid spacing is defined as gridspc(1:4), gridspc(3) should be equivalent to an
! uniform dz, since we are considering a non-uniform dz, gridspc(3) should be
! a 3D array in terrain following vertical coordinate. zz doesn't change with time.   HJ 6/15/2011 
! 
! Update: zz is 1D array, function of k only; jt3 is a function of i,j only based on its
! defination, both zz and jt2 are defined in vertical_z.f90. HJ 6/29/2011
!
    REAL,ALLOCATABLE :: ZZ(:)                 ! vertical grid spacing
    REAL,ALLOCATABLE :: jt1(:,:,:)            ! coordinate transformation in x
    REAL,ALLOCATABLE :: jt2(:,:,:)            ! coordinate transformation in y
    REAL,ALLOCATABLE :: jt3(:,:)              ! coordinate transformation in z, a fn of i,j only 
    REAL,ALLOCATABLE :: cor(:,:)              ! coriolis. HJ 6/20/2011
  END TYPE STMAS_bkgd

  TYPE STMAS_anal    ! STMAS analysis on current grid
    INTEGER :: numgrid(4)
    REAL,ALLOCATABLE :: analysis(:,:,:,:,:) ! HJ mod: changed from 4 to 5. 6/8/2011
  END TYPE STMAS_anal

  ! Global realization of defined type:
  TYPE(STMAS_obs)  :: observations    ! Observations
  TYPE(STMAS_gridded_obs) :: gridobs  ! Gridded observations
  TYPE(STMAS_bkgd) :: STMAS_mgrid     ! Current multigrid
  TYPE(STMAS_bkgd) :: STMAS_final     ! The finest grid
  TYPE(STMAS_bkgd) :: STMAS_tmgrid    ! temporary storage of mgrid
  
  
  !----------------
  ! Namelist:
  !----------------

! HJ added: STMASFC_varnames to read in sfc pressure. 6/22/2011
! STMASFC_varnames is no longer need for sigma_ht. 8/29/2011
!
  NAMELIST /STMAS_configure/STMAS_start_grdpts,STMAS_numlevels, &
                            STMAS_iters,STMAS_time_window,STMAS_debugging
  NAMELIST /STMAS_variables/STMAS_maxobs,STMAS_numvars,STMAS_varnames, &
                            STMASFC_varnames, &
                            STMAS_radius,STMAS_inc, &
                            STMAS_thresholds,STMAS_stddev,STMAS_penal,STMAS_smooth

  !----------------
  ! Routines:
  !----------------

CONTAINS

  SUBROUTINE read_namelist

  !doc==============================================================================================
  !doc  This routine reads in STMAS namelists.
  !doc
  !doc  History: April. 2010 by Yuanfu Xie.
  !doc==============================================================================================

    IMPLICIT NONE

    ! Local variables:
    CHARACTER*40, PARAMETER :: header = 'STMAS>read_namelist: '
    INTEGER :: i,istatus

    ! Read in configure namelist:
    OPEN(unit=11,file='stmas_conf.nl')
    READ(11,NML=STMAS_configure,IOSTAT=istatus)
    WRITE(*,1) header,STMAS_start_grdpts,STMAS_numlevels, STMAS_time_window
    CLOSE(11)
  1 FORMAT(A23,'Start levels: ',4I3,' Number of multigrid levels: ',i2,/,8x, &
          'Seconds before analysis time: ',i8,' Seconds after: ',i8)

! Read in analysis variable namelist:
    OPEN(unit=11,file='stmas_vars.nl')
    READ(11,NML=STMAS_variables,IOSTAT=istatus)
!
!
    WRITE(*,2) header,STMAS_numvars
  2 FORMAT(A23,'Number of analysis variables: ',i3)
! Check if too many analysis variables:
    IF (STMAS_numvars .GT. STMAS_maxvars) THEN
      PRINT*,'Too many analysis variables, increase STMAS_maxvars first'
      STOP
    ENDIF
!
    DO i=1,STMAS_numvars
      WRITE(*,3) STMAS_varnames(i),STMAS_radius(1:5,i),STMAS_thresholds(i)
    ENDDO
  3 FORMAT(8x,A4,': Influence radius in',/,5x,'x',4x,'|',4x,'y',4x,'|',4x,'t',4x,'|',4x,'z',4x,'|', &
           2x,'land',2x,'|',/,5e10.2, /,7x,' Threshold: ',e10.2,/)

  CLOSE(11)

  END SUBROUTINE read_namelist


  SUBROUTINE STMAS_memo_alloc
!doc==============================================================================================
!doc
!doc This routine allocates all needed memory for STMAS analysis:
!doc
!doc History: March 2010 by Yuanfu Xie.
!doc
!doc==============================================================================================

    IMPLICIT NONE

    ! Local variables:
    CHARACTER*41,PARAMETER :: header = 'STMAS>STMAS_anal>STMAS_memo_alloc: '

    INTEGER :: mm(4),num,istatus

    print*,'memo_alloc',STMAS_final%numgrid

    mm(1:4) = STMAS_final%numgrid(1:4)
    num = STMAS_final%numgrid(1)*STMAS_final%numgrid(2)*STMAS_final%numgrid(3)*STMAS_final%numgrid(4)

    ! Allocate memory:
! HJ add mm(4), 6/7/2011
! HJ added STMAS_final%ZZ(mm(3)), STMAS_mgrid%ZZ(mm(1),mm(2),mm(3)),
! STMAS_final%grdt, STMAS_mgrid%grdt, jt1, jt2, and jt3. 6/15/2011 
! add cor(mm(1),mm(2)), HJ 6/20/2011
! add bkgd2d(mm(1),mm(2),mm(4)), HJ 6/22/2011, no longer needed, removed. HJ 8/29/2011
! comment out allocation to STMAS_mgrid here, and will be allocated before STMAS_minimizer. 6/29/2011
print*,'Memory allocation: ',STMAS_maxobs,STMAS_numvars
    ALLOCATE(STMAS_final%lat (mm(1),mm(2)), &
             STMAS_final%lon (mm(1),mm(2)), &
             STMAS_final%topo(mm(1),mm(2)), &
             STMAS_final%land(mm(1),mm(2)), &
             STMAS_final%mapf(mm(1),mm(2)), &
             STMAS_final%anal(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
             STMAS_final%bkgd(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
             STMAS_final%grdt(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
             STMAS_final%zz(mm(3)), &
             STMAS_final%jt1(mm(1),mm(2),mm(3)), &
             STMAS_final%jt2(mm(1),mm(2),mm(3)), &
             STMAS_final%jt3(mm(1),mm(2)), &
             STMAS_final%cor(mm(1),mm(2)), &
             observations%stnames(STMAS_maxobs),observations%types(STMAS_maxobs), &
             observations%lat(STMAS_maxobs),observations%lon(STMAS_maxobs), &
             observations%time(STMAS_maxobs),observations%xyzt(4,STMAS_maxobs), &
             observations%value(STMAS_maxobs,STMAS_numvars), &
             observations%error(STMAS_maxobs,STMAS_numvars), &
             observations%bkgd(STMAS_maxobs,STMAS_numvars), &
             observations%topo(STMAS_maxobs),observations%land(STMAS_maxobs), &
             gridobs%ixyzt(4,num),gridobs%value(num,STMAS_numvars), &
             gridobs%error(num,STMAS_numvars), &
             STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*
      PRINT*,header,'Problem allocating memo'
      STOP
    ENDIF    

  END SUBROUTINE STMAS_memo_alloc


  SUBROUTINE STMAS_memo_release
  !doc------------------------------------------------------------------
  !doc This routine deallocates all needed memory for STMAS analysis:
  !doc
  !doc History: March 2010 by Yuanfu Xie.
  !doc------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables:
    CHARACTER*45,PARAMETER :: header = &
      'STMAS>STMAS_anal>STMAS_memo_release: '

    INTEGER :: istatus

! HJ added ZZ. 6/15/2011
! HJ added STMAS_final%zz, STMAS_mgrid%zz, STMAS_final%grdt, and STMAS_mgrid%grdt. 6/15/2011
    DEALLOCATE(STMAS_final%lat, STMAS_final%lon, STMAS_final%topo, &
               STMAS_final%land,STMAS_final%mapf,STMAS_final%anal, &
               STMAS_final%bkgd,STMAS_final%grdt,STMAS_final%zz, &
               STMAS_final%jt1, STMAS_final%jt2, STMAS_final%jt3, &
               STMAS_final%cor, &
               observations%stnames,observations%types, &
               observations%lat,observations%lon,observations%time, &
               observations%xyzt,observations%value,observations%error, &
               observations%bkgd,observations%topo,observations%land, &
               gridobs%ixyzt,gridobs%value,gridobs%error,&
               STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*,header,'Problem deallocating memo'
      STOP
    ENDIF

  END SUBROUTINE STMAS_memo_release

  SUBROUTINE STMAS_memo_alloc_mgrid(mg)
!doc==============================================================================================
!doc
!doc This routine is similar to STMAS_memo_alloc except for allocation of mgrid:
!doc
!doc History: 6/29/2011 by Hongli Jiang.
!doc
!doc==============================================================================================

    IMPLICIT NONE
    integer :: mg

    ! Local variables:
    CHARACTER*41,PARAMETER :: header = 'STMAS_anal>STMAS_memo_alloc_mgrid: '

    INTEGER :: mm(4),num,istatus,i

! The STMAS_mgrid%numgrid has to be initialized here before the mg loop in STMAS_analysis.f90. HJ 6/30/2011
! The coarset grid. 
    if(mg .eq. 1)then 
      STMAS_mgrid%numgrid = STMAS_start_grdpts
    else
      DO i=1,4
! Refining grid for the last maxlevel in each direction:
       IF (mg .GT. STMAS_numlevels-STMAS_maxlevels(i)+1) THEN
! No finer grid than the final analysis
         IF (2*(STMAS_mgrid%numgrid(i)-1)+1 .LE. STMAS_final%numgrid(i)) THEN
           STMAS_mgrid%numgrid(i) = 2*(STMAS_mgrid%numgrid(i)-1)+1
           STMAS_mgrid%incr(i) = STMAS_mgrid%incr(i)/2
         ENDIF
       ENDIF
      ENDDO
    endif

    mm(1:4) = STMAS_mgrid%numgrid(1:4)
    num = STMAS_mgrid%numgrid(1)*STMAS_mgrid%numgrid(2)*STMAS_mgrid%numgrid(3)*STMAS_mgrid%numgrid(4)

    print*,'mg, memo_alloc_mgrid',mg, STMAS_mgrid%numgrid

    ! Allocate memory for mgrid and tmgrid:
    ALLOCATE( STMAS_mgrid%lat (mm(1),mm(2)),STMAS_mgrid%lon (mm(1),mm(2)), &
              STMAS_mgrid%topo(mm(1),mm(2)),STMAS_mgrid%land(mm(1),mm(2)), &
              STMAS_mgrid%mapf(mm(1),mm(2)), &
              STMAS_mgrid%anal(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_mgrid%bkgd(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_mgrid%grdt(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_mgrid%zz(mm(3)), &
              STMAS_mgrid%jt1(mm(1),mm(2),mm(3)), &
              STMAS_mgrid%jt2(mm(1),mm(2),mm(3)), &
              STMAS_mgrid%jt3(mm(1),mm(2)), &
              STMAS_mgrid%cor(mm(1),mm(2)), STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*
      PRINT*,header,'Problem allocating memo'
      STOP
    ENDIF    

  END SUBROUTINE STMAS_memo_alloc_mgrid
!
  SUBROUTINE STMAS_memo_release_mgrid
  !doc------------------------------------------------------------------
  !doc Modified from STMAS_memo_release for multigrid. 
  !doc 6/30/2011 Hongli Jiang
  !doc------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables:
    CHARACTER*41,PARAMETER :: header = 'STMAS_anal>STMAS_memo_release_mgrid: '

    INTEGER :: istatus

    DEALLOCATE(STMAS_mgrid%lat, STMAS_mgrid%lon, STMAS_mgrid%topo, &
               STMAS_mgrid%land,STMAS_mgrid%mapf,STMAS_mgrid%anal, &
               STMAS_mgrid%bkgd,STMAS_mgrid%grdt,STMAS_mgrid%zz, &
               STMAS_mgrid%jt1, STMAS_mgrid%jt2, STMAS_mgrid%jt3, &
               STMAS_mgrid%cor, STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*,header,'Problem deallocating memo'
      STOP
    ENDIF

  END SUBROUTINE STMAS_memo_release_mgrid
!
  SUBROUTINE STMAS_memo_alloc_tmgrid(mg)
!doc==============================================================================================
!doc
!doc This routine is similar to STMAS_memo_alloc_mgrid except for allocation of tmgrid:
!doc
!doc History: 8/4/2011 by Hongli Jiang.
!doc
!doc==============================================================================================

    IMPLICIT NONE
    integer :: mg

    ! Local variables:
    CHARACTER*41,PARAMETER :: header = 'STMAS_anal>STMAS_memo_alloc_tmgrid: '

    INTEGER :: mm(4),num,istatus,i

! The STMAS_mgrid%numgrid has to be initialized here before the mg loop in STMAS_analysis.f90. HJ 6/30/2011
! The coarset grid. 
    if(mg .eq. 1)then 
      STMAS_tmgrid%numgrid = STMAS_start_grdpts
    else
      DO i=1,4
! Refining grid for the last maxlevel in each direction:
       IF (mg .GT. STMAS_numlevels-STMAS_maxlevels(i)+1) THEN
! No finer grid than the final analysis
         IF (2*(STMAS_tmgrid%numgrid(i)-1)+1 .LE. STMAS_final%numgrid(i)) THEN
           STMAS_tmgrid%numgrid(i) = 2*(STMAS_tmgrid%numgrid(i)-1)+1
           STMAS_tmgrid%incr(i) = STMAS_tmgrid%incr(i)/2
         ENDIF
       ENDIF
      ENDDO
    endif

    mm(1:4) = STMAS_tmgrid%numgrid(1:4)
    num = STMAS_tmgrid%numgrid(1)*STMAS_tmgrid%numgrid(2)*STMAS_tmgrid%numgrid(3)*STMAS_tmgrid%numgrid(4)

    ! Allocate memory for tmgrid:
    ALLOCATE( STMAS_tmgrid%lat (mm(1),mm(2)),STMAS_tmgrid%lon (mm(1),mm(2)), &
              STMAS_tmgrid%topo(mm(1),mm(2)),STMAS_tmgrid%land(mm(1),mm(2)), &
              STMAS_tmgrid%mapf(mm(1),mm(2)), &
              STMAS_tmgrid%anal(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_tmgrid%bkgd(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_tmgrid%grdt(mm(1),mm(2),mm(3),mm(4),STMAS_numvars), &
              STMAS_tmgrid%zz(mm(3)), &
              STMAS_tmgrid%jt1(mm(1),mm(2),mm(3)), &
              STMAS_tmgrid%jt2(mm(1),mm(2),mm(3)), &
              STMAS_tmgrid%jt3(mm(1),mm(2)), &
              STMAS_tmgrid%cor(mm(1),mm(2)), STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*
      PRINT*,header,'Problem allocating tmgrid'
      STOP
    ENDIF    

  END SUBROUTINE STMAS_memo_alloc_tmgrid
!
  SUBROUTINE STMAS_memo_release_tmgrid
  !doc------------------------------------------------------------------
  !doc Modified from STMAS_memo_release for tmgrid. 
  !doc 8/4/2011 Hongli Jiang
  !doc------------------------------------------------------------------

    IMPLICIT NONE

    ! Local variables:
    CHARACTER*41,PARAMETER :: header = 'STMAS_anal>STMAS_memo_release_tmgrid: '

    INTEGER :: istatus

    DEALLOCATE(STMAS_tmgrid%lat, STMAS_tmgrid%lon, STMAS_tmgrid%topo, &
               STMAS_tmgrid%land,STMAS_tmgrid%mapf,STMAS_tmgrid%anal, &
               STMAS_tmgrid%bkgd,STMAS_tmgrid%grdt,STMAS_tmgrid%zz, &
               STMAS_tmgrid%jt1, STMAS_tmgrid%jt2, STMAS_tmgrid%jt3, &
               STMAS_tmgrid%cor, STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT*,header,'Problem deallocating tmgrid'
      STOP
    ENDIF

  END SUBROUTINE STMAS_memo_release_tmgrid

  SUBROUTINE STMAS_copyto(grid_from,grid_to,grid_final)
  !doc------------------------------------------------------------------
  !doc  This routine copies a grid data type to another defined here.
  !doc
  !doc  History: Feb. 2011 by Yuanfu Xie.
  !doc
  !doc  Modified from STMASFC_CopyGrid to temorarily store STMAS_mgrid$bkgd and 
  !doc  STMAS_mgrid%anal before bi-linear interpolation in STMAS_analysis.f90
  !doc  7/13/2011 Hongli Jiang
  !doc------------------------------------------------------------------

    IMPLICIT NONE

    TYPE(STMAS_bkgd) :: grid_from,grid_to,grid_final

    ! Local variables:
    INTEGER :: i,j,k,t,iv,i2,j2,k2,t2

    print*,'copyto, from:',grid_from%numgrid
    print*,'copyto, to:',grid_to%numgrid
    print*,'copyto, final:',grid_final%numgrid

    ! Analysis & background:
    DO iv=1,STMAS_numvars
     DO t=1,grid_from%numgrid(4)
      DO k=1,grid_from%numgrid(3)
       DO j=1,grid_from%numgrid(2)
        DO i=1,grid_from%numgrid(1)
         i2=(i-1)*grid_from%incr(1)+1
         j2=(j-1)*grid_from%incr(2)+1
         k2=(k-1)*grid_from%incr(3)+1
         t2=(t-1)*grid_from%incr(4)+1
         grid_to%anal(i,j,k,t,iv) = grid_from%anal(i,j,k,t,iv) - grid_final%bkgd(i2,j2,k2,t2,iv)
         if(k .eq.1 .and. t .eq. 1 .and. iv .eq. 5 .and. &
            i .eq. grid_from%numgrid(1) .and. j.eq. grid_from%numgrid(2))then
           print*,'i2',i2,j2,grid_to%anal(i,j,k,t,iv),grid_final%bkgd(i2,j2,k2,t2,iv)
         endif
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO

  END SUBROUTINE STMAS_copyto

END MODULE STMAS

