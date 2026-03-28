SUBROUTINE STMAS_obsmapping(mg,obs_in,mgrid_in,missing_in,obsgrid_out)
!doc=======================================================================
!doc This routine maps observations on to grid using a given Gaussian
!doc function based on location, time, background/covariance and terrain 
!doc and land-factor.
!doc
!doc History: Dec. 2009 created by Yuanfu Xie.
!doc          Aug. 2010 modified by Yuanfu Xie.
!doc          Aug. 18 2011 modified from STMASFC_obsmapping by Hongli Jiang
!doc
!doc Algorithm: Map observation to gridpoints using a Gaussian with a 
!doc radii, topography, land-factor and covariance/background information.
!doc
!doc Assumption: All observations inside the domain and mgrid_in > 2.
!doc
!doc Input
!doc    obs_in:         Structured observations;
!doc	mgrid_in:	Multigrid structured background;
!doc    missing_in:     Value for missing data.
!doc
!doc Output
!doc    obsgrid_out:    Observations mapped onto gridpoints for analysis.
!doc=======================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: mg
  TYPE(STMAS_obs),  INTENT(IN) :: obs_in
  TYPE(STMAS_bkgd), INTENT(IN) :: mgrid_in
  REAL,    INTENT(IN) :: missing_in   ! Missing data

  TYPE(STMAS_gridded_obs), INTENT(INOUT) :: obsgrid_out

  ! Local variables:
  CHARACTER*20, PARAMETER :: header = &
    'STMASFC>obsmapping: '
  INTEGER, PARAMETER :: header_len = 18

  ! midx: index on current multigrid, mfine: on the finest multigrid
! HJ: changed 3 to 4. 8/18/2011
  INTEGER :: midx(-1:2,4),mfine(-1:2,4)        
  INTEGER :: i,j,k,t,io,iv,istatus

  REAL    :: bkg(STMAS_maxvars),xyzt(4)
  REAL    :: r2(5,STMAS_maxvars),s2(STMAS_maxvars)
  
  REAL,ALLOCATABLE :: weights(:,:,:,:,:)       ! Observation weight
  REAL,ALLOCATABLE :: wghtobs(:,:,:,:,:)       ! Weighted obs: sum w*obs
  REAL,ALLOCATABLE :: wghterr(:,:,:,:,:)       ! Weighted err: sum w*err
  REAL    :: dland,dtopo,gaussian            ! Save differences of land/topo
  REAL    :: obs_land,obs_topo               ! land and topo at obs location
real tem

  ! Calculate weight and weighted obs and errors:
! HJ added nfinest(4). 8/18/2011
  ALLOCATE(weights(mgrid_in%numgrid(1),mgrid_in%numgrid(2),mgrid_in%numgrid(3),mgrid_in%numgrid(4),obs_in%numvars), &
           wghtobs(mgrid_in%numgrid(1),mgrid_in%numgrid(2),mgrid_in%numgrid(3),mgrid_in%numgrid(4),obs_in%numvars), &
           wghterr(mgrid_in%numgrid(1),mgrid_in%numgrid(2),mgrid_in%numgrid(3),mgrid_in%numgrid(4),obs_in%numvars), &
           STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'Cannot allocate memory for weights'
    STOP
  ENDIF

  ! Initializing:
  weights = 0.0
  wghtobs = 0.0
  wghterr = 0.0

  ! 1/radius**2: influence radius in x, y, t, z and landwater
  DO iv=1,obs_in%numvars
    r2(1:5,iv) = 1.0/STMAS_radius(1:5,iv)**2
    s2(iv) = 1.0/STMAS_inc(iv)**2
  ENDDO

  DO io=1,obs_in%numobs ! * Loop all obs -- Map them to their nearest grid

    ! Multigrid indices containing the obs according the analysis grid:
! HJ changed i from 3 to 4, and xytz to xyzt 8/18/2011
    DO i=1,4
      midx(0,i) = (obs_in%xyzt(i,io)-STMAS_domain(1,i))/mgrid_in%gridspc(i)+1
      midx(1,i) = midx(0,i)+1
      midx(-1,i) = midx(0,i)-1
      midx(2,i) = midx(0,i)+2
      mfine(-1:2,i) = (midx(-1:2,i)-1)*mgrid_in%incr(i)+1
    ENDDO

    ! ** Loop of all corners of the multigrid box of x,y,t
    DO j=-1,2   ! y
      ! Out of domain:
      IF (midx(j,2) .GT. mgrid_in%numgrid(2) .OR. midx(j,2) .LT. 1) cycle

      ! Distance in y: from obs to grid:
      xyzt(2) = MOD(obs_in%xyzt(2,io)-STMAS_domain(1,2),mgrid_in%gridspc(2)) &
               -j*mgrid_in%gridspc(2)
      DO i=-1,2   ! x
        ! Out of domain:
        IF (midx(i,1) .GT. mgrid_in%numgrid(1) .OR. midx(i,1) .LT. 1) cycle

        ! Distance in x: from obs to grid:
        xyzt(1) = MOD(obs_in%xyzt(1,io)-STMAS_domain(1,1),mgrid_in%gridspc(1)) &
                 -i*mgrid_in%gridspc(1)

        ! Differences in landfactor and topography:
        ! Notce: mgrid topo and land are at the finest multigrid level:
        dland = (mgrid_in%land(mfine(i,1),mfine(j,2))-obs_in%land(io))**2
        dtopo = (mgrid_in%topo(mfine(i,1),mfine(j,2))-obs_in%topo(io))**2

! HJ: change from k to t; all 3 to 4. 
        DO t=0,1   ! t
          ! Out of domain:
          IF (midx(t,4) .GT. mgrid_in%numgrid(4)) cycle

          ! Distance in t: from obs to the current multigrid box:
          xyzt(4) = MOD(obs_in%xyzt(4,io)-STMAS_domain(1,4),mgrid_in%gridspc(4)) &
                   -t*mgrid_in%gridspc(4)

          DO iv=1,obs_in%numvars      ! *** Through all variables

            IF (obs_in%value(io,iv) .NE. missing_in) THEN ! Check if obs is valid
              ! Weight of Gaussian:
              ! The last term uses background pattern to map obs:
              ! where s2 is the variable increment scales**2

              gaussian = exp(-xyzt(1)*xyzt(1)*r2(1,iv) &
                             -xyzt(2)*xyzt(2)*r2(2,iv) &
                             -xyzt(3)*xyzt(3)*r2(3,iv) &
                             -xyzt(4)*xyzt(4)*r2(4,iv) &
              -(obs_in%bkgd(io,iv)-mgrid_in%bkgd(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv))**2*s2(iv) &
                             -dtopo*r2(4,iv)-dland*r2(5,iv))

              ! Save the sum of weights:
              weights(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv) = &
                weights(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv)+1.0 !gaussian

              ! Sum of weighted innovation: W*(obs-bkg)
              wghtobs(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv) = &
                wghtobs(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv)+gaussian* &
                (obs_in%value(io,iv)-obs_in%bkgd(io,iv))

              ! Sum of weighted reciprocal of obs errors: W/err
              wghterr(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv) = & 
                wghterr(mfine(i,1),mfine(j,2),mfine(k,3),mfine(t,4),iv)+gaussian/ &
                obs_in%error(io,iv)

            ENDIF                                       ! Check valid data

          ENDDO                      ! *** Through all variables

        ENDDO
      ENDDO
    ENDDO               ! ** Loop through all corners

  ENDDO                 ! * End of all obs loop -- Map them to their nearest grid


  ! * Save weighted grid observations:
  obsgrid_out%numvars = obs_in%numvars
  obsgrid_out%numobs = 0
! HJ add t loop
  DO t=1,mgrid_in%nfinest(4),mgrid_in%incr(4)
   DO k=1,mgrid_in%nfinest(3),mgrid_in%incr(3)
    DO j=1,mgrid_in%nfinest(2),mgrid_in%incr(2)
      DO i=1,mgrid_in%nfinest(1),mgrid_in%incr(1)

        ! Check if obs exists at this gridpoint:
        o = 0.0
        DO iv=1,obs_in%numvars
          o = o+weights(i,j,k,t,iv)
        ENDDO

        IF (o .GT. 0.0) THEN
          obsgrid_out%numobs = obsgrid_out%numobs+1
          obsgrid_out%value(obsgrid_out%numobs,1:obs_in%numvars) = missing_in
          obsgrid_out%ixyzt(1,obsgrid_out%numobs) = i
          obsgrid_out%ixyzt(2,obsgrid_out%numobs) = j
          obsgrid_out%ixyzt(3,obsgrid_out%numobs) = k
          obsgrid_out%ixyzt(4,obsgrid_out%numobs) = t
          DO iv=1,obs_in%numvars
            ! If the obs is valid: treat a weight of 1.0e-6 as no effect.
            IF (weights(i,j,k,t,iv) .GT. 0.0) THEN
              ! Weighted obs:
              obsgrid_out%value(obsgrid_out%numobs,iv) = &
                mgrid_in%bkgd(i,j,k,t,iv)+wghtobs(i,j,k,t,iv)/weights(i,j,k,t,iv)

              ! Note weighted error array saves the reciprocal of obs error:
              ! Thus: wghterr = sum w/e and the obs error is 
              ! 1/((sum w/e)/sum w) = sum w/sum w/e as calculated as follows:
              obsgrid_out%error(obsgrid_out%numobs,iv) = &
                1.0/weights(i,j,k,t,iv)/max(wghterr(i,j,k,t,iv),1.0e-6)
            ENDIF
          ENDDO
        ENDIF

      ENDDO
    ENDDO
  ENDDO     
  ENDDO        ! * End of saving grid observations

  DEALLOCATE(weights,wghtobs,wghterr)

END SUBROUTINE STMAS_obsmapping

