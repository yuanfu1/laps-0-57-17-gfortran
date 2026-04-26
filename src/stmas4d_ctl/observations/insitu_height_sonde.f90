SUBROUTINE insitu_height_sonde(nt,i4times,numgrid,lat,lon,topo,top, &
                               domain,verbal)

!doc==================================================================
!doc
!doc  This routine reads in observations observations and assign them 
!od   to the observation arrays following height vertical coordinate.
!doc
!doc  This routine uses LAPS ingest interface.
!doc
!doc  Input:
!doc        nt:       number of laps cycle times
!doc        i4times:  I4times for these cycle times
!doc        numgrid:  Numbers of gridpoints in all directions
!doc        lat/lon:  Lat and lon over the domain grid
!doc        topo:     Terrain
!doc        z:        Heights at all vertical levels
!doc        domain:   Analysis domain
!doc        verbal:   Option for printing debugging messages
!doc
!doc  Output:
!doc        observations:   output observations in STMAS data structure
!doc
!doc  NOTE: observation location: (1,numgrid) and vertical grid is uniform
!doc
!doc  History: June 2011 by Yuanfu Xie.
!doc==================================================================

  USE MEM_NAMELIST
  USE STMAS
  ! USE STMAS4d_obs

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt,i4times(nt),numgrid(4)
  REAL,    INTENT(IN) :: top,domain(2,4)
  REAL,    INTENT(IN) :: lat(numgrid(1),numgrid(2)), &
                         lon(numgrid(1),numgrid(2)), &
                         topo(numgrid(1),numgrid(2)),verbal

  ! Local variables:
  CHARACTER :: static*200,filename*200,ext*3
  CHARACTER,PARAMETER :: header*28 = 'STMAS4D>insitu_height_sonde: '
  INTEGER :: io,ip,il,isnd,it,length,istatus,k,ibeg,iend,observations_now,i,j
  INTEGER :: ntime_out,nspace_out,invalid
  REAL    :: x,y,ssh,sigma,vmissing

  ! Height to sigma conversion function:
  INCLUDE 'height2sigma.f90'

  ! For sonde:
  CHARACTER,ALLOCATABLE :: obs_type(:)*8,c5name(:)*5
  INTEGER :: n_snd,nobs
  INTEGER,ALLOCATABLE :: n_snd_lvl(:)    ! Number of profilers and levels
  INTEGER,ALLOCATABLE :: obs_tim(:,:)   ! Profiler obs times
  REAL,   ALLOCATABLE :: snd_lat(:,:),snd_lon(:,:),snd_elv(:) ! Obs lat/lon/elv
  REAL,   ALLOCATABLE :: obs_h(:,:),obs_u(:,:),obs_v(:,:)
  REAL,   ALLOCATABLE :: obs_t(:,:),obs_p(:,:),obs_d(:,:)
  REAL,   ALLOCATABLE :: press3d(:,:,:)

  ! Start and end time in integers:
  ibeg = domain(1,4)
  iend = domain(2,4)

  io = 12    ! I/O channel for reading profiler data

  ! Use LAPS missing value:
  CALL GET_R_MISSING_DATA(vmissing,istatus)

  CALL GET_DIRECTORY('static',static,length)
  filename = static(1:length)//'/wind.nl'
  CALL READ_NAMELIST_LAPS ('wind',filename) ! Get LAPS MAX_SND_GRID and LEVELS
  
  ! Allocate memory for profiler data:
  ALLOCATE(obs_type(MAX_SND_GRID),c5name(MAX_SND_GRID), &
           obs_tim(MAX_SND_GRID,MAX_SND_LEVELS),n_snd_lvl(MAX_SND_GRID), &
           snd_lat(MAX_SND_GRID,MAX_SND_LEVELS),snd_lon(MAX_SND_GRID,MAX_SND_LEVELS), &
           snd_elv(MAX_SND_GRID),obs_h(MAX_SND_GRID,MAX_SND_LEVELS), &
           obs_u(MAX_SND_GRID,MAX_SND_LEVELS),obs_v(MAX_SND_GRID,MAX_SND_LEVELS), &
           obs_t(MAX_SND_GRID,MAX_SND_LEVELS),obs_p(MAX_SND_GRID,MAX_SND_LEVELS), &
           obs_d(MAX_SND_GRID,MAX_SND_LEVELS), &
           press3d(numgrid(1),numgrid(2),numgrid(3)),STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*
    PRINT*,header,' Failed to allocate memory for reading profiler data'
    STOP
  ENDIF

  ! Record the number of observations data:
  observations_now = observations%numobs

  ! Check the index for background pressure field:
  DO ip=1,STMAS_numvars
    IF (STMAS_varnames(ip) .EQ. 'P3') exit
  ENDDO
  IF (ip .GT. STMAS_numvars) THEN
    print*,header,'No pressure background for reading in sonde data'
    return
  ENDIF    

  ! All time frames:
  ext = 'snd'
  DO it=1,nt         !%%%%%%%%%%%% through all time frames

print*,'RRR ',it
    ntime_out = 0 ! Count data out of time window
    nspace_out = 0 ! Count data out of space domain
    invalid = 0

    n_snd = 0     ! Initialized to zero
    n_snd_lvl = 0 ! Initialize to zero as LAPS may return non-zero value at invalid SND.
    CALL READ_SND_DATA2(io,i4times(it),ext,MAX_SND_GRID,MAX_SND_LEVELS, &
                       lat,lon,numgrid(1),numgrid(2),numgrid(3),STMAS_final%bkgd(1,1,1,it,ip), &
                       .true.,3,n_snd,snd_elv,n_snd_lvl,c5name, &
                       obs_type,obs_h,obs_p,obs_u,obs_v, &
                       obs_t,obs_d,snd_lat,snd_lon,obs_tim,istatus)
    ! Note: currently all profiler sfc obs are off;
    ! Surface obs:
    ! XXX

print*,'NNN: ',n_snd
    ! Assign the up air data to STMAS4D observation arrays:
    DO isnd=1,n_snd         !%%%%%%%%%%%% through all sondes
      
      IF (n_snd_lvl(isnd) .GT. 0) THEN   !%%%%%%%%%%%% if there are at least one level data
        ! Up air data is available:
        DO il=1,n_snd_lvl(isnd)          !%%%%%%%%%%%% through levels

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Check whether obs inside horizontal and time domain:
      
          ! Check if this profiler inside our horizontal domain:
          CALL LATLON_TO_RLAPSGRID(snd_lat(isnd,il),snd_lon(isnd,il),lat,lon,numgrid(1),numgrid(2), &
                                    x,y,istatus)
          IF (x .LT. 1 .OR. x .GT. numgrid(1) .OR. &
              y .LT. 1 .OR. y .GT. numgrid(2)) THEN
            nspace_out = nspace_out+1
            IF (verbal .EQ. 1) print*,header,' Station: ',c5name(isnd),' at level ',il,' Out of domain'
            CYCLE
          ENDIF

          ! Check if vertical location is available:
          IF (obs_h(isnd,il) .EQ. vmissing .AND. obs_p(isnd,il) .EQ. vmissing) THEN
            IF (verbal .EQ. 1) print*,header,' Station: ',c5name(isnd),' at level ',il,' No vertical info'
            CYCLE
          ENDIF
          ! Check if height is within analysis domain:
          IF (obs_h(isnd,il) .NE. vmissing) THEN
              ! Out of vertical domain:
              IF ((obs_h(isnd,il) .LT. topo(INT(x),INT(y))) .OR. obs_h(isnd,il) .GT. top) THEN
                IF (verbal .EQ. 1) print*,header,' Station: ',c5name(isnd),' at level ',il,' Out of height'
                CYCLE
              ENDIF
            ENDIF

          ! Check if this profiler inside our horizontal domain:
          IF (obs_tim(isnd,il) .LT. ibeg .OR. obs_tim(isnd,il) .GT. iend) THEN
            ntime_out = ntime_out+1
            IF (verbal .EQ. 1) print*,header,' Station: ',c5name(isnd),' at level ',il,' out of window'
            CYCLE
          ENDIF

          ! End of check
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          observations%numobs = observations%numobs+1
          observations%value(observations%numobs+1,1:5) = STMAS_invalid
          observations%error(observations%numobs+1,1:5) = STMAS_invalid
          nobs = 0

          ! Sonde has wind, t, dew, p/h: (see STMAS4d_obs for variable order):

          ! U and V obs: READ_SND_DATA2 returns both U and V in grid-UV and no conversion needed
          IF (obs_u(isnd,il) .NE. vmissing .AND. obs_v(isnd,il) .NE. vmissing) THEN
            nobs = nobs+2
            observations%value(observations%numobs,1) = obs_u(isnd,il)
            observations%error(observations%numobs,1) = 0.5             ! Hard coded error temporarily
            observations%value(observations%numobs,2) = obs_v(isnd,il)
            observations%error(observations%numobs,2) = 0.5             ! Hard coded error temporarily
          ENDIF
          ! P obs:
          IF (obs_p(isnd,il) .NE. vmissing) THEN
            nobs = nobs+1
            observations%value(observations%numobs,3) = obs_p(isnd,il)
            observations%error(observations%numobs,3) = 5.0             ! Hard coded error temporarily
          ENDIF
          ! T obs:
          IF (obs_t(isnd,il) .NE. vmissing) THEN
            nobs = nobs+1
            observations%value(observations%numobs,4) = obs_t(isnd,il)
            observations%error(observations%numobs,4) = 0.5             ! Hard coded error temporarily
          ENDIF
          ! Dew obs:
          IF (obs_d(isnd,il) .NE. vmissing .AND. obs_p(isnd,il) .NE. vmissing) THEN
            nobs = nobs+1
            observations%value(observations%numobs,5) = ssh(obs_p(isnd,il),obs_d(isnd,il))  ! G/Kg
            observations%error(observations%numobs,5) = 0.5             ! Hard coded error temporarily
          ENDIF

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Check if a valid obs:

          IF (nobs .EQ. 0) THEN
            observations%numobs = observations%numobs-1
            invalid = invalid+1
          ELSE
            observations%stnames(observations%numobs) = c5name(isnd)

            ! Save the horizontal location:
            observations%xyzt(1,observations%numobs) = x
            observations%xyzt(2,observations%numobs) = y

            ! For obs, height may be missing if pressure is observed:
            observations%xyzt(3,observations%numobs) = STMAS_invalid

            ! Find the vertical layer containing the height:
            sigma = height2sigma(obs_h(isnd,il),INT(x),INT(y))
            observations%xyzt(3,observations%numobs) = sigma/(domain(2,3)-domain(1,3))*(numgrid(3)-1)+1

            ! Save the obs time:
            observations%xyzt(4,observations%numobs) = &
              (obs_tim(isnd,il)-domain(1,4))/(domain(2,4)-domain(1,4))*(numgrid(4)-1)+1

          ENDIF

          ! End of check valid obs
          !+++++++++++++++++++++++++++++++++mrenew.o5777603+++++++++++++++++++++++++++++++++++++++++++++++++++

        ENDDO          !%%%%%%%%%%%% through levels

      ENDIF          !%%%%%%%%%%%% if there are at least one level data

    ENDDO         !%%%%%%%%%%%% through all profilers
    
  ENDDO        !%%%%%%%%%%%% through all time frames

  ! Deallocate memory for profiler data:
  DEALLOCATE(obs_type,c5name,obs_tim, &
             snd_lat,snd_lon,snd_elv,&
             obs_h,obs_u,obs_v,obs_t,obs_p,obs_d, &
             STAT=istatus)

  ! Print out sonde info:
  print*,''
  print*,'+----------------------------------------------------------+'
  WRITE(*,1) ntime_out,it
1 FORMAT(' | Number of sonde obs out of window: ',i4,' at time frame ',i2,' |')
  WRITE(*,2) nspace_out,it
2 FORMAT(' | Number of sonde obs out of domain: ',i4,' at time frame ',i2,' |')
  WRITE(*,3) invalid,it
3 FORMAT(' | Number of invalid sonde obs:       ',i4,' at time frame ',i2,' |')
  WRITE(*,4) ,observations%numobs-observations_now
4 FORMAT(' | Total sonde observations:        ',i6,'                  |')

!  print*,'| Total sonde observations: ',observations%numobs-observations_now,'                 |'
  print*,'+----------------------------------------------------------+'

  PRINT*
  PRINT*,'+----------------------------------+'
  PRINT*,'|   End of reading in sonde data   |'
  PRINT*,'+----------------------------------+'

END SUBROUTINE insitu_height_sonde

