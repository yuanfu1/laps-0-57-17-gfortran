SUBROUTINE insitu_height_profiler(nt,i4times,numgrid,lat,lon, &
                                  topo,top,domain)

!doc==================================================================
!doc
!doc  This routine reads in insitu observations and assign them to the
!doc  observation arrays following the height vertical coordinate.
!doc
!doc  This routine uses LAPS ingest interface.
!doc
!doc  Input:
!doc        nt:       Number of laps cycle times for ingesting data
!doc        i4times:  I4times for these cycle times
!doc        vmissing: Default value for missing data
!doc        numgrid:  Numbers of gridpoints in all directions
!doc        lat/lon:  Lat and lon grids
!doc        z:        Heights at each sigma grid
!doc        domain:   Analysis domain
!doc
!doc  Output:
!doc        observations:   A structured STMAS observation data type 
!doc                  for insitu obs
!doc
!doc  NOTE: obs location (1,numgrid) and vertical grid is uniform
!doc
!doc  History: June 2011 by Yuanfu Xie.
!doc==================================================================

  USE MEM_NAMELIST
  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nt,i4times(nt),numgrid(4)
  REAL,    INTENT(IN) :: top,domain(2,4)
  REAL,    INTENT(IN) :: lat(numgrid(1),numgrid(2)), &
                         lon(numgrid(1),numgrid(2)), &
                        topo(numgrid(1),numgrid(2))

  ! Local variables:
  CHARACTER :: static*200,filename*200
  CHARACTER,PARAMETER :: header*23 = 'STMAS4D>insitu_height: '
  INTEGER :: io,ip,il,it,length,istatus,k,ibeg,iend
  REAL    :: x,y,sigma,vmissing

  ! Height to sigma function:
  INCLUDE 'height2sigma.f90'

  ! For profilers:
  CHARACTER,ALLOCATABLE :: obs_type(:)*8,c5name(:)*5
  INTEGER :: n_pr
  INTEGER,ALLOCATABLE :: n_pr_lvl(:)    ! Number of profilers and levels
  INTEGER,ALLOCATABLE :: obs_time(:)    ! Profiler obs times
  REAL,   ALLOCATABLE :: pr_lat(:),pr_lon(:),pr_elv(:) ! Obs lat/lon/elv
  REAL,   ALLOCATABLE :: obs_h(:,:),obs_u(:,:),obs_v(:,:),obs_rms(:,:)
  REAL,   ALLOCATABLE :: obs_tsfc(:),obs_psfc(:),obs_rhsfc(:), &
                         obs_usfc(:),obs_vsfc(:)

  ! Start and end time in integers:
  ibeg = domain(1,4)
  iend = domain(2,4)

  io = 12    ! I/O channel for reading profiler data

  ! Use LAPS missing value:
  CALL GET_R_MISSING_DATA(vmissing,istatus)

  CALL GET_DIRECTORY('static',static,length)
  filename = static(1:length)//'/wind.nl'
  CALL READ_NAMELIST_LAPS ('wind',filename) ! Get LAPS MAX_PR and LEVELS
  
  ! Allocate memory for profiler data:
  ALLOCATE(obs_type(MAX_PR),c5name(MAX_PR), &
           obs_time(MAX_PR),n_pr_lvl(MAX_PR), &
           pr_lat(MAX_PR),pr_lon(MAX_PR),pr_elv(MAX_PR),&
           obs_h(MAX_PR,MAX_PR_LEVELS),obs_u(MAX_PR,MAX_PR_LEVELS), &
           obs_v(MAX_PR,MAX_PR_LEVELS),obs_rms(MAX_PR,MAX_PR_LEVELS), &
           obs_tsfc(MAX_PR),obs_psfc(MAX_PR),obs_rhsfc(MAX_PR), &
           obs_usfc(MAX_PR),obs_vsfc(MAX_PR),STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*
    PRINT*,header,' Failed to allocate memory for reading profiler data'
    STOP
  ENDIF

  ! All time frames:
  DO it=1,nt         !%%%%%%%%%%%% through all time frames

    CALL READ_PRO_DATA(io,i4times(it),'pro',MAX_PR,MAX_PR_LEVELS, &
                       n_pr,n_pr_lvl,pr_lat,pr_lon,pr_elv,c5name, &
                       obs_time,obs_type,obs_h,obs_u,obs_v, &
                       obs_rms,obs_tsfc,obs_psfc,obs_rhsfc, &
                       obs_usfc,obs_vsfc,istatus)
    ! Note: currently all profiler sfc obs are off;
    ! Surface obs:
    ! XXX

    ! Assign the up air data to STMAS4D observation arrays:
    DO ip=1,n_pr         !%%%%%%%%%%%% through all profilers

      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Check whether obs inside horizontal and time domain:
      
      ! Check if this profiler inside our horizontal domain:
      CALL LATLON_TO_RLAPSGRID(pr_lat(ip),pr_lon(ip), &
                               lat,lon,numgrid(1),numgrid(2), &
                               x,y,istatus)
      IF (x .LT. 1 .OR. x .GT. numgrid(1) .OR. &
          y .LT. 1 .OR. y .GT. numgrid(2)) CYCLE

      ! Check if this profiler inside our temporal domain:
      IF (obs_time(ip) .LT. ibeg .OR. obs_time(ip) .GT. iend) CYCLE
      
      ! End of check
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      IF (n_pr_lvl(ip) .GT. 0) THEN   !%%%%%%%%%%%% if there are at least one level data
        ! Up air data is available:
        DO il=1,n_pr_lvl(ip)          !%%%%%%%%%%%% through levels

          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! Check whether obs inside our vertical domain:
          IF (obs_h(ip,il) .LT. topo(INT(x),INT(y)) .OR. obs_h(ip,il) .GT. top) CYCLE
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          observations%numobs = observations%numobs+1
          observations%value(observations%numobs,:) = STMAS_invalid

          ! Profiler has wind only:

          ! U-V obs:
          IF (obs_u(ip,il) .NE. vmissing .AND. obs_v(ip,il) .NE. vmissing) THEN
            observations%value(observations%numobs,1) = obs_u(ip,il)
            observations%error(observations%numobs,1) = 0.5    ! Hard coded error temporarily
            observations%value(observations%numobs,2) = obs_v(ip,il)
            observations%error(observations%numobs,2) = 0.5    ! Hard coded error temporarily
            observations%stnames(observations%numobs) = c5name(ip) 

            ! Save the horizontal location:
            observations%xyzt(1,observations%numobs) = x
            observations%xyzt(2,observations%numobs) = y

            ! Save the vertical location:
            sigma = height2sigma(obs_h(ip,il),INT(x),INT(y))
            observations%xyzt(3,observations%numobs) = sigma/(domain(2,3)-domain(1,3))*(numgrid(3)-1)+1

            ! Save the temporal location:
            observations%xyzt(4,observations%numobs) = &
              (obs_time(ip)-domain(1,4))/(domain(2,4)-domain(1,4))*(numgrid(4)-1)+1

          ENDIF

          ! End of check valid obs
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ENDDO          !%%%%%%%%%%%% through levels

      ENDIF          !%%%%%%%%%%%% if there are at least one level data

    ENDDO         !%%%%%%%%%%%% through all profilers
    
  ENDDO        !%%%%%%%%%%%% through all time frames

  ! Deallocate memory for profiler data:
  DEALLOCATE(obs_type,c5name, &
             obs_time, &
             pr_lat,pr_lon,pr_elv,&
             obs_h,obs_u, &
             obs_v,obs_rms, &
             obs_tsfc,obs_psfc,obs_rhsfc, &
             obs_usfc,obs_vsfc,STAT=istatus)

  PRINT*
  PRINT*,'+----------------------------------+'
  PRINT*,'|  End of reading in profile data  |'
  PRINT*,'+----------------------------------+'

END SUBROUTINE insitu_height_profiler
