MODULE ingest_all

!doc==================================================================
!doc
!doc This module facilitates STMAS 4D analysis data ingest, including:
!doc 1. conventional obs;
!doc
!doc History: June 2011 by Yuanfu Xie (originated)
!doc
!doc==================================================================

  USE LAPS_ingest
  USE STMAS

  ! Local variables:
  CHARACTER :: vertical*10          ! What vertical coordinate used
  INTEGER   :: verbal

  CONTAINS

    ! Cover this while debugging as PGI debugger does not allow a stop in.
    ! INCLUDE 'insitu_height_profiler.f90'
    ! INCLUDE 'insitu_height_surfaces.f90'

    SUBROUTINE ingest(nt,i4times,numgrid,lat,lon,topo,top,domain,num_vars, &
                      varnames,max_obs)
    !doc================================================================
    !doc This routine ingest data at all requested time by observation
    !doc types according to the coordinate selected
    !doc 
    !doc Input:
    !doc       nt:       Number of laps cycle times for ingesting data
    !doc       i4times:  The i4times of these cycle times
    !doc       numgrid:  Numbers of gridpoints in all directions
    !doc       lat/lon:  Two-D lat-lon over the analysis domain
    !doc       topo:     Terrain
    !doc       z:        Vertical heights at each sigma level
    !doc       domain:   Analysis domain (1 - numgrid)
    !doc       num_vars: Number of analyzed variables
    !doc       varnames: Names of these variables
    !doc       max_obs:  Maximum number of observations
    !doc
    !doc History: June 2011 by Yuanfu Xie (originated)
    !doc================================================================

      IMPLICIT NONE

      CHARACTER, INTENT(IN) :: varnames(num_vars)*4
      INTEGER, INTENT(IN) :: nt,i4times(nt),numgrid(4),num_vars,max_obs
      REAL,    INTENT(IN) :: top,domain(2,4)  ! Heights and domain
      REAL,    INTENT(IN) :: lat(numgrid(1),numgrid(2)), &
                             lon(numgrid(1),numgrid(2)), &
                             topo(numgrid(1),numgrid(2))

      CHARACTER,PARAMETER :: header*16 = 'STMAS4D>ingest: '
      INTEGER :: istatus

      verbal = 1 ! Make it in the namelist later!!!!

      ! Initialize:
      observations%numobs = 0

      ! Ingest obs:
      SELECT CASE(vertical)
        CASE('HEIGHT')
          CALL insitu_height_surfaces(domain,num_vars,varnames,numgrid, &
                                      lat,lon,topo,top,max_obs,verbal)

          CALL insitu_height_profiler(nt,i4times,numgrid,lat,lon,topo,top,domain)

!          CALL insitu_height_sonde(nt,i4times,numgrid,lat,lon,topo,top, &
!                                   domain,verbal)
        CASE DEFAULT
          PRINT*
          PRINT*,header,'Vertical coordinate: ',vertical, &
            ' is not supported now; check with STMAS developers'
          STOP
      END SELECT

      ! Interpolate background fields to obs sites:
      CALL back2obs

      CALL quality_control

    END SUBROUTINE ingest

END MODULE ingest_all
