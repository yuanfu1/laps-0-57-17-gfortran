!doc================================================================== 
!doc This module is designed to use LAPS to ingest model backgrounds
!doc and observations for STMAS surface analysis.
!
!doc History: Dec. 2009 by Yuanfu Xie modifying the LAPSDATBKG.f90 for
!          stmas_mg under src/mesowave.
!doc
!doc Modified where it is appropriate to use in 4d. May 2011 by Hongli Jiang 
!doc==================================================================

MODULE LAPS_ingest

USE MEM_NAMELIST          ! USE LAPS NAMELIST, is in ../lib/modules/mem_namelist.f90. HJ 6/7/2011

INCLUDE 'lapsparms.cmn'  ! LAPS parameters: grid_spacing_m

PUBLIC get_LAPS_dimension,get_LAPS_config,get_LAPS_bkgd,get_LAPS_obs

! PRIVATE 
! Constants
REAL,PARAMETER :: knot2ms = 0.51444444
! A few parameters are re-defined in here to be free of STMAS_module.f90). HJ 6/29/2011
! re-define temp_00, was defined as temp_0 in STMAS_module.f90. HJ 6/13/2011
! re-define ztop as ztopp. HJ 6/29/2011
REAL,PARAMETER :: temp_00 = 273.16
REAL,PARAMETER :: ztopp=20000.

! Variables
! Note: i4time points to the second last time frame, i.e. numgrid(3)-1
! in order to use half cycle time data in the future, where background
! may not exist at numgrid(3):
INTEGER, ALLOCATABLE :: i4time_sequence(:)
INTEGER :: LAPS_cycle,max_sfc_stations
INTEGER :: levels ! Vertical levels for using snd surface data
!INTEGER :: max_profiles,max_profile_levels
REAL    :: missing,badsfc
! End of private variable declarations

!--------- Routines:

CONTAINS

SUBROUTINE get_LAPS_dimension(i4time,time_window_in,numgrid_out,domain_out, &
                                grid_spacing_out)

!doc------------------------------------------------------------------
!doc This routine reads in LAPS domain grid informaiton in x, y and t.
!doc
!doc History: March 2010 by Yuanfu Xie
!doc
!doc Parameters:
!doc
!doc   Input:
!doc     time_window_in:      (char*13) Analysis start and end time
!doc                          in form of YYYYMMDD_HHMM
!doc     i4time:              Current LAPS i4time
!doc
!doc   Output:
!doc     numgrid_out:         (int*4) Numbers of grid points in x,y,Z,t !HJ added Z, 6/7/2011
!doc     domain_out:          (real*2x4) Analysis domain in x,y,Z,t ! HJ added Z, 6/7/2011
!doc------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER*13, INTENT(IN) :: time_window_in(2)
  INTEGER, INTENT(IN) :: i4time

  INTEGER, INTENT(OUT) :: numgrid_out(4) !HJ changed from 3 to 4. 6/7/2011
  REAL,    INTENT(OUT) :: domain_out(2,4),grid_spacing_out(4) !HJ changed from 3 to 4. 6/7/2011

  ! Local variables:
  CHARACTER*9 :: time9(2),WFO_FNAME13_TO_FNAME9 ! LAPS function
  CHARACTER*42, PARAMETER :: header = 'STMAS>LAPS_ingest>get_LAPS_dimension: '
  CHARACTER*40 :: vert_grid
  INTEGER     :: begtime,endtime,istatus

  time9(1) = WFO_FNAME13_TO_FNAME9(time_window_in(1))
  time9(2) = WFO_FNAME13_TO_FNAME9(time_window_in(2))

  CALL CV_ASC_I4Time(time9(1),begtime,istatus)
  CALL CV_ASC_I4Time(time9(2),endtime,istatus)

  CALL GET_SYSTIME_I4(i4time,istatus)

  CALL GET_LAPS_CYCLE_TIME(LAPS_cycle,istatus)

  ! Require begin time coincident with a LAPS cycle:
  IF (MOD(i4time-begtime,LAPS_cycle) .NE. 0) THEN
    PRINT*
    PRINT*,header,'Start time is not at a LAPS cycle time,',&
                  ' adjust it and rerun please!'
    STOP
  ENDIF

  IF (endtime .GT. i4time+LAPS_cycle) THEN
    PRINT*
    PRINT*,header,'Warning!!! end time passes current ',&
                  'LAPS system time, data may not be available!'
  ENDIF

  ! Time frames: HJ change from 3 to 4. 6/7/2011
  numgrid_out(4) = (endtime-begtime)/LAPS_cycle+1

  ! I4 time sequence of the temporal domain:
  ALLOCATE(i4time_sequence(numgrid_out(4)))
  DO istatus=1,numgrid_out(4)
    i4time_sequence(istatus) = begtime+(istatus-1)*LAPS_cycle
  ENDDO

  ! Get numbers of horizontal gridpoints
  CALL GET_GRID_DIM_XY(numgrid_out(1),numgrid_out(2),istatus)
  CALL GET_LAPS_DIMENSIONS(numgrid_out(3),istatus)


  grid_spacing_out(1:2) = grid_spacing_m
!
! calculate vertical grid_spacing_out(3) here. HJ 6/29/2011
  grid_spacing_out(3) = ztopp/real(numgrid_out(3)-1.)
 
! HJ changed from 3 to 4, 4 for time and 3 for Z. 6/7/2011 
  grid_spacing_out(4) = FLOAT(LAPS_cycle)

  ! Grid spacing: horizontal center grid_spacing_m is defined in lapsparms.cmn
  PRINT*,'Horizontal gridspacing:',grid_spacing_out(1:2)
  PRINT*,'Time: ',FLOAT(LAPS_cycle)

! I guess zz should be defined at this point? or should zz be declared in STMAS_module.f90? HJ 6/10/2011
! Domain: HJ change from 3 to 4 6/7/2011
  domain_out(1,1:4) = 0.0
  domain_out(2,1:4) = FLOAT((numgrid_out(1:4)-1))*grid_spacing_out(1:4)

  ! Mis:
  ! Get number of vertical levels for using snd surface data:
  CALL GET_LAPS_DIMENSIONS(levels,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'Problem reading vertical levels!'
    STOP
  ENDIF

END SUBROUTINE get_LAPS_dimension

SUBROUTINE get_LAPS_config(ngrid_in,lat_out,lon_out, &
                           topo_out,land_out,mapfactor_out,rdplvl_out)

!doc==================================================================
!doc This routine reads in LAPS configuration information
!doc
! History: Dec. 2009 by Yuanfu Xie adapted from LAPSInfo under 
!                    src/mesowave/stmas_mg.
!doc
!doc    Input:
!doc	  ngrid_in:     (int*4) Number of gridpoints in x,y,Z,t, HJ changed from 3 to 4 6/7/2011
!doc
!doc	Output
!doc	  lat/lon_out:	(float ngrid_in(1)xngrid_in(2)
!doc                    Latitude/longitude of the domain.
!doc	  topo/land_out:(float ngrid_in(1)xngrid_in(2)
!doc                    Topography and landfactor over the domain.
!doc	  mapfactor_out:(float ngrid_in(1)xngrid_in(2)
!doc                    Projection mapping factor.
!doc==================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrid_in(4)

  REAL,INTENT(OUT) :: lat_out(ngrid_in(1),ngrid_in(2)), &
                      lon_out(ngrid_in(1),ngrid_in(2)), &
                      topo_out(ngrid_in(1),ngrid_in(2)), &
                      land_out(ngrid_in(1),ngrid_in(2)), &
                      mapfactor_out(ngrid_in(1),ngrid_in(2)), &
                      rdplvl_out

  ! Local variables:
  CHARACTER*39, PARAMETER :: header = &  ! Measure the string when it changes
    'STMAS>LAPS_ingest>get_LAPS_config: '
  INTEGER :: i,j
  INTEGER, PARAMETER :: len_header = 39             ! Match the header length

  CHARACTER :: c9time*9,filename*150,dir*150
  INTEGER :: istatus,length

  ! Invalide values from LAPS:
  CALL GET_R_MISSING_DATA(missing,istatus)
  CALL GET_SFC_BADFLAG(badsfc,istatus)
  PRINT*
  PRINT*,header,'Surface Missing/Bad flags: ',missing,badsfc

  ! Get lat/lon,topography and land-factor:
  CALL READ_STATIC_GRID(ngrid_in(1),ngrid_in(2),'LAT',lat_out,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'error in getting LAPS LAT'
    STOP
  ENDIF
  CALL READ_STATIC_GRID(ngrid_in(1),ngrid_in(2),'LON',lon_out,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'error in getting LAPS LON'
    STOP
  ENDIF
  CALL READ_STATIC_GRID(ngrid_in(1),ngrid_in(2),'AVG',topo_out,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'error in getting LAPS topo'
    STOP
  ENDIF
  CALL READ_STATIC_GRID(ngrid_in(1),ngrid_in(2),'LDF',land_out,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'error in getting LAPS land-factor'
    STOP
  ENDIF

  ! Use a point-wise routine from LAPS, Steve is developing an array version:
  DO j=1,ngrid_in(2)
    DO i=1,ngrid_in(1)
      CALL GET_SIGMA(lat_out(i,j),lon_out(i,j),mapfactor_out(i,j),istatus)
      IF (istatus .NE. 1) THEN
        PRINT*
        PRINT*,header,'error in getting projection map factor'
        STOP
      ENDIF
    ENDDO
  ENDDO

  ! Get maximum number of surface stations:
  CALL GET_MAXSTNS(max_sfc_stations,istatus)
  IF (istatus .NE. 1) THEN
    PRINT*
    PRINT*,header,'an error occurs accessing max_sfc_stations'
    STOP
  ENDIF
  PRINT*
  PRINT*,header,'Maximum number of surface stations allowed: ',max_sfc_stations

  ! Get LAPS wind parameters:
  CALL GET_DIRECTORY('static',dir,length)
  filename = dir(1:length)//'wind.nl'
  CALL READ_NAMELIST_LAPS('wind',filename)

  ! Get reduced pressure level:
  CALL GET_LAPS_REDP(rdplvl_out,istatus)

END SUBROUTINE get_LAPS_config

SUBROUTINE get_LAPS_bkgd(num_vars_in,var_names_in,ngrid_in,bkgnd_out,zz)
!doc==================================================================
!doc This routine reads in num_var LAPS background fields between 
!doc domain(1,3) and domain(2,3) with var_names into bkgnd array.
!doc
!doc History: Dec. 2009 by Yuanfu Xie
!doc
!doc  Purpose: returns background over the time window of 
!doc         domain(1,3) and domain(2,3). 
!doc  Modified. 6/8/2011 HJ
!doc         domain(1,4) and domain(2,4). 
!doc    
!doc  input
!doc      ngrid_in:         number of background time frames over
!doc                           the time window;
!doc      num_vars_in:      number of variables to be read;
!doc      var_names_in:     names of the variables;
!doc      var_names2d:      sfc pressure. HJ 6/22/2011
!doc    
!doc  output
!doc      bkgnd_out:        backgrounds between domain(1:2,4); changed from 3 to 4. HJ 6/10/2011
!doc      bkgnd2d:          backgrounds surface pressure. HJ 6/22/2011
!doc                        bkgnd2d is needed only in pressure coordinate. 
!doc                        For sigma_ht, bkgnd2d is no longer used. The file is saved in ./src_orig. 
!doc                        HJ 8/29/2011
!doc
!doc Modified to read in 3D fields. 6/8/2011 Hongli Jiang
!doc==================================================================

  IMPLICIT NONE

  INTEGER,      INTENT(IN) :: num_vars_in,ngrid_in(4) ! replace 3 with 4. HJ 6/8/2011 
  CHARACTER*4,  INTENT(IN) :: var_names_in(num_vars_in)

  REAL,        INTENT(OUT) :: bkgnd_out(ngrid_in(1),ngrid_in(2), &
                                        ngrid_in(3),ngrid_in(4),num_vars_in),zz(ngrid_in(3))

  ! Local variables:
  INTEGER :: icount,i,j,k,t,iv,istatus
  real :: omega_to_w

! Begin and end are seconds from the current time, 0:
! k was used to denote time, in ngrid_in(3), now it is ngrid_in(4). 
! the time loop index is replaced with t. HJ 6/8/2011
!
  DO t=1,ngrid_in(4)

! Read in all 3D variables:
    DO iv=1,num_vars_in                ! * Through variables

! LAPS GET_MODELFG_3D to read lga 3D data: HJ 6/8/2011
! NOTE: LAPS background pressure is in Pascal:
      CALL GET_MODELFG_3D(i4time_sequence(t),var_names_in(iv), &
                          ngrid_in(1),ngrid_in(2),ngrid_in(3), &
                          bkgnd_out(1,1,1,t,iv),istatus)

! HJ added: when there is no lga background, stop the run. Otherwise NaN will be produced 
! during the calulation of theta_v when using zero values of P, or T (C). 8/29/2011
      IF (istatus .NE. 1) THEN
        WRITE(*,1) var_names_in(iv),iv
1       FORMAT('STMAS>LAPS_ingest>get_LAPS_bkgd: No LGA background',&
               ' for: ',a4,i3,' stop the run')

        STOP
      ENDIF
      PRINT*,'LAPS background for variable ',var_names_in(iv),' is read in successfully!'

! print sample of bkgnd data. HJ 6/9/2011
      WRITE(*,'(A4,e13.5)') var_names_in(iv),bkgnd_out(2,2,15,t,iv)
    ENDDO                              ! * End through variables
  ENDDO
!
! get zz.
  if(var_names_in(4) .eq. 'P3')then
    call get_ht_1d(ngrid_in(3),zz,istatus)
    PRINT*,'XIE ZZ: ',MINVAL(zz), MAXVAL(zz)
    IF (istatus .NE. 1) THEN
      print*,'STMAS>LAPS_ingest>get_LAPS_bkgd: simga_ht is incorrect'
    endif
  endif
!
! OM=dp/dt, convert om to W by calling omega_to_w. convert_om2w is defined 
! in grid_interpolation.f90 HJ 8/5/2011
! when calling from /data/fab/stmas/stmas4d, var names and order: u3,v3,om,p3,t3,sh HJ 8/5/2011
  if(var_names_in(3) .eq. 'OM')then
    call convert_om2w(ngrid_in(1),ngrid_in(2),ngrid_in(3),ngrid_in(4),bkgnd_out(1,1,1,1,3),bkgnd_out(1,1,1,1,4),istatus)
    IF (istatus .NE. 1) THEN
      print*,'STMAS>LAPS_ingest>om2w is incorrect'
    endif
  endif
  
END SUBROUTINE get_LAPS_bkgd

SUBROUTINE get_LAPS_obs(domain_in,num_vars_in,var_names_in, &
                          invalid_in,ngrid_in,grid_lat_in, &
                          grid_lon_in,topo_in,max_obs_in, &
                          num_obs_out,obsvs_out,obserr_out, &
                          obsxytz_out,obslat_out,obslon_out, &
                          obsstn_out,obstyp_out,verbal,success)
!doc==================================================================
!doc This routine is to read in observation data through LAPS.
!doc
!  History: Jan. 2010 by Yuanfu Xie
!doc   
!doc   
!doc==================================================================
  
  IMPLICIT NONE

  INTEGER,     INTENT(IN) :: num_vars_in,ngrid_in(3),max_obs_in,verbal
  CHARACTER*4, INTENT(IN) :: var_names_in(num_vars_in)
  REAL,        INTENT(IN) :: grid_lat_in(ngrid_in(1),ngrid_in(2)), &
                             grid_lon_in(ngrid_in(1),ngrid_in(2)), &
                             topo_in(ngrid_in(1),ngrid_in(2))
  REAL,        INTENT(IN) :: invalid_in,domain_in(2,3)

  CHARACTER,INTENT(OUT) :: obsstn_out(max_obs_in)*10, &
                           obstyp_out(max_obs_in)*6
  INTEGER, INTENT(OUT) :: num_obs_out,success
  REAL,INTENT(OUT) :: obsvs_out(max_obs_in,num_vars_in), &
                      obserr_out(max_obs_in,num_vars_in), &
                      obsxytz_out(4,max_obs_in), &
                      obslat_out(max_obs_in),obslon_out(max_obs_in)

  ! Local variables:
  ! Lapse rates: see LAPS mdatlaps.f under sfc: lapses(1): temp; lapses(2): dewp
  REAL, PARAMETER :: lapses(2) = (/-0.01167, -0.007/)

  CHARACTER*36, PARAMETER :: header = &
    'STMAS>LAPS_ingest>get_LAPS_obs: '
  CHARACTER*40 :: atime    ! Use to print time for missing data
  INTEGER :: k,max_sfc_obs,n,istatus,iv,i4,out_xy,out_t,iuse_stn
  INTEGER :: totalobs(num_vars_in) ! Count total valid obs for each variable
  REAL    :: ALT_2_SFC_PRESS,x(2)  ! x holds the grid location of obs
  REAL    :: uerr(4),verr(4)       ! Save the u/v obs with obs errors
  REAL    :: tc,f_to_c,dwpt        ! Fahrenheit and dewpoint conversion

  ! LAPS LSO required arrays: ****
  CHARACTER,ALLOCATABLE :: station_names(:)*20, &
                           providers(:)*11,     & 
                           weathers(:)*25,      & ! Present weather
                           report_types(:)*6,   &
                           station_types(:)*6,  &
                           cloud_amount(:,:)*4

  CHARACTER :: time_obsfile*24

  INTEGER :: nobs_on_grid,nobs_on_box              ! At each frame
  INTEGER, ALLOCATABLE :: wmo_id(:),obs_time(:),num_cloud_layers(:)

  REAL,    ALLOCATABLE :: lat(:),lon(:),elevation(:), &
                          temp_2m(:),temp_err(:), &
                          dew(:),dew_err(:), &
                          rh(:),rh_err(:), &
                          wind_dir(:),dir_err(:), &
                          wind_spd(:),spd_err(:), &
                          gust_dir(:),gust_spd(:), &
                          altimeter(:),alt_err(:), &
                          stn_prs(:),prs_err(:), &  ! Station pressure
                          msl_prs(:), & 
                          press3_change_character(:), & ! 3 hour pressure 
                          press3_cc_err(:), &           ! change character
                          press3_change(:), &  ! 3 hour pressure change
                          press3_c_err(:), &   ! 3 hour pressure change error
                          visibility(:), visibility_err(:), &
                          solar(:),solar_err(:), &
                          soil_temp(:),soil_t_err(:), &
                          soil_moist(:),soil_m_err(:), &
                          precip1(:),precip_err(:), & ! 1 hour precipitation
                          precip3(:),precip6(:), & ! 3,6 hour precipitation
                          precip24(:), & ! 24 hour precipitation
                          snow_depth(:),snow_err(:), &
                          temp_max(:),temp_min(:), &
                          cloud_layer_height(:,:)                         

  ! End LAPS LSO array declarations ****
                      

  ! Total time frames of data to read:
  max_sfc_obs = max_sfc_stations*ngrid_in(3)

  ! Allocate memory for reading LSO files:
  ALLOCATE(station_names(max_sfc_obs), &
           providers(max_sfc_obs), &
           weathers(max_sfc_obs), &
           report_types(max_sfc_obs), &
           station_types(max_sfc_obs), &
           cloud_amount(max_sfc_obs,5), STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'error for allocating memory for chars'
    success = 0
    STOP
  ENDIF

  ALLOCATE(wmo_id(max_sfc_obs), &
           obs_time(max_sfc_obs), &
           num_cloud_layers(max_sfc_obs), STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'error for allocating memory for ints'
    success = 0
    STOP
  ENDIF

  ALLOCATE(lat(max_sfc_obs), &
           lon(max_sfc_obs), &
           elevation(max_sfc_obs), &
           temp_2m(max_sfc_obs), &
           temp_err(max_sfc_obs), &
           dew(max_sfc_obs), &
           dew_err(max_sfc_obs), &
           rh(max_sfc_obs), &
           rh_err(max_sfc_obs), &
           wind_dir(max_sfc_obs), &
           wind_spd(max_sfc_obs), &
           dir_err(max_sfc_obs), &
           spd_err(max_sfc_obs), &
           gust_dir(max_sfc_obs), &
           gust_spd(max_sfc_obs), &
           altimeter(max_sfc_obs), &
           alt_err(max_sfc_obs), &
           stn_prs(max_sfc_obs), &
           prs_err(max_sfc_obs), &
           msl_prs(max_sfc_obs), &
           press3_change_character(max_sfc_obs), &
           press3_cc_err(max_sfc_obs), &
           press3_change(max_sfc_obs), &
           press3_c_err(max_sfc_obs), &
           visibility(max_sfc_obs), &
           visibility_err(max_sfc_obs), &
           solar(max_sfc_obs), &
           solar_err(max_sfc_obs), &
           soil_temp(max_sfc_obs), &
           soil_t_err(max_sfc_obs), &
           soil_moist(max_sfc_obs), &
           soil_m_err(max_sfc_obs), &
           precip1(max_sfc_obs), &
           precip_err(max_sfc_obs), &
           precip3(max_sfc_obs), &
           precip6(max_sfc_obs), &
           precip24(max_sfc_obs), &
           snow_depth(max_sfc_obs), &
           snow_err(max_sfc_obs), &
           temp_max(max_sfc_obs), &
           temp_min(max_sfc_obs), &
           cloud_layer_height(max_sfc_obs,5), STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'error for allocating memory for reals'
    success = 0
    STOP
  ENDIF

  ! Read in LSO data: **********
  n = 1
  DO k=1,ngrid_in(3)
    time_obsfile = ' '

    ! Read in surface observation data frame by frame:
    ! Note: nobs_on_box is initialized by reading lso file:
    CALL READ_SURFACE_DATA(i4time_sequence(k), &   ! Time frame
         time_obsfile,nobs_on_grid,nobs_on_box, &
         obs_time(n),wmo_id(n),station_names(n), &
         providers(n),weathers(n),report_types(n), &
         station_types(n),lat(n),lon(n),elevation(n), &
         temp_2m(n),dew(n),rh(n),wind_dir(n),wind_spd(n), &
         gust_dir(n),gust_spd(n),altimeter(n),stn_prs(n), &
         msl_prs(n),press3_change_character(n), &
         press3_change(n),visibility(n),solar(n), &
         soil_temp(n),soil_moist(n),precip1(n),precip3(n), &
         precip6(n),precip24(n),snow_depth(n), &
         num_cloud_layers(n),temp_max(n),temp_min(n), &
         temp_err(n),dew_err(n),rh_err(n),dir_err(n), &
         spd_err(n),alt_err(n),prs_err(n),visibility_err(n), &
         solar_err(n),soil_t_err(n),soil_m_err(n), &
         precip_err(n),snow_err(n),cloud_amount(n,1), &
         cloud_layer_height(n,1),max_sfc_stations,istatus)

    IF (istatus .NE. 1) THEN
      CALL CV_I4TIM_ASC_LP(i4time_sequence(k),atime,istatus)
      PRINT*,header,'No sfc data at time: ',atime
    ELSE
      PRINT*,header,'Num_SFC_obs on grid/box: ', &
        nobs_on_grid,nobs_on_box,' at: ',time_obsfile

      ! Convert observation time to i4time:
      DO iv=0,nobs_on_box-1
        ! LAPS obs_time is in format of HHMM:
        IF (obs_time(n+iv) .GE. 0 .AND. obs_time(n+iv) .LE. 2400) THEN
          CALL GET_SFC_OBTIME(obs_time(n+iv),i4time_sequence(k),i4,istatus)
          obs_time(n+iv) = i4-i4time_sequence(1)
        ELSE
          obs_time(n+iv) = -86400*365 ! Unused by setting it past year
        ENDIF
      ENDDO
      n = n+nobs_on_box     ! Count sfc obs

      IF (n .GT. max_sfc_obs) THEN
        PRINT*,header,'too many sfc obs: ',n,max_sfc_obs
        success = 0
        STOP
      ENDIF
    ENDIF


    ! End of read of LSO data  **********


    ! Read in SND SFC data: **********

    ! Read in SND surface observation data frame by frame:
    ! Must initial nobs_on_box before calling:
    nobs_on_grid = 0; nobs_on_box = 0
    time_obsfile = ''
    CALL READ_SFC_SND(i4time_sequence(k), &   ! Time frame
         time_obsfile,nobs_on_grid,nobs_on_box, &
         obs_time(n),wmo_id(n),station_names(n), &
         providers(n),weathers(n),report_types(n), &
         station_types(n),lat(n),lon(n),elevation(n), &
         temp_2m(n),dew(n),rh(n),wind_dir(n),wind_spd(n), &
         gust_dir(n),gust_spd(n),altimeter(n),stn_prs(n), &
         msl_prs(n),press3_change_character(n), &
         press3_change(n),visibility(n),solar(n), &
         soil_temp(n),soil_moist(n),precip1(n),precip3(n), &
         precip6(n),precip24(n),snow_depth(n), &
         num_cloud_layers(n),temp_max(n),temp_min(n), &
         temp_err(n),dew_err(n),rh_err(n),dir_err(n), &
         spd_err(n),alt_err(n),prs_err(n),visibility_err(n), &
         solar_err(n),soil_t_err(n),soil_m_err(n), &
         precip_err(n),snow_err(n),cloud_amount(n,1), &
         cloud_layer_height(n,1),max_sfc_stations, &
         grid_lat_in,grid_lon_in,ngrid_in(1),ngrid_in(2), & ! SND SFC additions
         levels,max_pr,max_pr_levels,topo_in,istatus)
    IF (istatus .NE. 1) THEN
      PRINT*,header,'No sonde sfc data at: ',k
    ELSE
      PRINT*,header,'Num_SND_obs on grid/box: ', &
        nobs_on_grid,nobs_on_box,' at: ',time_obsfile

      ! Convert observation time to i4time:
      DO iv=1,nobs_on_box
        IF (obs_time(n+iv) .GE. 0 .AND. obs_time(n+iv) .LE. 2400) THEN
          CALL GET_SFC_OBTIME(obs_time(n+iv),i4time_sequence(k),i4,istatus)
          obs_time(n+iv) = i4-i4time_sequence(1)
        ELSE
          obs_time(n+iv) = -86400*365 ! Unused by setting it past year
        ENDIF
      ENDDO
      n = n+nobs_on_box       ! Count SND sfc data

      IF (n .GT. max_sfc_obs) THEN
        PRINT*,header,'Sfc obs > max_sfc_obs: ',n,max_sfc_obs
        success = 0
        STOP
      ENDIF
    ENDIF

    ! End of read of SND SFC data  **********

  ENDDO

  n = n-1 ! As n indicates the next record
  PRINT*
  PRINT*,'                      ++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT*,header,'| Total sfc obs read:    ',n
  PRINT*,'                      ++++++++++++++++++++++++++++++++++++++++++++++++++++'

  ! Check if the number of obs is greater than the one STMAS namelist allows
  IF (n .GT. max_obs_in) THEN
    PRINT*
    PRINT*,header,'Number of obs is greater than STMAS namelist allows,',&
                  ' Change the parameter in namelist and rerun!'
    success = 0
    STOP
  ENDIF

  ! Pass the data to the output arrays:
  !++++++++++++++++++++++++++++++++++++
  ! Pass locations:
  num_obs_out = 0
  out_xy = 0
  out_t = 0
  totalobs = 0
  DO k=1,n

    ! Initialize iuse_stn: 0 no valid obs at this station
    iuse_stn = 0

    ! X/Y: Note: latlon_to_rlapsgrid returns value 1-ngrid_in:
    IF (lat(k) .NE. missing .AND. lat(k) .NE. badsfc .AND. &
        lon(k) .NE. missing .AND. lon(k) .NE. badsfc .AND. &
        elevation(k) .NE. missing .AND. elevation(k) .NE. badsfc .AND. &
        obs_time(k) .GE. domain_in(1,3) .AND. &
        obs_time(k) .LE. domain_in(2,3)) THEN

      ! Valid location:
      CALL LATLON_TO_RLAPSGRID(lat(k),lon(k),grid_lat_in,grid_lon_in, &
                               ngrid_in(1),ngrid_in(2), &
                               x(1),x(2),istatus)
      ! Checking horizontal domain:
      IF (x(1) .GE. 1 .AND. x(1) .LE. ngrid_in(1) .AND. &
          x(2) .GE. 1 .AND. x(2) .LE. ngrid_in(2) ) THEN

        ! Found obs inside domain:
        num_obs_out = num_obs_out+1
        obsxytz_out(1:2,num_obs_out) = grid_spacing_m*(x(1:2)-1.0)
        obsxytz_out(3,num_obs_out) = obs_time(k)
        obsxytz_out(4,num_obs_out) = elevation(k)

        ! Pass station names and report types:
        obsstn_out(num_obs_out)(1:10) = station_names(k)(1:10)
        obslat_out(num_obs_out) = lat(k)
        obslon_out(num_obs_out) = lon(k)
        obstyp_out(num_obs_out) = report_types(k)

        ! Pass observation values:
        DO iv=1,num_vars_in  ! ** through all analysis variables
          ! Check through all variables
          SELECT CASE(var_names_in(iv))
            CASE('USF')
              IF (var_names_in(iv+1) .NE. 'VSF') THEN
                PRINT*,header,'Please use VSF following USF as variable name!'
                STOP
              ENDIF

              ! Save wind fields:
              IF (wind_dir(k) .NE. missing .AND. &
                  wind_dir(k) .NE. badsfc .AND. &
                  wind_spd(k) .NE. missing .AND. &
                  wind_spd(k) .NE. badsfc .AND. &
                  dir_err(k) .NE. missing .AND. &
                  dir_err(k) .NE. badsfc .AND. &
                  spd_err(k) .NE. missing .AND. &
                  spd_err(k) .NE. badsfc) THEN

                ! Find valid wind obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                totalobs(iv+1) = totalobs(iv+1)+1
                CALL disptrue_to_uvgrid(wind_dir(k),wind_spd(k), &
                                        obsvs_out(num_obs_out,iv), &
                                        obsvs_out(num_obs_out,iv+1),lon(k))

                ! Calculate the obs errors for u and v:
                !**************************************
                IF (dir_err(k) .LE. 1.0e-10 .AND. spd_err(k) .LE. 1.0e-10) THEN
                  PRINT*,header,'Invalid obs errors from LAPS:',dir_err(k),spd_err(k),k
                  STOP
                ENDIF

                ! First add or subtract the error from the obs and compute the u/v with max error:
                CALL disptrue_to_uvgrid(wind_dir(k)+dir_err(k),wind_spd(k)+spd_err(k), &
                                        uerr(1),verr(1),lon(k))
                CALL disptrue_to_uvgrid(wind_dir(k)-dir_err(k),wind_spd(k)+spd_err(k), &
                                        uerr(2),verr(2),lon(k))
                CALL disptrue_to_uvgrid(wind_dir(k)+dir_err(k),wind_spd(k)-spd_err(k), &
                                        uerr(3),verr(3),lon(k))
                CALL disptrue_to_uvgrid(wind_dir(k)-dir_err(k),wind_spd(k)-spd_err(k), &
                                        uerr(4),verr(4),lon(k))
                ! Second choose the maximum as the obs error for u/v:
                obserr_out(num_obs_out,iv  ) = MAXVAL(ABS(uerr(1:4)-obsvs_out(num_obs_out,iv  )))
                obserr_out(num_obs_out,iv+1) = MAXVAL(ABS(verr(1:4)-obsvs_out(num_obs_out,iv+1)))
                !**************************************

                ! Convert obs in knots to m/s
                obsvs_out(num_obs_out,iv:iv+1) = &
                  obsvs_out(num_obs_out,iv:iv+1)*knot2ms
                ! Convert obs errors in knots to absolute m/s:
                obserr_out(num_obs_out,iv:iv+1) = &
                  obserr_out(num_obs_out,iv:iv+1)*knot2ms
              ENDIF
            CASE('VSF')
              IF (var_names_in(iv-1) .NE. 'USF') THEN
                PRINT*,'Please use USF before VSF as variable name!'
                STOP
              ENDIF
            CASE('DSF')
              ! LAPS dewpint is in Fahrenheit:
              obsvs_out(num_obs_out,iv) = invalid_in  ! Default
              IF (dew(k) .NE. missing .AND. dew(k) .NE. badsfc) THEN
                ! Valid temperature obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                obsvs_out(num_obs_out,iv) = (dew(k)-32.0)*5.0/9.0+temp_00
              ELSEIF (rh(k) .NE. missing .AND. rh(k) .NE. badsfc .AND. &
                      temp_2m(k) .NE. missing .AND. temp_2m(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                tc = (temp_2m(k)-32.0)*5.0/9.0
                obsvs_out(num_obs_out,iv) = dwpt(tc,rh(k))+temp_00
              ENDIF

              obserr_out(num_obs_out,iv) = invalid_in  ! Default
              IF (dew_err(k) .NE. missing .AND. dew_err(k) .NE. badsfc) THEN
                obserr_out(num_obs_out,iv) = dew_err(k)*5.0/9.0
              ELSEIF (rh_err(k) .NE. missing .AND. rh_err(k) .NE. badsfc .AND. &
                      temp_err(k) .NE. missing .AND. temp_err(k) .NE. badsfc) THEN
                tc = (temp_err(k)-32.0)*5.0/9.0
                obserr_out(num_obs_out,iv) = dwpt(tc,rh(k)) ! Error in C is the same as in K
              ENDIF
            CASE('TSF')
              ! LAPS temp is in Fahrenheit:
              obsvs_out(num_obs_out,iv) = invalid_in  ! Default
              IF (temp_2m(k) .NE. missing .AND. temp_2m(k) .NE. badsfc) THEN
                ! Valid temperature obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                obsvs_out(num_obs_out,iv) = (temp_2m(k)-32.0)*5.0/9.0+temp_00
              ENDIF

              obserr_out(num_obs_out,iv) = invalid_in  ! Default
              IF (temp_err(k) .NE. missing .AND. temp_err(k) .NE. badsfc) &
                obserr_out(num_obs_out,iv) = temp_err(k)*5.0/9.0
            CASE('PSF')
              IF (stn_prs(k) .NE. missing .AND. stn_prs(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                obsvs_out(num_obs_out,iv) = stn_prs(k)  ! LAPS in mb
              ELSEIF (altimeter(k) .NE. missing .AND. &
                      altimeter(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                obsvs_out(num_obs_out,iv) = &
                  ALT_2_SFC_PRESS(altimeter(k),elevation(k))
              ELSEIF (msl_prs(k) .NE. missing .AND. &
                      msl_prs(k) .NE. badsfc  .AND. &
                      temp_2m(k) .NE. missing .AND. &
                      temp_2m(k) .NE. badsfc  .AND. &
                      dew(k) .NE. missing .AND. dew(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                CALL REDUCE_P(1.8*(temp_2m(k)-temp_00)+32.0, &
                              1.8*(dew(k)-temp_00)+32.0, &
                              msl_prs(k),0.0,lapses(1),lapses(2), &
                              obsvs_out(num_obs_out,iv),elevation(k),invalid_in)
              ELSE
                obsvs_out(num_obs_out,iv) = invalid_in
              ENDIF            
          
              IF (prs_err(k) .NE. missing .AND. prs_err(k) .NE. badsfc) THEN
                obserr_out(num_obs_out,iv) = prs_err(k)  ! LAPS in mb
              ELSEIF (alt_err(k) .NE. missing .AND. &
                      alt_err(k) .NE. badsfc) THEN
                obserr_out(num_obs_out,iv) = &
                  ALT_2_SFC_PRESS(alt_err(k),elevation(k))
              ELSE
                obserr_out(k,iv) = invalid_in
              ENDIF
            CASE DEFAULT
              PRINT*
              PRINT*,header,'No such observations: ',var_names_in(iv)
              success = 0
              STOP
          END SELECT
        ENDDO ! ** through all analysis variables

        ! If no valid obs for the analysis variables, remove this record
        IF (iuse_stn .EQ. 0) THEN
          num_obs_out = num_obs_out-1
          IF (verbal .NE. 0) &
            print*,'No valid obs at this station: ',station_names(k)(1:10)
        ENDIF

      ELSE
        out_xy = out_xy+1
        IF (verbal .GE. 1) THEN
          PRINT*
          PRINT*,header,'station: ',station_names(k),' out of x-y domain', &
            x(1),x(2)
          PRINT*,' (',domain_in(1:2,1:2)/grid_spacing_m,')'
        ENDIF
      ENDIF   ! End checking horizontal domain
    ELSE
      IF (verbal .GE. 1) THEN
        PRINT*,header,'Warning! station: ',station_names(k)
        IF (lat(k) .EQ. missing .OR. lat(k) .EQ. badsfc .OR. &
            lon(k) .EQ. missing .OR. lon(k) .EQ. badsfc .OR. &
            elevation(k) .EQ. missing .OR. elevation(k) .EQ. badsfc) THEN
          PRINT*,' missing lat/lon/elev: ', &
          lat(k),lon(k),elevation(k),report_types(k)
        ELSE IF (obs_time(k) .LT. domain_in(1,3) .OR. &
                 obs_time(k) .GT. domain_in(2,3)) THEN
          out_t = out_t+1
          PRINT*,' out of time window: ',obs_time(k),domain_in(1:2,3)
        ELSE
          PRINT*,header,'Unknown reason!'
        ENDIF
      ENDIF
    ENDIF     ! End of checking lat/lon/elev/time window

  ENDDO ! *** End of loop through all obs

  ! Print information about the observations:
  PRINT*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT*,header,'| Total valid obs: '
  DO iv=1,num_vars_in
    PRINT*,var_names_in(iv),' total obs: ',totalobs(iv)
  ENDDO
  PRINT*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT*
  PRINT*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT*,header,'| Total number of obs in domain',num_obs_out
  PRINT*,'++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  IF (verbal .GE. 1) &
    PRINT*,header,'| Number out of xy domain: ',out_xy,' of time domain: ',out_t

  !++++++++++++++++++++++++++++++++++++
  
  ! Deallocate the temporary memory:
  DEALLOCATE(station_names,providers,weathers,report_types, &
             station_types,cloud_amount,STAT=istatus)

  DEALLOCATE(wmo_id,obs_time,num_cloud_layers)
  DEALLOCATE(temp_2m,temp_err,dew,dew_err,rh,rh_err, &
             wind_dir,wind_spd,gust_dir,gust_spd, &
             altimeter,alt_err,stn_prs,prs_err,msl_prs, &
             press3_change_character,press3_cc_err, &
             press3_change,press3_c_err,visibility,visibility_err, &
             solar,solar_err,soil_temp,soil_t_err, &
             soil_moist,soil_m_err,precip1,precip_err, &
             precip3,precip6,precip24,snow_depth,snow_err, &
             temp_max,temp_min,cloud_layer_height)

  ! Also deallocating i4time_sequence: Assumption last use is in this routine:
  DEALLOCATE(i4time_sequence)

END SUBROUTINE get_LAPS_obs

END MODULE LAPS_ingest
