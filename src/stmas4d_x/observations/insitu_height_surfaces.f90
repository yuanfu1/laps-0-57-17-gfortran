SUBROUTINE insitu_height_surfaces(domain_in,num_vars_in,var_names_in, &
                                  ngrid_in,grid_lat_in,grid_lon_in, &
                                  topo,top,max_obs_in,verbal)
!doc==================================================================
!doc This routine is to read in observation data through LAPS.
!doc
!  History: Jan. 2010 by Yuanfu Xie
!doc   
!doc   
!doc==================================================================

  USE LAPS_ingest
  USE STMAS
  
  IMPLICIT NONE

  INTEGER,     INTENT(IN) :: num_vars_in,ngrid_in(4),max_obs_in,verbal
  CHARACTER*4, INTENT(IN) :: var_names_in(num_vars_in)
  REAL,        INTENT(IN) :: grid_lat_in(ngrid_in(1),ngrid_in(2)), &
                             grid_lon_in(ngrid_in(1),ngrid_in(2)), &
                             topo(ngrid_in(1),ngrid_in(2)),top
  REAL,        INTENT(IN) :: domain_in(2,4)

  ! Height to sigma conversion function:
  INCLUDE 'height2sigma.f90'

  ! Local variables:
  ! Lapse rates: see LAPS mdatlaps.f under sfc: lapses(1): temp; lapses(2): dewp
  REAL, PARAMETER :: t_liquid = 0.0, t_ice = -132.0  ! in C

  CHARACTER*32, PARAMETER :: header = &
    'STMAS4D>insitu_height_surfaces: '
  CHARACTER*40 :: atime    ! Use to print time for missing data
  INTEGER :: k,max_sfc_obs,n,istatus,iv,i4,out_xy,out_t,iuse_stn,nbefore
  INTEGER :: totalobs(num_vars_in) ! Count total valid obs for each variable
  REAL    :: ALT_2_SFC_PRESS,x(2)  ! x holds the grid location of obs
  REAL    :: uerr(4),verr(4)       ! Save the u/v obs with obs errors
  REAL    :: tc,dc,f_to_c,dwpt     ! Fahrenheit and dewpoint conversion
  REAL    :: sigma                 ! Sigma value

  ! Real functions for specific humidity
  REAL    :: ssh2,make_ssh

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
  max_sfc_obs = max_sfc_stations*ngrid_in(4)

  ! Allocate memory for reading LSO files:
  ALLOCATE(station_names(max_sfc_obs), &
           providers(max_sfc_obs), &
           weathers(max_sfc_obs), &
           report_types(max_sfc_obs), &
           station_types(max_sfc_obs), &
           cloud_amount(max_sfc_obs,5), STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'error for allocating memory for chars'
    STOP
  ENDIF

  ALLOCATE(wmo_id(max_sfc_obs), &
           obs_time(max_sfc_obs), &
           num_cloud_layers(max_sfc_obs), STAT=istatus)
  IF (istatus .NE. 0) THEN
    PRINT*,header,'error for allocating memory for ints'
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
    STOP
  ENDIF

  ! Read in LSO data: **********
  n = 1
  DO k=1,ngrid_in(4)
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
        STOP
      ENDIF
    ENDIF

    ! End of read of LSO data  **********


    ! Read in SND SFC data: **********
goto 111 ! Temporarily skip sonde surface obs
    PRINT*,'Start reading sonde surface data ...'

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
         levels,max_pr,max_pr_levels,topo,istatus)

    IF (istatus .NE. 1) THEN
      PRINT*,header,'No sonde sfc data at time frame: ',k
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
        STOP
      ENDIF
    ENDIF

    ! End of read of SND SFC data  **********
111 continue

  ENDDO

  n = n-1 ! As n indicates the next record
  PRINT*
  PRINT*,'       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  PRINT*,header,'| Total sfc obs read:    ',n
  PRINT*,'       +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'

  ! Check if the number of obs is greater than the one STMAS namelist allows
  IF (n .GT. max_obs_in) THEN
    PRINT*
    PRINT*,header,'Number of obs is greater than STMAS namelist allows,',&
                  ' Change the parameter in namelist and rerun!'
    print*,'Max allowed: ',max_obs_in,' actual ingested by far: ',n
    STOP
  ENDIF

  ! Pass the data to the output arrays:
  !++++++++++++++++++++++++++++++++++++
  ! Pass locations:
  nbefore = observations%numobs 
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
        obs_time(k) .GE. domain_in(1,4) .AND. &
        obs_time(k) .LE. domain_in(2,4)) THEN

      ! Valid location:
      CALL LATLON_TO_RLAPSGRID(lat(k),lon(k),grid_lat_in,grid_lon_in, &
                               ngrid_in(1),ngrid_in(2), &
                               x(1),x(2),istatus)
      ! Checking horizontal domain:
      IF (x(1) .GE. 1 .AND. x(1) .LE. ngrid_in(1) .AND. &
          x(2) .GE. 1 .AND. x(2) .LE. ngrid_in(2) ) THEN

        ! Found obs inside domain:
        sigma = height2sigma(elevation(k),INT(x(1)),INT(x(2)))
        IF (sigma .LT. domain_in(1,3) .OR. sigma .GT. domain_in(2,3)) THEN
          IF (verbal .GT. 0) PRINT*,'Out of vertical bounds',k,station_names(k)
          CYCLE
        ENDIF
        observations%numobs = observations%numobs+1
        observations%xyzt(1:2,observations%numobs) = x(1:2)
        observations%xyzt(3,observations%numobs) = &
          (sigma-domain_in(1,3))/(domain_in(2,3)-domain_in(1,3))*(ngrid_in(3)-1)+1
        observations%xyzt(4,observations%numobs) = &
          (obs_time(k)-domain_in(1,4))/(domain_in(2,4)-domain_in(1,4))*(ngrid_in(4)-1)+1

        ! Pass station names and report types:
        observations%stnames(observations%numobs)(1:5) = station_names(k)(1:5)

        ! Pass observation values:
        DO iv=1,num_vars_in  ! ** through all analysis variables
          ! Check through all variables
          SELECT CASE(var_names_in(iv)(1:2))
            CASE('U3')
              IF (var_names_in(iv+1) .NE. 'V3') THEN
                PRINT*,header,'Please use V3 following U3 as variable name!'
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
                                        observations%value(observations%numobs,iv), &
                                        observations%value(observations%numobs,iv+1),lon(k))

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
                observations%error(observations%numobs,iv  ) = &
                  MAXVAL(ABS(uerr(1:4)-observations%value(observations%numobs,iv  )))
                observations%error(observations%numobs,iv+1) = &
                  MAXVAL(ABS(verr(1:4)-observations%value(observations%numobs,iv+1)))
                !**************************************

                ! Convert obs in knots to m/s
                observations%value(observations%numobs,iv:iv+1) = &
                  observations%value(observations%numobs,iv:iv+1)*knot2ms
                ! Convert obs errors in knots to absolute m/s:
                observations%error(observations%numobs,iv:iv+1) = &
                  observations%error(observations%numobs,iv:iv+1)*knot2ms
              ENDIF
            CASE('V3')
              IF (var_names_in(iv-1) .NE. 'U3') THEN
                PRINT*,'Please use U3 before V3 as variable name!'
                STOP
              ENDIF
            CASE('W3')
              ! No vertical wind insitu observations
            CASE('SH')
              ! LAPS dewpint is in Fahrenheit:
              observations%value(observations%numobs,iv) = STMAS_invalid  ! Default
              IF (dew(k) .NE. missing .AND. dew(k) .NE. badsfc .AND. &
                  temp_2m(k) .NE. missing .AND. temp_2m(k) .NE. badsfc .AND. &
                  stn_prs(k) .NE. missing .AND. stn_prs(k) .NE. badsfc) THEN
                ! From valid temperature, dewpoint and pressure obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                ! SSH2 take mb, c, c units:
                tc = (temp_2m(k)-32.0)*5.0/9.0
                dc = (dew(k)-32.0)*5.0/9.0
                ! Output SH in kg/kg as LAPS ssh2 in g/kg:
                observations%value(observations%numobs,iv) = ssh2(stn_prs(k),tc,dc,t_liquid)*0.001
              ELSEIF (rh(k) .NE. missing .AND. rh(k) .NE. badsfc .AND. &
                      temp_2m(k) .NE. missing .AND. temp_2m(k) .NE. badsfc .AND. &
                      stn_prs(k) .NE. missing .AND. stn_prs(k) .NE. badsfc) THEN
                ! From valid RH, temperature and pressure obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                ! make_ssh takes mb, c,fraction:
                tc = (temp_2m(k)-32.0)*5.0/9.0
                ! Output SH in kg/kg as LAPS make_ssh in g/kg:
                observations%value(observations%numobs,iv) = make_ssh(stn_prs(k),tc,rh(k),t_liquid)*0.001
print*,'SH from RH: ',observations%value(observations%numobs,iv),rh(k),temp_2m(k),stn_prs(k)
stop
              ENDIF

              observations%error(observations%numobs,iv) = STMAS_invalid  ! Default
              IF (dew_err(k) .NE. missing .AND. dew_err(k) .NE. badsfc) THEN
                observations%error(observations%numobs,iv) = dew_err(k)*5.0/9.0
              ELSEIF (rh_err(k) .NE. missing .AND. rh_err(k) .NE. badsfc .AND. &
                      temp_err(k) .NE. missing .AND. temp_err(k) .NE. badsfc) THEN
                tc = (temp_err(k)-32.0)*5.0/9.0
                observations%error(observations%numobs,iv) = dwpt(tc,rh(k)) ! Error in C is the same as in K
              ENDIF
            CASE('T3')
              ! LAPS temp is in Fahrenheit:
              observations%value(observations%numobs,iv) = STMAS_invalid  ! Default
              IF (temp_2m(k) .NE. missing .AND. temp_2m(k) .NE. badsfc) THEN
                ! Valid temperature obs:
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                observations%value(observations%numobs,iv) = (temp_2m(k)-32.0)*5.0/9.0+temp_00
If (observations%numobs .eq. 1) print*,'Surface T: ',observations%value(observations%numobs,iv),iv
              ENDIF

              observations%error(observations%numobs,iv) = STMAS_invalid  ! Default
              IF (temp_err(k) .NE. missing .AND. temp_err(k) .NE. badsfc) &
                observations%error(observations%numobs,iv) = temp_err(k)*5.0/9.0
            CASE('P3')
              IF (stn_prs(k) .NE. missing .AND. stn_prs(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                observations%value(observations%numobs,iv) = stn_prs(k)*100.0  ! LAPS in mb STMAS uses Pascal
              ELSEIF (altimeter(k) .NE. missing .AND. &
                      altimeter(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                observations%value(observations%numobs,iv) = &
                  ALT_2_SFC_PRESS(altimeter(k),elevation(k))
                observations%value(observations%numobs,iv) = &
                  observations%value(observations%numobs,iv)*100.0  ! Convert LAPS mb to Pascal
              ELSEIF (msl_prs(k) .NE. missing .AND. &
                      msl_prs(k) .NE. badsfc  .AND. &
                      temp_2m(k) .NE. missing .AND. &
                      temp_2m(k) .NE. badsfc  .AND. &
                      dew(k) .NE. missing .AND. dew(k) .NE. badsfc) THEN
                iuse_stn = iuse_stn+1
                totalobs(iv) = totalobs(iv)+1
                CALL REDUCE_P(temp_2m(k),dew(k), &
                              msl_prs(k),0.0,lapses(1),lapses(2), &
                              observations%value(observations%numobs,iv),elevation(k),STMAS_invalid)
                observations%value(observations%numobs,iv) = &
                  observations%value(observations%numobs,iv)*100.0  ! Convert LAPS mb to Pascal
              ELSE
                observations%value(observations%numobs,iv) = STMAS_invalid
              ENDIF            
          
              IF (prs_err(k) .NE. missing .AND. prs_err(k) .NE. badsfc) THEN
                observations%error(observations%numobs,iv) = prs_err(k)*100.0  ! Convert LAPS mb to Pascal
              ELSEIF (alt_err(k) .NE. missing .AND. &
                      alt_err(k) .NE. badsfc) THEN
                observations%error(observations%numobs,iv) = &
                  ALT_2_SFC_PRESS(alt_err(k),elevation(k))
                observations%error(observations%numobs,iv) = &
                  observations%error(observations%numobs,iv)*100.0   ! Convert LAPS mb to Pascal
              ELSE
                observations%error(k,iv) = STMAS_invalid
              ENDIF
            CASE DEFAULT
              PRINT*
              PRINT*,header,'No such observations: ',var_names_in(iv)
              STOP
          END SELECT
        ENDDO ! ** through all analysis variables

        ! If no valid obs for the analysis variables, remove this record
        IF (iuse_stn .EQ. 0) THEN
          observations%numobs = observations%numobs-1
          IF (verbal .NE. 0) &
            print*,'No valid obs at this station: ',station_names(k)(1:10)
        ENDIF

      ELSE
        out_xy = out_xy+1
        IF (verbal .GE. 1) THEN
          PRINT*
          PRINT*,header,'station: ',station_names(k),' out of x-y domain', &
            x(1),x(2)
          PRINT*,' (',domain_in(1:2,1:2)/grid_spacing_m+1,')'
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
        ELSE IF (obs_time(k) .LT. domain_in(1,4) .OR. &
                 obs_time(k) .GT. domain_in(2,4)) THEN
          out_t = out_t+1
          PRINT*,' out of time window: ',obs_time(k),domain_in(1:2,4)
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
  PRINT*,header,'| Total number of surface obs in domain',observations%numobs-nbefore
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

  PRINT*
  PRINT*,'+----------------------------------+'
  PRINT*,'|  End of reading in surface data  |'
  PRINT*,'+----------------------------------+'


END SUBROUTINE insitu_height_surfaces
