!doc==================================================================
!doc
!doc This routine reads in background, obs, etc into STMAS structure
!doc using LAPS libraries and data root.
!doc
!doc History: April 2010 by Yuanfu Xie.
!doc
!doc This is a copy from STMASFC_LAPS_input.f90, and modified for 4d
!doc June 2011 by Hongli Jiang
!doc==================================================================

SUBROUTINE STMAS_LAPS_input

!  USE ingest_all
  USE LAPS_ingest
  USE STMAS

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*13 :: ctime_window(2)
  CHARACTER*20 :: header = 'STMAS_LAPS_input: '
  INTEGER :: i,j,n

  ! Analysis time window:
  CALL get_analysis_window(LAPS_i4time,ctime_window, &
         STMAS_time_window(1), STMAS_time_window(2))

  WRITE(*,1) header,ctime_window(1:2)
1 FORMAT(/,A20,'Requested time window: from ',A13,' to ',A13/)
  print*,'laps_i4time',laps_i4time

  ! Analysis grid dimensions:
  CALL get_LAPS_dimension(LAPS_i4time,ctime_window,STMAS_final%numgrid, &
                                      STMAS_domain,STMAS_final%gridspc)

  WRITE(*,2) header,STMAS_final%numgrid(1:4) !HJ change from 1:3 to 1:4 6/7/2011
2 FORMAT(A20,'Number of gridpoints of the analysis grid: ',4i5)

!  DO i=1,STMAS_maxdim
!    WRITE(*,3) i,STMAS_domain(1:2,i)
!  ENDDO
!3 FORMAT(3X,'Domain in dimention ',i2,' from ',e10.4,' to ',e10.4)

!  WRITE(*,4) header,STMAS_final%gridspc(1:4)
!4 FORMAT(A20,'Gridspacing of the analysis grid: ',/,4e12.4)

  ! Compute the maximum number of mulgrid levels in all dimensions:
  DO i=1,STMAS_maxdim
    n = STMAS_start_grdpts(i)
    STMAS_maxlevels(i) = 1
    ! Count the start grid as one level and then keep counting for the rest levels:
    DO j=2,STMAS_numlevels
      IF (2*(n-1)+1 .LE. STMAS_final%numgrid(i)) THEN
        n = 2*(n-1)+1
        STMAS_maxlevels(i) = STMAS_maxlevels(i)+1
      ENDIF
    ENDDO
  ENDDO

  PRINT*
  PRINT*,'------------------------------------------------------------'
  PRINT*,'|                                                          |'
  PRINT*,'|    STMAS_LAPS_input: array allocation                    |'
  PRINT*,'|                                                          |'
  PRINT*,'------------------------------------------------------------'

  CALL STMAS_memo_alloc

  PRINT*
  PRINT*,'------------------------------------------------------------'
  PRINT*,'|                                                          |'
  PRINT*,'|    STMAS_LAPS_input: get_LAPS_config:                    |'
  PRINT*,'|                      Lat/Lon/Topography/Land/Map-factor  |'
  PRINT*,'|                                                          |'
  PRINT*,'------------------------------------------------------------'
  CALL get_LAPS_config(STMAS_final%numgrid,STMAS_final%lat, &
                       STMAS_final%lon,STMAS_final%topo, & 
                       STMAS_final%land,STMAS_final%mapf, &
                       LAPS_rdplvl)

  PRINT*
  PRINT*,'------------------------------------------------------------'
  PRINT*,'|                                                          |'
  PRINT*,'|    STMAS_LAPS_input: get_LAPS_back                       |'
  PRINT*,'|                                                          |'
  PRINT*,'------------------------------------------------------------'
! Get background fields: and readin sigma_ht. Otherwise, remove the last argument. HJ 7/29/2011
Print*,'YUAN: get_laps_bkgd',STMAS_varnames(4)

  CALL get_LAPS_bkgd(STMAS_numvars,STMAS_varnames, &
                     STMAS_final%numgrid,STMAS_final%bkgd, &
                     STMAS_final%zz)
Print*,'YUAN end: get_laps_bkgd',STMAS_final%bkgd(1,1,1,1,1:6)
!
! For sigma_ht. HJ 8/29/2011
  PRINT*,'------------------------------------------------------------'
  PRINT*,'|                                                          |'
  PRINT*,'|    STMAS_LAPS_input: calculating Jacobian from sigma_ht  |'
  PRINT*,'|                                                          |'
  PRINT*,'------------------------------------------------------------'
   CALL jacobian_cal(STMAS_final%numgrid,STMAS_final%gridspc,STMAS_final%topo, &
                  STMAS_final%jt1,STMAS_final%jt2,STMAS_final%jt3,STMAS_final%zz, &
                  STMAS_final%cor,STMAS_final%bkgd,STMAS_numvars)
Print*,'YUAN End: get_laps_bkgd',STMAS_final%bkgd(1,1,1,1,1:6)

! Define the grid info for the final grid:
  STMAS_final%nfinest = STMAS_final%numgrid
  STMAS_final%incr = 1

  ! PRINT*
  ! PRINT*,'---------------------------------------------------------'
  ! PRINT*,'|                                                       |'
  ! PRINT*,'|           STMAS: Reading observation data             |'
  ! PRINT*,'|                                                       |'
  ! PRINT*,'---------------------------------------------------------'

  ! ! Pass the namelist variable values to observations:
  ! observations%numvars = STMAS_numvars

  ! CALL get_LAPS_obs(STMAS_domain,STMAS_numvars,STMAS_varnames, &
  !                   STMAS_invalid,STMAS_final%numgrid, &
  !                   STMAS_final%lat,STMAS_final%lon, &
  !                   STMAS_final%topo,STMAS_maxobs, &
  !                   observations%numobs,observations%value, &
  !                   observations%error,observations%xyzt, &
  !                   observations%lat,observations%lon, &
  !                   observations%stnames,observations%types, &
  !                   STMAS_debugging,STMAS_success)

!hj  CALL background_at_obs(STMAS_final,observations)

END SUBROUTINE STMAS_LAPS_input

SUBROUTINE background_at_obs(bkg,obs)
!doc==================================================================
!doc
!doc This routine interpolates the background fields to observation's
!doc location and time and saves the values in the observation data
!doc structure. In the mean time, a preliminary QC threshold value 
!doc check is applied.
!doc
!doc History: October 2010 by Yuanfu Xie.
!doc
!doc==================================================================

  USE STMAS

  IMPLICIT NONE

  TYPE(STMAS_bkgd), INTENT(IN) :: bkg
  TYPE(STMAS_obs), INTENT(INOUT) :: obs

  ! Local variables:
  CHARACTER*10 :: sttname(STMAS_numvars)
  CHARACTER*41, PARAMETER :: header = 'STMAS_LAPS_input: background_at_obs: '
  DOUBLE PRECISION :: std_deviation(STMAS_numvars)
  INTEGER :: i,j,k,t,TEMP,DEWP,PRES
! HJ 4 was 3. 6/10/2011
  INTEGER :: iobs,ivar,indx(2,4),i_o(STMAS_numvars),nvalid(STMAS_numvars)
  REAL    :: coef(2,4),b_o(STMAS_numvars),dff,temperature,dewpoint

  ! Find the index for obs variables:
  TEMP = 0
  DEWP = 0
  PRES = 0
  DO ivar=1,obs%numvars
    IF (STMAS_varnames(ivar) .EQ. 'TSF') TEMP = ivar
    IF (STMAS_varnames(ivar) .EQ. 'DSF') DEWP = ivar
    IF (STMAS_varnames(ivar) .EQ. 'PSF') PRES = ivar
  ENDDO

  ! For all observations:
  b_o = 0.0
  std_deviation = 0.0d0
  nvalid = 0
  DO iobs=1,obs%numobs

    ! Compute the interpolation coefficients:
    indx(1,1:4) = INT((obs%xyzt(1:4,iobs)-STMAS_domain(1,1:4))/STMAS_final%gridspc(1:4))+1
    coef(2,1:4) = (obs%xyzt(1:4,iobs)-(indx(1,1:4)-1)*STMAS_final%gridspc(1:4)-STMAS_domain(1,1:4))/ &
                  STMAS_final%gridspc(1:4)
    indx(2,1:4) = MIN(indx(1,1:4)+1,bkg%numgrid(1:4))
    coef(1,1:4) = 1.0-coef(2,1:4)            ! Assuming grid distance 1

    ! For all observation variables, compute the background values at the obs site:
    DO ivar=1,obs%numvars
      obs%bkgd(iobs,ivar) = 0.0
! HJ add t loop 8/18/2011
     DO t=1,2
      DO k=1,2
        DO j=1,2
          DO i=1,2
            obs%bkgd(iobs,ivar) = obs%bkgd(iobs,ivar)+ &
              bkg%bkgd(indx(i,1),indx(j,2),indx(k,3),indx(t,4),ivar)* &
              coef(i,1)*coef(j,2)*coef(k,3)*coef(t,4)
          ENDDO
        ENDDO
       ENDDO
      ENDDO
    ENDDO

    ! For all terrain and land at observation location as they will be used for adjust obs:
    obs%land(iobs)=0.0    ! land
    DO j=1,2
      DO i=1,2
        obs%land(iobs) = obs%land(iobs)+ &
          bkg%land(indx(i,1),indx(j,2))* &
          coef(i,1)*coef(j,2)
      ENDDO
    ENDDO
    obs%topo(iobs)=0.0    ! topo
    DO j=1,2
      DO i=1,2
        obs%topo(iobs) = obs%topo(iobs)+ &
          bkg%topo(indx(i,1),indx(j,2))* &
          coef(i,1)*coef(j,2)
      ENDDO
    ENDDO

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  Adjust obs according the obs elevation and domain topography height at obs
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Use LAPS lapse rates for temperature and dew point: 
    IF (TEMP .NE. 0 .AND. obs%value(iobs,TEMP) .NE. STMAS_invalid) &
      obs%value(iobs,TEMP) = obs%value(iobs,TEMP)- &
        lapses(1)*(obs%xyzt(4,iobs)-obs%topo(iobs))
    IF (DEWP .NE. 0 .AND. obs%value(iobs,DEWP) .NE. STMAS_invalid) &
      obs%value(iobs,DEWP) = obs%value(iobs,DEWP)- &
        lapses(2)*(obs%xyzt(4,iobs)-obs%topo(iobs))
    ! For saturated obs:
    IF (obs%value(iobs,DEWP) .NE. STMAS_invalid .AND. &
        obs%value(iobs,DEWP) .GE. obs%value(iobs,TEMP)) &
      obs%value(iobs,DEWP)=obs%value(iobs,TEMP)

    ! Reduce pressure to analysis terrain: 
    ! using LAPS reduce_p routine and obs of t, dt, p:
    IF (TEMP .NE. 0 .AND. obs%value(iobs,TEMP) .NE. STMAS_invalid .AND. &
        DEWP .NE. 0 .AND. obs%value(iobs,DEWP) .NE. STMAS_invalid .AND. &
        PRES .NE. 0 .AND. obs%value(iobs,PRES) .NE. STMAS_invalid) THEN
      CALL REDUCE_P(1.8*(obs%value(iobs,TEMP)-T00)+32.0, &
                    1.8*(obs%value(iobs,DEWP)-T00)+32.0, &
                    obs%value(iobs,PRES),obs%xyzt(4,iobs),lapses(1),lapses(2), &
                    obs%value(iobs,PRES),obs%topo(iobs),STMAS_invalid)
    !ELSEIF (PRES .NE. 0 .AND. obs%value(iobs,PRES) .NE. STMAS_invalid) THEN
      ! Using LAPS reduce_p routine and obs and background to fill t or dt:
      ! IF (TEMP .EQ. 0) THEN
      !   temperature = STMAS_invalid
      ! ELSE
      !   temperature = obs%value(iobs,TEMP)
      ! ENDIF
      ! IF (DEWP .EQ. 0) THEN
      !   dewpoint = STMAS_invalid
      ! ELSE
      !   dewpoint = obs%value(iobs,DEWP)
      ! ENDIF

      ! Notice using '+' sign to adjust: from analysis terrain to obs terrain:
      ! IF (temperature .EQ. STMAS_invalid) temperature = obs%bkgd(iobs,TEMP)+ &
      !   lapses(1)*(obs%xytz(4,iobs)-obs%topo(iobs))
      ! IF (dewpoint .EQ. STMAS_invalid) dewpoint = obs%bkgd(iobs,DEWP)+ &
      !   lapses(2)*(obs%xytz(4,iobs)-obs%topo(iobs))
      ! Saturation:
      ! IF (dewpoint .GE. temperature) dewpoint = temperature

      ! Reduce pressure to analysis terrain using obs and background
      ! CALL REDUCE_P(1.8*(temperature-T00)+32.0, &
      !               1.8*(dewpoint-T00)+32.0, &
      !               obs%value(iobs,PRES),obs%xytz(4,iobs),lapses(1),lapses(2), &
      !               obs%value(iobs,PRES),obs%topo(iobs),STMAS_invalid)

    ELSEIF (ABS(obs%xyzt(4,iobs)-obs%topo(iobs)) .GE. 50.0) THEN
      ! Too large height difference but cannot reduce:
      obs%value(iobs,PRES) = STMAS_invalid
    ENDIF

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  End of adjustment of these obs
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! For all observed fields:
    DO ivar=1,obs%numvars

      ! Skip invalid data:
      IF (obs%value(iobs,ivar) .EQ. STMAS_invalid) cycle

      ! Standard deviation from the background computation:
      dff = obs%bkgd(iobs,ivar)-obs%value(iobs,ivar)
      std_deviation(ivar) = std_deviation(ivar) + dff*dff
      nvalid(ivar) = nvalid(ivar)+1

      ! Plot innovation distribution:
      if (STMAS_debugging .EQ. ivar) write(30,*) nvalid(STMAS_debugging),dff

    ENDDO

  ENDDO

  ! Compute standard deviation:
  DO ivar=1,obs%numvars
    ! Standard deviation:
    IF (nvalid(ivar) .GT. 1) &
      std_deviation(ivar) = SQRT(std_deviation(ivar)/nvalid(ivar))
  ENDDO

  ! Standard deviation QC:
  PRINT*,header
  PRINT*,'+----------------------------------------------------+'
  PRINT*,'|       Standard deviation QC in progress...         |'
  PRINT*,'+----------------------------------------------------+'
  PRINT*
  DO iobs=1,obs%numobs
    ! For all observed fields:
    DO ivar=1,obs%numvars
      ! If background and observation differs greater than the factor*std_dev, remove:
      IF (ABS(obs%bkgd(iobs,ivar)-obs%value(iobs,ivar)) .GT. std_deviation(ivar)*STMAS_stddev(ivar)) THEN
        obs%value(iobs,ivar) = STMAS_invalid
      ELSE
        ! Debugging: find out the largest innovation for each var
        IF (STMAS_debugging .GT. 0) THEN
          IF (b_o(ivar) .LT. ABS(obs%bkgd(iobs,ivar)-obs%value(iobs,ivar))) THEN
            b_o(ivar) =  ABS(obs%bkgd(iobs,ivar)-obs%value(iobs,ivar))
            sttname(ivar) = obs%stnames(iobs)
            i_o(ivar) = iobs
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! Print debugging message:
  IF (STMAS_debugging .GT. 0) THEN
    DO ivar=1,obs%numvars
      PRINT*,header
      WRITE(*,1) ivar,b_o(ivar),sttname(ivar),i_o(ivar)
1     FORMAT('    Max innovation in obs: ',i2,' is, ',f10.2,' at: ',A10,' at obs: ',i8)

      PRINT*,header
      WRITE(*,2) nvalid(ivar),std_deviation(ivar),STMAS_varnames(ivar),STMAS_stddev(ivar)
2     FORMAT('    Number of valid obs: ',i8,/, &
             ' Standard Deviation: ',f12.4,' for ',a8,' Factors of deviation: ',f6.1)
    ENDDO
  ENDIF

END SUBROUTINE background_at_obs

