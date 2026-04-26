SUBROUTINE get_analysis_window(i4time,time_window_out,before_in,after_in)

!doc==================================================================
!doc This routine is provided by a user to determine the current time
!doc for analysis and use the information of before/after in seconds
!doc to return correct analysis time window in char 13 format.
!doc
!doc History: March 2010 by Yuanfu Xie.
!doc
!doc Parameters:
!doc   Input:
!doc     before_in:           (int) Seconds before the current time (<0)
!doc     after_in:            (int) Seconds after the current time (>=0)
!doc
!doc   Output:
!doc     time_window:         (char*13)x2 Time-window in char 13.
!doc==================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: before_in,after_in
  CHARACTER*13, INTENT(OUT) :: time_window_out(2)
  INTEGER, INTENT(OUT) :: i4time

  ! Local variables:
  INTEGER :: istatus,seconds,year,month,day,hour,minute

  ! No zero time interval allowed:
  IF (after_in-before_in .LE. 0) THEN
    PRINT*,'In subroutine get_analysis_window: no zero time window allowed!', &
      ' check before and after in stmasfc.nl and rerun, please!'
    STOP
  ENDIF

  ! Use LAPS system time:
  CALL GET_SYSTIME_I4(i4time,istatus)

  ! Start time:
  CALL CV_I4TIM_INT_LP(i4time+before_in,year,month,day,hour,minute, &
                        seconds,istatus)
  IF (seconds .GT. 0) THEN
    PRINT*
    PRINT*,'In subroutine get_analysis_window: Cannot handle seconds', &
      ' check parameter BEFORE in namelist and rerun please!',before_in
    STOP
  ENDIF

  ! Year is the years from 1900:
  year = year+1900

  WRITE(time_window_out(1)(1:13),1) year,month,day,hour,minute
1 FORMAT(i4,i2,i2,'_',i2,i2)
  IF (month  .LE. 9) time_window_out(1)(5:5) = '0'
  IF (day    .LE. 9) time_window_out(1)(7:7) = '0'
  IF (hour   .LE. 9) time_window_out(1)(10:10) = '0'
  IF (minute .LE. 9) time_window_out(1)(12:12) = '0'

  ! End time:
  CALL CV_I4TIM_INT_LP(i4time+after_in,year,month,day,hour,minute, &
                        seconds,istatus)
  IF (seconds .GT. 0) THEN
    PRINT*
    PRINT*,'STMAS>get_analysis_window: Cannot handle seconds', &
      ' check parameter AFTER in namelist and rerun please!',after_in
    STOP
  ENDIF

  ! Year is the years from 1900:
  year = year+1900

  WRITE(time_window_out(2)(1:13),1) year,month,day,hour,minute
  IF (month  .LE. 9) time_window_out(2)(5:5) = '0'
  IF (day    .LE. 9) time_window_out(2)(7:7) = '0'
  IF (hour   .LE. 9) time_window_out(2)(10:10) = '0'
  IF (minute .LE. 9) time_window_out(2)(12:12) = '0'

END SUBROUTINE get_analysis_window
