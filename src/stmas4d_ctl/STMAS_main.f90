PROGRAM STMAS_main
!doc==================================================================
!doc This is a main driver program for STMAS surface analysis.
!doc
! History: Dec. 2009 by Yuanfu Xie.
!doc		Recoding of the surface analysis: 
!doc		1. Multivariate analysis P&T reduced to surface; U&V;
!doc		2. Constraints of normal velocity to zero
!doc==================================================================

  ! USE LAPS_ingest        ! LAPS data ingest module
  USE STMAS                ! STMAS analysis module

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*13 :: time_window(2) 
  INTEGER :: i,j

  ! Initializing:
  STMAS_success = 1

  ! Read STMAS namelist:
  PRINT*
  PRINT*,'*********************************************************'
  PRINT*,'*                                                       *'
  PRINT*,'*             STMAS: analysis begins                    *'
  PRINT*,'*                                                       *'
  PRINT*,'*********************************************************'
  
  PRINT*
  PRINT*,'---------------------------------------------------------'
  PRINT*,'|                                                       |'
  PRINT*,'|             STMAS: Reading namelist                   |'
  PRINT*,'|                                                       |'
  PRINT*,'---------------------------------------------------------'
  CALL read_namelist

  PRINT*
  PRINT*,'---------------------------------------------------------'
  PRINT*,'|                                                       |'
  PRINT*,'|            STMAS: Calling STMAS_LAPS_INPUT            |'
  PRINT*,'|                                                       |'
  PRINT*,'---------------------------------------------------------'
  CALL STMAS_LAPS_input

  PRINT*
  PRINT*,'---------------------------------------------------------'
  PRINT*,'|                                                       |'
  PRINT*,'|            STMAS: Calling STMAS_Analyzing...          |'
  PRINT*,'|                                                       |'
  PRINT*,'---------------------------------------------------------'
  CALL STMAS_analysis

!  PRINT*
!  PRINT*,'---------------------------------------------------------'
!  PRINT*,'|                                                       |'
!  PRINT*,'|         STMAS: Output analysis                        |'
!  PRINT*,'|                                                       |'
!  PRINT*,'---------------------------------------------------------'
  CALL STMAS_output

  PRINT*
  PRINT*,'---------------------------------------------------------'
  PRINT*,'|                                                       |'
  PRINT*,'|             STMAS: Deallocating memory                |'
  PRINT*,'|                                                       |'
  PRINT*,'---------------------------------------------------------'
  ! STMAS memory deallocation:
  CALL STMAS_memo_release

  IF (STMAS_success .EQ. 1) THEN
    ! Success of analysis:
    PRINT*
    PRINT*,'*********************************************************'
    PRINT*,'*                                                       *'
    PRINT*,'*            STMAS analysis succeeds                    *'
    PRINT*,'*                                                       *'
    PRINT*,'*********************************************************'
  ELSE
    PRINT*
    PRINT*,'*********************************************************'
    PRINT*,'*                                                       *'
    PRINT*,'*             STMAS analysis fails                      *'
    PRINT*,'*                                                       *'
    PRINT*,'*********************************************************'
  ENDIF

END PROGRAM STMAS_main

