SUBROUTINE quality_control
!doc====================================================================
!doc  This routine calculates the innovations of the observations and 
!doc  determines the weights of these observations in the cost function.
!doc
!doc  History: Yuanfu Xie May 2012
!doc====================================================================

  USE STMAS

  IMPLICIT NONE

  ! Local variables:
  INTEGER :: iv,io,nc(20),nchannel(20)

  ! Output bias channels:
  DO iv=1,20
    nchannel(iv) = 20+iv
  ENDDO

  ! Observation counter for each variable:
  nc = 0
  DO io=1,observations%numobs

    DO iv=1,observations%numvars

      IF (observations%value(io,iv) .EQ. STMAS_invalid) CYCLE

      nc(iv) = nc(iv)+1

      WRITE(nchannel(iv),1) nc(iv),observations%value(io,iv)-observations%bkgd(io,iv)
1     FORMAT(i6,e15.4)

    ENDDO

  ENDDO

  ! Observation weights in the cost function:
  ! Temporarily set as constants:
  observations%weights(1) = 1.0  ! U
  observations%weights(2) = 1.0  ! V
  observations%weights(3) = 1.0  ! W
  observations%weights(4) = 1.0  ! P
  observations%weights(5) = 1.0  ! T
  observations%weights(6) = 1.0  ! SH

END SUBROUTINE quality_control
