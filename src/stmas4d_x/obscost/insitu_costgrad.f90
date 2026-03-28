SUBROUTINE insitu_costgrad
!====================================================================
!  This routine adds insitu observation terms to the cost function.
!
!  History: Yuanfu Xie in Dec 2011
!====================================================================

  USE STMAS

  IMPLICIT NONE

  ! Local variables:
  INTEGER :: io,iv
  DOUBLE PRECISION :: cost ! Ensure summation accuracy

  cost = 0.0D0

  DO iv=1,gridobs%numvars
    DO io=1,gridobs%numobs
      IF (gridobs%value(io,iv) .NE. STMAS_invalid) THEN
        cost = cost + (STMAS_mgrid%anal(gridobs%ixyzt(1,io), &
                                        gridobs%ixyzt(2,io), &
                                        gridobs%ixyzt(3,io), &
                                        gridobs%ixyzt(4,io),iv) &
                      -gridobs%value(io,iv))**2
write(*,1) gridobs%ixyzt(1:4,io),io,iv
1 format('idx: ',6i6)
  print*,'Ana vs obs: ',STMAS_mgrid%anal(gridobs%ixyzt(1,io), &
                                         gridobs%ixyzt(2,io), &
                                         gridobs%ixyzt(3,io), &
                                         gridobs%ixyzt(4,io),iv),gridobs%value(io,iv),cost
if (iv .eq. 4 .and. io .eq. 6) then
print*,'Cost jump'
endif
      ENDIF
    ENDDO
  ENDDO

print*,'Cost: ',cost
END SUBROUTINE insitu_costgrad
