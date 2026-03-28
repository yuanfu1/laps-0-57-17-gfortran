SUBROUTINE COST_OBS(func,grad,num_controls,mgrid)
!==============================================================================================
!  This routine evaluates the observation term for the STMASFC cost function and its gradients.
!
!  History: Feb. 2011 by Yuanfu Xie at NOAA/ESRL/GSD.
!
!  Description:
!    This code assumes that the actual observations mapped to a gridded observation data 
!    structure, STMASFC_gridded_obs. As discussed by a presentation LAPS_STMAS2011Feb10.pptx
!    under ~xiey/work/prjcts/STMAS/ on jinx, the observation term is simply to compute the 
!    differences between observations and analysis at the grid points where the gridded obs 
!    resides.
!==============================================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(OUT) :: func,grad(num_controls)
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  INTEGER :: nobs,nv,io,i,j,k
  REAL    :: term(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3)) ! Save each obs term
  DOUBLE PRECISION :: f                           ! Save function value reducing round off

  ! Pass through all variables:
  DO nv=1,gridobs%numvars

    ! Re-initialize for each variables:
    f = 0.0d0
    nobs = 0
    term = 0.0

    ! Pass through all observations:
    DO io=1,gridobs%numobs

      IF (gridobs%value(io,nv) .NE. STMAS_invalid) THEN

        nobs = nobs+1

        ! Save indirect addressing:
        i = (gridobs%ixyzt(1,io)-1)/mgrid%incr(1)+1
        j = (gridobs%ixyzt(2,io)-1)/mgrid%incr(2)+1
        k = (gridobs%ixyzt(3,io)-1)/mgrid%incr(3)+1

        ! Each obs term:
        term(i,j,k) = (mgrid%anal(gridobs%ixyzt(1,io), &
                                  gridobs%ixyzt(2,io), &
                                  gridobs%ixyzt(3,io),nv)- &
                       gridobs%value(io,nv))/gridobs%error(io,nv)

        ! Sum of obs terms:
        f = f + term(i,j,k)*term(i,j,k)

        ! Save term as a gradient before scaling:
        term(i,j,k) = term(i,j,k)/gridobs%error(io,nv)
      ENDIF

    ENDDO

    IF (nobs .GT. 0) THEN

      ! Function:
      func = func+0.5*f/nobs
      ! Scaling gradient:
      DO k=1,mgrid%numgrid(3)
        DO j=1,mgrid%numgrid(2)
          DO i=1,mgrid%numgrid(1)
            grad(i  +mgrid%numgrid(1)*( &
                 j-1+mgrid%numgrid(2)*( &
                 k-1+mgrid%numgrid(3)*(nv-1)))) = &
              term(i,j,k)/nobs
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  ENDDO

END SUBROUTINE COST_OBS
