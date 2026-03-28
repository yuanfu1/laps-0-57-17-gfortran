SUBROUTINE back2obs
!doc====================================================================
!doc  This routine interpolates background fields to their values at the
!doc  observation sites.
!doc
!doc  History: Yuanfu Xie Dec 2011
!doc====================================================================

  USE STMAS

  IMPLICIT NONE

  ! Local variables:
  INTEGER :: iv,io,idx(2,4),i,j,k,l
  REAL    :: coe(2,4)

  ! Every obs:
  DO io=1,observations%numobs

    ! Indices:
    idx(1,1:4) = observations%xyzt(1:4,io)
    idx(2,1:4) = idx(1,1:4)
    DO iv=1,4
      IF (idx(2,iv) .GT. STMAS_final%numgrid(iv)) &
        idx(2,iv) = STMAS_final%numgrid(iv)
    ENDDO

    ! Coefficients: assuming obs location on grid scale
    coe(2,1:4) = observations%xyzt(1:4,io)-idx(1,1:4)
    coe(1,1:4) = 1.0-coe(2,1:4)

    ! Land-water and topo at obs:
    observations%land(io) = 0.0
    observations%topo(io) = 0.0
    DO j=1,2
    DO i=1,2
      observations%land(io) = observations%land(io)+ &
        coe(i,1)*coe(j,2)*STMAS_final%land(idx(i,1),idx(j,2))
      observations%topo(io) = observations%topo(io)+ &
        coe(i,1)*coe(j,2)*STMAS_final%topo(idx(i,1),idx(j,2))
    ENDDO
    ENDDO

    ! Every states:
    DO iv=1,observations%numvars
      observations%bkgd(io,iv) = 0.0
      DO l=1,2
      DO k=1,2
      DO j=1,2
      DO i=1,2
        observations%bkgd(io,iv) = observations%bkgd(io,iv)+ &
          coe(i,1)*coe(j,2)*coe(k,3)*coe(l,4)* &
          STMAS_final%bkgd(idx(i,1),idx(j,2),idx(k,3),idx(l,4),iv)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE back2obs
