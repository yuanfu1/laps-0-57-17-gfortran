SUBROUTINE qv_costgrad(func,num_controls,mgrid)

!doc=============================================================================================
!doc This routine computes the cost function and its gradient of moisture eqn for STMAS4D
!doc
!doc History: 6/24/2011 by Hongli Jiang
!doc=============================================================================================

  USE STMAS, ONLY: STMAS_bkgd,STMAS_penal

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(OUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  INTEGER :: i,j,k,t,i1,i2,j1,j2,k1,k2,t1
  REAL    :: qx,qy,qz,j13,j23,dzk,tdx,tdy,penal,CE
  REAL    :: amu(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amv(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amw(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amq(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  DOUBLE PRECISION :: f                     ! Save function value reducing round off

  tdx=2.*mgrid%gridspc(1)
  tdy=2.*mgrid%gridspc(2)
! add the factor 1.e+9 here, and reduce penal(6) in the file stmas_vars.nl
  penal=1.e+9*STMAS_penal(6)

  DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
      DO j=1,mgrid%numgrid(2)
        DO i=1,mgrid%numgrid(1)
          i1=i-1
          j1=j-1
          k1=k-1
          t1=t-1
          amu(i,j,k) = mgrid%anal(i,j,k,t,1)
          amv(i,j,k) = mgrid%anal(i,j,k,t,2)
          amw(i,j,k) = mgrid%anal(i,j,k,t,3)
          amq(i,j,k) = mgrid%anal(i,j,k,t,6)
        ENDDO
      ENDDO
    ENDDO

    f = 0.0d0
    CE = 0.0

    DO k=2,mgrid%numgrid(3)-1
      DO j=2,mgrid%numgrid(2)-1
        DO i=2,mgrid%numgrid(1)-1
          j13=mgrid%jt1(i,j,k)/mgrid%jt3(i,j)
          j23=mgrid%jt2(i,j,k)/mgrid%jt3(i,j)
          dzk=mgrid%zz(k+1)-mgrid%zz(k-1)
          i1=i-1
          i2=i+1
          j1=j-1
          j2=j+1
          k1=k-1
          k2=k+1
          qx = (amq(i2,j,k)-amq(i1,j,k))/tdx
          qy = (amq(i,j2,k)-amq(i,j1,k))/tdy
          qz = (amq(i,j,k2)-amq(i,j,k1))/dzk

          CE = amu(i,j,k)*qx+amv(i,j,k)*qy+(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))*qz

          f = f + penal*CE*CE

! nv=1 for uu
          mgrid%grdt(i,j,k,t,1) =  mgrid%grdt(i,j,k,t,1)+2.*penal*CE*(qx - j13*qz)
! nv=2 for vv
          mgrid%grdt(i,j,k,t,2) =  mgrid%grdt(i,j,k,t,2)+2.*penal*CE*(qy - j23*qz)
! nv=3 for ww
          mgrid%grdt(i,j,k,t,3) =  mgrid%grdt(i,j,k,t,3)+2.*penal*CE*qz
! nv=6 for qv
          mgrid%grdt(i2,j,k,t,6) = mgrid%grdt(i2,j,k,t,6)+2.*penal*CE*amu(i,j,k)/tdx
          mgrid%grdt(i1,j,k,t,6) = mgrid%grdt(i1,j,k,t,6)-2.*penal*CE*amu(i,j,k)/tdx
          mgrid%grdt(i,j2,k,t,6) = mgrid%grdt(i,j2,k,t,6)+2.*penal*CE*amv(i,j,k)/tdy
          mgrid%grdt(i,j1,k,t,6) = mgrid%grdt(i,j1,k,t,6)-2.*penal*CE*amv(i,j,k)/tdy
          mgrid%grdt(i,j,k2,t,6) = mgrid%grdt(i,j,k2,t,6)+2.*penal*CE*(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))/dzk
          mgrid%grdt(i,j,k1,t,6) = mgrid%grdt(i,j,k1,t,6)-2.*penal*CE*(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))/dzk
        ENDDO
      ENDDO
    ENDDO
! cost function. HJ 6/27/2011
    func=func+f
  ENDDO
  write(*,'(A9,e13.6)')'func, qv:', func

END SUBROUTINE qv_costgrad
