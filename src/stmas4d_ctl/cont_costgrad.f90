SUBROUTINE cont_costgrad(func,num_controls,mgrid)

!doc=============================================================================================
!doc This routine computes the cost function and its gradient of continuity eqn for STMAS4D
!doc
!doc History: 6/24/2011 by Hongli Jiang
!doc=============================================================================================

  USE STMAS, ONLY: STMAS_bkgd,kapp1,STMAS_penal

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(OUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  INTEGER :: i,j,k,t,i1,i2,j1,j2,k1,k2,t1
  REAL    :: ux,uz,vy,vz,wz,px,py,pz,tx,ty,tz,j1z,j2z,j13,j23,dzk,tdx,tdy,penal,CE
  REAL    :: amu(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amv(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amw(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amp(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amt(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt1(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt2(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt3(mgrid%numgrid(1),mgrid%numgrid(2))
  DOUBLE PRECISION :: f             ! Save function value reducing round off

  tdx=2.*mgrid%gridspc(1)
  tdy=2.*mgrid%gridspc(2)
  penal=STMAS_penal(4)

  DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
      DO j=1,mgrid%numgrid(2)
        DO i=1,mgrid%numgrid(1)
          amu(i,j,k) = mgrid%anal(i,j,k,t,1)
          amv(i,j,k) = mgrid%anal(i,j,k,t,2)
          amw(i,j,k) = mgrid%anal(i,j,k,t,3)
          amp(i,j,k) = mgrid%anal(i,j,k,t,4)
          amt(i,j,k) = mgrid%anal(i,j,k,t,5)
          jt1(i,j,k) = mgrid%jt1(i,j,k)
          jt2(i,j,k) = mgrid%jt2(i,j,k)
          jt3(i,j) = mgrid%jt3(i,j)
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
          ux = (amu(i2,j,k)-amu(i1,j,k))/tdx
          uz = (amu(i,j,k2)-amu(i,j,k1))/dzk
          vy = (amv(i,j2,k)-amv(i,j1,k))/tdy
          vz = (amv(i,j,k2)-amv(i,j,k1))/dzk
          wz = (amw(i,j,k2)-amw(i,j,k1))/dzk
          px = (amp(i2,j,k)-amp(i1,j,k))/tdx
          py = (amp(i,j2,k)-amp(i,j1,k))/tdy
          pz = (amp(i,j,k2)-amp(i,j,k1))/dzk
          tx = (amt(i2,j,k)-amt(i1,j,k))/tdx
          ty = (amt(i,j2,k)-amt(i,j1,k))/tdy
          tz = (amt(i,j,k2)-amt(i,j,k1))/dzk
          j1z =(jt1(i,j,k2)-jt1(i,j,k1))/dzk
          j2z =(jt2(i,j,k2)-jt2(i,j,k1))/dzk

          CE = -kapp1*(amu(i,j,k)*px+amv(i,j,k)*py+amw(i,j,k)*pz)/amp(i,j,k)- &
                            (amu(i,j,k)*tx+amv(i,j,k)*ty+amw(i,j,k)*tz)/amt(i,j,k)+ &
                             ux+vy+wz-j13*uz-j23*vz + &
                            (-amu(i,j,k)*j1z-amv(i,j,k)*j2z)/jt3(i,j)

          f = f + penal*CE*CE

! nv=1 for uu
          mgrid%grdt(i,j,k,t,1) =  mgrid%grdt(i,j,k,t,1) -2.*penal*CE*(kapp1*px/amp(i,j,k) + tx/amt(i,j,k) + j1z/jt3(i,j))
          mgrid%grdt(i2,j,k,t,1) = mgrid%grdt(i2,j,k,t,1)+2.*penal*CE/tdx
          mgrid%grdt(i1,j,k,t,1) = mgrid%grdt(i1,j,k,t,1)-2.*penal*CE/tdx
          mgrid%grdt(i,j,k2,t,1) = mgrid%grdt(i,j,k2,t,1)-2.*penal*CE*j13/dzk
          mgrid%grdt(i,j,k1,t,1) = mgrid%grdt(i,j,k1,t,1)+2.*penal*CE*j13/dzk
! nv=2 for vv
          mgrid%grdt(i,j,k,t,2) =  mgrid%grdt(i,j,k,t,2) -2.*penal*CE*(kapp1*py/amp(i,j,k) + ty/amt(i,j,k) + j2z/jt3(i,j))
          mgrid%grdt(i,j2,k,t,2) = mgrid%grdt(i,j2,k,t,2)+2.*penal*CE/tdy
          mgrid%grdt(i,j1,k,t,2) = mgrid%grdt(i,j1,k,t,2)-2.*penal*CE/tdy
          mgrid%grdt(i,j,k2,t,2) = mgrid%grdt(i,j,k2,t,2)-2.*penal*CE*j23/dzk
          mgrid%grdt(i,j,k1,t,2) = mgrid%grdt(i,j,k1,t,2)+2.*penal*CE*j23/dzk
! nv=3 for ww
          mgrid%grdt(i,j,k,t,3) =  mgrid%grdt(i,j,k,t,3) -2.*penal*CE*(kapp1*pz/amp(i,j,k) + tz/amt(i,j,k))
          mgrid%grdt(i,j,k2,t,3) = mgrid%grdt(i,j,k2,t,3)+2.*penal*CE/dzk
          mgrid%grdt(i,j,k1,t,3) = mgrid%grdt(i,j,k1,t,3)-2.*penal*CE/dzk
! nv=4 for pp
          mgrid%grdt(i,j,k,t,4) = mgrid%grdt(i,j,k,t,4)+2.*penal*CE*kapp1/(amp(i,j,k)**2)*(amu(i,j,k)*px+amv(i,j,k)*py+amw(i,j,k)*pz)
          mgrid%grdt(i2,j,k,t,4) = mgrid%grdt(i2,j,k,t,4)-2.*penal*CE*kapp1*amu(i,j,k)/(amp(i,j,k)*tdx)
          mgrid%grdt(i1,j,k,t,4) = mgrid%grdt(i1,j,k,t,4)+2.*penal*CE*kapp1*amu(i,j,k)/(amp(i,j,k)*tdx)
          mgrid%grdt(i,j2,k,t,4) = mgrid%grdt(i,j2,k,t,4)-2.*penal*CE*kapp1*amv(i,j,k)/(amp(i,j,k)*tdy)
          mgrid%grdt(i,j1,k,t,4) = mgrid%grdt(i,j1,k,t,4)+2.*penal*CE*kapp1*amv(i,j,k)/(amp(i,j,k)*tdy)
          mgrid%grdt(i,j,k2,t,4) = mgrid%grdt(i,j,k2,t,4)-2.*penal*CE*kapp1*amw(i,j,k)/(amp(i,j,k)*dzk)
          mgrid%grdt(i,j,k1,t,4) = mgrid%grdt(i,j,k1,t,4)+2.*penal*CE*kapp1*amw(i,j,k)/(amp(i,j,k)*dzk)
! nv=5 for tt
          mgrid%grdt(i,j,k,t,5) = mgrid%grdt(i,j,k,t,5)+2.*penal*CE*(amu(i,j,k)*tx+amv(i,j,k)*ty+amw(i,j,k)*tz)/(amt(i,j,k)**2)
          mgrid%grdt(i2,j,k,t,5) = mgrid%grdt(i2,j,k,t,5)-2.*penal*CE*amu(i,j,k)/(amt(i,j,k)*tdx)
          mgrid%grdt(i1,j,k,t,5) = mgrid%grdt(i1,j,k,t,5)+2.*penal*CE*amu(i,j,k)/(amt(i,j,k)*tdx)
          mgrid%grdt(i,j2,k,t,5) = mgrid%grdt(i,j2,k,t,5)-2.*penal*CE*amv(i,j,k)/(amt(i,j,k)*tdy)
          mgrid%grdt(i,j1,k,t,5) = mgrid%grdt(i,j1,k,t,5)+2.*penal*CE*amv(i,j,k)/(amt(i,j,k)*tdy)
          mgrid%grdt(i,j,k2,t,5) = mgrid%grdt(i,j,k2,t,5)-2.*penal*CE*amw(i,j,k)/(amt(i,j,k)*dzk)
          mgrid%grdt(i,j,k1,t,5) = mgrid%grdt(i,j,k1,t,5)+2.*penal*CE*amw(i,j,k)/(amt(i,j,k)*dzk)
        ENDDO
      ENDDO
    ENDDO
! cost function. HJ 6/20/2011
    func=func+f
  ENDDO
  write(*,'(A9,e13.6)')'func, ct:', func

END SUBROUTINE cont_costgrad
