SUBROUTINE mz_costgrad(func,num_controls,mgrid)

!doc=============================================================================================
!doc This routine computes the cost function and its gradient of z-momentum for STMAS4D
!doc
!doc History: 6/24/2011 by Hongli Jiang
!doc=============================================================================================

  USE STMAS, ONLY: STMAS_bkgd,rpk,kapp1,kapp2,G,STMAS_penal

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(INOUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  INTEGER :: i,j,k,t,i1,i2,j1,j2,k1,k2,t1
  REAL    :: ux,uy,uz,vx,vy,vz,wx,wy,wz,pz,j13,j23,dzk,tdx,tdy,penal,CE
  REAL    :: amu(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: j1u(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amv(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: j2v(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amw(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: j3w(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amp(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amt(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt1(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt2(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: jt3(mgrid%numgrid(1),mgrid%numgrid(2))
  DOUBLE PRECISION :: f                     ! Save function value reducing round off

  tdx=2.*mgrid%gridspc(1)
  tdy=2.*mgrid%gridspc(2)
  penal=STMAS_penal(3)

  DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
      DO j=1,mgrid%numgrid(2)
        DO i=1,mgrid%numgrid(1)
          amu(i,j,k) = mgrid%anal(i,j,k,t,1)
          j1u(i,j,k) = mgrid%anal(i,j,k,t,1)*mgrid%jt1(i,j,k)
          amv(i,j,k) = mgrid%anal(i,j,k,t,2)
          j2v(i,j,k) = mgrid%anal(i,j,k,t,2)*mgrid%jt2(i,j,k)
          amw(i,j,k) = mgrid%anal(i,j,k,t,3)
          j3w(i,j,k) = mgrid%anal(i,j,k,t,3)*mgrid%jt3(i,j)
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
          ux = (j1u(i2,j,k)-j1u(i1,j,k))/tdx
          uy = (j1u(i,j2,k)-j1u(i,j1,k))/tdy
          uz = (j1u(i,j,k2)-j1u(i,j,k1))/dzk
          vx = (j2v(i2,j,k)-j2v(i1,j,k))/tdx
          vy = (j2v(i,j2,k)-j2v(i,j1,k))/tdy
          vz = (j2v(i,j,k2)-j2v(i,j,k1))/dzk
          wx = (j3w(i2,j,k)-j3w(i1,j,k))/tdx
          wy = (j3w(i,j2,k)-j3w(i,j1,k))/tdy
          wz = (j3w(i,j,k2)-j3w(i,j,k1))/dzk
          pz = (amp(i,j,k2)-amp(i,j,k1))/dzk

          CE = rpk*amt(i,j,k)*amp(i,j,k)**kapp1*pz/jt3(i,j)+ &
                     (amu(i,j,k)*wx+amv(i,j,k)*wy+amw(i,j,k)*wz)- &
                     (amu(i,j,k)*ux+amv(i,j,k)*uy+amw(i,j,k)*uz)- &
                     (amu(i,j,k)*vx+amv(i,j,k)*vy+amw(i,j,k)*vz)- &
                     (j13*amu(i,j,k)+j23*amv(i,j,k))*(-uz-vz+wz)+G

          f = f + penal*CE*CE
!      if(i.eq.2 .and. j.eq.2 .and. k.eq.10)then
!      write(*,'(A3,2e13.5)')'mz',f, CE
!      endif
! nv=1 for uu
          mgrid%grdt(i,j,k,t,1) =  mgrid%grdt(i,j,k,t,1) +2.*penal*CE*(wx-ux-vx-j13*(-uz-vz+wz))
          mgrid%grdt(i2,j,k,t,1) = mgrid%grdt(i2,j,k,t,1)-2.*penal*CE*amu(i,j,k)*jt1(i2,j,k)/tdx
          mgrid%grdt(i1,j,k,t,1) = mgrid%grdt(i1,j,k,t,1)+2.*penal*CE*amu(i,j,k)*jt1(i1,j,k)/tdx
          mgrid%grdt(i,j2,k,t,1) = mgrid%grdt(i,j2,k,t,1)-2.*penal*CE*amv(i,j,k)*jt1(i,j2,k)/tdy
          mgrid%grdt(i,j1,k,t,1) = mgrid%grdt(i,j1,k,t,1)+2.*penal*CE*amv(i,j,k)*jt1(i,j2,k)/tdy
          mgrid%grdt(i,j,k2,t,1) = mgrid%grdt(i,j,k2,t,1)-2.*penal*CE*amw(i,j,k)*jt1(i,j,k2)/dzk
          mgrid%grdt(i,j,k1,t,1) = mgrid%grdt(i,j,k1,t,1)+2.*penal*CE*amw(i,j,k)*jt1(i,j,k1)/dzk
! nv=2 for vv
          mgrid%grdt(i,j,k,t,2) =  mgrid%grdt(i,j,k,t,2) +2.*penal*CE*(wy-uy-vy-j23*(-uz-vz+wz))
          mgrid%grdt(i2,j,k,t,2) = mgrid%grdt(i2,j,k,t,2)-2.*penal*CE*amu(i,j,k)*jt2(i2,j,k)/tdx
          mgrid%grdt(i1,j,k,t,2) = mgrid%grdt(i1,j,k,t,2)+2.*penal*CE*amu(i,j,k)*jt2(i1,j,k)/tdx
          mgrid%grdt(i,j2,k,t,2) = mgrid%grdt(i,j2,k,t,2)-2.*penal*CE*amv(i,j,k)*jt2(i,j2,k)/tdy
          mgrid%grdt(i,j1,k,t,2) = mgrid%grdt(i,j1,k,t,2)+2.*penal*CE*amv(i,j,k)*jt2(i,j1,k)/tdy
          mgrid%grdt(i,j,k2,t,2) = mgrid%grdt(i,j,k2,t,2)-2.*penal*CE*amw(i,j,k)*jt2(i,j,k2)/dzk
          mgrid%grdt(i,j,k1,t,2) = mgrid%grdt(i,j,k1,t,2)+2.*penal*CE*amw(i,j,k)*jt2(i,j,k1)/dzk
! nv=3 for ww
          mgrid%grdt(i,j,k,t,3) = mgrid%grdt(i,j,k,t,3)+2.*penal*CE*(wz-uz-vz)
          mgrid%grdt(i2,j,k,t,3) = mgrid%grdt(i2,j,k,t,3)+2.*penal*CE*amu(i,j,k)*jt3(i2,j)/tdx
          mgrid%grdt(i1,j,k,t,3) = mgrid%grdt(i1,j,k,t,3)-2.*penal*CE*amu(i,j,k)*jt3(i1,j)/tdx
          mgrid%grdt(i,j2,k,t,3) = mgrid%grdt(i,j2,k,t,3)+2.*penal*CE*amv(i,j,k)*jt3(i,j2)/tdy
          mgrid%grdt(i,j1,k,t,3) = mgrid%grdt(i,j1,k,t,3)-2.*penal*CE*amv(i,j,k)*jt3(i,j1)/tdy
          mgrid%grdt(i,j,k2,t,3) = mgrid%grdt(i,j,k2,t,3)+2.*penal*CE*(amw(i,j,k)-j13*amu(i,j,k)-j23*amv(i,j,k))*jt3(i,j)/dzk
          mgrid%grdt(i,j,k1,t,3) = mgrid%grdt(i,j,k1,t,3)-2.*penal*CE*(amw(i,j,k)-j13*amu(i,j,k)-j23*amv(i,j,k))*jt3(i,j)/dzk
! nv=4 for pp
          mgrid%grdt(i,j,k,t,4)  = mgrid%grdt(i,j,k,t,4)+2.*penal*CE*rpk*(amp(i,j,k)**kapp2)*kapp1*amt(i,j,k)*pz/jt3(i,j)
          mgrid%grdt(i,j,k2,t,4) = mgrid%grdt(i,j,k2,t,4)+2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)/(jt3(i,j)*dzk)
          mgrid%grdt(i,j,k1,t,4) = mgrid%grdt(i,j,k1,t,4)-2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)/(jt3(i,j)*dzk)
! nv=5 for tt
          mgrid%grdt(i,j,k,t,5) = mgrid%grdt(i,j,k,t,5)+2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*pz/jt3(i,j)
        ENDDO
      ENDDO
    ENDDO
! cost function. HJ 6/20/2011
    func=func+f
  ENDDO
  write(*,'(A9,e13.6)')'func, mz:', func

END SUBROUTINE mz_costgrad
