SUBROUTINE my_costgrad(func,num_controls,mgrid)

!doc=============================================================================================
!doc This routine computes the cost function and its gradient of y-momentum for STMAS4D
!doc
!doc History: May. 2011 by Hongli Jiang
!doc=============================================================================================

  USE STMAS, ONLY: STMAS_bkgd,rpk,kapp1,kapp2,STMAS_penal

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(OUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  INTEGER :: i,j,k,t,i1,i2,j1,j2,k1,k2,t1
  REAL    :: vx,vy,vz,py,pz,j13,j23,dzk,tdx,tdy,penal,CE
  REAL    :: amu(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amv(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amw(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amp(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  REAL    :: amt(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  DOUBLE PRECISION :: f                     ! Save function value reducing round off

  tdx=2.*mgrid%gridspc(1)
  tdy=2.*mgrid%gridspc(2)
  penal=STMAS_penal(2)

!  print*,'my', mgrid%gridspc(1),mgrid%gridspc(2)

  DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
      DO j=1,mgrid%numgrid(2)
        DO i=1,mgrid%numgrid(1)
          amu(i,j,k) = mgrid%anal(i,j,k,t,1)
          amv(i,j,k) = mgrid%anal(i,j,k,t,2)
          amw(i,j,k) = mgrid%anal(i,j,k,t,3)
          amp(i,j,k) = mgrid%anal(i,j,k,t,4)
          amt(i,j,k) = mgrid%anal(i,j,k,t,5)
        ENDDO
      ENDDO
    ENDDO

    f = 0.0d0
    CE = 0.0

    DO k=2,mgrid%numgrid(3)-1
      DO j=2,mgrid%numgrid(2)-1
        DO i=2,mgrid%numgrid(1)-1
          i1=i-1
          i2=i+1
          j1=j-1
          j2=j+1
          k1=k-1
          k2=k+1
          j13=mgrid%jt1(i,j,k)/mgrid%jt3(i,j)
          j23=mgrid%jt2(i,j,k)/mgrid%jt3(i,j)
          dzk=mgrid%zz(k+1)-mgrid%zz(k-1)
          vx = (amv(i2,j,k)-amv(i1,j,k))/tdx
          vy = (amv(i,j2,k)-amv(i,j1,k))/tdy
          vz = (amv(i,j,k2)-amv(i,j,k1))/dzk
          py = (amp(i,j2,k)-amp(i,j1,k))/tdy
          pz = (amp(i,j,k2)-amp(i,j,k1))/dzk

          CE = rpk*amt(i,j,k)*amp(i,j,k)**kapp1*(py+j23*pz) &
                    +amu(i,j,k)*vx+amv(i,j,k)*vy+(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))*vz + mgrid%cor(i,j)*amu(i,j,k)
          f = f + penal*CE*CE
!      if(i.eq.2 .and. j.eq.2 .and. k.eq.10)then
!        write(*,'(A3,2e13.5)')'my',f, CE
!      endif
! nv=1 for uu
          mgrid%grdt(i,j,k,t,1) = mgrid%grdt(i,j,k,t,1)+2.*penal*CE*(vx-j13*vz+mgrid%cor(i,j))
! nv=2 for vv
          mgrid%grdt(i, j,k,t,2) = mgrid%grdt(i, j,k,t,2)+2.*penal*CE*(vy-j23*vz)
          mgrid%grdt(i2,j,k,t,2) = mgrid%grdt(i2,j,k,t,2)+2.*penal*CE*amu(i,j,k)/tdx
          mgrid%grdt(i1,j,k,t,2) = mgrid%grdt(i1,j,k,t,2)-2.*penal*CE*amu(i,j,k)/tdx
          mgrid%grdt(i,j2,k,t,2) = mgrid%grdt(i,j2,k,t,2)+2.*penal*CE*amv(i,j,k)/tdy
          mgrid%grdt(i,j1,k,t,2) = mgrid%grdt(i,j1,k,t,2)-2.*penal*CE*amv(i,j,k)/tdy
          mgrid%grdt(i,j,k2,t,2) = mgrid%grdt(i,j,k2,t,2)+2.*penal*CE*(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))/dzk
          mgrid%grdt(i,j,k1,t,2) = mgrid%grdt(i,j,k2,t,2)-2.*penal*CE*(-j13*amu(i,j,k)-j23*amv(i,j,k)+amw(i,j,k))/dzk
! nv=3 for ww
          mgrid%grdt(i,j,k,t,3) = mgrid%grdt(i,j,k,t,3)+2.*penal*CE*vz
! nv=4 for pp
          mgrid%grdt(i, j,k,t,4) = mgrid%grdt(i, j,k,t,4)+2.*penal*CE*rpk*(amp(i,j,k)**kapp2)*kapp1*amt(i,j,k)*(py+j23*pz)
          mgrid%grdt(i,j2,k,t,4) = mgrid%grdt(i,j2,k,t,4)+2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)/tdy
          mgrid%grdt(i,j1,k,t,4) = mgrid%grdt(i,j1,k,t,4)-2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)/tdy
          mgrid%grdt(i,j,k2,t,4) = mgrid%grdt(i,j,k2,t,4)+2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)*j23/dzk
          mgrid%grdt(i,j,k1,t,4) = mgrid%grdt(i,j,k1,t,4)-2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*amt(i,j,k)*j23/dzk
! nv=5 for tt
          mgrid%grdt(i,j,k,t,5) = mgrid%grdt(i,j,k,t,5)+2.*penal*CE*rpk*(amp(i,j,k)**kapp1)*(py+j23*pz)
        ENDDO
      ENDDO
    ENDDO
! cost function. HJ 6/20/2011
    func=func+f
  ENDDO
  write(*,'(A9,e13.6)'),'func, my:', func

END SUBROUTINE my_costgrad
