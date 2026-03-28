SUBROUTINE cost_smooth(func,num_controls,mgrid)
!==============================================================================================
!  This routine evaluates the smoothing terms for the STMASFC cost function and its gradients.
!
!  History: May. 2011 by Yuanfu Xie at NOAA/ESRL/GSD.
!
!  Modified from STMASFC_smooth.f90 for 4d. 8/18/2011 By Hongli Jiang
!
!  Description:
!    This code adds a smoothing term of the increments to the cost function.
!==============================================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: num_controls
  REAL,    INTENT(OUT) :: func
  TYPE(STMAS_bkgd)   :: mgrid

  ! Local variables:
  INTEGER :: iv,i,j,k,t
  REAL    :: dx2,dy2,dx4,dy4,vxx,vyy,vzz,z1,z2,z3,az,bz,cz
  REAL    :: d(mgrid%numgrid(1),mgrid%numgrid(2),mgrid%numgrid(3))
  DOUBLE PRECISION :: f                           ! Save function value reducing round off

  ! Pass through all variables at interior points:
  dx2 = 1.0/(mgrid%gridspc(1)*mgrid%gridspc(1))
  dy2 = 1.0/(mgrid%gridspc(2)*mgrid%gridspc(2))
  dx4=dx2*dx2
  dy4=dy2*dy2

  DO t=1,mgrid%numgrid(4)
    DO iv=1,STMAS_numvars
      f = 0.0d0
      d = 0.0
      DO k=2,mgrid%numgrid(3)-1
        z1=mgrid%zz(k-1)
        z2=mgrid%zz(k)
        z3=mgrid%zz(k+1)
        vxx = 0.
        vyy = 0.
        vzz = 0.
        DO j=2,mgrid%numgrid(2)-1
          DO i=2,mgrid%numgrid(1)-1
            call g2orderit(z1,z2,z3,az,bz,cz)
            az=az*(z2-z1)*(z3-z2)
            bz=bz*(z2-z1)*(z3-z2)
            cz=cz*(z2-z1)*(z3-z2)
!      
            vxx = mgrid%anal(i+1,j,k,t,iv)-2.0*mgrid%anal(i,j,k,t,iv)+mgrid%anal(i-1,j,k,t,iv)
            vyy = mgrid%anal(i,j+1,k,t,iv)-2.0*mgrid%anal(i,j,k,t,iv)+mgrid%anal(i,j-1,k,t,iv)
            vzz = az*mgrid%anal(i,j,k-1,t,iv)+bz*mgrid%anal(i,j,k,t,iv)+cz*mgrid%anal(i,j,k+1,t,iv)

            d(i,j,k) = (vxx*dx2)**2+(vyy*dy2)**2+(vzz)**2
            f = f + d(i,j,k)

!            if(i.eq.2 .and. j.eq.2 .and. k.eq.15)then
!             write(*,'(A4,I4,8e13.5)')'sm',t,f,d(i,j,k),mgrid%anal(i,j,k-1,t,iv), &
!                         mgrid%anal(i,j,k,t,iv),mgrid%anal(i,j,k+1,t,iv),az,bz,cz
!            endif

! iv=1 for uu, for mgrig%grdt(i,j,k,t,iv) term, the factor is 2*2 (one is from the -2, and
! the other is from dJ^2)/du=2*J*dJ/du. HJ 8/24/2011
            mgrid%grdt(i,j,k,t,iv)  =mgrid%grdt(i,j,k,t,iv)  -2.*STMAS_smooth(iv)*(2*vxx*dx4+2*vyy*dy4-vzz*bz)
            mgrid%grdt(i-1,j,k,t,iv)=mgrid%grdt(i-1,j,k,t,iv)+2.*STMAS_smooth(iv)*vxx*dx4
            mgrid%grdt(i+1,j,k,t,iv)=mgrid%grdt(i+1,j,k,t,iv)+2.*STMAS_smooth(iv)*vxx*dx4
            mgrid%grdt(i,j-1,k,t,iv)=mgrid%grdt(i,j-1,k,t,iv)+2.*STMAS_smooth(iv)*vyy*dy4
            mgrid%grdt(i,j+1,k,t,iv)=mgrid%grdt(i,j+1,k,t,iv)+2.*STMAS_smooth(iv)*vyy*dy4
            mgrid%grdt(i,j,k-1,t,iv)=mgrid%grdt(i,j,k-1,t,iv)+2.*STMAS_smooth(iv)*vzz*az
            mgrid%grdt(i,j,k+1,t,iv)=mgrid%grdt(i,j,k+1,t,iv)+2.*STMAS_smooth(iv)*vzz*cz
          ENDDO
        ENDDO
      ENDDO
! Add smooth cost to the cost function:
    func = func+0.5*STMAS_smooth(iv)*f
  enddo ! iv-loop
  ENDDO !t-llop
  print*,'func:smooth ', func

END SUBROUTINE cost_smooth
!
!
SUBROUTINE G2ORDERIT(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL 2-ORDER DERIVATIVE OF INTERIOR (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
    REAL X,Y
    REAL Z1,Z2,Z3
    REAL A,B,C
    X=Z2-Z1
    Y=Z3-Z2
    A=2.0/(X*X+X*Y)
    B=-2.0/(X*Y)
    C=2.0/(X*Y+Y*Y)

END SUBROUTINE G2ORDERIT
