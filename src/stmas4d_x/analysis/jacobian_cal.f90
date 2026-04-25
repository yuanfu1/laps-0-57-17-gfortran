SUBROUTINE jacobian_cal(ngrid,grid_spacing_out,topo_out,jt1_out,jt2_out,jt3_out, &
                       zz_out,cor_out,bkgnd,nv)
!
!doc=========================================================================
!doc The vertical grid is defined here. Figure out later whether this is 
!doc the best place to do so.  
!doc Hongli Jiang 6/15/2011
!doc=========================================================================

USE STMAS, ONLY: ZTOP

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ngrid(4),nv
  real, INTENT(OUT) :: grid_spacing_out(4), &
                       topo_out(ngrid(1),ngrid(2)), &
                       jt1_out(ngrid(1),ngrid(2),ngrid(3)), &
                       jt2_out(ngrid(1),ngrid(2),ngrid(3)), &
                       jt3_out(ngrid(1),ngrid(2)), &
                       bkgnd(ngrid(1),ngrid(2),ngrid(3),ngrid(4),nv),&
                       zz_out(ngrid(3)),cor_out(ngrid(1),ngrid(2))
       
! local variables
  character*30, PARAMETER :: header = 'STMAS_LAPS_INPUT>vertical_z: '
  integer I,J,K,t,i1,i2,j1,j2
  real dzoz
!
  do j=1,ngrid(2)
    do i=1,ngrid(1)
      cor_out(i,j)=1.e-4
    enddo
  enddo

! coriolis parameter. set it to be a constant for now. HJ 6/20/2011

  do k=1,ngrid(3)
    do j=1,ngrid(2)
      do i=1,ngrid(1)
        i2=min0(i+1,ngrid(1))
        i1=max0(i-1,1)
        j2=min0(j+1,ngrid(2))
        j1=max0(j-1,1)
        DZOZ=(zz_out(k)-ZTOP)/ZTOP
        jt1_out(i,j,k)=DZOZ*(topo_out(i2,j)-topo_out(i1,j))/(2.*grid_spacing_out(1))
        jt2_out(i,j,k)=DZOZ*(topo_out(i,j2)-topo_out(i,j1))/(2.*grid_spacing_out(2))
        jt3_out(i,j)=(ZTOP-topo_out(i,j))/ZTOP
      enddo
    enddo
  enddo
!
! converting T to thv.  HJ 6/22/2011
  call thv_cal(ngrid(1),ngrid(2),ngrid(3),ngrid(4),bkgnd(1,1,1,1,4), &
                  bkgnd(1,1,1,1,6),bkgnd(1,1,1,1,5))
!
  RETURN
  END
!
!=========thv_cal
! thv=Tv*(P00/p)**kappa, tv=t*(1.+0.61*rv*1.e-3) if rv in g/kg
   subroutine thv_cal(nx,ny,nz,nt,pres,rv,temp)

use STMAS, ONLY: kappa,P00

   integer i,j,k,t,nx,ny,nz,nt
   real tv,pres(nx,ny,nz,nt),temp(nx,ny,nz,nt),rv(nx,ny,nz,nt)

  PRINT*,'Debugging: in thv_cal: kappa/P00: ', kappa, P00
   do t=1,nt
     do k=1,nz
       do j=1,ny
         do i=1,nx
           tv=temp(i,j,k,t)*(1.+0.61*rv(i,j,k,t))
           temp(i,j,k,t)=tv*(P00/pres(i,j,k,t))**kappa
         enddo
       enddo
     enddo
   enddo
   return
   end
