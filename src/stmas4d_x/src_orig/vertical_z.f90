SUBROUTINE vertical_z(ngrid_in,grid_spacing_out,topo_out,jt1_out,jt2_out,jt3_out,zz_out,cor_out)
!
!doc=========================================================================
!doc The vertical grid is defined here. Figure out later whether this is 
!doc the best place to do so.  
!doc Hongli Jiang 6/15/2011
!doc=========================================================================

USE STMAS, ONLY: ZTOP

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ngrid_in(4)
  real, INTENT(OUT) :: grid_spacing_out(4), &
                       topo_out(ngrid_in(1),ngrid_in(2)), &
                       jt1_out(ngrid_in(1),ngrid_in(2),ngrid_in(3)), &
                       jt2_out(ngrid_in(1),ngrid_in(2),ngrid_in(3)), &
                       jt3_out(ngrid_in(1),ngrid_in(2)), &
                       zz_out(ngrid_in(3)),cor_out(ngrid_in(1),ngrid_in(2))
       
! local variables
  character*30, PARAMETER :: header = &
    'STMAS_LAPS_INPUT>vertical_z: '
  integer I,J,K,i1,i2,j1,j2
  real dzoz
!
  do j=1,ngrid_in(2)
    do i=1,ngrid_in(1)
      cor_out(i,j)=1.e-4
    enddo
  enddo
!
! Try to use z information from height.nl to define zz similar to sigma_ht. 
! HJ 8/11/2011
  do k=1,ngrid_in(3)
    zz_out(k)=real(k-1)*ZTOP/real(ngrid_in(3)-1)
  enddo
! coriolis parameter. set it to be a constant for now. HJ 6/20/2011

  do k=1,ngrid_in(3)
    do j=1,ngrid_in(2)
      do i=1,ngrid_in(1)
        i2=min0(i+1,ngrid_in(1))
        i1=max0(i-1,1)
        j2=min0(j+1,ngrid_in(2))
        j1=max0(j-1,1)
        DZOZ=(zz_out(k)-ZTOP)/ZTOP
        jt1_out(i,j,k)=DZOZ*(topo_out(i2,j)-topo_out(i1,j))/(2.*grid_spacing_out(1))
        jt2_out(i,j,k)=DZOZ*(topo_out(i,j2)-topo_out(i,j1))/(2.*grid_spacing_out(2))
        jt3_out(i,j)=(ZTOP-topo_out(i,j))/ZTOP
      enddo
    enddo
  enddo

  RETURN
  END
