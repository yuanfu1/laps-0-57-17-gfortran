SUBROUTINE vinterp_zz(ngrid,zz,bkgnd,nv,ter,jt3,bkgd2d)
!
!doc=========================================================================
!doc call interpolate_3dfield to vertically interpret all background variables  
!doc from LAPS to sigma_ht (defined in vertical_z) levels. This function 
!doc may not be needed once the analysis from sigma_ht is working.
!doc 
!doc Hongli Jiang 6/22/2011
!doc=========================================================================

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ngrid(4),nv
  real, INTENT(OUT) :: zz(ngrid(3)), &
                       ter(ngrid(1),ngrid(2)), &
                       jt3(ngrid(1),ngrid(2)), &
                       bkgnd(ngrid(1),ngrid(2),ngrid(3),ngrid(4),nv),&
                       bkgd2d(ngrid(1),ngrid(2),ngrid(4))
       
! local variables
  integer I,J,K,t,iv,i1,i2,j1,j2
  character*30, PARAMETER :: header = 'STMAS_LAPS_INPUT>vinterp_zz: '
! define temporarily arrays for storing variables. 
  real :: uu(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: vv(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: ww(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: tt(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: qv(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: pp(ngrid(1),ngrid(2),ngrid(3),ngrid(4))
  real :: h2(ngrid(1),ngrid(2),ngrid(3))
!
! call h2_call first since zz is actually eta, h2 and eta are related
! by eta=ztop*(h2-ter_in)/(ztop-ter_in). HJ 6/22/2011
!
  call h2_cal(ngrid(1),ngrid(2),ngrid(3),ter,zz,jt3,h2)
!
! call vtint_4dht
! u
  call vtint_4dht(bkgnd(1,1,1,1,1),ngrid(1),ngrid(2),ngrid(3),ngrid(4), &
                  bkgnd(1,1,1,1,4),ngrid(3),uu,h2)
! v
  call vtint_4dht(bkgnd(1,1,1,1,2),ngrid(1),ngrid(2),ngrid(3),ngrid(4), &
                  bkgnd(1,1,1,1,4),ngrid(3),vv,h2)
! w
  call vtint_4dht(bkgnd(1,1,1,1,3),ngrid(1),ngrid(2),ngrid(3),ngrid(4), &
                  bkgnd(1,1,1,1,4),ngrid(3),ww,h2)
! t
  call vtint_4dht(bkgnd(1,1,1,1,5),ngrid(1),ngrid(2),ngrid(3),ngrid(4), &
                  bkgnd(1,1,1,1,4),ngrid(3),tt,h2)
! qv
  call vtint_4dht(bkgnd(1,1,1,1,6),ngrid(1),ngrid(2),ngrid(3),ngrid(4), &
                  bkgnd(1,1,1,1,4),ngrid(3),qv,h2)
!
! calculate pressure based on hydrostatic balance. HJ 6/22/2011
  call pp_eta(ngrid(1),ngrid(2),ngrid(3),ngrid(4),bkgnd(1,1,1,1,5), &
                  bkgnd(1,1,1,1,6),zz,jt3,bkgd2d(1,1,1),pp)

! converting T to thetav.  HJ 6/22/2011
  call thetav_cal(ngrid(1),ngrid(2),ngrid(3),ngrid(4),pp, &
                  bkgnd(1,1,1,1,6),tt)
!
! assign interpolated field back to bkgnd. pp should be in pascal. 
!
  do t=1,ngrid(4)
    do k=1,ngrid(3)
      do j=1,ngrid(2)
        do i=1,ngrid(1)
          bkgnd(i,j,k,t,1)=uu(i,j,k,t)
          bkgnd(i,j,k,t,2)=vv(i,j,k,t)
          bkgnd(i,j,k,t,3)=ww(i,j,k,t)
          bkgnd(i,j,k,t,4)=100.*pp(i,j,k,t)
          bkgnd(i,j,k,t,5)=tt(i,j,k,t)
          bkgnd(i,j,k,t,6)=qv(i,j,k,t)
!           if(i.eq.1 .and. j.eq.1)then
!             write(*,'(A4,2I4,7e13.5)')'bk',t,k,bkgnd(i,j,k,t,1),bkgnd(i,j,k,t,2),&
!             bkgnd(i,j,k,t,3),bkgnd(i,j,k,t,4),bkgnd(i,j,k,t,5),bkgnd(i,j,k,t,6),zz(k)
!           endif
        enddo
      enddo
    enddo
  enddo

  RETURN
  END
!
!==========vtint_4dht
   subroutine vtint_4dht(v1,nx,ny,nz1,nt,h1,nz2,v2,h2)
   implicit none
   integer :: nz1,nz2,nt,nx,ny
   real :: v1(nx,ny,nz1,nt),v2(nx,ny,nz2,nt)
   real :: h1(nx,ny,nz1,nt),h2(nx,ny,nz2)

! when using bkgnd(1,1,1,1,4), it is 4d
   integer :: i,j,l,k,it,kk
   real :: wt

    do it=1,nt
    do j=1,ny
    do i=1,nx
      l=1
      do 20 k=1,nz2
 30   continue
      if(h2(i,j,k).lt.h1(i,j,1,it))go to 35
      if(h2(i,j,k).ge.h1(i,j,l,it).and.h2(i,j,k).le.h1(i,j,l+1,it))go to 35
      if(h2(i,j,k).gt.h1(i,j,nz1,it))go to 36
      l=l+1
      if(l.eq.nz1) then
        print *,'vtint_3dht:nz1',nz1
        stop 'htint'
      endif
      go to 30
 35   continue
      wt=(h2(i,j,k)-h1(i,j,l,it))/(h1(i,j,l+1,it)-h1(i,j,l,it))
      v2(i,j,k,it)=v1(i,j,l,it)+(v1(i,j,l+1,it)-v1(i,j,l,it))*wt
      go to 20
 36   continue
      wt=(h2(i,j,k)-h1(i,j,nz1,it))/(h1(i,j,nz1-1,it)-h1(i,j,nz1,it))
      v2(i,j,k,it)=v1(i,j,nz1,it)+(v1(i,j,nz1-1,it)-v1(i,j,nz1,it))*wt
 20   continue
     enddo !i-loop
     enddo !j-loop
     enddo !it-loop

     return
     end
!
!=================
! h2=jt3*eta+ter by following the eta defination. HJ 6/22/2011
   subroutine h2_cal(nx,ny,nz,ter,zfine,jt3,h2)

   integer i,j,k,nx,ny,nz
   real ter(nx,ny),zfine(nz),jt3(nx,ny),h2(nx,ny,nz)
     

   call azero(nx*ny*nz,h2)
   do k=1,nz
     do j=1,ny
       do i=1,nx
         h2(i,j,k)=jt3(i,j)*zfine(k)+ter(i,j)
       enddo
     enddo
   enddo
   return
   end
!
!=========thetav_cal
! thetav=Tv*(P00/p)**kappa, tv=t*(1.+0.61*rv*1.e-3) if rv in g/kg
   subroutine thetav_cal(nx,ny,nz,nt,pres,rv,temp)

use STMAS, ONLY: kappa,P00

   integer i,j,k,t,nx,ny,nz,nt
   real tv,pres(nx,ny,nz,nt),temp(nx,ny,nz,nt),rv(nx,ny,nz,nt)

   do t=1,nt
     do k=1,nz
       do j=1,ny
         do i=1,nx
           tv=temp(i,j,k,t)*(1.+0.61*rv(i,j,k,t)*1.e-3)
           temp(i,j,k,t)=tv*(P00/(100.*pres(i,j,k,t)))**kappa
         enddo
       enddo
     enddo
   enddo
   return
   end
!
!=========pp_eta
! compute pressure from height levels based on hydrostatic eq. HJ 6/22/2011
! dp/dz=-rho*g=(1/j3)dp/deta, P=rho*Rd*T --> rho=P/Rd*T --> dlog(P)=-(j3*g/(Rd*T))*deta
      subroutine pp_eta(nx,ny,nz,nt,temp,rv,zfine,jt3,psf,pres)

use STMAS, ONLY: G, gascnt

      integer i,j,k,t,nx,ny,nz,nt
      real tavg,pres(nx,ny,nz,nt),temp(nx,ny,nz,nt),zfine(nz)
      real rv(nx,ny,nz,nt),jt3(nx,ny),psf(nx,ny,nt)
      double precision dz

      call azero(nx*ny*nz*nt,pres)
      call a2b(nx,ny,nz,nt,1,pres,psf)
      do t=1,nt
      do k=2,nz
        do j=1,ny
          do i=1,nx
            tavg=0.5*(temp(i,j,k  ,t)*(1.+0.61*rv(i,j,k,t)*1.e-3)+ &
                      temp(i,j,k-1,t)*(1.+0.61*rv(i,j,k-1,t)*1.e-3))
            dz=(zfine(k)-zfine(k-1))
            pres(i,j,k,t)=pres(i,j,k-1,t) &
                        *exp(-jt3(i,j)*G*dz/(gascnt*tavg))
          enddo
        enddo
      enddo
      enddo
      return
      end
!
!===azero
    subroutine azero(N,varin)
!
    integer k,N
    real varin(N)

    do k=1,N
      varin(k) = 0.
    enddo
    return
    end
!
!===a2b
    subroutine a2b(nx,ny,nz,nt,k,varin,var2d)
!
    integer i,j,k,t,nx,ny,nz,nt
    real varin(nx,ny,nz,nt),var2d(nx,ny,nt)

    do t=1,nt
    do j=1,ny
      do i=1,nx
        varin(i,j,k,t) = var2d(i,j,t)
      enddo
    enddo
    enddo
    return
    end
