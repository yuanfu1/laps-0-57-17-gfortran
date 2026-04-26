SUBROUTINE grid1_interpolation(from_in,nfrom_in,to_out,nto_in)
!doc==============================================================================================
!doc  This routine interpolates a 1-D in k grid function to another.
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc
!doc  modified from grid5 to interpolate z only
!doc
!doc  7/7/2011 Hongli Jiang
!doc==============================================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nfrom_in(3),nto_in(3)
  REAL,    INTENT(IN) :: from_in(nfrom_in(3))
  REAL,    INTENT(OUT) :: to_out(nto_in(3))

  ! Local variables:
  INTEGER :: i,j,k,iv,idx(2,3),ii,jj,kk
  REAL    :: x,y,z,coe(2,3)

  ! Interpolation:
  DO k=1,nto_in(3)
    z = (k-1)*(nfrom_in(3)-1)/FLOAT(nto_in(3)-1)+1
    idx(1,3) = INT(z)  ! index for the third dimension
    idx(2,3) = MIN(nfrom_in(3),idx(1,3)+1)
    coe(2,3) = z-idx(1,3)
    coe(1,3) = 1.0-coe(2,3)
    ! Interpolates:
    to_out(k) = 0.0
    DO kk=1,2
     to_out(k) = to_out(k)+from_in(idx(kk,3))*coe(kk,3)
    ENDDO
  ENDDO

END SUBROUTINE grid1_interpolation

SUBROUTINE grid2_interpolation(from_in,nfrom_in,to_out,nto_in)
!doc==============================================================================================
!doc  This routine interpolates a 2-D grid function to another.
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc 
!doc  Removed nvar_in. 7/7/2011 Hongli Jiang
!doc==============================================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nfrom_in(2),nto_in(2)
  REAL,    INTENT(IN) :: from_in(nfrom_in(1),nfrom_in(2))
  REAL,    INTENT(OUT) :: to_out(nto_in(1),nto_in(2))

  ! Local variables:
  INTEGER :: i,j,idx(2,2),ii,jj
  REAL    :: x,y,coe(2,2)

  ! Interpolation:
  DO j=1,nto_in(2)
    y = (j-1)*(nfrom_in(2)-1)/FLOAT(nto_in(2)-1)+1 
    idx(1,2) = INT(y)  ! index for the second dimension
    idx(2,2) = MIN(nfrom_in(2),idx(1,2)+1)
    coe(2,2) = y-idx(1,2)
    coe(1,2) = 1.0-coe(2,2)
    DO i=1,nto_in(1)
      x = (i-1)*(nfrom_in(1)-1)/FLOAT(nto_in(1)-1)+1 
      idx(1,1) = INT(x)   ! index for the first dimension
      idx(2,1) = MIN(nfrom_in(1),idx(1,1)+1)
      coe(2,1) = x-idx(1,1)
      coe(1,1) = 1.0-coe(2,1)

      ! Interpolates:
      to_out(i,j) = 0.0
      DO jj=1,2
        DO ii=1,2
          to_out(i,j) = to_out(i,j)+ from_in(idx(ii,1),idx(jj,2))*coe(ii,1)*coe(jj,2)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE grid2_interpolation

!
SUBROUTINE grid4_interpolation(from_in,nfrom_in,to_out,nto_in,nvars_in)
!doc==============================================================================================
!doc  This routine interpolates a 3-D grid function to another.
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc
!doc  Modified from grid3_interpolation. 6/10/2011 Hongli Jiang.
!doc==============================================================================================

  IMPLICIT NONE

! HJ mod: 4 was 3, nfrom_in(4) and nto_in(4) are added. 6/10/2011
  INTEGER, INTENT(IN) :: nfrom_in(4),nto_in(4),nvars_in
  REAL,    INTENT(IN) :: from_in(nfrom_in(1),nfrom_in(2),nfrom_in(3),nfrom_in(4),nvars_in)
  REAL,    INTENT(OUT) :: to_out(nto_in(1),nto_in(2),nto_in(3),nto_in(4),nvars_in)

  ! Local variables:
! HJ mod: 4 was 3. 6/10/2011
  INTEGER :: i,j,k,it,iv,idx(2,4),ii,jj,kk,tt
  REAL    :: x,y,z,t,coe(2,4)

! HJ mod: the third dimension (old) should be t, and is 4th dimension now. By following how it was done
! in grid3_interpolation. Add t-loop. 6/10/2011
  ! Interpolation:
DO it=1,nto_in(4)
  t = (it-1)*(nfrom_in(4)-1)/FLOAT(nto_in(4)-1)+1
  idx(1,4) = INT(t)  ! index for the 4th dimension
  idx(2,4) = MIN(nfrom_in(4),idx(1,4)+1)
  coe(2,4) = t-idx(1,4)
  coe(1,4) = 1.0-coe(2,4)
  DO k=1,nto_in(3)
    z = (k-1)*(nfrom_in(3)-1)/FLOAT(nto_in(3)-1)+1
    idx(1,3) = INT(z)  ! index for the third dimension
    idx(2,3) = MIN(nfrom_in(3),idx(1,3)+1)
    coe(2,3) = z-idx(1,3)
    coe(1,3) = 1.0-coe(2,3)
    DO j=1,nto_in(2)
      y = (j-1)*(nfrom_in(2)-1)/FLOAT(nto_in(2)-1)+1 
      idx(1,2) = INT(y)  ! index for the second dimension
      idx(2,2) = MIN(nfrom_in(2),idx(1,2)+1)
      coe(2,2) = y-idx(1,2)
      coe(1,2) = 1.0-coe(2,2)
      DO i=1,nto_in(1)
        x = (i-1)*(nfrom_in(1)-1)/FLOAT(nto_in(1)-1)+1 
        idx(1,1) = INT(x)   ! index for the first dimension
        idx(2,1) = MIN(nfrom_in(1),idx(1,1)+1)
        coe(2,1) = x-idx(1,1)
        coe(1,1) = 1.0-coe(2,1)

        ! All variables:
        DO iv=1,nvars_in

          ! Interpolates:
          to_out(i,j,k,it,iv) = 0.0
          DO tt=1,2
           DO kk=1,2
            DO jj=1,2
              DO ii=1,2
                to_out(i,j,k,it,iv) = to_out(i,j,k,it,iv)+ &
                  from_in(idx(ii,1),idx(jj,2),idx(kk,3),idx(tt,4),iv)* &
                  coe(ii,1)*coe(jj,2)*coe(kk,3)*coe(tt,4)
              ENDDO
            ENDDO
           ENDDO
          ENDDO

        ENDDO
        if(it .eq. 1 .and. k .eq. 10 .and. i .eq. nto_in(1) .and. j.eq. nto_in(2))then
          print*,'grid4',i,j,to_out(i,j,k,it,5)
        endif
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE grid4_interpolation

SUBROUTINE grid5_interpolation(from_in,nfrom_in,to_out,nto_in)
!doc==============================================================================================
!doc  This routine interpolates a 3-D grid function to another.
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc
!doc  This subroutine is used to be grid3 for x,y,t,nvars. Modified to interpolate x,y,z only
!doc
!doc  6/30/2011 Hongli Jiang
!doc==============================================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nfrom_in(3),nto_in(3)
  REAL,    INTENT(IN) :: from_in(nfrom_in(1),nfrom_in(2),nfrom_in(3))
  REAL,    INTENT(OUT) :: to_out(nto_in(1),nto_in(2),nto_in(3))

  ! Local variables:
  INTEGER :: i,j,k,iv,idx(2,3),ii,jj,kk
  REAL    :: x,y,z,coe(2,3)

  ! Interpolation:
  DO k=1,nto_in(3)
    z = (k-1)*(nfrom_in(3)-1)/FLOAT(nto_in(3)-1)+1
    idx(1,3) = INT(z)  ! index for the third dimension
    idx(2,3) = MIN(nfrom_in(3),idx(1,3)+1)
    coe(2,3) = z-idx(1,3)
    coe(1,3) = 1.0-coe(2,3)
    DO j=1,nto_in(2)
      y = (j-1)*(nfrom_in(2)-1)/FLOAT(nto_in(2)-1)+1 
      idx(1,2) = INT(y)  ! index for the second dimension
      idx(2,2) = MIN(nfrom_in(2),idx(1,2)+1)
      coe(2,2) = y-idx(1,2)
      coe(1,2) = 1.0-coe(2,2)
      DO i=1,nto_in(1)
        x = (i-1)*(nfrom_in(1)-1)/FLOAT(nto_in(1)-1)+1 
        idx(1,1) = INT(x)   ! index for the first dimension
        idx(2,1) = MIN(nfrom_in(1),idx(1,1)+1)
        coe(2,1) = x-idx(1,1)
        coe(1,1) = 1.0-coe(2,1)

        ! Interpolates:
        to_out(i,j,k) = 0.0
        DO kk=1,2
          DO jj=1,2
            DO ii=1,2
              to_out(i,j,k) = to_out(i,j,k)+ &
              from_in(idx(ii,1),idx(jj,2),idx(kk,3))* &
              coe(ii,1)*coe(jj,2)*coe(kk,3)
            ENDDO
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE grid5_interpolation
!
SUBROUTINE convert_om2w(n1,n2,n3,n4,w,pres)
!doc==============================================================================================
!doc  converting omega to w by calling omega_to_w
!doc
!doc  8/5/2011 Hongli Jiang
!doc==============================================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n1,n2,n3,n4
  REAL,    INTENT(OUT) :: w(n1,n2,n3,n4),pres(n1,n2,n3,n4)

  ! Local variables:
  INTEGER :: i,j,k,t
  real :: omega_to_w

  do t=1,n4
    do k=1,n3
      do j=1,n2
        do i=1,n1
!          if(t.eq.1 .and. i.eq.1 .and. j.eq.1) print*,'bef',k,w(i,j,k,t)
         
          w(i,j,k,t) =omega_to_w(w(i,j,k,t),pres(i,j,k,t))

!          if(t.eq.1 .and. i.eq.1 .and. j.eq.1) print*,'aft',k,w(i,j,k,t)
        enddo
      enddo
    enddo
  enddo

END SUBROUTINE convert_om2w
