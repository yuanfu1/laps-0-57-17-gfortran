SUBROUTINE mapping(nob_in,nvs_in,ngd_in,fgd_in,obs_in,err_in,pos_in, &
                   bkg_in,b2o_in,rd2_in,opt_in,frc_in,rdc_in,mob_in, &
                   lwo_in,vod_in,nob_out,loc_out,obs_out,err_out)
!doc
!doc===================================================================
!doc This routine maps an observation dateset from observation location
!doc to grid location using height, land-water, influence radius and
!doc background flow structure.
!doc
!doc Input:
!doc       nob_in: number of observations;
!doc       nvs_in: number of state variables;
!doc       ngd_in: numbers of gridpoints of the background;
!doc       fgd_in: numbers of finest gridpoint;
!doc       obs_in: observations;
!doc       err_in: observation errors;
!doc       pos_in: observation locations (1 -- numgrid);
!doc       bkg_in: background;
!doc       b2o_in: background values at obs locations;
!doc       rd2_in: influence radius square, x,y,z,t,flow,land-water;
!doc       opt_in: options - 0 no background flow dependent;
!doc                         1 with background flow dependent;
!doc       frc_in: fraction of land-water;
!doc       rdc_in: a vertical profile reducing land-water influence;
!doc       mob_in: the maximum obs defined for obs_in and b2o_in, 
!doc               and output obs allowed;
!doc       lwo_in: land-water factor at obs sites;
!doc       vod_in: void obs value;
!doc Output:
!doc       nob_out
!doc       obs_out: mapped observations at grid locations;
!doc       loc_out: gridded obs locations;
!doc
!doc History: Yuanfu Xie Nov. 2011
!doc===================================================================
!doc

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nob_in           ! Number of obs
  INTEGER, INTENT(IN) :: nvs_in           ! Number of state variables
  INTEGER, INTENT(IN) :: ngd_in(4)        ! Numbers of gridpoints
  INTEGER, INTENT(IN) :: fgd_in(4)        ! Numbers of finest gridpoints
  INTEGER, INTENT(IN) :: opt_in           ! Option to use background
  INTEGER, INTENT(IN) :: mob_in           ! Maximum obs out allowed
  REAL,    INTENT(IN) :: obs_in(mob_in,nvs_in)   ! Observed values
  REAL,    INTENT(IN) :: err_in(mob_in,nvs_in)   ! Observation errors
  REAL,    INTENT(IN) :: pos_in(4,nob_in) ! Observation positions in 
                                          ! finest grid scale
  REAL,    INTENT(IN) :: rdc_in(fgd_in(3),nvs_in)! land-water influence reduction
                                          ! as height at finest grid

  ! Background:
  REAL,    INTENT(IN) :: bkg_in(ngd_in(1),ngd_in(2),ngd_in(3),ngd_in(4),nvs_in)
  REAL,    INTENT(IN) :: b2o_in(mob_in,nvs_in)  ! Background values at obs loc
  ! Fraction of land-water:
  REAL,    INTENT(IN) :: frc_in(ngd_in(1),ngd_in(2)) 
  REAL,    INTENT(IN) :: rd2_in(6,nvs_in) ! Influence radius square
                                          ! (r*r) in grid scale 
                                          ! for x,y,z,t,flow and land-water
  REAL,    INTENT(IN) :: lwo_in(nob_in)   ! fraction of land-water at obs
  REAL,    INTENT(IN) :: vod_in           ! invalid obs value

  INTEGER, INTENT(OUT) :: nob_out
  INTEGER, INTENT(OUT) :: loc_out(4,fgd_in(1)*fgd_in(2)*fgd_in(3)*fgd_in(4))
  REAL,    INTENT(OUT) :: obs_out(fgd_in(1)*fgd_in(2)*fgd_in(3)*fgd_in(4),nvs_in)
  REAL,    INTENT(OUT) :: err_out(fgd_in(1)*fgd_in(2)*fgd_in(3)*fgd_in(4),nvs_in)

  ! Local variables:
  INTEGER :: io,i,j,k,l,iv                ! looping integers
  INTEGER :: ix,iy,iz,it,ih               ! indices
  REAL    :: ncount(4,4,4,4,nob_in,nvs_in)       ! count of number of obs
  REAL    :: wobsvs(4,4,4,4,nob_in,nvs_in)       ! weighted obs
  REAL    :: werror(4,4,4,4,nob_in,nvs_in)       ! weighted err
  ! REAL    :: gb,gw,gx,gy,gz,gt            ! Gaussian weightings
  ! Gaussian weights:
  REAL    :: gb(nvs_in),gw(nvs_in),gx(nvs_in),gy(nvs_in),gz(nvs_in),gt(nvs_in)

  ! Initialize:
  nob_out = 0
  ncount = 0
  wobsvs = 0.0
  werror = 0.0

  DO io=1,nob_in ! for each obs

    ! Indices of box corner:

    ! Use a scheme spreading obs to 4 adjacent gridpoints:
    DO l=0,3  ! Time
      it = ((pos_in(4,io)-1)/(fgd_in(4)-1))*(ngd_in(4)-1)+l
      ! Check gridpoint is in the analysis domain:
      IF (it .LT. 1 .OR. it .GT. ngd_in(4)) cycle

      ! dt:
      gt(1) = ((pos_in(4,io)-1)/(fgd_in(4)-1))*(ngd_in(4)-1)-it
      ! e**(-dt*dt/rr)
      gt = EXP(-gt(1)*gt(1)/rd2_in(4,1:nvs_in))

      DO k=0,3  ! Z
        ih = pos_in(3,io)  ! determine boundary layer for applying land-water

        iz = ((pos_in(3,io)-1)/(fgd_in(3)-1))*(ngd_in(3)-1)+k
        ! Check gridpoint is in the analysis domain:
        IF (iz .LT. 1 .OR. iz .GT. ngd_in(3)) cycle

        ! dz:
        gz(1) = ((pos_in(3,io)-1)/(fgd_in(3)-1))*(ngd_in(3)-1)-iz
        ! e**(-dz*dz/rr)
        gz = EXP(-gz(1)*gz(1)/rd2_in(3,1:nvs_in))

        DO j=0,3  ! Y
          iy = ((pos_in(2,io)-1)/(fgd_in(2)-1))*(ngd_in(2)-1)+j
          ! Check gridpoint is in the analysis domain:
          IF (iy .LT. 1 .OR. iy .GT. ngd_in(2)) cycle

          ! dy:
          gy(1) = ((pos_in(2,io)-1)/(fgd_in(2)-1))*(ngd_in(2)-1)-iy
          ! e**(-dy*dy/rr)
          gy = EXP(-gy(1)*gy(1)/rd2_in(2,1:nvs_in))

          DO i=0,3  ! X
            ix = ((pos_in(1,io)-1)/(fgd_in(1)-1))*(ngd_in(1)-1)+i
            ! Check gridpoint is in the analysis domain:
            IF (ix .LT. 1 .OR. ix .GT. ngd_in(1)) cycle

            ! dx:
            gx(1) = ((pos_in(1,io)-1)/(fgd_in(1)-1))*(ngd_in(1)-1)-ix
            ! e**(-dx*dx/rr)
            gx = EXP(-gx(1)*gx(1)/rd2_in(1,1:nvs_in))

            ! background:
            gb = 1.0
            IF (opt_in .EQ. 1) &
              gb = EXP(-(b2o_in(io,1:nvs_in)-bkg_in(ix,iy,iz,it,1:nvs_in))**2/rd2_in(5,1:nvs_in))

            ! land-water:
            gw = EXP(-(lwo_in(io)-frc_in(ix,iy))**2/rd2_in(6,1:nvs_in)*rdc_in(ih,1:nvs_in))

            ! Count each gridpoint for all variables:
            nob_out = nob_out+1

if (ix .eq. 2 .and. iy .eq. 57 .and. iz .eq. 1 .and. it .eq. 1 .and. nob_out .eq. 6) then
print*,'Found xie'
endif
            loc_out(1,nob_out) = ix
            loc_out(2,nob_out) = iy
            loc_out(3,nob_out) = iz
            loc_out(4,nob_out) = it

            ! For all variables counted at the gridpoint:
            DO iv=1,nvs_in

              ! Skip void obs:
              IF (obs_in(io,iv) .EQ. vod_in) CYCLE

              ncount(i+1,j+1,k+1,l+1,io,iv) = ncount(i+1,j+1,k+1,l+1,io,iv)+1

              ! Weighted obs:
              wobsvs(i+1,j+1,k+1,l+1,io,iv) = wobsvs(i+1,j+1,k+1,l+1,io,iv)+ &
                gb(iv)*gw(iv)*gx(iv)*gy(iv)*gz(iv)*gt(iv)*(obs_in(io,iv)-b2o_in(io,iv))

              ! Weighted err:
              werror(i+1,j+1,k+1,l+1,io,iv) = werror(i+1,j+1,k+1,l+1,io,iv)+ &
                gb(iv)*gw(iv)*gx(iv)*gy(iv)*gz(iv)*gt(iv)*(err_in(io,iv))

            ENDDO     ! for each variable

            ! No valid obs at this location:
            IF (SUM(ncount(1:4,1:4,1:4,1:4,io,1:nvs_in)) .EQ. 0) nob_out = nob_out-1

          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDDO          ! for each obs

  ! Expected observations at grid locations:
  nob_out = 0

  ! Testing for 1 obs
  DO io=1,1 !nob_in
    DO l=1,4
    DO k=1,4
    DO j=1,4
    DO i=1,4
  
      IF (SUM(ncount(i,j,k,l,io,1:nvs_in)) .GT. 0) THEN
        nob_out = nob_out+1
        obs_out(nob_out,1:nvs_in) = vod_in

        DO iv=1,nvs_in

          IF (ncount(i,j,k,l,io,iv) .GT. 0) THEN
            ! Expected obs value:

            obs_out(nob_out,iv) = bkg_in(loc_out(1,nob_out), &
                                         loc_out(2,nob_out), &
                                         loc_out(3,nob_out), &
                                         loc_out(4,nob_out),iv) + &
                             wobsvs(i,j,k,l,io,iv)/ncount(i,j,k,l,io,iv)

            ! Expected obs error:
            err_out(nob_out,iv) = werror(i,j,k,l,io,iv)/ncount(i,j,k,l,io,iv)
          ENDIF
        ENDDO
      ENDIF

    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

END SUBROUTINE mapping
