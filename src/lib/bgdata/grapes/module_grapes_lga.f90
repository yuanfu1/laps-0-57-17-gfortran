MODULE grapes_lga

! A module containing necessary items to convert a GRAPES forecast
! file to LGA for the current LAPS domain

  USE map_utils
  USE wrf_netcdf
  USE time_utils
  USE constants
  USE horiz_interp
  IMPLICIT NONE

  PRIVATE
  integer                :: icentw, jcentw
  integer                :: icentl, jcentl
  REAL                   :: rmissingflag
  ! LAPS Pressure Levels 
  REAL, ALLOCATABLE      :: pr_laps(:)

  ! LGA variables
  REAL, ALLOCATABLE      :: ht(:,:,:)
  REAL, ALLOCATABLE      :: t3(:,:,:)
  REAL, ALLOCATABLE      :: sh(:,:,:)
  REAL, ALLOCATABLE      :: u3(:,:,:)
  REAL, ALLOCATABLE      :: v3(:,:,:)
  REAL, ALLOCATABLE      :: om(:,:,:)
  ! LGB Variables
  REAL, ALLOCATABLE      :: usf(:,:)
  REAL, ALLOCATABLE      :: vsf(:,:)
  REAL, ALLOCATABLE      :: tsf(:,:)
  REAL, ALLOCATABLE      :: tsk(:,:) ! surface skin temp. (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: dsf(:,:)
  REAL, ALLOCATABLE      :: slp(:,:)
  REAL, ALLOCATABLE      :: psf(:,:)
  REAL, ALLOCATABLE      :: rsf(:,:)
  REAL, ALLOCATABLE      :: p(:,:)
  REAL, ALLOCATABLE      :: pcp(:,:) ! RAINNC+RAINC (added by Wei-Ting 130312)
  ! LAPS static variables
  REAL, ALLOCATABLE      :: topo_laps(:,:)
  REAL, ALLOCATABLE      :: lat(:,:)
  REAL, ALLOCATABLE      :: lon(:,:)
  INTEGER                :: nxl, nyl, nzl
  CHARACTER(LEN=200)     :: laps_data_root
  CHARACTER(LEN=10)      :: laps_domain_name  
  REAL                   :: redp_lvl
  ! WRF on pressure levels
  REAL, ALLOCATABLE      :: ht_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: t3_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: sh_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: u3_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: v3_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: om_grapesp(:,:,:)
  REAL, ALLOCATABLE      :: pr_grapesp(:,:,:) ! 3D pressure on sigma-height levels (SIGMA_HT)
  REAL, ALLOCATABLE      :: prgd(:,:,:)        ! 3D pressure on LAPS grid (SIGMA_HT)
  REAL, ALLOCATABLE      :: usf_grapes(:,:)
  REAL, ALLOCATABLE      :: vsf_grapes(:,:)
  REAL, ALLOCATABLE      :: tsf_grapes(:,:)
  REAL, ALLOCATABLE      :: tsk_grapes(:,:) ! surface skin temp. (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: dsf_grapes(:,:)
  REAL, ALLOCATABLE      :: slp_grapes(:,:)
  REAL, ALLOCATABLE      :: psf_grapes(:,:)
  REAL, ALLOCATABLE      :: rsf_grapes(:,:)
  REAL, ALLOCATABLE      :: p_grapes(:,:)
  REAL, ALLOCATABLE      :: pcp_grapes(:,:) ! RAINNC+RAINC (added by Wei-Ting 130312)
  REAL, ALLOCATABLE      :: tvb_grapes(:,:)  ! Mean virtual temperature in lowest 60mb
  ! WRF on native variables
  REAL, ALLOCATABLE      :: pr_grapess(:,:,:)
  REAL, ALLOCATABLE      :: ht_grapess(:,:,:)
  REAL, ALLOCATABLE      :: t3_grapess(:,:,:)
  REAL, ALLOCATABLE      :: sh_grapess(:,:,:)
  REAL, ALLOCATABLE      :: u3_grapess(:,:,:)
  REAL, ALLOCATABLE      :: v3_grapess(:,:,:)
  REAL, ALLOCATABLE      :: om_grapess(:,:,:)
  REAL, ALLOCATABLE      :: rho_grapess(:,:,:) ! Density
  REAL, ALLOCATABLE      :: mr_grapess(:,:,:) ! Mixing Ratio
  ! WRF static variables 
  INTEGER                :: cdf,cdp ! added cdp by Wei-Ting (130312)
  TYPE(proj_info)        :: grapesgrid
  REAL, ALLOCATABLE      :: topo_grapes(:,:)
  CHARACTER(LEN=19)      :: reftime
  INTEGER                :: tau_hr, tau_min,tau_sec
  INTEGER                :: itimestep,projcode
  INTEGER                :: nxw,nyw,nzw
  REAL                   :: dx_m, dy_m,dt
  REAL                   :: lat1_grapes, lon1_grapes
  REAL                   :: truelat1_grapes, truelat2_grapes, stdlon_grapes

  PUBLIC grapes2lga_griddata, vinterp_grapes2p_standalone, hinterp_grapes2lga
CONTAINS

  SUBROUTINE grapes2lga_griddata(filenames, files_i4time, grapes_files, grid, cmodel, i4time_now, istatus)
    ! Convert GRAPES forecast data (already read into grid) to LAPS LGA/LGB format.
    ! Mirrors the logic of wrf2lga in module_wrf_lga.f90.
    USE module_domain
    IMPLICIT NONE

    ! Input parameters
    CHARACTER(LEN=256), INTENT(IN) :: filenames(:)
    INTEGER, INTENT(IN)            :: files_i4time(2)
    INTEGER, INTENT(IN)            :: grapes_files
    TYPE(domain_t), INTENT(IN)     :: grid(:)
    CHARACTER(LEN=12), INTENT(IN)  :: cmodel
    INTEGER, INTENT(IN)            :: i4time_now
    INTEGER, INTENT(OUT)           :: istatus

    ! Local variables
    INTEGER :: i, j, k, ix, jy, kz
    INTEGER :: ids, ide, kds, kde, jds, jde
    INTEGER :: i4reftime, bg_valid, idx1, idx2
    REAL    :: rh, time_alpha, time_beta
    REAL    :: tvbar, tvbar_nlevs
    REAL, PARAMETER :: tvbar_thick = 6000.  ! Pa, depth for mean Tv computation
    REAL, EXTERNAL  :: relhum, dewpt
    CHARACTER(LEN=40) :: v_g

    istatus = 1

    ! ------------------------------------------------------------------
    ! Step 1: LAPS domain setup
    ! ------------------------------------------------------------------
    CALL find_domain_name(laps_data_root, laps_domain_name, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Could not get LAPS domain info"
      RETURN
    ENDIF
    CALL get_grid_dim_xy(nxl, nyl, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Could not get LAPS xy dims"
      RETURN
    ENDIF
    CALL get_laps_dimensions(nzl, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Could not get LAPS z dim"
      RETURN
    ENDIF
    PRINT *, "LAPS dims: nxl,nyl,nzl = ", nxl, nyl, nzl
    CALL get_r_missing_data(rmissingflag, istatus)
    CALL get_laps_redp(redp_lvl, istatus)

    IF (.NOT. ALLOCATED(pr_laps))   ALLOCATE(pr_laps(nzl))
    IF (.NOT. ALLOCATED(lat))       ALLOCATE(lat(nxl, nyl))
    IF (.NOT. ALLOCATED(lon))       ALLOCATE(lon(nxl, nyl))
    IF (.NOT. ALLOCATED(topo_laps)) ALLOCATE(topo_laps(nxl, nyl))

    CALL get_laps_domain(nxl, nyl, laps_domain_name, lat, lon, topo_laps, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Error reading LAPS static info."
      RETURN
    ENDIF
    CALL get_vertical_grid(v_g, istatus)
    CALL upcase(v_g, v_g)
    PRINT *, "grapes2lga_griddata: VERTICAL_GRID = ", TRIM(v_g)
    IF (TRIM(v_g) .EQ. 'SIGMA_HT') THEN
      CALL get_ht_1d(nzl, pr_laps, istatus)
      IF (istatus .NE. 1) THEN
        PRINT *, "grapes2lga_griddata: Error reading LAPS sigma-height levels"
        RETURN
      ENDIF
    ELSE
      CALL get_pres_1d(i4time_now, nzl, pr_laps, istatus)
      IF (istatus .NE. 1) THEN
        PRINT *, "grapes2lga_griddata: Error reading LAPS pressure levels"
        RETURN
      ENDIF
    END IF

    ! ------------------------------------------------------------------
    ! Step 2: Extract GRAPES grid dimensions and set up projection
    ! NOTE: GRAPES arrays are indexed (ids:ide, kds:kde, jds:jde) = (x,z,y).
    !       Interior vertical levels run from kds+1 to kde-1.
    ! ------------------------------------------------------------------
    ids = grid(1)%ids
    ide = grid(1)%ide
    kds = grid(1)%kds
    kde = grid(1)%kde
    jds = grid(1)%jds
    jde = grid(1)%jde

    nxw = ide - ids + 1
    nyw = jde - jds + 1
    nzw = kde - kds - 1   ! interior levels only

    PRINT *, "GRAPES grid dimensions: nxw,nyw,nzw = ", nxw, nyw, nzw

    ! GRAPES uses a regular lat/lon grid (PROJ_LATLON = 0).
    ! stdlon stores the longitude increment; truelat1 stores the latitude increment.
    CALL map_set(PROJ_LATLON, REAL(grid(1)%config%ys_sn), REAL(grid(1)%config%xs_we), &
                 0.0, REAL(grid(1)%config%xd), REAL(grid(1)%config%yd), 0.0, &
                 nxw, nyw, grapesgrid)

    icentw = nxw / 2
    jcentw = nyw / 2
    icentl = nxl / 2
    jcentl = nyl / 2

    ! ------------------------------------------------------------------
    ! Step 3: Time bookkeeping
    ! For a GRAPES analysis file the reference and valid times are the same.
    ! ------------------------------------------------------------------
    IF (i4time_now .EQ. files_i4time(1)) THEN
      time_alpha = 1.0D0
      time_beta  = 0.0D0
      idx1 = 1; idx2 = 1
    ELSE IF (i4time_now .EQ. files_i4time(2)) THEN
      time_alpha = 0.0D0
      time_beta  = 1.0D0
      idx1 = 2; idx2 = 2
    ELSE
      IF (i4time_now .LT. files_i4time(1) .AND. i4time_now .GT. files_i4time(2)) THEN
        time_alpha = DBLE(files_i4time(1)-i4time_now)/DBLE(files_i4time(1)-files_i4time(2))
        time_beta  = DBLE(i4time_now-files_i4time(2))/DBLE(files_i4time(1)-files_i4time(2))
        idx1 = 1; idx2 = 2
      ELSE
        PRINT*,'i4time_now in module_grapes_lga.f90 is not in the files windown! ', &
          ' files_i4time: ',files_i4time, ' i4time_now: ',i4time_now
        RETURN
      END IF
    END IF
    i4reftime = files_i4time(1)
    bg_valid  = i4time_now

    ! ------------------------------------------------------------------
    ! Step 4: Allocate GRAPES native (model-level) work arrays
    ! ------------------------------------------------------------------
    IF (.NOT. ALLOCATED(pr_grapess))  ALLOCATE(pr_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(ht_grapess))  ALLOCATE(ht_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(t3_grapess))  ALLOCATE(t3_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(sh_grapess))  ALLOCATE(sh_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(u3_grapess))  ALLOCATE(u3_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(v3_grapess))  ALLOCATE(v3_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(om_grapess))  ALLOCATE(om_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(mr_grapess))  ALLOCATE(mr_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(rho_grapess)) ALLOCATE(rho_grapess(nxw, nyw, nzw))
    IF (.NOT. ALLOCATED(topo_grapes)) ALLOCATE(topo_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(usf_grapes))  ALLOCATE(usf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(vsf_grapes))  ALLOCATE(vsf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(tsf_grapes))  ALLOCATE(tsf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(tsk_grapes))  ALLOCATE(tsk_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(rsf_grapes))  ALLOCATE(rsf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(dsf_grapes))  ALLOCATE(dsf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(slp_grapes))  ALLOCATE(slp_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(psf_grapes))  ALLOCATE(psf_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(p_grapes))    ALLOCATE(p_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(pcp_grapes))  ALLOCATE(pcp_grapes(nxw, nyw))
    IF (.NOT. ALLOCATED(tvb_grapes))  ALLOCATE(tvb_grapes(nxw, nyw))

    ! ------------------------------------------------------------------
    ! Step 5: Extract and convert GRAPES variables
    ! Physics:
    !   Pressure      : p = P0 * pi^(cp/R)   (CPOR = cp/R from constants)
    !   Temperature   : T = (th + thp) * pi   (th = base-state theta,
    !                                           thp = perturbation theta,
    !                                           pi  = total Exner function)
    !   Height        : zz is geometric height in metres
    !   Mixing ratio  : q  (mois_2)
    !   Specific hum. : sh = q / (1+q)
    !   Winds         : u, v already on mass points in m/s
    !   Omega         : set to 0 (wzet not read from file)
    ! ------------------------------------------------------------------
    PRINT *, "Extracting GRAPES model-level variables...", UBOUND(grid,1)
    DO j = 1, nyw
      jy = jds + (j - 1)
      DO i = 1, nxw
        ix = ids + (i - 1)
        DO k = 1, nzw
          kz = kds + k   ! GRAPES file is bottom-to-top: k=1->ground(kds+1), k=nzw->top(kds+nzw)

          pr_grapess(i,j,k) = P0 * ((time_alpha*grid(idx1)%pi(ix, kz, jy)+ &
                                     time_beta* grid(idx2)%pi(ix, kz, jy))**CPOR)
          t3_grapess(i,j,k) = time_alpha*(grid(idx1)%th(ix, kz, jy) + grid(idx1)%thp(ix, kz, jy)) * &
                                          grid(idx1)%pi(ix, kz, jy) + &
                              time_beta *(grid(idx2)%th(ix, kz, jy) + grid(idx2)%thp(ix, kz, jy)) * &
                                          grid(idx2)%pi(ix, kz, jy)
          ht_grapess(i,j,k) = time_alpha*grid(idx1)%zz(ix, kz, jy) + time_beta*grid(idx2)%zz(ix, kz, jy)
          mr_grapess(i,j,k) = time_alpha*MAX(0.0, grid(idx1)%q(ix, kz, jy)) + time_beta*MAX(0.0, grid(idx2)%q(ix, kz, jy))
          sh_grapess(i,j,k) = time_alpha*mr_grapess(i,j,k) / (1.0 + mr_grapess(i,j,k)) + &
                              time_beta *mr_grapess(i,j,k) / (1.0 + mr_grapess(i,j,k))
          u3_grapess(i,j,k) = time_alpha*grid(idx1)%u(ix, kz, jy) + time_beta*grid(idx2)%u(ix, kz, jy)
          v3_grapess(i,j,k) = time_alpha*grid(idx1)%v(ix, kz, jy) + time_beta*grid(idx2)%v(ix, kz, jy)
          om_grapess(i,j,k) = 0.0   ! vertical velocity not available

        ENDDO
      ENDDO
    ENDDO
    PRINT*,'Checking T3: ', MINVAL(t3_grapess), MAXVAL(t3_grapess)

    ! Surface fields (use lowest model level where 2m/10m fields unavailable)
    PRINT *, "Extracting GRAPES surface fields..."
    DO j = 1, nyw
      jy = jds + (j - 1)
      DO i = 1, nxw
        ix = ids + (i - 1)
        topo_grapes(i,j) = time_alpha*grid(idx1)%ht(ix, jy) + time_beta*grid(idx2)%ht(ix, jy)
        tsk_grapes(i,j)  = time_alpha*grid(idx1)%tsk(ix, jy) + time_beta*grid(idx2)%tsk(ix, jy)
        usf_grapes(i,j)  = u3_grapess(i,j,1)
        vsf_grapes(i,j)  = v3_grapess(i,j,1)
        tsf_grapes(i,j)  = t3_grapess(i,j,1)
        rsf_grapes(i,j)  = sh_grapess(i,j,1)
        psf_grapes(i,j)  = pr_grapess(i,j,1)
        pcp_grapes(i,j)  = 0.0   ! accumulated precipitation not available
      ENDDO
    ENDDO

    ! Mean virtual temperature in lowest tvbar_thick Pa (for extrapolation)
    DO j = 1, nyw
      DO i = 1, nxw
        tvbar       = 0.0
        tvbar_nlevs = 0.0
        DO k = 1, nzw
          IF ((psf_grapes(i,j) - pr_grapess(i,j,k)) .LE. tvbar_thick) THEN
            tvbar       = tvbar + t3_grapess(i,j,k) * (1.0 + 0.61*sh_grapess(i,j,k))
            tvbar_nlevs = tvbar_nlevs + 1.0
          ELSE
            EXIT
          ENDIF
        ENDDO
        IF (tvbar_nlevs .GT. 0.0) THEN
          tvb_grapes(i,j) = tvbar / tvbar_nlevs
        ELSE
          tvb_grapes(i,j) = t3_grapess(i,j,1) * (1.0 + 0.61*sh_grapess(i,j,1))
        ENDIF
      ENDDO
    ENDDO
    CALL smooth2(nxw, nyw, 4, tvb_grapes)

    ! Surface dewpoint
    DO j = 1, nyw
      DO i = 1, nxw
        rh = relhum(tsf_grapes(i,j), rsf_grapes(i,j), psf_grapes(i,j))
        dsf_grapes(i,j) = dewpt(tsf_grapes(i,j), rh)
      ENDDO
    ENDDO

    PRINT *, "Min/Max GRAPES 3D Pressure:", MINVAL(pr_grapess), MAXVAL(pr_grapess)
    PRINT *, "GRAPES model-level data at domain centre"
    PRINT *, "K   PRESS     HEIGHT   TEMP   SH       U       V      OM"
    PRINT *, "--- --------  -------  -----  -------  ------  ------ -----------"
    DO k = 1, nzw
      PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
        k, pr_grapess(icentw,jcentw,k), ht_grapess(icentw,jcentw,k), &
        t3_grapess(icentw,jcentw,k), sh_grapess(icentw,jcentw,k), &
        u3_grapess(icentw,jcentw,k), v3_grapess(icentw,jcentw,k), om_grapess(icentw,jcentw,k)
    ENDDO

    ! Debugging:
    WRITE(*,11) MINVAL(t3_grapess), MAXVAL(t3_grapess), MINVAL(v3_grapess), MAXVAL(v3_grapess)
11  FORMAT('Min t3_grapess: ', F10.3, ' Max t3_grapess: ', F10.3, &
           ' Min v3_grapess: ', F10.3, ' Max v3_grapess: ', F10.3)

    ! ------------------------------------------------------------------
    ! Step 6: Vertical interpolation to LAPS pressure levels
    ! ------------------------------------------------------------------
    PRINT *, "Allocating arrays for GRAPES on pressure levels"
    IF (.NOT. ALLOCATED(ht_grapesp)) ALLOCATE(ht_grapesp(nxw, nyw, nzl))
    IF (.NOT. ALLOCATED(t3_grapesp)) ALLOCATE(t3_grapesp(nxw, nyw, nzl))
    IF (.NOT. ALLOCATED(sh_grapesp)) ALLOCATE(sh_grapesp(nxw, nyw, nzl))
    IF (.NOT. ALLOCATED(u3_grapesp)) ALLOCATE(u3_grapesp(nxw, nyw, nzl))
    IF (.NOT. ALLOCATED(v3_grapesp)) ALLOCATE(v3_grapesp(nxw, nyw, nzl))
    IF (.NOT. ALLOCATED(om_grapesp)) ALLOCATE(om_grapesp(nxw, nyw, nzl))
    IF (TRIM(v_g) .EQ. 'SIGMA_HT') THEN
      IF (.NOT. ALLOCATED(pr_grapesp)) ALLOCATE(pr_grapesp(nxw, nyw, nzl))
      pr_grapesp = rmissingflag
    END IF
    ht_grapesp = rmissingflag
    t3_grapesp = rmissingflag
    sh_grapesp = rmissingflag
    u3_grapesp = rmissingflag
    v3_grapesp = rmissingflag
    om_grapesp = rmissingflag

    IF (TRIM(v_g) .EQ. 'SIGMA_HT') THEN
      PRINT *, "Calling vinterp_grapes2ht (height-based vertical interpolation)"
      CALL vinterp_grapes2ht(istatus)
    ELSE
      PRINT *, "Calling vinterp_grapes2p (pressure-based vertical interpolation)"
      CALL vinterp_grapes2p(istatus)
    END IF
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Problem vertically interpolating GRAPES data"
      RETURN
    ENDIF

    ! Debugging:
    WRITE(*,13) MINVAL(u3_grapesp), MAXVAL(u3_grapesp), MINVAL(v3_grapesp), &
      MAXVAL(v3_grapesp), MINVAL(t3_grapesp), MAXVAL(t3_grapesp), &
      MINVAL(ht_grapesp), MAXVAL(ht_grapesp)
13  FORMAT('Min u3_grapesp: ', F10.3, ' Max u3_grapesp: ', F10.3, &
           ' Min v3_grapesp: ', F10.3, ' Max v3_grapesp: ', F10.3, &
           ' Min t3_grapesp: ', F10.3, ' Max t3_grapesp: ', F10.3, &
           ' Min ht_grapesp: ', F10.3, ' Max ht_grapesp: ', F10.3)

    PRINT *, "Deallocating GRAPES model-level vars"
    IF (ALLOCATED(pr_grapess))  DEALLOCATE(pr_grapess)
    IF (ALLOCATED(ht_grapess))  DEALLOCATE(ht_grapess)
    IF (ALLOCATED(t3_grapess))  DEALLOCATE(t3_grapess)
    IF (ALLOCATED(u3_grapess))  DEALLOCATE(u3_grapess)
    IF (ALLOCATED(v3_grapess))  DEALLOCATE(v3_grapess)
    IF (ALLOCATED(om_grapess))  DEALLOCATE(om_grapess)
    IF (ALLOCATED(mr_grapess))  DEALLOCATE(mr_grapess)
    IF (ALLOCATED(rho_grapess)) DEALLOCATE(rho_grapess)

    ! ------------------------------------------------------------------
    ! Step 7: Allocate LAPS-grid (LGA/LGB) arrays and horizontal interpolation
    ! ------------------------------------------------------------------
    PRINT *, "Allocating LGA variables"
    IF (.NOT. ALLOCATED(ht)) ALLOCATE(ht(nxl, nyl, nzl))
    IF (.NOT. ALLOCATED(sh)) ALLOCATE(sh(nxl, nyl, nzl))
    IF (.NOT. ALLOCATED(t3)) ALLOCATE(t3(nxl, nyl, nzl))
    IF (.NOT. ALLOCATED(u3)) ALLOCATE(u3(nxl, nyl, nzl))
    IF (.NOT. ALLOCATED(v3)) ALLOCATE(v3(nxl, nyl, nzl))
    IF (.NOT. ALLOCATED(om)) ALLOCATE(om(nxl, nyl, nzl))
    PRINT *, "Allocating LGB variables"
    IF (.NOT. ALLOCATED(usf)) ALLOCATE(usf(nxl, nyl))
    IF (.NOT. ALLOCATED(vsf)) ALLOCATE(vsf(nxl, nyl))
    IF (.NOT. ALLOCATED(tsf)) ALLOCATE(tsf(nxl, nyl))
    IF (.NOT. ALLOCATED(tsk)) ALLOCATE(tsk(nxl, nyl))
    IF (.NOT. ALLOCATED(rsf)) ALLOCATE(rsf(nxl, nyl))
    IF (.NOT. ALLOCATED(dsf)) ALLOCATE(dsf(nxl, nyl))
    IF (.NOT. ALLOCATED(slp)) ALLOCATE(slp(nxl, nyl))
    IF (.NOT. ALLOCATED(psf)) ALLOCATE(psf(nxl, nyl))
    IF (.NOT. ALLOCATED(p))   ALLOCATE(p  (nxl, nyl))
    IF (.NOT. ALLOCATED(pcp)) ALLOCATE(pcp(nxl, nyl))
    ht  = rmissingflag;  sh  = rmissingflag;  t3  = rmissingflag
    u3  = rmissingflag;  v3  = rmissingflag;  om  = rmissingflag
    usf = rmissingflag;  vsf = rmissingflag;  tsf = rmissingflag
    tsk = rmissingflag;  rsf = rmissingflag;  dsf = rmissingflag
    slp = rmissingflag;  psf = rmissingflag;  p   = rmissingflag
    pcp = rmissingflag
    IF (TRIM(v_g) .EQ. 'SIGMA_HT') THEN
      IF (.NOT. ALLOCATED(prgd))  ALLOCATE(prgd(nxl, nyl, nzl))
      prgd = rmissingflag
    END IF

    PRINT *, "Calling hinterp_grapes2lga", MINVAL(t3_grapesp), MAXVAL(t3_grapesp)
    CALL hinterp_grapes2lga(istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Problem in hinterp_grapes2lga"
      istatus = 0
      RETURN
    ENDIF

    PRINT *, "Deallocating GRAPES pressure-level vars"
    DEALLOCATE(ht_grapesp, t3_grapesp, sh_grapesp, u3_grapesp, v3_grapesp, om_grapesp)

    PRINT *, "Deallocating GRAPES surface vars"
    IF (ALLOCATED(usf_grapes))  DEALLOCATE(usf_grapes)
    IF (ALLOCATED(vsf_grapes))  DEALLOCATE(vsf_grapes)
    IF (ALLOCATED(tsf_grapes))  DEALLOCATE(tsf_grapes)
    IF (ALLOCATED(tsk_grapes))  DEALLOCATE(tsk_grapes)
    IF (ALLOCATED(rsf_grapes))  DEALLOCATE(rsf_grapes)
    IF (ALLOCATED(dsf_grapes))  DEALLOCATE(dsf_grapes)
    IF (ALLOCATED(slp_grapes))  DEALLOCATE(slp_grapes)
    IF (ALLOCATED(psf_grapes))  DEALLOCATE(psf_grapes)
    IF (ALLOCATED(p_grapes))    DEALLOCATE(p_grapes)
    IF (ALLOCATED(pcp_grapes))  DEALLOCATE(pcp_grapes)
    IF (ALLOCATED(tvb_grapes))  DEALLOCATE(tvb_grapes)

    PRINT *, "Deallocating LAPS static vars"
    IF (ALLOCATED(lat))       DEALLOCATE(lat)
    IF (ALLOCATED(lon))       DEALLOCATE(lon)
    IF (ALLOCATED(topo_laps)) DEALLOCATE(topo_laps)
    IF (ALLOCATED(topo_grapes))  DEALLOCATE(topo_grapes)

    ! ------------------------------------------------------------------
    ! Step 8: Compute derived pressures (MSLP, reduced pressure)
    ! ------------------------------------------------------------------
    PRINT *, "Creating derived pressure arrays"
    CALL make_derived_pressures
    PRINT *, "slp/psf/p/tsf/tsk/dsf/rsf/usf/vsf:", &
      slp(icentl,jcentl), psf(icentl,jcentl), p(icentl,jcentl), &
      tsf(icentl,jcentl), tsk(icentl,jcentl), dsf(icentl,jcentl), &
      rsf(icentl,jcentl), usf(icentl,jcentl), vsf(icentl,jcentl)

    ! Debugging:
    WRITE(*,12) MINVAL(u3), MAXVAL(u3), MINVAL(v3), MAXVAL(v3), MINVAL(t3), MAXVAL(t3)
12  FORMAT('XIE Min u3: ', F10.3, ' Max u3: ', F10.3, &
           ' Min v3: ', F10.3, ' Max v3: ', F10.3, &
           ' Min t3: ', F10.3, ' Max t3: ', F10.3)

    ! ------------------------------------------------------------------
    ! Step 9: Write LGA and LGB output files
    ! ------------------------------------------------------------------
    PRINT *, "Writing LGA at time: ", i4time_now, bg_valid, ' File times: ', files_i4time(idx2), files_i4time(idx1)
    IF (TRIM(v_g) .EQ. 'SIGMA_HT') THEN
      ! pr_laps holds 1D sigma-height levels (m, from get_ht_1d)
      ! prgd holds 3D pressure field (Pa) on those levels
      CALL write_lgap(nxl, nyl, nzl, i4time_now, bg_valid, cmodel, rmissingflag, &
                     pr_laps, prgd, t3, sh, u3, v3, om, istatus)
    ELSE
      ! pr_laps holds 1D pressure levels (Pa); convert to mb for write_lga
      CALL write_lga(nxl, nyl, nzl, i4time_now, bg_valid, cmodel, rmissingflag, &
                    pr_laps*0.01, ht, t3, sh, u3, v3, om, istatus)
    END IF
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Error writing LGA"
      RETURN
    ENDIF

    PRINT *, "Writing LGB"
    CALL write_lgb(nxl, nyl, i4time_now, bg_valid, cmodel, rmissingflag, &
                   usf, vsf, tsf, tsk, rsf, psf, slp, dsf, p, pcp, istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, "grapes2lga_griddata: Error writing LGB"
      RETURN
    ENDIF

    ! ------------------------------------------------------------------
    ! Step 10: Final cleanup
    ! ------------------------------------------------------------------
    PRINT *, "Deallocating LGA vars"
    IF (ALLOCATED(ht))        DEALLOCATE(ht)
    IF (ALLOCATED(t3))        DEALLOCATE(t3)
    IF (ALLOCATED(sh))        DEALLOCATE(sh)
    IF (ALLOCATED(u3))        DEALLOCATE(u3)
    IF (ALLOCATED(v3))        DEALLOCATE(v3)
    IF (ALLOCATED(om))        DEALLOCATE(om)
    IF (ALLOCATED(pr_laps))   DEALLOCATE(pr_laps)
    IF (ALLOCATED(pr_grapesp)) DEALLOCATE(pr_grapesp)
    IF (ALLOCATED(prgd))      DEALLOCATE(prgd)

    PRINT *, "Deallocating LGB vars"
    IF (ALLOCATED(usf)) DEALLOCATE(usf)
    IF (ALLOCATED(vsf)) DEALLOCATE(vsf)
    IF (ALLOCATED(tsf)) DEALLOCATE(tsf)
    IF (ALLOCATED(tsk)) DEALLOCATE(tsk)
    IF (ALLOCATED(rsf)) DEALLOCATE(rsf)
    IF (ALLOCATED(dsf)) DEALLOCATE(dsf)
    IF (ALLOCATED(slp)) DEALLOCATE(slp)
    IF (ALLOCATED(psf)) DEALLOCATE(psf)
    IF (ALLOCATED(p))   DEALLOCATE(p)
    IF (ALLOCATED(pcp)) DEALLOCATE(pcp)

    PRINT *, "grapes2lga_griddata: Successful processing of GRAPES data"
    RETURN

  END SUBROUTINE grapes2lga_griddata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fill_wrfs(istatus)
 
    IMPLICIT NONE
    INTEGER, INTENT(OUT)  :: istatus
    INTEGER               :: status,i,j,k
    REAL, ALLOCATABLE     :: dum3d(:,:,:)
    REAL, ALLOCATABLE     :: dum3df(:,:,:)
    REAL, ALLOCATABLE     :: dum3df2(:,:,:)
    REAL, ALLOCATABLE     :: dum2dt1(:,:) ! added by Wei-Ting (130312) to get RAINNC
    REAL, ALLOCATABLE     :: dum2dt2(:,:) ! added by Wei-Ting (130312) to get RAINNC
    REAL, EXTERNAL        :: mixsat, relhum, dewpt
    REAL                  :: rh 
    REAL                  :: tvbar, tvbar_nlevs
    REAL, PARAMETER       :: tvbar_thick = 6000.
    ! Varialbles have already been allocated by our driver routine
    ! so just start getting them
    istatus = 1
    PRINT *, " Allocating arrays"
    IF (.NOT. ALLOCATED(dum3d)) ALLOCATE(dum3d(nxw,nyw,nzw))
    IF (.NOT. ALLOCATED(dum3df))  ALLOCATE(dum3df(nxw,nyw,nzw+1))
    IF (.NOT. ALLOCATED(dum3df2)) ALLOCATE(dum3df2(nxw,nyw,nzw+1))
    IF (.NOT. ALLOCATED(dum2dt1)) ALLOCATE(dum2dt1(nxw,nyw))
    IF (.NOT. ALLOCATED(dum2dt2)) ALLOCATE(dum2dt2(nxw,nyw))

    ! Get 3D pressure array
    ! Get pressures
    PRINT *, "Getting PB"
    CALL get_wrfnc_3d(cdf,"PB","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status .GT. 0) THEN
      PRINT *, "+++CRITICAL:  Could not get base pressure!"
      istatus = 0
      RETURN
    ENDIF
    pr_grapess = dum3d
                        
    PRINT *, "Getting P"                                                                    
    CALL get_wrfnc_3d(cdf,"P","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status .GT. 0) THEN
      PRINT *, "+++CRITICAL:  Could not get pert pressure!"
      istatus = 0
      RETURN
    ENDIF
    pr_grapess = pr_grapess + dum3d
    print *, "Min/Max WRF 3D Pressure: ",minval(pr_grapess),maxval(pr_grapess)
  
    ! Get heights
    print *, "Getting PHB" 
    CALL get_wrfnc_3d(cdf,"PHB","T",nxw,nyw,nzw+1,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF base-state geopotential.'
      istatus = 0
      RETURN
    ENDIF
    dum3df2 = dum3df
  
    PRINT *, "Getting PH"
    CALL get_wrfnc_3d(cdf,"PH","T",nxw,nyw,nzw+1,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF geopotential.'
      istatus = 0
      RETURN
    ENDIF
    PRINT *,"Destaggering (vertically) heights"
    dum3df2 = (dum3df2 + dum3df) / grav
    DO k = 1,nzw
      ht_grapess(:,:,k) = 0.5 * (dum3df2(:,:,k) + dum3df2(:,:,k+1))
    ENDDO

    ! Get theta and convert to temperature
    PRINT *, "Getting Theta"
    CALL get_wrfnc_3d(cdf, "T","T",nxw,nyw,nzw,1,dum3d,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF perturbation theta.'
      istatus = 0
      RETURN
    ENDIF

    PRINT *, "Computing temp"
    dum3d = dum3d + 300.
    DO k = 1, nzw
      DO j = 1, nyw
        DO i = 1,nxw
          t3_grapess(i,j,k) = dum3d(i,j,k)/ ((100000./pr_grapess(i,j,k))**kappa)
        ENDDO
      ENDDO
    ENDDO

    ! Get Q on sigma
    PRINT *, "Getting Q"
    CALL get_wrfnc_3d(cdf, "QVAPOR","T",nxw,nyw,nzw,1,mr_grapess,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF mixing ratio.'
      istatus = 0
      RETURN
    ENDIF

    ! Derive specific humidity
    PRINT *, "Computing SH"
    sh_grapess = mr_grapess / (1. + mr_grapess)

    ! Get U on sigma   
    PRINT *, "Getting U"
    CALL get_wrfnc_3d(cdf, "U","T",nxw,nyw,nzw,1,u3_grapess,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF U-comp.'
      istatus = 0
      RETURN
    ENDIF

    ! Get V on sigma
    PRINT *, "Getting V"
    CALL get_wrfnc_3d(cdf, "V","T",nxw,nyw,nzw,1,v3_grapess,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF V-comp.'
      istatus = 0
      RETURN
    ENDIF

    ! Get W on sigma
    PRINT *, "Getting W"
    CALL get_wrfnc_3d(cdf, "W","T",nxw,nyw,nzw+1,1,dum3df,status)
    IF (status.NE.0) THEN
      PRINT *, 'Could not properly obtain WRF W-comp'
      istatus = 0
      RETURN 
    ENDIF
    PRINT*, "Destaggering (vertically) w"
    DO k = 1,nzw
      om_grapess(:,:,k) = 0.5*(dum3df(:,:,k)+dum3df(:,:,k+1))
    ENDDO

    ! Now, derive density, virtual potential temp, and omega
    PRINT *, "Computing Rho, Theta-V, and Omega"
    DO k = 1,nzw
      DO j=1,nyw
        DO i = 1,nxw
          rho_grapess(i,j,k) = pr_grapess(i,j,k) / ( r * t3_grapess(i,j,k)*(1.+0.61*sh_grapess(i,j,k)))
          om_grapess(i,j,k) = -1. * rho_grapess(i,j,k) * grav * om_grapess(i,j,k)
        ENDDO
      ENDDO
    ENDDO

    ! Get surface fields
    
    PRINT *, "Getting WRF TOPO"
    CALL get_wrfnc_2d(cdf, "HGT","T",nxw,nyw,1,topo_grapes,status)
    IF (status .NE. 0) THEN
      PRINT *, "Could not get topo from WRF"
      istatus = 0 
      RETURN
    ENDIF
    PRINT *, "Getting USF"
    usf_grapes = u3_grapess(:,:,1)
    PRINT *, "Getting VSF"
    vsf_grapes = v3_grapess(:,:,1)


    PRINT *, "Getting PSF"
    CALL get_wrfnc_2d(cdf, "PSFC","T",nxw,nyw,1,psf_grapes,status)
    IF ((status .NE. 0).OR.(MAXVAL(psf_grapes) .LT. 10000.))THEN
      PRINT *, "Could not get PSFC, using lowest sigma level"
      psf_grapes = pr_grapess(:,:,1)
    ENDIF
 
    PRINT *, "Getting T2" 
    CALL get_wrfnc_2d(cdf, "T2","T",nxw,nyw,1,tsf_grapes,status)
    IF ((status .NE. 0).OR.(MAXVAL(tsf_grapes) .LT. 100.))THEN
      PRINT *, "Could not get T2, using lowest sigma level"
      tsf_grapes = t3_grapess(:,:,1)
    ENDIF
    
    ! added tsk(skin temp.) by Wei-Ting (130312)
    PRINT *, "Getting TSK" 
    CALL get_wrfnc_2d(cdf, "TSK","T",nxw,nyw,1,tsk_grapes,status)
    IF ((status .NE. 0).OR.(MAXVAL(tsk_grapes) .LT. 100.))THEN
      PRINT *, "Could not get TSK, using lowest sigma level"
      tsk_grapes = t3_grapess(:,:,1)
    ENDIF
    
    ! added PCP (RAINNC+RAINC) by Wei-Ting (130312) & Modified (130326)
    PRINT *, "Getting Precipitation"
    PRINT *, "!!!!! This Precipitaion is an accumulation per N hours. !!!!!"
    PRINT *, "!!!!! N depends on the time difference of each WRFOUT.  !!!!!"
    PRINT *, "   Getting RAINNC(t)"
    CALL get_wrfnc_2d(cdf,"RAINNC","T",nxw,nyw,1,dum2dt2,status)
    IF (status .NE. 0) THEN
      PRINT *, "   Could not get RAINNC(t), setting the value = 0"
      dum2dt2 = 0
    ENDIF
    pcp_grapes = dum2dt2
    PRINT *, "   Getting RAINC(t)"
    CALL get_wrfnc_2d(cdf,"RAINC","T",nxw,nyw,1,dum2dt2,status)
    IF (status .NE. 0) THEN
      PRINT *, "   Could not get RAINC(t), setting the value = 0"
      dum2dt2 = 0
    ENDIF
    pcp_grapes = pcp_grapes+dum2dt2

    PRINT *, "   Getting RAINNC(t-1)"
    IF (cdp .LT. 0 .AND. cdf .GT. 0) THEN
      PRINT *, "   Could not get RAINNC(t-1), maybe result from t = initial time!"
      PRINT *, "   Set RAINNC(t-1) = 0"
      dum2dt1 = 0
    ELSE
      CALL get_wrfnc_2d(cdp,"RAINNC","T",nxw,nyw,1,dum2dt1,status)
      IF (status .NE. 0) THEN
         PRINT *, "   Could not get RAINNC(t-1), setting the value = 0"
         dum2dt1 = 0
      ENDIF
    ENDIF
    pcp_grapes = pcp_grapes-dum2dt1
    PRINT *, "   Getting RAINC(t-1)"
    IF (cdp .LT. 0 .AND. cdf .GT. 0) THEN
      PRINT *, "   Could not get RAINC(t-1), maybe result from t = initial time!"
      PRINT *, "   Set RAINC(t-1) = 0"
      dum2dt1 = 0
    ELSE
      CALL get_wrfnc_2d(cdp,"RAINC","T",nxw,nyw,1,dum2dt1,status)
      IF (status .NE. 0) THEN
         PRINT *, "   Could not get RAINC(t-1), setting the value = 0"
         dum2dt1 = 0
      ENDIF
    ENDIF
    pcp_grapes = pcp_grapes-dum2dt1
    where ( pcp_grapes < 0 ) ; pcp_grapes = 0 ; endwhere ! keep pcp >= 0
    print *, "Min/Max WRF Precipitation : ",minval(pcp_grapes),maxval(pcp_grapes)
    ! End of reading RAINNC+RAINC

    ! qvapor at 2m
    PRINT *, "Getting Q2"
    CALL get_wrfnc_2d(cdf, "Q2","T",nxw,nyw,1,rsf_grapes,status)
    IF ((status .NE. 0).OR.(MAXVAL(rsf_grapes) .LT. 0.0001))THEN
      PRINT *, "Could not get Q2, using lowest sigma level"
      rsf_grapes = sh_grapess(:,:,1)
    ELSE 

      ! Because 2m qv and T are derived from the PBL scheme and
      ! the WRF is apparently not checking for saturation, clean this
      ! up now
      PRINT *, "Checking Q2 for supersaturation"
      DO j = 1, nyw
        DO i= 1, nxw
          rsf_grapes(i,j) = MIN(rsf_grapes(i,j),mixsat(rsf_grapes(i,j),psf_grapes(i,j)))
          ! Compute dewpoint
          rh = relhum(tsf_grapes(i,j),rsf_grapes(i,j),psf_grapes(i,j))
          dsf_grapes(i,j) = dewpt(tsf_grapes(i,j),rh)
          ! Compute tvbar
          tvbar = 0.
          tvbar_nlevs = 0.
          comptvb:  DO k = 1,nzw
            IF ((psf_grapes(i,j)-pr_grapess(i,j,k)) .LE. tvbar_thick) THEN
               tvbar = tvbar + (t3_grapess(i,j,k)*(1.+0.61*sh_grapess(i,j,k)) )
               tvbar_nlevs = tvbar_nlevs + 1.
            ELSE
              exit comptvb
            ENDIF
          ENDDO comptvb
          tvb_grapes(i,j) = tvbar / tvbar_nlevs
        ENDDO
      ENDDO
    ENDIF 
    ! Smooth the tvb_grapes field
    CALL smooth2(nxw,nyw,4,tvb_grapes)
    ! Convert mr to sh
    rsf_grapes(:,:) = rsf_grapes(:,:)/(1. + rsf_grapes(:,:))

    ! Diagnostics
    PRINT *, "WRF Sigma data from center of WRF domain"
    print *, "K   PRESS     HEIGHT   TEMP   SH       U       V      OM"
    print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    DO k = 1,nzw
      PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_grapess(icentw,jcentw,k),ht_grapess(icentw,jcentw,k), t3_grapess(icentw,jcentw,k), &
          sh_grapess(icentw,jcentw,k),u3_grapess(icentw,jcentw,k),v3_grapess(icentw,jcentw,k),&
          om_grapess(icentw,jcentw,k)
    ENDDO
    PRINT *, "SFC:"
    PRINT ('(F8.1,2x,F7.0,2x,F5.1,2x,F5.1,2x,F5.1,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1)'), &
        psf_grapes(icentw,jcentw),topo_grapes(icentw,jcentw),tsf_grapes(icentw,jcentw), &
        tsk_grapes(icentw,jcentw),dsf_grapes(icentw,jcentw),tvb_grapes(icentw,jcentw), &
        rsf_grapes(icentw,jcentw), usf_grapes(icentw,jcentw), vsf_grapes(icentw,jcentw)
        ! added tsk_grapes by Wei-Ting (130312)
    
    PRINT *, "Deallocating arrays"
    IF ( ALLOCATED(dum3d)) DEALLOCATE(dum3d)
    IF ( ALLOCATED(dum3df)) DEALLOCATE(dum3df)
    IF ( ALLOCATED(dum3df2)) DEALLOCATE(dum3df2)
    IF ( ALLOCATED(dum2dt1)) DEALLOCATE(dum2dt1)
    IF ( ALLOCATED(dum2dt2)) DEALLOCATE(dum2dt2)
    RETURN
  END SUBROUTINE fill_wrfs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_grapes2p_standalone(nxw_in, nyw_in, nzw_in, nzl_in, &
       pr_grapess_in, pr_laps_in, ht_grapess_in, t3_grapess_in, sh_grapess_in, &
       u3_grapess_in, v3_grapess_in, om_grapess_in, mr_grapess_in, tvb_grapes_in, &
       ht_grapesp_out, t3_grapesp_out, sh_grapesp_out, u3_grapesp_out, v3_grapesp_out, &
       om_grapesp_out, icentw_in, jcentw_in, istatus)

    IMPLICIT NONE

    ! Input dimensions
    INTEGER, INTENT(IN)    :: nxw_in, nyw_in, nzw_in, nzl_in
    INTEGER, INTENT(IN)    :: icentw_in, jcentw_in
    
    ! Input arrays (on sigma/model levels)
    REAL, INTENT(IN)       :: pr_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: pr_laps_in(nzl_in)
    REAL, INTENT(IN)       :: ht_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: t3_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: sh_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: u3_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: v3_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: om_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: mr_grapess_in(nxw_in,nyw_in,nzw_in)
    REAL, INTENT(IN)       :: tvb_grapes_in(nxw_in,nyw_in)
    
    ! Output arrays (on pressure levels)
    REAL, INTENT(OUT)      :: ht_grapesp_out(nxw_in,nyw_in,nzl_in)
    REAL, INTENT(OUT)      :: t3_grapesp_out(nxw_in,nyw_in,nzl_in)
    REAL, INTENT(OUT)      :: sh_grapesp_out(nxw_in,nyw_in,nzl_in)
    REAL, INTENT(OUT)      :: u3_grapesp_out(nxw_in,nyw_in,nzl_in)
    REAL, INTENT(OUT)      :: v3_grapesp_out(nxw_in,nyw_in,nzl_in)
    REAL, INTENT(OUT)      :: om_grapesp_out(nxw_in,nyw_in,nzl_in)
    INTEGER, INTENT(OUT)   :: istatus
    
    ! Local variables
    INTEGER  :: i,j,k,ks,kp, ksb,kst
    REAL     :: lpb,lpt,lp, wgtb, wgtt
    REAL,PARAMETER     :: dTdlnPBase = 50.0
    REAL               :: tvbot, tvbar, deltalnp
    REAL               :: dz
    
    istatus = 1

    ! Loop over horizontal domain, interpolating vertically for each column
    DO j = 1, nyw_in
      DO i = 1, nxw_in   

        pressloop: DO kp = 1,nzl_in
  
          ! Initialize kst and ksb, which will hold the vertical sigma
          ! index values of the top and bottom bounding layers
          kst = 0
          ksb = 0 

          ! Find bounding levels in raw data
          sigmaloop: DO ks = 1,nzw_in
            IF (pr_grapess_in(i,j,ks) .LE. pr_laps_in(kp)) THEN   

              kst = ks
              ksb = kst - 1
              EXIT sigmaloop
            ENDIF
          ENDDO sigmaloop

          IF (kst .GT. 1) THEN ! Interpolate between two bounding points
            lp = ALOG(pr_laps_in(kp))
            lpt = ALOG(pr_grapess_in(i,j,kst))
            lpb = ALOG(pr_grapess_in(i,j,ksb))
            wgtb = (lpt - lp) / (lpt - lpb)
            wgtt = 1.0 - wgtb

            ! Height
            ht_grapesp_out(i,j,kp) = wgtb * ht_grapess_in(i,j,ksb) + &
                              wgtt * ht_grapess_in(i,j,kst)

            ! Temp
            t3_grapesp_out(i,j,kp) = wgtb * t3_grapess_in(i,j,ksb) + &
                              wgtt * t3_grapess_in(i,j,kst)

            ! SH
            sh_grapesp_out(i,j,kp) = wgtb * sh_grapess_in(i,j,ksb) + &
                              wgtt * sh_grapess_in(i,j,kst)

            ! U3
            u3_grapesp_out(i,j,kp) = wgtb * u3_grapess_in(i,j,ksb) + &
                              wgtt * u3_grapess_in(i,j,kst)

            ! V3
            v3_grapesp_out(i,j,kp) = wgtb * v3_grapess_in(i,j,ksb) + &
                              wgtt * v3_grapess_in(i,j,kst)
 
            ! OM
            om_grapesp_out(i,j,kp) = wgtb * om_grapess_in(i,j,ksb) + &
                              wgtt * om_grapess_in(i,j,kst)


          ELSEIF (kst .EQ. 1) THEN ! Extrapolate downward
            lpt = ALOG(pr_grapess_in(i,j,kst))
            lpb = ALOG(pr_laps_in(kp))
            deltalnp = lpb - lpt
            tvbot = tvb_grapes_in(i,j) + deltalnp*dTdlnPBase
            tvbar = 0.5*(tvb_grapes_in(i,j) + tvbot)
            ! Height
 
            IF ((pr_laps_in(kp) - pr_grapess_in(i,j,1)).LT. 500.) THEN
              ! Very small difference in pressures, so
              ! assume 10 m per 100 Pa, because
              ! hypsometric eq breaks down in these cases
              dz =  0.1 * (pr_laps_in(kp) - pr_grapess_in(i,j,1))
            ELSE
              dz = tvbar * rog * ALOG(pr_laps_in(kp)/pr_grapess_in(i,j,1))
            ENDIF
            ht_grapesp_out(i,j,kp) = ht_grapess_in(i,j,1) - dz
            ! Temp
            t3_grapesp_out(i,j,kp) = tvbot/(1.+0.61*mr_grapess_in(i,j,1))
            ! SH
            sh_grapesp_out(i,j,kp) = sh_grapess_in(i,j,1)
            ! U3
            u3_grapesp_out(i,j,kp) = u3_grapess_in(i,j,1)

            ! V3
            v3_grapesp_out(i,j,kp) = v3_grapess_in(i,j,1)
           ! OM
            om_grapesp_out(i,j,kp) = 0.
          ELSE ! kst never got set .. extrapolate upward
            ! Assume isothermal (above tropopause)
            t3_grapesp_out(i,j,kp) = t3_grapess_in(i,j,nzw_in)
            u3_grapesp_out(i,j,kp) = u3_grapess_in(i,j,nzw_in)
            v3_grapesp_out(i,j,kp) = v3_grapess_in(i,j,nzw_in)
            om_grapesp_out(i,j,kp) = 0.
            dz = t3_grapess_in(i,j,nzw_in) * rog * ALOG(pr_grapess_in(i,j,nzw_in)/pr_laps_in(kp))
            ht_grapesp_out(i,j,kp) = ht_grapess_in(i,j,nzw_in) + dz 
            ! Reduce moisture toward zero
            IF ( pr_laps_in(kp-1) .GT. pr_grapess_in(i,j,nzw_in) ) THEN
              sh_grapesp_out(i,j,kp) = 0.5*sh_grapess_in(i,j,nzw_in)
            ELSE
              sh_grapesp_out(i,j,kp) = 0.5*sh_grapesp_out(i,j,kp-1)
            ENDIF
          ENDIF

        ENDDO pressloop
      ENDDO
    ENDDO

    ! If we do smoothing, do it here
    
    ! Print some diagnostics
    PRINT *, "GRAPES Press data from center of domain"
    print *, "KP  PRESS     HEIGHT   TEMP   SH       U       V      OM"
    print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    DO k = 1,nzl_in
      PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_laps_in(k),ht_grapesp_out(icentw_in,jcentw_in,k), t3_grapesp_out(icentw_in,jcentw_in,k), &
          sh_grapesp_out(icentw_in,jcentw_in,k),u3_grapesp_out(icentw_in,jcentw_in,k),v3_grapesp_out(icentw_in,jcentw_in,k),&
          om_grapesp_out(icentw_in,jcentw_in,k)
    ENDDO
    
    RETURN
  END SUBROUTINE vinterp_grapes2p_standalone
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_grapes2p(istatus)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)  :: istatus
    INTEGER  :: i,j,k,ks,kp, ksb,kst
    REAL     :: lpb,lpt,lp, wgtb, wgtt
    REAL,PARAMETER     :: dTdlnPBase = 50.0
    REAL               :: tvbot, tvbar, deltalnp
    REAL               :: dz
    istatus = 1

    ! Loop over horizontal domain, interpolating vertically for each
    ! column
    DO j = 1, nyw
      DO i = 1, nxw   

        pressloop: DO kp = 1,nzl
  
          ! Initialize kst and ksb, which will hold the vertical sigma
          ! index values of the top and bottom bounding layers
          kst = 0
          ksb = 0 

          ! Find bounding levels in raw data
          sigmaloop: DO ks = 1,nzw
            IF (pr_grapess(i,j,ks) .LE. pr_laps(kp)) THEN   

              kst = ks
              ksb = kst - 1
              EXIT sigmaloop
            ENDIF
          ENDDO sigmaloop

          IF (kst .GT. 1) THEN ! Interpolate between two bounding points
            lp = ALOG(pr_laps(kp))
            lpt = ALOG(pr_grapess(i,j,kst))
            lpb = ALOG(pr_grapess(i,j,ksb))
            wgtb = (lpt - lp) / (lpt - lpb)
            wgtt = 1.0 - wgtb

            ! Height
            ht_grapesp(i,j,kp) = wgtb * ht_grapess(i,j,ksb) + &
                              wgtt * ht_grapess(i,j,kst)

            ! Temp
            t3_grapesp(i,j,kp) = wgtb * t3_grapess(i,j,ksb) + &
                              wgtt * t3_grapess(i,j,kst)

            ! SH
            sh_grapesp(i,j,kp) = wgtb * sh_grapess(i,j,ksb) + &
                              wgtt * sh_grapess(i,j,kst)

            ! U3
            u3_grapesp(i,j,kp) = wgtb * u3_grapess(i,j,ksb) + &
                              wgtt * u3_grapess(i,j,kst)

            ! V3
            v3_grapesp(i,j,kp) = wgtb * v3_grapess(i,j,ksb) + &
                              wgtt * v3_grapess(i,j,kst)
 
            ! OM
            om_grapesp(i,j,kp) = wgtb * om_grapess(i,j,ksb) + &
                              wgtt * om_grapess(i,j,kst)


          ELSEIF (kst .EQ. 1) THEN ! Extrapolate downward
            lpt = ALOG(pr_grapess(i,j,kst))
            lpb = ALOG(pr_laps(kp))
            deltalnp = lpb - lpt
            tvbot = tvb_grapes(i,j) + deltalnp*dTdlnPBase
            tvbar = 0.5*(tvb_grapes(i,j) + tvbot)
            ! Height
 
            IF ((pr_laps(kp) - pr_grapess(i,j,1)).LT. 500.) THEN
              ! Very small difference in pressures, so
              ! assume 10 m per 100 Pa, because
              ! hypsometric eq breaks down in these cases
              dz =  0.1 * (pr_laps(kp) - pr_grapess(i,j,1))
            ELSE
              dz = tvbar * rog * ALOG(pr_laps(kp)/pr_grapess(i,j,1))
            ENDIF
            ht_grapesp(i,j,kp) = ht_grapess(i,j,1) - dz
            ! Temp
            t3_grapesp(i,j,kp) = tvbot/(1.+0.61*mr_grapess(i,j,1))
            ! SH
            sh_grapesp(i,j,kp) = sh_grapess(i,j,1)
            ! U3
            u3_grapesp(i,j,kp) = u3_grapess(i,j,1)

            ! V3
            v3_grapesp(i,j,kp) = v3_grapess(i,j,1)
           ! OM
            om_grapesp(i,j,kp) = 0.
          ELSE ! kst never got set .. extrapolate upward
            ! Assume isothermal (above tropopause)
            t3_grapesp(i,j,kp) = t3_grapess(i,j,nzw)
            u3_grapesp(i,j,kp) = u3_grapess(i,j,nzw)
            v3_grapesp(i,j,kp) = v3_grapess(i,j,nzw)
            om_grapesp(i,j,kp) = 0.
            dz = t3_grapess(i,j,nzw) * rog * ALOG(pr_grapess(i,j,nzw)/pr_laps(kp))
            ht_grapesp(i,j,kp) = ht_grapess(i,j,nzw) + dz 
            ! Reduce moisture toward zero
            IF ( pr_laps(kp-1) .GT. pr_grapess(i,j,nzw) ) THEN
              sh_grapesp(i,j,kp) = 0.5*sh_grapess(i,j,nzw)
            ELSE
              sh_grapesp(i,j,kp) = 0.5*sh_grapesp(i,j,kp-1)
            ENDIF
          ENDIF

        ENDDO pressloop
      ENDDO
    ENDDO

    ! Diagnostic: find any unreasonable t3_grapesp values
    DO kp = 1, nzl
      DO j = 1, nyw
        DO i = 1, nxw
          IF (t3_grapesp(i,j,kp) .GT. 350.0 .OR. t3_grapesp(i,j,kp) .LT. 150.0) THEN
            PRINT *, 'BAD t3_grapesp: i,j,kp,T,P_laps,P_grapes_sfc,P_grapes_top=', &
              i, j, kp, t3_grapesp(i,j,kp), pr_laps(kp), pr_grapess(i,j,1), pr_grapess(i,j,nzw)
          ENDIF
          IF (ht_grapesp(i,j,kp) .LT. -900.0) THEN
            PRINT *, 'BAD ht_grapesp: i,j,kp,ht,P_laps,P_grapes_sfc,ht_grapess_sfc=', &
              i, j, kp, ht_grapesp(i,j,kp), pr_laps(kp), pr_grapess(i,j,1), ht_grapess(i,j,1)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ! If we do smoothing, do it here
    
    ! Print some diagnostics
    PRINT *, "GRAPES Press data from center of GRAPES domain"
    print *, "KP  PRESS     HEIGHT   TEMP   SH       U       V      OM"
    print *, "--- --------  -------  -----  -------  ------  ------ -----------"
    DO k = 1,nzl
      PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_laps(k),ht_grapesp(icentw,jcentw,k), t3_grapesp(icentw,jcentw,k), &
          sh_grapesp(icentw,jcentw,k),u3_grapesp(icentw,jcentw,k),v3_grapesp(icentw,jcentw,k),&
          om_grapesp(icentw,jcentw,k)
    ENDDO
    
    RETURN
  END SUBROUTINE vinterp_grapes2p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE vinterp_grapes2ht(istatus)
  ! Vertical interpolation from GRAPES native levels to LAPS SIGMA_HT levels.
  ! For SIGMA_HT grids pr_laps holds the target height levels (m) loaded via
  ! get_ht_1d.  GRAPES native heights ht_grapess are used as the vertical
  ! coordinate.  Outputs: pr_grapesp (3D pressure at target heights, Pa),
  ! ht_grapesp (diagnostic, equals target heights), t3_grapesp, sh_grapesp,
  ! u3_grapesp, v3_grapesp, om_grapesp.

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus
    INTEGER :: i, j, k, ks, kp, ksb, kst
    REAL    :: htb, htt, ht_tgt, wgtb, wgtt
    istatus = 1

    DO j = 1, nyw
      DO i = 1, nxw
        htloop: DO kp = 1, nzl
          ht_tgt = pr_laps(kp)  ! pr_laps holds sigma-height levels (m) for SIGMA_HT
          kst = 0
          ksb = 0
          ! GRAPES native levels: k=1 is bottom (lowest height), k=nzw is top
          sigloop: DO ks = 1, nzw
            IF (ht_grapess(i,j,ks) .GE. ht_tgt) THEN
              kst = ks
              ksb = kst - 1
              EXIT sigloop
            ENDIF
          ENDDO sigloop

          IF (kst .GT. 1) THEN  ! interpolate between ksb and kst
            htt  = ht_grapess(i,j,kst)
            htb  = ht_grapess(i,j,ksb)
            wgtb = (htt - ht_tgt) / (htt - htb)
            wgtt = 1.0 - wgtb
            pr_grapesp(i,j,kp) = wgtb*pr_grapess(i,j,ksb) + wgtt*pr_grapess(i,j,kst)
            ht_grapesp(i,j,kp) = ht_tgt
            t3_grapesp(i,j,kp) = wgtb*t3_grapess(i,j,ksb) + wgtt*t3_grapess(i,j,kst)
            sh_grapesp(i,j,kp) = wgtb*sh_grapess(i,j,ksb) + wgtt*sh_grapess(i,j,kst)
            u3_grapesp(i,j,kp) = wgtb*u3_grapess(i,j,ksb) + wgtt*u3_grapess(i,j,kst)
            v3_grapesp(i,j,kp) = wgtb*v3_grapess(i,j,ksb) + wgtt*v3_grapess(i,j,kst)
            om_grapesp(i,j,kp) = wgtb*om_grapess(i,j,ksb) + wgtt*om_grapess(i,j,kst)

          ELSEIF (kst .EQ. 1) THEN  ! target below lowest GRAPES level: use bottom
            pr_grapesp(i,j,kp) = pr_grapess(i,j,1)
            ht_grapesp(i,j,kp) = ht_tgt
            t3_grapesp(i,j,kp) = t3_grapess(i,j,1)
            sh_grapesp(i,j,kp) = sh_grapess(i,j,1)
            u3_grapesp(i,j,kp) = u3_grapess(i,j,1)
            v3_grapesp(i,j,kp) = v3_grapess(i,j,1)
            om_grapesp(i,j,kp) = 0.0

          ELSE  ! kst == 0: target above highest GRAPES level: use top
            pr_grapesp(i,j,kp) = pr_grapess(i,j,nzw)
            ht_grapesp(i,j,kp) = ht_tgt
            t3_grapesp(i,j,kp) = t3_grapess(i,j,nzw)
            sh_grapesp(i,j,kp) = MAX(0.0, 0.5*sh_grapess(i,j,nzw))
            u3_grapesp(i,j,kp) = u3_grapess(i,j,nzw)
            v3_grapesp(i,j,kp) = v3_grapess(i,j,nzw)
            om_grapesp(i,j,kp) = 0.0
          ENDIF
        ENDDO htloop
      ENDDO
    ENDDO

    PRINT *, "GRAPES sigma-height data from center of GRAPES domain"
    PRINT *, "KP  HT_TGT(m)  PRESS(Pa)  TEMP(K)  SH       U(m/s)  V(m/s)"
    PRINT *, "--- ----------  ---------  -------  ------   ------  ------"
    DO k = 1, nzl
      PRINT ('(I3,1x,F10.1,2x,F9.1,2x,F7.2,2x,F6.4,2x,F6.1,2x,F6.1)'), &
        k, pr_laps(k), pr_grapesp(icentw,jcentw,k), t3_grapesp(icentw,jcentw,k), &
        sh_grapesp(icentw,jcentw,k), u3_grapesp(icentw,jcentw,k), v3_grapesp(icentw,jcentw,k)
    ENDDO

    RETURN
  END SUBROUTINE vinterp_grapes2ht

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE hinterp_grapes2lga(istatus)

  IMPLICIT NONE
  INTEGER, INTENT(OUT)  :: istatus
  INTEGER               :: i,j,k
  INTEGER               :: n_outside
  REAL, ALLOCATABLE     :: xloc(:,:), yloc(:,:),dum2d(:,:)
  REAL, ALLOCATABLE     :: topo_grapesl(:,:)  ! WRF Topo interpolated to LAPS grid
  REAL                  :: ri,rj,dtopo,dz, wgt1, wgt2,lp
  REAL                  :: pct_outside
  istatus = 1 

  IF (.NOT. ALLOCATED(dum2d))  ALLOCATE(dum2d(nxl,nyl))
  ! First, generate xloc/yloc locations
  IF (.NOT. ALLOCATED(xloc))  ALLOCATE(xloc(nxl,nyl))
  IF (.NOT. ALLOCATED(yloc))  ALLOCATE(yloc(nxl,nyl))
  n_outside = 0
  DO j = 1, nyl
    DO i = 1, nxl
      CALL latlon_to_ij(grapesgrid,lat(i,j),lon(i,j),ri,rj)
      IF (ri < 1.0 .OR. ri > REAL(nxw) .OR. rj < 1.0 .OR. rj > REAL(nyw)) THEN
        n_outside = n_outside + 1
      ENDIF
      ! Clamp to valid GRAPES grid bounds. LAPS domain may extend outside
      ! the GRAPES domain; unbound indices cause out-of-bounds array access
      ! in interpolate_standard (which only guards the upper limit).
      xloc(i,j) = MAX(1.0, MIN(REAL(nxw), ri))
      yloc(i,j) = MAX(1.0, MIN(REAL(nyw), rj))
    ENDDO
  ENDDO

  ! Check how many LAPS points fall outside the GRAPES domain
  IF (n_outside > 0) THEN
    pct_outside = 100.0 * REAL(n_outside) / REAL(nxl*nyl)
    WRITE(6,'(A,I0,A,F5.1,A)') &
      'WARNING: hinterp_grapes2lga: ', n_outside, ' LAPS grid points (', &
      pct_outside, '%) fall outside the GRAPES domain.'
    WRITE(6,'(A,I0,A,I0,A)') &
      '  GRAPES domain: ', nxw, ' x ', nyw, ' grid points'
    WRITE(6,'(A,I0,A,I0,A)') &
      '  LAPS   domain: ', nxl, ' x ', nyl, ' grid points'
    WRITE(6,'(A)') &
      '  Points outside the GRAPES boundary are filled by nearest-edge extrapolation.'
    IF (pct_outside > 50.0) THEN
      WRITE(6,'(A)') &
        'ERROR: More than 50% of the LAPS domain lies outside the GRAPES domain.'
      WRITE(6,'(A)') &
        '       Please configure a LAPS domain that is covered by your GRAPES data.'
      istatus = 0
      RETURN
    ENDIF
  ENDIF

  ! Now, horizontally interpolate level by level
  DO k=1,nzl
    ! Height
    CALL interpolate_standard(nxw,nyw,ht_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    ht(:,:,k) = dum2d


    ! Temp
    CALL interpolate_standard(nxw,nyw,t3_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    t3(:,:,k) = dum2d

    ! SH
    CALL interpolate_standard(nxw,nyw,sh_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    sh(:,:,k) = dum2d

    ! U3
    CALL interpolate_standard(nxw,nyw,u3_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    u3(:,:,k) = dum2d

    ! v3
    CALL interpolate_standard(nxw,nyw,v3_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    v3(:,:,k) = dum2d

    ! om
    CALL interpolate_standard(nxw,nyw,om_grapesp(:,:,k), &
                              nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                              dum2d)
    om(:,:,k) = dum2d

    ! 3D pressure (SIGMA_HT only)
    IF (ALLOCATED(pr_grapesp) .AND. ALLOCATED(prgd)) THEN
      CALL interpolate_standard(nxw,nyw,pr_grapesp(:,:,k), &
                                nxl,nyl,xloc,yloc,METHOD_LINEAR, &
                                dum2d)
      prgd(:,:,k) = dum2d
    END IF

  ENDDO

  PRINT *, "hinterp_grapes2lga: LGA Press data from center of LAPS domain"
  print *, "KP  PRESS     HEIGHT   TEMP   SH       U       V      OM"
  print *, "--- --------  -------  -----  -------  ------  ------ -----------"
  DO k = 1,nzl
    PRINT ('(I3,1x,F8.1,2x,F7.0,2x,F5.1,2x,F7.5,2x,F6.1,2x,F6.1,2x,F11.8)'), &
          k,pr_laps(k),ht(icentl,jcentl,k), t3(icentl,jcentl,k), &
          sh(icentl,jcentl,k),u3(icentl,jcentl,k),v3(icentl,jcentl,k),&
          om(icentl,jcentl,k)
  ENDDO

  IF (ALLOCATED(dum2d))  DEALLOCATE(dum2d)

  ALLOCATE(topo_grapesl(nxl,nyl))
  ! Do the surface variables
  CALL interpolate_standard(nxw,nyw,usf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,usf)
  CALL interpolate_standard(nxw,nyw,vsf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,vsf)
  CALL interpolate_standard(nxw,nyw,tsf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,tsf)
  CALL interpolate_standard(nxw,nyw,tsk_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,tsk) ! added tsk by Wei-Ting (130312)
  CALL interpolate_standard(nxw,nyw,dsf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,dsf)
  CALL interpolate_standard(nxw,nyw,rsf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,rsf)
  CALL interpolate_standard(nxw,nyw,psf_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,psf)
  CALL interpolate_standard(nxw,nyw,topo_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,topo_grapesl)
  CALL interpolate_standard(nxw,nyw,pcp_grapes, nxl,nyl,xloc,yloc, &
      METHOD_LINEAR,pcp) ! added pcp by Wei-Ting (130312)
  where ( pcp < 0 ) ; pcp = 0 ; endwhere ! keep pcp >= 0 added by Wei-Ting (130312)

  ! Adjust for terrain differences between WRF and LAPS
  PRINT *, "Adjusting WRF surface to LAPS surface"
  DO j = 1 , nyl
    DO i = 1 , nxl
      dtopo = topo_grapesl(i,j) - topo_laps(i,j)
      IF (ABS(dtopo) .GT. 10.) THEN
        IF (dtopo .GT. 0) THEN ! Move downward to LAPS level
          downloop:  DO k = nzl , 1, -1
            IF (ht(i,j,k) .LE. topo_laps(i,j)) THEN
              dz = topo_grapesl(i,j) - ht(i,j,k)
              wgt1 = (topo_grapesl(i,j) - topo_laps(i,j)) / dz
              wgt2 = 1.0 - wgt1
              IF (ALLOCATED(prgd)) THEN
                lp = wgt1*ALOG(prgd(i,j,k)) + wgt2*ALOG(psf(i,j))
              ELSE
                lp = wgt1*ALOG(pr_laps(k)) + wgt2*ALOG(psf(i,j))
              END IF
              psf(i,j) = EXP(lp)
              tsf(i,j) = wgt1*t3(i,j,k)  + wgt2*tsf(i,j)
              tsk(i,j) = wgt1*t3(i,j,k)  + wgt2*tsk(i,j) ! added tsk by Wei-Ting (130312)
              rsf(i,j) = wgt1*sh(i,j,k)  + wgt2*rsf(i,j)
              usf(i,j) = wgt1*u3(i,j,k)  + wgt2*usf(i,j)
              vsf(i,j) = wgt1*v3(i,j,k)  + wgt2*vsf(i,j)

              EXIT downloop
            ENDIF
          ENDDO downloop
        ELSE  ! Move upward to LAPS level
          uploop:  DO k = 1, nzl
            IF (ht(i,j,k) .GE. topo_laps(i,j)) THEN
              dz = ht(i,j,k) - topo_grapesl(i,j)
              wgt2 = (topo_laps(i,j) - topo_grapesl(i,j)) / dz
              wgt1 = 1.0 - wgt2
              IF (ALLOCATED(prgd)) THEN
                psf(i,j) = wgt1 * psf(i,j) + wgt2*prgd(i,j,k)
              ELSE
                psf(i,j) = wgt1 * psf(i,j) + wgt2*pr_laps(k)
              END IF
              tsf(i,j) = wgt1 * tsf(i,j) + wgt2*t3(i,j,k)
              tsk(i,j) = wgt1 * tsk(i,j) + wgt2*t3(i,j,k) ! added tsk by Wei-Ting (130312)
              rsf(i,j) = wgt1 * rsf(i,j) + wgt2*sh(i,j,k)
              usf(i,j) = wgt1 * usf(i,j) + wgt2*u3(i,j,k)
              vsf(i,j) = wgt1 * vsf(i,j) + wgt2*v3(i,j,k) 
 
              EXIT uploop
            ENDIF
          ENDDO uploop
        ENDIF
      ENDIF
   
    ENDDO
  ENDDO
  IF (ALLOCATED(xloc))  DEALLOCATE(xloc)
  IF (ALLOCATED(yloc))  DEALLOCATE(yloc)
  IF (ALLOCATED(topo_grapesl))  DEALLOCATE(topo_grapesl)

  END SUBROUTINE hinterp_grapes2lga

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_derived_pressures

    IMPLICIT NONE
    INTEGER :: i,j,k
    REAL    :: wgt1, wgt2, dz, lp, lp1, lp2
    LOGICAL :: sigma_ht_mode

    sigma_ht_mode = ALLOCATED(prgd)  ! prgd only allocated for SIGMA_HT

    DO j = 1, nyl
      DO i = 1, nxl

        ! Sea-level pressure: find lowest level at or above z=0
        ! Start at k=2 so that k-1 is always valid
        slploop: DO k = 2, nzl
          IF (ht(i,j,k) .GE. 0.0) THEN
            IF (sigma_ht_mode) THEN
              lp1 = ALOG(prgd(i,j,k-1))
              lp2 = ALOG(prgd(i,j,k))
            ELSE
              lp1 = ALOG(pr_laps(k-1))
              lp2 = ALOG(pr_laps(k))
            END IF
            dz   = ht(i,j,k) - ht(i,j,k-1)
            IF (dz .NE. 0.0) THEN
              wgt1 = ht(i,j,k) / dz
              wgt2 = 1.0 - wgt1
              lp   = wgt1 * lp1 + wgt2 * lp2
            ELSE
              lp = lp2
            END IF
            slp(i,j) = EXP(lp)
            EXIT slploop
          ENDIF
        ENDDO slploop

        ! LAPS reduced pressure: interpolate to redp_lvl
        ! Start at k=2 so that k-1 is always valid
        redploop: DO k = 2, nzl
          IF (ht(i,j,k) .GE. redp_lvl) THEN
            IF (sigma_ht_mode) THEN
              lp1 = ALOG(prgd(i,j,k-1))
              lp2 = ALOG(prgd(i,j,k))
            ELSE
              lp1 = ALOG(pr_laps(k-1))
              lp2 = ALOG(pr_laps(k))
            END IF
            dz   = ht(i,j,k) - ht(i,j,k-1)
            IF (dz .NE. 0.0) THEN
              wgt1 = (ht(i,j,k) - redp_lvl) / dz
              wgt2 = 1.0 - wgt1
              lp   = wgt1 * lp1 + wgt2 * lp2
            ELSE
              lp = lp2
            END IF
            p(i,j) = EXP(lp)
            EXIT redploop
          ENDIF
        ENDDO redploop
      ENDDO
    ENDDO
    RETURN
  END SUBROUTINE make_derived_pressures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE grapes_lga
