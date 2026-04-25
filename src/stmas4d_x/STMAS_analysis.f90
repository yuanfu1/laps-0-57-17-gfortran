SUBROUTINE STMAS_analysis

!doc==============================================================================================
!doc  This routine performs STMAS surface analysis using a single multigrid cycle from coarse to
!doc  finest grids.
!doc
!doc  History: June. 2010 by Yuanfu Xie.
!doc
!doc  Modified for STMAS 4D analysis. June 10, 2011 Hongli Jiang 
!doc==============================================================================================

  USE STMAS

  IMPLICIT NONE
 
  ! Local variables:
  
  INTEGER :: mg,num_controls
  CHARACTER*17, PARAMETER :: header = 'STMAS_analysis: '
  integer :: i,j,k,t

  !-----------------------------------------------------------------------------------------------
  ! Multigrid cycle: Start a multigrid:
  !-----------------------------------------------------------------------------------------------
  DO mg=1,STMAS_numlevels

! HJ add: allocate memory for mgrid. 6/29/2011
    call STMAS_memo_alloc_mgrid(mg)

    !---------------------------------------------------------------------------------------------
    ! Initial guess:
    !---------------------------------------------------------------------------------------------
    IF (mg .eq. 1) THEN
      CALL get_initial       ! coarsest grid use background
      call STMAS_memo_alloc_tmgrid(mg)
    ELSE
! Project previous analysis onto a finer grid as an initial guess:
      CALL projection(mg,STMAS_tmgrid)
! release memory of tmgrid from mg-1 grid
      call STMAS_memo_release_tmgrid
! allocate memory of tmgrid for current mg
      call STMAS_memo_alloc_tmgrid(mg)
    ENDIF

    !---------------------------------------------------------------------------------------------
    ! Map observations onto grid:
    !---------------------------------------------------------------------------------------------
!    CALL STMAS_obsmapping(mg,observations,STMAS_mgrid,STMAS_invalid,gridobs)

    !---------------------------------------------------------------------------------------------
    ! Minimization at the current multgrid level:
    !---------------------------------------------------------------------------------------------
! HJ add STMAS_mgrid%numgrid(4). 6/10/2011
    num_controls = STMAS_mgrid%numgrid(1)*STMAS_mgrid%numgrid(2)* &
                   STMAS_mgrid%numgrid(3)*STMAS_mgrid%numgrid(4)*STMAS_numvars

    CALL STMAS_minimizer(num_controls,mg,STMAS_mgrid)

! Copy STMAS_mgrid%anal only to the temporary array STMAS_tmgrid by calling STMAS_copyto,
! then call STMAS_memo_release_mgrid to deallocate memory. HJ 8/9/2011
! Passing in STMAS_final (use STMAS_final%bkgd only) to store only the increment, i.e. anal-bkgd. 
! this can be done in the subroutine projection as well.  HJ 8/17/2011
    CALL STMAS_copyto(STMAS_mgrid,STMAS_tmgrid,STMAS_final)

! HJ add: deallocate memory for mgrid. 6/30/2011
    call STMAS_memo_release_mgrid
  ENDDO
 
  !-----------------------------------------------------------------------------------------------
  ! Minimization at the final analysis level:
  !-----------------------------------------------------------------------------------------------

! HJ add STMAS_mgrid%numgrid(4). 6/10/2011
  IF (STMAS_final%numgrid(1) .NE. STMAS_mgrid%numgrid(1) .OR. &
      STMAS_final%numgrid(2) .NE. STMAS_mgrid%numgrid(2) .OR. &
      STMAS_final%numgrid(3) .NE. STMAS_mgrid%numgrid(3) .OR. &
      STMAS_final%numgrid(4) .NE. STMAS_mgrid%numgrid(4) ) THEN

    PRINT*,'+------------------------------------------------------------------------+'
    PRINT*,header,' Additional analysis at the final grid'
    PRINT*,'+------------------------------------------------------------------------+'


! HJ mod: call to grid4_interpolation. 6/10/2011
    CALL grid4_interpolation(STMAS_tmgrid%anal(1:STMAS_tmgrid%numgrid(1), &
                                               1:STMAS_tmgrid%numgrid(2), &
                                               1:STMAS_tmgrid%numgrid(3), &
                                               1:STMAS_tmgrid%numgrid(4), &
                                               1:STMAS_numvars),STMAS_tmgrid%numgrid, &
                             STMAS_final%anal(1:STMAS_final%numgrid(1), &
                                              1:STMAS_final%numgrid(2), &
                                              1:STMAS_final%numgrid(3), &
                                              1:STMAS_final%numgrid(4), &
                                              1:STMAS_numvars),STMAS_final%numgrid, &
                                              STMAS_numvars)

    print*,'final-tmgrid',STMAS_tmgrid%numgrid,' final-tgrid pressure: ', &
      MINVAL(STMAS_tmgrid%anal(:,:,:,:,4)),MAXVAL(STMAS_tmgrid%anal(:,:,:,:,4))
    print*,'final-num',STMAS_final%numgrid,' final pressure: ', &
      MINVAL(STMAS_final%anal(:,:,:,:,4)),MAXVAL(STMAS_final%anal(:,:,:,:,4))

! Add the increment to the background as an initial guess:
    STMAS_final%anal = STMAS_final%anal+STMAS_final%bkgd
    print*,'final-tot',STMAS_final%anal(2,2,10,1,5)

    ! Map observation onto the finest grid:
! comment out. HJ 6/10/2011
!    CALL STMAS_obsmapping(STMAS_numlevels,observations,STMAS_final,STMAS_invalid,gridobs)

    ! Run final minimization on the finest grid:
! HJ add STMAS_final%numgrid(4). 6/10/2011
    num_controls = STMAS_final%numgrid(1)*STMAS_final%numgrid(2)* &
                   STMAS_final%numgrid(3)*STMAS_final%numgrid(4)*STMAS_numvars
    CALL STMAS_minimizer(num_controls,STMAS_numlevels,STMAS_final)
  ENDIF

END SUBROUTINE STMAS_analysis


SUBROUTINE get_initial
!doc==============================================================================================
!doc  This routine initializes a coarse grid using the first guess
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc==============================================================================================

  USE STMAS

  IMPLICIT NONE

  ! Local variables:
  INTEGER :: i,j
 
  !-----------------------------------------------------------------------------------------------
  ! Initial grid points:
  !-----------------------------------------------------------------------------------------------

!  STMAS_mgrid%numgrid = STMAS_start_grdpts
! define STMAS_final%incr for use with projection of background during bi-linear interpolation 
  ! Calculate finest multigrid:
  STMAS_mgrid%nfinest = STMAS_start_grdpts
  STMAS_mgrid%incr = 1  ! Default increment
  DO i=1,STMAS_numlevels
    DO j=1,STMAS_maxdim
      IF (2*(STMAS_mgrid%nfinest(j)-1)+1 .LE. STMAS_final%numgrid(j)) THEN
        STMAS_mgrid%nfinest(j) = 2*(STMAS_mgrid%nfinest(j)-1)+1
        STMAS_mgrid%incr(j) = STMAS_mgrid%incr(j)*2
        STMAS_final%incr(j) = STMAS_final%incr(j)*2
      ENDIF
    ENDDO
  ENDDO
!  print*,'nfinest',STMAS_mgrid%nfinest
!  print*,'final%inc',STMAS_final%incr

  !-----------------------------------------------------------------------------------------------
  ! Gridspacing:
  !-----------------------------------------------------------------------------------------------
  STMAS_mgrid%gridspc(1:4) = STMAS_final%gridspc(1:4)* &
    (STMAS_final%numgrid(1:4)-1)/FLOAT(STMAS_mgrid%numgrid(1:4)-1)

  !-----------------------------------------------------------------------------------------------
  ! Grid interpolation:
  !-----------------------------------------------------------------------------------------------
  ! background: interpolating final grid background to the finest multigrid level
! HJ mod: was grid3. 6/10/2011
! HJ: I think the grid interpolation should be from STMAS_final to STMAS_mgrid(current level). 6/30/2011
! since bkgd is not used during minimization process. 
! But, STMAS_final%bkgd should be interpolated to STMAS_mgrid%bkgd(1:STMAS_final%nfinest) for later use, 
! the largest that near to the final grid.  8/8/2011

  !-----------------------------------------------------------------------------------------------
  ! Use background as initial guesses:
  !-----------------------------------------------------------------------------------------------

  CALL grid4_interpolation(STMAS_final%bkgd(1:STMAS_final%numgrid(1), &
                                            1:STMAS_final%numgrid(2), &
                                            1:STMAS_final%numgrid(3), &
                                            1:STMAS_final%numgrid(4), &
                                            1:STMAS_numvars),STMAS_final%numgrid, &
                           STMAS_mgrid%anal(1:STMAS_mgrid%numgrid(1), &
                                            1:STMAS_mgrid%numgrid(2), &
                                            1:STMAS_mgrid%numgrid(3), &
                                            1:STMAS_mgrid%numgrid(4), &
                                            1:STMAS_numvars),STMAS_mgrid%numgrid, &
                                            STMAS_numvars)

  !-----------------------------------------------------------------------------------------------
  ! topography: interpolates to the finest multigrid level:
  !-----------------------------------------------------------------------------------------------
  CALL grid2_interpolation(STMAS_final%topo(1:STMAS_final%numgrid(1), &
                                            1:STMAS_final%numgrid(2)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%topo(1:STMAS_mgrid%numgrid(1), &
                                            1:STMAS_mgrid%numgrid(2)), &
                           STMAS_mgrid%numgrid)

  ! land-water factor: interpolates to the finest multigrid level:
  CALL grid2_interpolation(STMAS_final%land(1:STMAS_final%numgrid(1), &
                                            1:STMAS_final%numgrid(2)), & 
                           STMAS_final%numgrid, &
                           STMAS_mgrid%land(1:STMAS_mgrid%numgrid(1), &
                                            1:STMAS_mgrid%numgrid(2)), &
                           STMAS_mgrid%numgrid)

  ! mapping factor: interpolates to the finest multigrid level:
  CALL grid2_interpolation(STMAS_final%mapf(1:STMAS_final%numgrid(1), &
                                            1:STMAS_final%numgrid(2)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%mapf(1:STMAS_mgrid%numgrid(1), &
                                            1:STMAS_mgrid%numgrid(2)), &
                           STMAS_mgrid%numgrid)
  ! Jacobin: jt1,jt2,jt3, grid5 is modified from grid3 to exclude nvars. HJ 6/30/2011
  CALL grid5_interpolation(STMAS_final%jt1(1:STMAS_final%numgrid(1), &
                                           1:STMAS_final%numgrid(2), &
                                           1:STMAS_final%numgrid(3)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%jt1(1:STMAS_mgrid%numgrid(1), &
                                           1:STMAS_mgrid%numgrid(2), &
                                           1:STMAS_mgrid%numgrid(3)), &
                           STMAS_mgrid%numgrid)
  CALL grid5_interpolation(STMAS_final%jt2(1:STMAS_final%numgrid(1), &
                                           1:STMAS_final%numgrid(2), &
                                           1:STMAS_final%numgrid(3)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%jt2(1:STMAS_mgrid%numgrid(1), &
                                           1:STMAS_mgrid%numgrid(2), &
                                           1:STMAS_mgrid%numgrid(3)), &
                           STMAS_mgrid%numgrid)
  CALL grid2_interpolation(STMAS_final%jt3(1:STMAS_final%numgrid(1), &
                                           1:STMAS_final%numgrid(2)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%jt3(1:STMAS_mgrid%numgrid(1), &
                                           1:STMAS_mgrid%numgrid(2)), &
                           STMAS_mgrid%numgrid)
  CALL grid2_interpolation(STMAS_final%cor(1:STMAS_final%numgrid(1), &
                                           1:STMAS_final%numgrid(2)), &
                           STMAS_final%numgrid, &
                           STMAS_mgrid%cor(1:STMAS_mgrid%numgrid(1), &
                                           1:STMAS_mgrid%numgrid(2)), &
                           STMAS_mgrid%numgrid)
! zz
  CALL grid1_interpolation(STMAS_final%zz(1:STMAS_final%numgrid(3)),STMAS_final%numgrid, &
                           STMAS_mgrid%zz(1:STMAS_mgrid%numgrid(3)),STMAS_mgrid%numgrid)

END SUBROUTINE get_initial


SUBROUTINE projection(mg,grid_from)
!doc==============================================================================================
!doc  This routine projects a multigrid to a finer grid.
!doc
!doc  History: October 2010 by Yuanfu Xie
!doc
!doc  Modified from the original projection with additional argument grid_from, and assign
!doc  grid_from to STMAS_mgrid first.  7/15/2011 Hongli Jiang
!doc==============================================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mg                      ! Multigrid level projected to
  TYPE(STMAS_bkgd) :: grid_from

  ! Local variables:
  CHARACTER*29 :: header= 'STMAS_analysis>projection: '
  INTEGER :: i,j,k,t,inc(4),i1,j1,k1,t1,i2,j2,k2,t2,im1,ip1,jm1,jp1,km1,kp1,tm1,tp1

  !-----------------------------------------------------------------------------------------------
  ! Finer grid number to be used:
  !-----------------------------------------------------------------------------------------------
  inc = 1      ! Increase or fix grid increment: Default not increase
! Memory is allocated in STMAS_module.f90/memo_alloc_mgrid. Only inc is defined here. HJ 7/8/2011 
  DO i=1,4
    IF (mg .GT. STMAS_numlevels-STMAS_maxlevels(i)+1) THEN
      IF (2*(STMAS_tmgrid%numgrid(i)-1)+1 .LE. STMAS_final%numgrid(i)) THEN  
        inc(i) = 2
      ENDIF
    ENDIF
  ENDDO

  print*,header
  print*,'inc',inc
  print*,'mgrid%incr',STMAS_mgrid%incr

  !-----------------------------------------------------------------------------------------------
  ! Gridspacing:
  !-----------------------------------------------------------------------------------------------
! for cwb domain nx=153, ny=149, projected to 129x129 resulting different gridspacing at multigrid level.
!
  STMAS_mgrid%gridspc(1:4) = STMAS_final%gridspc(1:4)* &
    (STMAS_final%numgrid(1:4)-1)/FLOAT(STMAS_mgrid%numgrid(1:4)-1)

  !-----------------------------------------------------------------------------------------------
  ! Distribute coarser grid increment function to finer one:
  !-----------------------------------------------------------------------------------------------
! Project coast grid onto finer grid first. HJ 7/15/2011
! At the coarsest level, e.g at 9x9, i1 (j1,k1,t1) index changes from 1,2,3.. 
! while the index of i,(j,k,t) for STMAS_mgrid%anal changes at inc spacing. HJ 8/8/2011  
! copy and subtract bkgd to produce increment only. Define i2,j2,k2,and t2 to find bkgd values
! at the coarse grid point. HJ 8/9/2011
  DO t=1,STMAS_mgrid%numgrid(4),inc(4)
    t1=(t-1)/inc(4)+1
    DO k=1,STMAS_mgrid%numgrid(3),inc(3)
      k1=(k-1)/inc(3)+1
      DO j=1,STMAS_mgrid%numgrid(2),inc(2)
        j1=(j-1)/inc(2)+1
        DO i=1,STMAS_mgrid%numgrid(1),inc(1)
          i1=(i-1)/inc(1)+1
! i2,j2,k2,t2 are defined and used for 
          i2=(i-1)*STMAS_mgrid%incr(1)+1
          j2=(j-1)*STMAS_mgrid%incr(2)+1
          k2=(k-1)*STMAS_mgrid%incr(3)+1
          t2=(t-1)*STMAS_mgrid%incr(4)+1
          STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) = grid_from%anal(i1,j1,k1,t1,1:STMAS_numvars)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
!
!  
! Bi-linear interpolation: Only the increment --- fill in 2,4,6,8... for anal after coping above.
! Bi-linear first, subtract bkgd before exit this subroutine. HJ 8/8/2011
!X:
  IF (inc(1) .EQ. 2) THEN
    DO t=1,STMAS_mgrid%numgrid(4),inc(4)
      DO k=1,STMAS_mgrid%numgrid(3),inc(3)
        DO j=1,STMAS_mgrid%numgrid(2),inc(2)
          DO i=2,STMAS_mgrid%numgrid(1),inc(1)
           im1 = i-1
           ip1 = i+1
           STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) =  &
                 0.5*(STMAS_mgrid%anal(im1,j,k,t,1:STMAS_numvars) + &
                      STMAS_mgrid%anal(ip1,j,k,t,1:STMAS_numvars))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Y:
  IF (inc(2) .EQ. 2) THEN
    DO t=1,STMAS_mgrid%numgrid(4),inc(4)
      DO k=1,STMAS_mgrid%numgrid(3),inc(3)
        DO j=2,STMAS_mgrid%numgrid(2),inc(2)
          DO i=1,STMAS_mgrid%numgrid(1)
            jm1 = j-1
            jp1 = j+1
            STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) = &
                  0.5*(STMAS_mgrid%anal(i,jm1,k,t,1:STMAS_numvars) + &
                       STMAS_mgrid%anal(i,jp1,k,t,1:STMAS_numvars))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Z:
  IF (inc(3) .EQ. 2) THEN
    DO t=1,STMAS_mgrid%numgrid(4),inc(4)
      DO k=2,STMAS_mgrid%numgrid(3),inc(3)
        km1 = k-1
        kp1 = k+1
        DO j=1,STMAS_mgrid%numgrid(2)
          DO i=1,STMAS_mgrid%numgrid(1)
            STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) = &  
                  0.5*(STMAS_mgrid%anal(i,j,km1,t,1:STMAS_numvars) + &
                       STMAS_mgrid%anal(i,j,kp1,t,1:STMAS_numvars))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF
! HJ add t. 7/12/2011
  ! t:
  IF (inc(4) .EQ. 2) THEN
! Fill in finer grid with coarse increment: incr*2
! HJ add t-loop. 6/10/2011
    DO t=2,STMAS_mgrid%numgrid(4),inc(4)
      tm1 = t-1
      tp1 = t+1
      DO k=1,STMAS_mgrid%numgrid(3)
        DO j=1,STMAS_mgrid%numgrid(2)
          DO i=1,STMAS_mgrid%numgrid(1)
            STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) = &  
                  0.5*(STMAS_mgrid%anal(i,j,k,tm1,1:STMAS_numvars) - &
                       STMAS_mgrid%anal(i,j,k,tp1,1:STMAS_numvars))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDIF
! 
! subtract bkgd from anal, and update lat, long, mapf, land, topo, and etc. HJ 8/8/2011
    DO t=1,STMAS_mgrid%numgrid(4)
      DO k=1,STMAS_mgrid%numgrid(3)
        DO j=1,STMAS_mgrid%numgrid(2)
          DO i=1,STMAS_mgrid%numgrid(1)
            i1=(i-1)*STMAS_mgrid%incr(1)+1
            j1=(j-1)*STMAS_mgrid%incr(2)+1
            k1=(k-1)*STMAS_mgrid%incr(3)+1
            t1=(t-1)*STMAS_mgrid%incr(4)+1
            STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) = &
                  STMAS_mgrid%anal(i,j,k,t,1:STMAS_numvars) + &
                  STMAS_final%bkgd(i1,j1,k1,t1,1:STMAS_numvars)
! 3d-arrays, doing it at t=1 only. HJ 8/4/2011
           if(t.eq.1)then
             STMAS_mgrid%jt1(i,j,k)  = STMAS_final%jt1(i1,j1,k1)
             STMAS_mgrid%jt2(i,j,k)  = STMAS_final%jt2(i1,j1,k1)
! 2d arrays doint it at k=1 only.
             if(k .eq. 1)then
               STMAS_mgrid%lat(i,j)  = STMAS_final%lat(i1,j1)
               STMAS_mgrid%lon(i,j)  = STMAS_final%lon(i1,j1)
               STMAS_mgrid%topo(i,j) = STMAS_final%topo(i1,j1)
               STMAS_mgrid%land(i,j) = STMAS_final%land(i1,j1)
               STMAS_mgrid%mapf(i,j) = STMAS_final%mapf(i1,j1)
               STMAS_mgrid%jt3(i,j)  = STMAS_final%jt3(i1,j1)
               STMAS_mgrid%cor(i,j)  = STMAS_final%cor(i1,j1)
             endif
             if(i.eq.1 .and. j.eq.1)then
               STMAS_mgrid%zz(k)  = STMAS_final%zz(k1)
             endif
           endif
          ENDDO
        ENDDO
      ENDDO
    ENDDO

END SUBROUTINE projection

