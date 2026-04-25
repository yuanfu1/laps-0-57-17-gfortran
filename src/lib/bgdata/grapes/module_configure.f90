!!--------------------------------------------------------------------------------------------------
! PROJECT           : GRAPES IO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Zhao Liu (liuzhao@nsccsz.cn), 2021/3/18, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # GRAPES IO Module
!!
!!  *This module defines data structures for GRAPES input namelist*
!! @author Sanshan Tu
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
MODULE module_configure
  USE NMLRead_m
  TYPE grid_config_rec_type
    INTEGER :: s_we
    INTEGER :: e_we
    INTEGER :: s_sn
    INTEGER :: e_sn
    INTEGER :: s_vert
    INTEGER :: e_vert
    LOGICAL :: global_opt
    integer :: spec_bdy_width
    INTEGER :: nh
    INTEGER :: iotype
    REAL*8:: xs_we
    REAL*8:: ys_sn
    REAL*8:: xd
    REAL*8:: yd
    REAL*4 :: cen_lat
    INTEGER :: num_soil_layers
    INTEGER :: start_year
    INTEGER :: start_month
    INTEGER :: start_day
    INTEGER :: start_hour
    INTEGER :: start_minute
    INTEGER :: start_second
  END TYPE grid_config_rec_type

CONTAINS

  SUBROUTINE initial_config(nmlst_fn, config_flags)

    IMPLICIT NONE
    TYPE(grid_config_rec_type) :: config_flags
    CHARACTER(*)        :: nmlst_fn

    call namelist_read(nmlst_fn, 's_we', config_flags%s_we)
    call namelist_read(nmlst_fn, 'e_we', config_flags%e_we)
    call namelist_read(nmlst_fn, 's_sn', config_flags%s_sn)
    call namelist_read(nmlst_fn, 'e_sn', config_flags%e_sn)
    call namelist_read(nmlst_fn, 's_vert', config_flags%s_vert)
    call namelist_read(nmlst_fn, 'e_vert', config_flags%e_vert)
    call namelist_read(nmlst_fn, 'spec_bdy_width', config_flags%spec_bdy_width)
    CALL namelist_read(nmlst_fn, 'xd', config_flags%xd)
    CALL namelist_read(nmlst_fn, 'yd', config_flags%yd)
    CALL namelist_read(nmlst_fn, 'cen_lat', config_flags%cen_lat)
    CALL namelist_read(nmlst_fn, 'start_year', config_flags%start_year)
    CALL namelist_read(nmlst_fn, 'start_month', config_flags%start_month)
    CALL namelist_read(nmlst_fn, 'start_day', config_flags%start_day)
    CALL namelist_read(nmlst_fn, 'start_hour', config_flags%start_hour)
    CALL namelist_read(nmlst_fn, 'num_soil_layers', config_flags%num_soil_layers)
    CALL namelist_read(nmlst_fn, 'xs_we', config_flags%xs_we)
    CALL namelist_read(nmlst_fn, 'ys_sn', config_flags%ys_sn)
    CALL namelist_read(nmlst_fn, 'nh', config_flags%nh)

  END SUBROUTINE initial_config

END MODULE module_configure
