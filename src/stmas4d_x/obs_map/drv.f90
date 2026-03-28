PROGRAM drv
!==========
! Testing
!==========

  USE mem_grid
  USE STMAS
  USE LAPS_ingest

  ! REAL :: rdlvl

  ! CALL get_LAPS_config(STMAS_final%numgrid,STMAS_final%lat,STMAS_final%lon, &
  !                     STMAS_final%topo,STMAS_final%land,STMAS_final%mapf,rdlvl)

  CALL insitu_map

END PROGRAM drv
