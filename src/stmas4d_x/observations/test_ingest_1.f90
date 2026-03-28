PROGRAM test_ingest

  USE mem_grid
  USE STMAS
  USE LAPS_ingest
  USE ingest_all

  IMPLICIT NONE

  CHARACTER :: header*20
  INTEGER :: k,istatus,nvars,max_sfcobs
  REAL    :: top,domain(2,4),spacings(4),rdlvl

  STMAS_time_window(1) = -7200
  STMAS_time_window(2) = 0

  top = 10000.0

  call read_namelist

  ! LAPS ingest:
  CALL STMAS_LAPS_input

  vertical = 'HEIGHT'
print*,'NUmber var: ',STMAS_numvars,STMAS_varnames
print*,'Background: ',STMAS_final%bkgd(1,1,1,1,1:6)

  ! Max surface obs in multiple frames:
  max_sfcobs = STMAS_maxobs*STMAS_final%numgrid(4)

  CALL ingest(1,i4time_sequence,STMAS_final%numgrid, &
              STMAS_final%lat,STMAS_final%lon, &
              STMAS_final%topo,top,STMAS_domain, &
              STMAS_numvars,STMAS_varnames,max_sfcobs)

  print*,'Number obs: ',observations%numobs
  print*,observations%xyzt(1:4,1:10)
  do k=1,10
  print*,'Station names: ',observations%stnames(k)
  enddo
  print*,'Numgrid: ',STMAS_final%numgrid,STMAS_maxobs

  STMAS_mgrid%numgrid = STMAS_final%numgrid
  allocate(STMAS_mgrid%land(STMAS_mgrid%numgrid(1),STMAS_mgrid%numgrid(2)))
  allocate(STMAS_mgrid%anal(STMAS_mgrid%numgrid(1),STMAS_mgrid%numgrid(2), &
                            STMAS_mgrid%numgrid(3),STMAS_mgrid%numgrid(4), &
                            STMAS_numvars))
  allocate(STMAS_mgrid%bkgd(STMAS_mgrid%numgrid(1),STMAS_mgrid%numgrid(2), &
                            STMAS_mgrid%numgrid(3),STMAS_mgrid%numgrid(4), &
                            STMAS_numvars))
  STMAS_mgrid%land = STMAS_final%land 
   STMAS_mgrid%bkgd = STMAS_final%bkgd 
   STMAS_mgrid%anal = STMAS_final%bkgd 

  do k=1,observations%numobs
    if (observations%xyzt(3,k) .le. 0) then
      print*,'Zero?'
      stop
    endif
  enddo

  CALL insitu_map

  CALL obscost

END PROGRAM test_ingest
