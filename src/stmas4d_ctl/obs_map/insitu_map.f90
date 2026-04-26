SUBROUTINE insitu_map
!doc
!doc==================================================================
!doc This routine maps insitu observations to a given multigrid level
!doc using a terrain, land-water and flow dependent scheme.
!doc 
!doc Input: structured background and insitu observation -
!doc        background: type(STMAS_bkgd)
!doc                    STMAS_mgrid (at the current multigrid level)
!doc                    STMAS_final (at the final grid);
!doc        observations: type(STMAS4d_insitu)
!doc                    insitu
!doc
!doc Output: a structured observation grid -
!doc         gridded obs: type(STMAS_gridded_obs)
!doc                      gridobs
!doc
!doc History: Yuanfu Xie Nov. 2011
!doc==================================================================
!doc

  USE STMAS

  IMPLICIT NONE

  INTEGER :: i,iv
  REAL    :: rdc(STMAS_final%numgrid(3),6)

  ! Land water influence reduction in vertical:  ! Temporary;
  DO i=1,STMAS_final%numgrid(3)
    rdc(i,1:6) = i
  ENDDO

  gridobs%numvars = observations%numvars
print*,'OBS: ',observations%value(1,4)

  CALL mapping(observations%numobs,observations%numvars, &
               STMAS_mgrid%numgrid,STMAS_final%numgrid, &
               observations%value,observations%error,observations%xyzt, &
               STMAS_final%bkgd,observations%bkgd,STMAS_radius, &
               1,STMAS_mgrid%land,rdc,STMAS_maxobs,observations%land, &
               STMAS_invalid, &
               gridobs%numobs,gridobs%ixyzt,gridobs%value,gridobs%error)

END SUBROUTINE insitu_map
