!====================================================================
! This routine reads in all observation datasets.
!
! History: 
!  Creation: Dec. 2009 by Yuanfu Xie
!====================================================================

module observations

! Process for reading observations:
public read_observation

! Variables for passing observation datasets:
public convention_obs, radial_wind, reflectivity

real, allocatable :: convention_obs(:,:), radial_wind(:,:), reflectivity(:,:)

include 'read_insitu.f90'

end module observations
