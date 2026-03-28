MODULE STMAS4d_obs

!doc==================================================================
!doc This module defines STMAS 4D analysis data structures
!doc
!doc History: June 2011 by Yuanfu Xie (originated)
!doc==================================================================

  INTEGER,PARAMETER :: MAX_insitu=10000

  REAL,PARAMETER :: STMAS4d_missing = 10.0e10

  TYPE STMAS4D_insitu
    INTEGER   :: num_obs
    CHARACTER :: station(MAX_insitu)*5
    REAL      :: obs_val(MAX_insitu,5)  ! U,V,P/H,T,SH
    REAL      :: obs_err(MAX_insitu,5)  ! U,V,P/H,T,SH
    REAL      :: xyzt(4,MAX_insitu)
  END TYPE STMAS4D_insitu

  ! Instantiation:
  TYPE(STMAS4d_insitu) :: insitu

END MODULE STMAS4d_obs
