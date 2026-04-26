SUBROUTINE OBSCOST
!====================================================================
!  This subroutine adds observation terms into the STMAS cost function.
!
! History: Yuanfu Xie Dec 2011
!====================================================================

  USE STMAS

  print*,'Number of gridded obs: ',gridobs%numobs

  CALL insitu_cost

END SUBROUTINE

