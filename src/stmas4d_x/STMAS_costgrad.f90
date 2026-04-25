SUBROUTINE costfunc(mg,func,num_controls,mgrid)

!doc=============================================================================================
!doc This routine computes the cost function value for STMAS surface analysis.
!doc
!doc Cost function description:
!doc   It contains terms of observations and physical balance as weak constraints.
!doc   Observation: the mapped observation data structure.
!doc   balance:     horizontal momentum equations with terrain normal vector's 
!doc                vertical component as penalty weights.
!doc
!doc History: Jan. 2011 by Yuanfu Xie
!doc
!doc This routine is updated and computes the cost function value for STMAS 4D (space-time).
!doc
!doc   balance: updated. momentum equations (u,v,w) with terrain following vertical coordinate, 
!doc            continuity, thermodynamic, and moisture (water vapor) equation as constraints.
!doc
!doc   call added: cost_bal(func,num_controls,mgrid), done in STMAS_balcost.f90
!doc
!doc History: May. 2011, 6/14/2011,  by Hongli Jiang
!doc=============================================================================================

  USE STMAS, ONLY: STMAS_numlevels, STMAS_bkgd

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls,mg
  REAL, INTENT(INOUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

  func = 0.0d0
  mgrid%grdt = 0.0d0

  !--------------
  ! Observations:
  !--------------
! HJ comment out. 6/10/2011
!  CALL cost_obs(func,num_controls,mgrid)

  !-------------
  ! Backgrounds:
  !-------------

  !-------------
  ! Constraints:
  !-------------
! HJ mod: 5/12/2011
  CALL cost_bal(func,num_controls,mgrid)

  !------------
  ! Smoothings:
  !------------
! HJ mod: 8/24/2011
  if (mg .eq. STMAS_numlevels) CALL cost_smooth(func,num_controls,mgrid)

END SUBROUTINE costfunc
