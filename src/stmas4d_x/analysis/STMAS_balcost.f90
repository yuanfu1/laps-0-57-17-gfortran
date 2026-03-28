SUBROUTINE COST_BAL(func,num_controls,mgrid)
!==============================================================================================
!  This routine calculates the constrains applied to the STMAS4D cost function and its gradients.
!
!  History: May. 2011, 6/14/2011 by Hongli Jiang at NOAA/ESRL/GSD.
!
!  Description:
!    Six equations include both horizontal and veritical momentum equations, continuity equation
!    thermodynamic equation, and moisture equation. The moisture equation is for Qv (water vapor) only.
!    Eventually the Qv equation will be replaced with or in addition to Qc (liquid water) equation
!    with source term (P_cond: condensation rate).  
!==============================================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls
  REAL,    INTENT(OUT) :: func
  TYPE(STMAS_bkgd) :: mgrid

! calculate x-momentum-grad
  call mx_costgrad(func,num_controls,mgrid)

! calculate y-momentum-grad
  call my_costgrad(func,num_controls,mgrid)

! calculate z-momentum-grad
  call mz_costgrad(func,num_controls,mgrid)

! calculate continuity eqn
  call cont_costgrad(func,num_controls,mgrid)

! calculate thermodynamic eqn
  call tv_costgrad(func,num_controls,mgrid)

! calculate moisture eqn (water vapor, qv)
  call qv_costgrad(func,num_controls,mgrid)

END SUBROUTINE COST_BAL
