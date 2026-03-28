SUBROUTINE STMAS_output

!doc================================================================================================
!doc  This routine writes out the analysis.
!doc
!doc  Combined from STMASFC_output and output_anal and modified for use in stmas4d. Hongli Jiang 8/24/2011
!doc================================================================================================

  USE STMAS

  IMPLICIT NONE

  ! Write out to bal files:
  CALL write_bal(STMAS_final%numgrid,STMAS_maxvars, &
                 STMAS_numvars,STMAS_varnames,STMAS_final%anal, &
                 STMAS_final%zz,LAPS_i4time,INT(STMAS_final%gridspc(4)), &
                 KAPPA,P00)

END SUBROUTINE STMAS_output
!
SUBROUTINE write_bal(numgrid,max_vars,numvar,varnam,analys,zz,i4time,lapsdt,kappa,p00)
!
!doc================================================================================================
!doc  This routine writes out the analyses into lw3 gridded
!doc  data in NetCDF format.
!doc
!doc  modified from write_lsx. Hongli Jiang 8/30/2011
!doc  Input:
!doc        numgrid:         Number of gridpoints in each direction
!doc        max_vars:        Maximum number of analyzed variables
!doc        numvar:          Number of analysis variables
!doc        i4time:          Current time frame in i4 format
!doc        lapsdt:          LAPS cycle time
!doc        varnam:          Analysis variable names
!doc        analys:          Analysis
!doc================================================================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: numgrid(4),max_vars,numvar,i4time,lapsdt
  real, INTENT(IN) :: kappa,p00
  CHARACTER*4, INTENT(IN) :: varnam(max_vars)
  REAL, INTENT(IN) :: analys(numgrid(1),numgrid(2),numgrid(3),numgrid(4),max_vars)
  REAL, INTENT(IN) :: zz(numgrid(3))

!doc============================
! Local variables:
  INTEGER  :: I,J,K,T,itm,istatus,i4t
  REAL :: dat(numgrid(1),numgrid(2),numgrid(3),max_vars)

! PARAMETERS FOR CONVERTING Q, P AND T TO RH (COPIED FROM lib/degrib/rrpr.F90:
  real, parameter :: svp1=611.2
  real, parameter :: svp2=17.67
  real, parameter :: svp3=29.65
  real, parameter :: svpt0=273.15
  real, parameter :: eps = 0.622
  real            :: tmp,sph,pval,ssh2,w_to_omega

  REAL          :: HT(numgrid(1),numgrid(2),numgrid(3)) ! 3D-Height
  REAL          :: RH(numgrid(1),numgrid(2),numgrid(3)) ! RELATIVE HUMIDITY
  CHARACTER*125 :: cmt(2)=(/'Temperature','Pressure   '/) ! lt1 COMMENTS
! --------------------
! USE LAPS WRITE BALANCED FIELD:

! Time frame to write out:
  DO itm = numgrid(4),max0(1,numgrid(4)-2),-1     ! Time frame
    i4t = i4time-(numgrid(4)-itm)*lapsdt         ! i4time corresponding to itm

!U3
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),1) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,1)
!V3
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),2) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,2)
!W3
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),3) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,3)
!P3
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),4) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,4)
!T3
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),5) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,5)
!SH
    dat(1:numgrid(1),1:numgrid(2),1:numgrid(3),6) = &
         analys(1:numgrid(1),1:numgrid(2),1:numgrid(3),itm,6)

    DO K=1,numgrid(3)
      DO J=1,numgrid(2)
        DO I=1,numgrid(1)
! Calculate RH FROM P, Q T:
          pval = dat(i,j,k,4)
          tmp = dat(i,j,k,5)
          sph = dat(i,j,k,6)
          RH(i,j,k) = 1.E2 * (pval*sph/(sph*(1.-eps) + eps))/(svp1*exp(svp2*(tmp-svpt0)/(tmp-svp3)))
!
! specific humidity is adjusted by Q_r (rain content) using ssh2 routine of LAPS:
! ../lib/ssh2.f/ssh2(p,td,t,t_ref). The output of SSH2 is in g/kg. HJ 8/31/2011
! Assume saturation: TD = T:
! The following call puts dat(i,j,k,6) to 1.0. HJ 11/2/2011
!          if (dat(i,j,k,6) .GT. 0.0) then
!            dat(i,j,k,6) = 0.001*SSH2(dat(i,j,k,4),dat(i,j,k,5),dat(i,j,k,5),-132.0)
!          endif
! converting potential Temp back to Temp in Kelvin. HJ 11/2/2011
          dat(i,j,k,5)=dat(i,j,k,5)*(dat(i,j,k,4)/p00)**kappa/(1.0+0.61*dat(i,j,k,6))
! converting w3 back to om in pa/s. Use the function w_to_omega for now. 
! Will replace it with omega=-w*g*rho. HJ 11/9/2011
!          dat(i,j,k,3)=w_to_omega(dat(i,j,k,3),dat(i,j,k,4))
        ENDDO
      ENDDO
    ENDDO
! HT
    HT(1:numgrid(1),1:numgrid(2),1:numgrid(3)) = 0.
    DO K=2,numgrid(3)
      DO J=1,numgrid(2)
        DO I=1,numgrid(1)
          HT(i,j,k) = HT(i,j,k-1) - 287.0/9.8*0.5*(dat(i,j,k,5)+dat(i,j,k-1,5))* &
                      (log(dat(i,j,k,4)) - log(dat(i,j,k-1,4)))
          if(i.eq.1 .and. j.eq.1)then
           print*,'k',k,dat(i,j,k,4),dat(i,j,k,5)
          endif
        enddo
      enddo
    enddo

!====================
!  Write to bal file:
!====================

!  CALL write_temp_anal_ht(i4t,numgrid(1),numgrid(2),numgrid(3),dat(1,1,1,5), &
!                     dat(1,1,1,4),cmt,istatus)

! HJ: ../balance/writeballaps/
! write_bal_laps_ht(i4time,ht,u,v,temp,om,rh,sh,imax,jmax,kmax,pres,istatus) 9/30/2011
!  CALL WRITE_BAL_LAPS_HT(i4t,dat(1,1,1,4),dat(1,1,1,1),dat(1,1,1,2),dat(1,1,1,5), &
!                     dat(1,1,1,3),RH,dat(1,1,1,6), &
!                     numgrid(1),numgrid(2),numgrid(3),zz,istatus)

  ENDDO

END SUBROUTINE write_bal
