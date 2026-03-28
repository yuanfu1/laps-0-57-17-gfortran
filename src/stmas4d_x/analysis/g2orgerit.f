SUBROUTINE G2ORDERIT(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL 2-ORDER DERIVATIVE OF INTERIOR (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
    REAL X,Y
    REAL Z1,Z2,Z3
    REAL A,B,C
    X=Z2-Z1
    Y=Z3-Z2
    A=2.0/(X*X+X*Y)
    B=-2.0/(X*Y)
    C=2.0/(X*Y+Y*Y)
    RETURN
    END

