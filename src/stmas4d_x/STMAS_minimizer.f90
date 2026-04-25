SUBROUTINE STMAS_minimizer(num_controls,mg,mgrid)

!doc=======================================================================
!doc This routine minimizes a STMASFC cost function by using LBFGS_B method
!doc combining observations and constraints into its cost function.
!doc
!doc History: Dec. 2010 by Yuanfu Xie
!doc
!doc Modified to minimize a STMAS 4D cost function. June 2011, Hongli Jiang
!doc=======================================================================

  USE STMAS

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: num_controls,mg
  TYPE(STMAS_bkgd) :: mgrid

  ! Local variables:
  CHARACTER*19, PARAMETER :: header = 'STMAS_minimizer: '
  INTEGER :: i,j,k,nv,ic,t                            ! Looping variables
  INTEGER :: iteration                                ! Number of iterations
  REAL,ALLOCATABLE :: new_sol(:)                      ! Save the newest solution

  !** lbfgs_b variables:
  REAL :: fctn
  REAL,ALLOCATABLE :: controls(:),grad(:)
  INTEGER, PARAMETER :: msave=7                       ! max iterations saved
  
  CHARACTER*60 :: ctask,csave                         ! evaluation flag
  
  REAL,ALLOCATABLE :: wkspc(:)                        ! working space
  REAL,ALLOCATABLE :: bdlow(:)                        ! lower bounds
  REAL,ALLOCATABLE :: bdupp(:)                        ! upper bounds
  REAL :: factr,pgtol,dsave(29)

  integer :: iprint,isbmn,isave(44)
  integer,ALLOCATABLE :: nbund(:)                     ! Bound flags
  INTEGER,ALLOCATABLE :: iwkspc(:)

  LOGICAL :: lsave(4)

  !** End of LBFGS_B declarations.

  ! Guard against 32-bit integer overflow: wkspc size = nc*(2*msave+4)+O(msave^2)
  ! For large domains this can overflow INT32 and cause a silent bad ALLOCATE -> SIGSEGV.
  IF (INT(num_controls,8)*(2*msave+4) > INT(HUGE(num_controls),8)) THEN
    PRINT*,'FATAL: num_controls=',num_controls,' is too large for 32-bit LBFGSB workspace.'
    PRINT*,'Reduce domain size (nx,ny,nz,nt) before running STMAS on this grid.'
    STOP
  ENDIF

  ! Allocate the space:
  ALLOCATE(new_sol(num_controls),controls(num_controls),grad(num_controls), &
           bdlow(num_controls),bdupp(num_controls),nbund(num_controls), &
           wkspc(num_controls*(2*msave+4)+12*msave*msave+12*msave), &
           iwkspc(3*num_controls))

  PRINT*,''
  PRINT*,'minimization: ',mgrid%numgrid

  !------------------------------------------------------------------------
  ! Minimization:
  !------------------------------------------------------------------------
  ! LBFGS_B parameter setup:
  ctask = 'START'
  factr = 1.0
  pgtol = 1.0e-4
  iprint = 1
  isbmn = 1
  grad = 0.0

  ! Total number of control variables:
!HJ added numgrid(4). 6/10/2011
  IF (num_controls .NE. mgrid%numgrid(1)*mgrid%numgrid(2)* &
                        mgrid%numgrid(3)*mgrid%numgrid(4)*STMAS_numvars) THEN
    PRINT*,header,'Number of controls does not match. Check!'
    STOP
  ENDIF

  ! Simple bound constraints for controls (only low bound used here):
  nbund = 0 ! To be set!!! bund
  if (1 .EQ. 1) bdlow = 0.0
  
  ! Start of minimization:
  iteration = 0

  ! Grid to controls for initial guess:
  ic = 0

  DO nv=1,STMAS_numvars
!HJ add t-loop. 6/10/2011
    DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
    DO j=1,mgrid%numgrid(2)
    DO i=1,mgrid%numgrid(1)
      ic = ic+1
      controls(ic) = mgrid%anal(i,j,k,t,nv)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

1 CONTINUE     ! Minimization loop

  ! Before minimization, the initial guess is the newest solution:
  IF (iteration .EQ. 0) new_sol = controls
  !#########################################################################
  CALL SETULB(num_controls,msave,controls,bdlow,bdupp,nbund,fctn,grad, &
              factr,pgtol,wkspc,iwkspc,ctask,iprint,csave,lsave,isave,dsave)
  !#########################################################################
  ! Controls to grid:
  ic = 0
!HJ add t-loop. 6/10/2011
  DO nv=1,STMAS_numvars
    DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
    DO j=1,mgrid%numgrid(2)
    DO i=1,mgrid%numgrid(1)
      ic = ic+1
      mgrid%anal(i,j,k,t,nv) = controls(ic)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

! Exit if succeeds:
  IF (ctask(1:11) .EQ. 'CONVERGENCE') GOTO 2

! Function and gradient values are needed:
  IF (ctask(1:2) .EQ. 'FG') THEN
! Function and gradient values:
    CALL costfunc(mg,fctn,num_controls,mgrid)
!
!HJ mgrid%grdt to grad. 6/24/2011
    ic = 0
    DO nv=1,STMAS_numvars
    DO t=1,mgrid%numgrid(4)
    DO k=1,mgrid%numgrid(3)
    DO j=1,mgrid%numgrid(2)
    DO i=1,mgrid%numgrid(1)
      ic = ic+1
      grad(ic)=mgrid%grdt(i,j,k,t,nv)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
!
    GOTO 1
  ENDIF


! Exit if irregularity is encountered:
  IF ((ctask(1:2) .NE. 'FG') .AND. (ctask(1:5) .NE. 'NEW_X')) THEN
    WRITE(*,*) header,'Irregularity termination of LBFGSB'
    GOTO 2
  ENDIF

! A new approximation is found:
  IF (ctask(1:5) .EQ. 'NEW_X') THEN
    ! Save the newest solution:
    new_sol = controls
    iteration = iteration+1
    IF (iteration .LT. STMAS_iters) GOTO 1      ! Loop of the minimization
  ENDIF

! Exit of LBGGSB iteration:
2 WRITE(*,*) header,'Exit status: ',ctask(1:10)

! Save the newest solution to mgrid as the analysis of this level:
  IF (ctask(1:4) .NE. 'CONV') THEN
    ic = 0
    DO nv=1,STMAS_numvars
!HJ add t-loop. 6/10/2011
      DO t=1,mgrid%numgrid(4)
      DO k=1,mgrid%numgrid(3)
      DO j=1,mgrid%numgrid(2)
      DO i=1,mgrid%numgrid(1)
        ic = ic+1
        mgrid%anal(i,j,k,t,nv) = new_sol(ic)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDIF

  ! Deallocate the local LBFGS arrays:
  DEALLOCATE(new_sol,controls,grad, &
           bdlow,bdupp,nbund, &
           wkspc, &
           iwkspc)

END SUBROUTINE STMAS_minimizer
