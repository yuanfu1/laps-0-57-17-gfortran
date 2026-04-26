!**********************************************************************************************!
! MODULE normal_subs                                                                           !
!**********************************************************************************************!
!                                                                                              !
!---   Purpose - Holds subroutines written for normal calculations                             !
!==============================================================================================!

MODULE normal_subs
IMPLICIT NONE
CONTAINS

!==============================================================================================!
! SUBROUTINE | MakeNorm(InField, DimX, DimY, Reach, WeightSwitch,                              !
!                       NormSwitch, BoundSwitch, OutField)                                     !
!==============================================================================================!
!---   Purpose - Takes in a field of height values on an xy grid and                           !
!                outputs a field of normal vectors of the same size.                           !
!                                                                                              !
!---   Inputs  - InField:      The input field of height values.                               !
!                Dim_x:        int_X dimension of input vector.                                !
!                dim_y:        int_Y Dimension of input vector.                                !
!                Dx:           real_Grid size in X (to dimensionalize).                        !
!                Dy:           real_Grid size in Y (to dimensionalize).                        !
!                Reach:        int_How many points away to consider height values.             !
!                WeightSwitch: Weighting strategy to use, 0:None  1: (1/distance^2)            !
!                NormSwitch:   Whether to output normalized vectors 0:No 1:Yes                 !
!                BoundSwitch:  Strategy to handle boundaries                                   !
!                              0:Zero Them  1:Calculate  2:Nearest Neighbor Fill 3:Flat Terrain!
!                                                                                              !
!---   Output  - OutField:     The output field of normal vectors (same size as InField)       !
!                                                                                              !
!---   Example - Call MakeNorm(Data2D, 256, 128, 1.0, 1.0, 1, 0, 1, 0, Out2D)                  !
!----------------------------------------------------------------------------------------------!
SUBROUTINE MakeNorm(InField, dim_X, dim_y, Dx, Dy, &
                    Reach, WeightSwitch, NormSwitch, BoundSwitch, OutField)
  Implicit None 
  Integer, Intent(In):: dim_x, dim_y, Reach, WeightSwitch, NormSwitch, BoundSwitch
  Integer:: XCount, YCount, DimCount, curx, cury, dim_r, oob
  Real, Intent(In):: Dx,Dy
  Real,Dimension(3):: NormVector, HoldVec
  Real,Dimension(dim_x, dim_y), Intent(In):: InField
  Real,Dimension(dim_x, dim_y, 3), Intent(Out):: OutField
  Real,Allocatable,Dimension(:):: InDatax, InDatay, InDataz, InWeight, newInWeight
  

  !-------------------------------------------------------------------------------------------!
  !-=A=- Set the InWeight Vector based on method choice.                                      !
  !-------------------------------------------------------------------------------------------!
  dim_r = ((Reach*2)+1)**2   !--- Find number of data points to be inputted from Reach
  Allocate(InDatax(dim_r))   !--- Set Up matrices used in norm calculation
  Allocate(InDatay(dim_r))   
  Allocate(InDataz(dim_r))   
  Allocate(InWeight(dim_r))  
  Allocate(newInWeight(dim_r))  

  If (WeightSwitch.eq.0) InWeight(:) = 1.0  !--- Uniform Weights 

  If (WeightSwitch.eq.1) then !--- Use inverse distance squared as weights
    Do DimCount = 1,dim_r
      curx = mod((DimCount-1),(Reach*2+1)) - Reach                       
      cury = mod((DimCount-1)/((Reach*2)+1),dim_x) - Reach    
      if (DimCount.ne.(dim_r/2)+1) InWeight(DimCount)=1.0/sqrt(float(curx*curx)+float(cury*cury))
      if (DimCount.eq.(dim_r/2)+1) InWeight(DimCount)=1.0 ! set the zero point to weight 1.0
    EndDo
  EndIf

  !-------------------------------------------------------------------------------------------!
  !-=B=- Loop through the matrix and calculate all normal vectors except for the boudaries.   !
  !-------------------------------------------------------------------------------------------!
  Do XCount = 1,dim_x  
  Do YCount = 1,dim_y  !--- X&YCount will mark the center point in the analysis
    newInWeight = InWeight
    oob = 0
    Do DimCount = 1,dim_r    !--- Now set the analysis point vectors used in FindNorm
      curx = (XCount-Reach) + mod((DimCount-1),(Reach*2+1))                       
      cury = (YCount-Reach) + mod((DimCount-1)/((Reach*2)+1),dim_x)     
      if (((curx.le.0).or.(cury.le.0)).or.((curx.gt.dim_x).or.(cury.gt.dim_y))) then
        oob = 1  !--- Out of Bounds!
        newInWeight(DimCount)=0.0
        InDatax(DimCount)= curx*dx
        InDatay(DimCount)= cury*dy
        InDataz(DimCount)= 0.0
      else
        InDatax(DimCount)= curx*dx
        InDatay(DimCount)= cury*dy
        InDataz(DimCount)= InField(curx,cury)
      endIf
    EndDo
    
    Call FindNorm(dim_r, InDatax, InDatay, InDataz, newInWeight, NormSwitch, NormVector)
    OutField(XCount,YCount,:) = NormVector(:) !--- Found the norm now save it.
    If ((BoundSwitch.eq.0).and.(oob.eq.1)) OutField(XCount,YCount,:) = 0.0

  EndDo
  EndDo
  
  Deallocate(InDatax)  !--- Tidy Up
  Deallocate(InDatay)
  Deallocate(InDataz)
  Deallocate(InWeight)  
  Deallocate(newInWeight)  

  !-------------------------------------------------------------------------------------------!
  !-=C=- Set boundaries according to the method chosen.                                       !
  !-------------------------------------------------------------------------------------------!

  If (BoundSwitch.eq.2) then  !--- Use nearest neighbor to fill the edges.
    Do YCount = 1+Reach, dim_y-Reach
    HoldVec(:) = OutField(Reach+1, YCount,:)
    Do XCount = 1, Reach                        !--- Left_Edge
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    Do YCount = 1+Reach, dim_y-Reach
    HoldVec(:) = OutField(dim_x-Reach, YCount,:)
    Do XCount = dim_x-Reach+1, dim_x              !--- Right_Edge
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    Do XCount = 1+Reach, dim_x-Reach             !--- Bottom Edge
    HoldVec(:) = OutField(XCount, Reach+1,:)
    Do YCount = 1, Reach
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    Do XCount = Reach+1, dim_x-Reach             !--- Top Edge
    HoldVec(:) = OutField(XCount, dim_y-Reach,:)
    Do YCount = dim_y-Reach+1, dim_y
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    HoldVec(:) = Outfield(Reach+1, dim_y-Reach, :)
    Do XCount = 1, Reach ; Do YCount = dim_y-Reach+1, dim_y           !--- Top Left Corner
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    HoldVec(:) = Outfield(dim_x-Reach, dim_y-Reach, :)
    Do XCount = dim_x-Reach+1, dim_x ; Do YCount = dim_y-Reach+1, dim_y !---  Top Right Corner
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    HoldVec(:) = Outfield(Reach+1, Reach+1, :)
    Do XCount = 1, Reach ; Do YCount = 1, Reach                     !--- Bottom Left Corner
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
    HoldVec(:) = Outfield(dim_x-Reach, Reach+1, :)
    Do XCount = dim_x-Reach+1, dim_x ; Do YCount = 1, Reach           !--- Bottom Right Corner
      OutField(XCount,YCount,:) = HoldVec(:)
    EndDo;EndDo    
  EndIf

  If (BoundSwitch.eq.3) then  !--- Use flat boundary conditions (normal = <0,0,1>)
    Do YCount = 1+Reach, dim_y-Reach ; Do XCount = 1, Reach           !--- Left_Edge
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do YCount = 1+Reach, dim_y-Reach ; Do XCount = dim_x-Reach+1, dim_x !--- Right_Edge
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = 1+Reach, dim_x-Reach ; Do YCount = 1, Reach           !--- Bottom Edge
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = Reach+1, dim_x-Reach ; Do YCount = dim_y-Reach+1, dim_y !--- Top Edge
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = 1, Reach ; Do YCount = dim_y-Reach+1, dim_y            !--- Top Left Corner
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = dim_x-Reach+1, dim_x ; Do YCount = dim_y-Reach+1, dim_y  !---  Top Right Corner
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = 1, Reach ; Do YCount = 1, Reach                      !--- Bottom Left Corner
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
    Do XCount = dim_x-Reach+1, dim_x ; Do YCount = 1, Reach            !--- Bottom Right Corner
      OutField(XCount,YCount,:) = [0.0,0.0,1.0]
    EndDo;EndDo    
  EndIf
 
  Return
End SUBROUTINE MakeNorm


!==============================================================================================!
! SUBROUTINE  |  FindNorm(numin, InDatax, InDatay, InDataz, InWeight, NormSwitch, OutNormal)   !
!==============================================================================================!
!---   Purpose - Finds the optimal 3D normal vector given a set of <numin> data                !
!                points, and optionally their weights.                                         !
!                                                                                              !
!---   Inputs  - numin:      Dimension of input vector (Number of data points considered)      !
!                InDatax:    x data point indexes (Has size numin)                             !
!                InDatay:    y data point indexes (Has size numin)                             !
!                InDataz:    z data point indexes (Has size numin)                             !
!                InWeight:   Weights for each data point (Has size numin)                      !
!                NormSwitch: Switch to turn on normalization                                   !
!                                                                                              !
!---   Output  - OutNormal:  The 3 dimensional normal vector                                   !
!                                                                                              !
!---   Example - Call FindNorm(3, [1,2,3], [1,2,1], [0,1,0], [1,1,1], 1, OutNormal)            !
!----------------------------------------------------------------------------------------------!
SUBROUTINE FindNorm(numin, InDatax, InDatay, InDataz, InWeight, NormSwitch, OutNormal)
  USE LAPACK_subs, ONLY : sgesv
  Implicit None
  Integer:: numin, NormSwitch, Loopcount, INFO
  Integer,Dimension(numin):: IPIV
  Real:: norm
  Real,Dimension(numin):: InDatax, InDatay, InDataz, InWeight
  Real,Dimension(3,3):: LeftSum
  Real,Dimension(3)::   OutNormal, RightSum   !--- Normal vector <A,B,C>

    LeftSum(:,:)=0.0
    RightSum(:) =0.0
    !--- Create Left and Right Matrices
    Do LoopCount=1,numin
      LeftSum(1,1) = LeftSum(1,1) + InWeight(Loopcount)*(InDatax(Loopcount)*InDatax(Loopcount))
      LeftSum(1,2) = LeftSum(1,2) + InWeight(Loopcount)*(InDatax(Loopcount)*InDatay(Loopcount))
      LeftSum(1,3) = LeftSum(1,3) + InWeight(Loopcount)*(InDatax(Loopcount))

      LeftSum(2,2) = LeftSum(2,2) + InWeight(Loopcount)*(InDatay(Loopcount)*InDatay(Loopcount))
      LeftSum(2,3) = LeftSum(2,3) + InWeight(Loopcount)*(InDatay(Loopcount))          
      
      LeftSum(3,3) = LeftSum(3,3) + InWeight(Loopcount)

      RightSum(1) = RightSum(1) + InWeight(Loopcount)*(InDatax(Loopcount)*InDataz(Loopcount))
      RightSum(2) = RightSum(2) + InWeight(Loopcount)*(InDatay(Loopcount)*InDataz(Loopcount))
      RightSum(3) = RightSum(3) + InWeight(Loopcount)*(InDataz(Loopcount))
    EndDo
    LeftSum(2,1) = LeftSum(1,2)
    LeftSum(3,1) = LeftSum(1,3)
    LeftSum(3,2) = LeftSum(2,3)

    Call sgesv(3,1,LeftSum,3,IPIV,RightSum,3,INFO)

    If (INFO.ne.0) then 
      write(6,*) 'Problem While Solving Equations... Info: ',INFO
      OutNormal=[0.0,0.0,0.0]
      Return
    EndIf

    OutNormal(1) = -RightSum(1)
    OutNormal(2) = -RightSum(2)
    OutNormal(3) =  1.0

    If (NormSwitch.eq.1) then
      norm = sqrt((OutNormal(1)*OutNormal(1)) + (OutNormal(2)*OutNormal(2)) + 1.0)
      OutNormal = OutNormal / norm
    EndIf

    Return
End SUBROUTINE FindNorm


END MODULE normal_subs
