
PROGRAM TESTNORM
  USE lapack_subs
  USE normal_subs
  IMPLICIT NONE     
  Real :: Data(9) = (/ (/1.0, 2.0, 1.0/), (/2.0, 4.0, 2.0/), (/1.0, 2.0, 1.0/) /)
  Real, Dimension(3,3) :: Data2d
  Real, Dimension(3,3,3) :: Out2d
  Character(len=80) :: fstr

  Data2d = reshape(Data, (/ 3, 3 /))
  write(6,*) '';write(6,*) ''
  
  Call MakeNorm(Data2D, 3, 3, 1.0, 1.0, 1, 0, 1, 1, Out2D)
  
  fstr = '(a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a,f4.2,a)'

  write(6,fstr) '<',Out2D(1,1,1),' ',Out2D(1,1,2),' ',Out2D(1,1,3),'>, <', &
                    Out2D(2,1,1),' ',Out2D(2,1,2),' ',Out2D(2,1,3),'>, <', &
                    Out2D(3,1,1),' ',Out2D(3,1,2),' ',Out2D(3,1,3),'>'

  write(6,fstr) '<',Out2D(1,2,1),' ',Out2D(1,2,2),' ',Out2D(1,2,3),'>, <', &
                    Out2D(2,2,1),' ',Out2D(2,2,2),' ',Out2D(2,2,3),'>, <', &
                    Out2D(3,2,1),' ',Out2D(3,2,2),' ',Out2D(3,2,3),'>'

  write(6,fstr) '<',Out2D(1,3,1),' ',Out2D(1,3,2),' ',Out2D(1,3,3),'>, <', &
                    Out2D(2,3,1),' ',Out2D(2,3,2),' ',Out2D(2,3,3),'>, <', &
                    Out2D(3,3,1),' ',Out2D(3,3,2),' ',Out2D(3,3,3),'>'


  write(6,'(a)') ' '
  write(6,'(a)') 'Successful completion of TESTNORM!'
  write(6,'(a)') '-------------------------------------------------------------------------'
  write(6,'(a)') ' '
  
END PROGRAM TESTNORM
