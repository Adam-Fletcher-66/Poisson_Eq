! Name: Adam_Fletcher
! Date: 4/12/2022
! Purpose: Solve the Poisson equation for a 2-D steady state                  &
!          electromagnetics problem. The electrical potential of cell (25,50) &
!          is 1.73006161e-8 volts. It took 80 iterations to converged to this &
!          number. We could obtain a more precise value at this point if we   &
!          used the double precision function "d" instead of "e". We could    &
!          also have increased the parameters of our max change loop to run   &
!          more iterations which would converge at a more precise answer.     &

program poisson_eq
  implicit none


                               ! Variable Dictionary
  
  integer, parameter :: X=100                  ! Domain of x dimension
  integer, parameter :: Y=100                  ! Domain of y dimension
  real, parameter :: H=0.01                    ! Height of each cell
  real, parameter :: PI=3.141592654            ! Constant value of pi
  
  real :: U(0:Y+1, 0:X+1)                      ! Old elec potential array
  real :: Unew(0:Y+1, 0:X+1)                   ! New elec potential array
  real :: q(0:Y+1, 0:X+1)                      ! Charge density array
  real :: max_change=1.0                       ! Var of max change
  real :: w,z                                  ! Var of columns in .dat file 
  
  integer :: num_iter=0                        ! Iterations to converge to new EP
  integer :: i,j                               ! Index var
  integer :: lun1                              ! Log unit number var

  U = 0.0                                      ! Initialize EP old and new
  Unew = U

  q = 0.0                                      ! Initialize charge density at 2 &
  q(25,25) = -4.0                              ! points and at boundary cells
  q(75,75) = -4.0

  num_iter_loop: do while(max_change > 1.0e-5) 

     do i=1,Y
        do j=1,X                               ! Loop over interior cells
           Unew(i,j) = (U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1)-4.0*PI*(H**2)       &
                *q(i,j))/4.0
        enddo
     enddo

     max_change = maxval(abs(U-Unew))
    
     U = Unew                                  ! Make EP old the EP new

     num_iter = num_iter + 1                   ! Count number of iterations

  enddo num_iter_loop


  write(*,*) "Number of iterations:", num_iter          

  write(*,*) "Electrical potential at (25,50)=",U(25,50)
                                               ! Write the EP at the given cell                                        
 
  
  !==========================================================================

  open(newunit=lun1,file='poisson.dat',status='REPLACE')
                                               ! Open new file to write to
  
  z = -0.5*H                                   ! Initialize vertical coordinates
  do j=0,X+1                                   ! Loop over each column
     w = -0.5*H                                ! Initialize horizontal coordinate
     do i=0,Y+1                                ! Loop over each row
        write(lun1,*) w,z,Unew(i,j)
        w = w+h                                ! Increment horizontal coord
     enddo
     write(lun1,*) " "
     z = z+h                                   ! Increment vertical coord
  enddo
  close(unit=lun1)                             ! Close file

  
  stop 0                                       ! Terminate program execution
end program poisson_eq
