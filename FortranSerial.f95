Program Poisson_Equation
Implicit none
real, dimension(:,:), allocatable    :: Grid, GridNew,error
real,dimension (:), allocatable                  :: X,Y ! ---- Coordinate
integer                   :: l,i,j, it
real                      :: h,k,B !------- h (x increment), k (y increament), B (B = h/k)
real                      :: C ! ---------- C (Coefficient)
integer                   :: iMax,ConElents ! ------- The Max Number of Iterations
real                      :: eps!, ier  ! - eps (Stoping ceriateria)
real                      :: StartTime,EndTime
real, parameter           :: pi = 3.1416
logical                   :: file_exists, Converged


INQUIRE(FILE="Grid.dat", EXIST=file_exists) ! ------------ Transiant state file data
if (file_exists) then
    open(3, file = "Grid.dat" , status = 'old') 
else 
    open(3, file = "Grid.dat" , status = 'new')
end if

l = 100
allocate (Grid(l,l))
allocate (GridNew(l,l))
allocate (error(l,l))
allocate (X(l))
allocate (Y(l))

h = 3/(real(l)-1); k = 3/(real(l)-1); 

print *, h

it = 1  ! -- The main loop iterations 
iMax = 100000
eps  = 10.0**(-5)

do j=1,l
   do i=1,l 
       Grid(i,j) = 0
   end do 
   X(j) = (j-1)*h     ! (i) is the row number    >> indection of y position
   Y(j) = (j-1)*k     ! (j) is the colmun number >> indection of x position
end do 

! Boundary Condition
do i =1,l-1
    Grid(i,l) = sin(pi*y(i)) * exp(2.0)
    Grid(i,1) = sin(pi*y(i))
end do



! ===== Initialize coefficients ======
GridNew = Grid
B = (h/k)**2
C = 2 * (1 + B) 
Converged = .False.


call CPU_Time(StartTime)
! The main Loop for iterations

do while ((it <= iMax) .and. (Converged .eqv. .false.))
    do i = 2 , l-1
        do j = 2 , l-1
            GridNew(i,j) = (Grid(i+1,j) + Grid(i-1,j) + &
                           (h**2)*((pi**2) - 1)*exp(y(i))*sin(y(j)*pi) + &
                           Grid(i,j+1)+Grid(i,j-1))/4
        end do
    End do 
    ConElents = 0
    do i = 1 , l
      do j = 1 , l
         !Grid(i,j) = exp((j-1)*h) * sin((i-1)*k*pi)
         error (i,j)  = abs(GridNew(i,j)- Grid(i,j))
         if (error (i,j) <= eps) then
             ConElents = ConElents + 1
         End If
      end do
    end do
    If (ConElents == l**2) then
      Converged = .True.
    else
      Converged = .False.
    end If
    Grid = GridNew
    IF(MOD(it,5).EQ.0) WRITE(6,41,ADVANCE='NO') "Loading: ", (real(it)/real(iMax))*100, ' %'//CHAR(13)
      41  FORMAT(A, f 6.2, A)
    it = it + 1
end do

do i = lbound(Grid,1), ubound(Grid,1)
   write(3,*) (Grid(i,j), j = lbound(Grid,2), ubound(Grid,2))
   print *, (Grid(i,j), j = lbound(Grid,2), ubound(Grid,2))
end do

call system('gnuplot -p Grid.gnu')

 
call CPU_Time(EndTime)

If (Converged) then
    print *, "Simulation Time:",EndTime-StartTime, "Seconds", "       Convergence achieved after:",it, "iterations"
else
    print *, "Simulation Time:",EndTime-StartTime, "Seconds", "       No convergence after",it, "iterations"
end If

deallocate (Grid)
deallocate (GridNew)
deallocate (error)

END Program Poisson_Equation


