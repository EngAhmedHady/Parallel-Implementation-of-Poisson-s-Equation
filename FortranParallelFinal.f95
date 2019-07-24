Program Poisson_Equation
USE MPI
Implicit none

real, dimension (:,:), allocatable  :: Grid
real, dimension (:)  , allocatable  :: Y,A
real, parameter                     :: pi = 3.1416
integer   :: ierr,R,S,P ! ----------------- Processors information
integer   :: Status(MPI_STATUS_SIZE)
integer   :: l,m          ! --------------- Grid Size (Number of nods)
integer   :: i,j,n,W,o,D   ! -------------- Problem Size
real      :: h              ! ------------- Delta x (dx) or Delta y (dy)
real      :: StartTime,EndTime
logical   :: file_exists

CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,R,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,S,ierr)
P = S - 1

IF(S < 3)THEN
   PRINT *,"The minimum valid number of processors is 3"
   PRINT *,"Abort the program ..."
   CALL MPI_FINALIZE(ierr)
   STOP
END IF

IF (R == 0) THEN
    call CPU_Time(StartTime)
    INQUIRE(FILE="GridParallel.dat", EXIST=file_exists) ! ------------ Transiant state file data
    if (file_exists) then
        open(3, file = "GridParallel.dat" , status = 'old') 
    else 
        open(3, file = "GridParallel.dat" , status = 'new')
    end if
    
    PRINT *, "Master Processor, Insert Grid size"
    READ  *, l
    allocate (Grid(l,l))
    allocate (Y(l))
    m = (l/P)
    D = 3                                       ! ------------ Domain Size
    h = D/(real(l) - 1); 
    
    if ((abs((real(l)/real(P)) - m)) > 0) print *, "Unbalanced (m = (l / (np-1)) is not integer) ",abs((real(l)/real(P)) - m)
    If ((m > 1) .and. (abs((real(l)/real(P)) - m) == 0)) then
        print *, "Allocated,  m :",m
        do i = 1, l
           Y(i) = (i-1) * h     ! (j) is the colmun nimber >> indection of x position
        end do 

        DO  o= 1 , P
            W = l * (m + 1)
            If (o > 2) W = l * (m + 2)
            
            allocate (A(W)) 
            
            A = Initiate(l,m,Y,o,W)
            
            CALL MPI_SEND(m,1,MPI_real,o,3,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(W,1,MPI_integer,o,0,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(l,1,MPI_integer,o,2,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(A,W,MPI_real,o,1,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(Y,l,MPI_real,o,4,MPI_COMM_WORLD,ierr)
            
            deallocate (A)
        END DO
        
        DO  o = 1 , P
        
            W = l * (m + 1)
            If (o > 2) W = l * (m + 2)
            
            allocate (A(W))

            CALL MPI_RECV(A,W,MPI_real,o,8,MPI_COMM_WORLD,Status,ierr)
            if (o == 1) then
                n = 1
                
                do i = 1, l
                    do j = 1, m
                        Grid(i,j) = A(n)
                        n = n + 1
                    end do
                    n = n + 1
                end do
            else if (o ==2) then
                n = 2
                do i = 1, l
                    do j = l-(m-1), l
                        Grid(i,j) = A(n)
                        n = n + 1
                    end do
                    n = n + 1
                end do
            else
                n = 2
                do i = 1, l
                    do j = (o-2)*m + 1, (o-1)*m
                        Grid(i,j) = A(n)
                        n = n + 1
                    end do
                    n = n + 2
                end do
            End if
            deallocate (A) 
        END DO
        
        do i = lbound(Grid,1), ubound(Grid,1)
            write(3,*) (Grid(i,j), j = lbound(Grid,2), ubound(Grid,2))
        end do
        call CPU_Time(EndTime)
        print *, "Simulation Time:",EndTime-StartTime, "Seconds"
        
        call system('gnuplot -p GridParallel.gnu')
    Else 
        DO  o= 1 , P
            CALL MPI_SEND(m,1,MPI_real,o,3,MPI_COMM_WORLD,ierr)
            CALL MPI_SEND(l,1,MPI_integer,o,2,MPI_COMM_WORLD,ierr)
        End Do
        PRINT *,"Insufficient Grid"
        PRINT *,"Abort the program ..."
        CALL MPI_FINALIZE(ierr)
        STOP
    End If

END IF

IF (R /= 0) THEN
    call CPU_Time(StartTime)
    CALL MPI_RECV(l,1,MPI_integer,0,2,MPI_COMM_WORLD,Status,ierr)
    CALL MPI_RECV(m,1,MPI_real,0,3,MPI_COMM_WORLD,Status,ierr)
    If ((m > 1) .and. (abs((real(l)/real(P)) - m) == 0)) then
        CALL MPI_RECV(W,1,MPI_integer,0,0,MPI_COMM_WORLD,Status,ierr)
        allocate (A(W))
        allocate (Y(l))
        CALL MPI_RECV(A,W,MPI_real,0,1,MPI_COMM_WORLD,Status,ierr)
        CALL MPI_RECV(Y,l,MPI_real,0,4,MPI_COMM_WORLD,Status,ierr)
        CALL Solution(l,P,R,A,m)
        CALL MPI_SEND(A,W,MPI_real,0,8,MPI_COMM_WORLD,ierr)
        call CPU_Time(EndTime)
        print *, R ,"Simulation Time:",EndTime-StartTime, "Seconds"
    Else 
        CALL MPI_FINALIZE(ierr)
        STOP
    End If
    
END IF

CALL MPI_FINALIZE(ierr)


contains
function Initiate(l,m,Y,o,W) 
    integer                    :: W,l,m                !                Matrix Size, Matrix Rows, Matrix Columns
    real, dimension (W)        :: Initiate !                            Main Matrix
    real, dimension (l)        :: Y !                                   coordinates matrix
    integer                    :: o                                !    Processor index
   
    n = 1
    If ( o == 1) then
        do i = 1, l
            do j = 1 , m + 1
                if ((j == 1) .and. (i < l)) then
                    Initiate(n) = sin(pi*y(i))
                else
                    Initiate(n) = 0
                end if
                n = n + 1
            end do 
        end do
    Else If ( o == 2) then
        do i = 1, l
            do j = 1 , m + 1
                if ((j == m + 1) .and. (i < l)) then
                    Initiate(n) = sin(pi*y(i)) * exp(2.0)
                else
                    Initiate(n) = 0
                end if
                n = n + 1
            end do 
        end do
    Else
        n = 1
        do i = 1, W
            Initiate(i) = 0
        end do
    End IF

end function Initiate


subroutine Solution(l,P,R,A,m)
Implicit none

real, dimension (W) :: A
real, dimension (l) :: FirstC, LastC, MidC1, MidC2!                       First column In The Matrix, First column In The Matrix
integer             :: l,m,P,R,o,x,oo
integer             :: iMax, CountParallel,CountParallel2,Count
real                :: eps, D
logical             :: LocalConvergence,GlobalConvergence

D = 3                                       ! ------------ Domain Size
h = D/(real(l) - 1); 

do x = 1, l
    LastC(x)  = 0
    FirstC(x) = 0
end do

iMax = 100000
eps  = 10.0**(-5)
Count = 0

If (R == 1) then
    print *, "solution on ",R
    LocalConvergence = .False.
    GlobalConvergence = .False.
    
    do while (GlobalConvergence .eqv. .false.)
        CountParallel = 0

        A = Solver(A,l,m,W,iMax,eps,R)
    
        do x = 0, l
            LastC(x+1) = A((x+1)*(m+1)-1)
        end do
    
        IF (P>2) THEN
            CALL MPI_SEND(LastC,l,MPI_real,3,6,MPI_COMM_WORLD,ierr)
            CALL MPI_RECV(FirstC,l,MPI_real,3,6,MPI_COMM_WORLD,Status,ierr)
        ELSE
            CALL MPI_SEND(LastC,l,MPI_real,2,6,MPI_COMM_WORLD,ierr)
            CALL MPI_RECV(FirstC,l,MPI_real,2,6,MPI_COMM_WORLD,Status,ierr)
        END IF
    
        LastC = FirstC

 
        do x = 1, l
            If (abs(LastC(x) - A(x*(m+1))) > 10.0**(-8)) Then
                A(x*(m+1)) = LastC(x)
            Else 
                CountParallel = CountParallel + 1
           End If
        end do
        
        Count = Count + 1
        
        If (CountParallel >= l) then
            LocalConvergence = .True.
        else
            LocalConvergence = .False.
        end If
        
        GlobalConvergence = LocalConvergence

        do oo = 2, P
            CALL MPI_SEND(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,ierr)
            CALL MPI_RECV(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,Status,ierr)
            If (GlobalConvergence .neqv. LocalConvergence) then
                GlobalConvergence = .False.
                LocalConvergence = .False.
            end if
        end do 
        
    End Do

Else if (R == 2) then
    print *, "solution on ",R
    LocalConvergence = .False.
    GlobalConvergence = .False.
    do while (GlobalConvergence .eqv. .false.)
        CountParallel = 0
        A = Solver(A,l,m,W,iMax,eps,R)
    
        do x = 0, l
            FirstC(x+1) = A(2+(x*(m+1)))
        end do
    
        IF (P>2) THEN
            CALL MPI_SEND(FirstC,l,MPI_real,P,6,MPI_COMM_WORLD,ierr)
            CALL MPI_RECV(LastC,l,MPI_real,P,6,MPI_COMM_WORLD,Status,ierr)
        ELSE
            CALL MPI_SEND(FirstC,l,MPI_real,1,6,MPI_COMM_WORLD,ierr)
            CALL MPI_RECV(LastC,l,MPI_real,1,6,MPI_COMM_WORLD,Status,ierr)
        END IF
        
        FirstC = LastC
        
        do x = 0, l
            If (abs(FirstC(x+1) - A((x*(m+1))+1)) > 10.0**(-8)) Then
                A((x*(m+1))+1) = LastC(x+1)
            Else 
                CountParallel = CountParallel + 1
           End If
        end do

        Count = Count + 1
        
        If (CountParallel >= l) then
            LocalConvergence = .True.
        else
            LocalConvergence = .False.
        end If
        
        GlobalConvergence = LocalConvergence
 
        do oo = 1, P
            If (oo /= 2) then
                CALL MPI_SEND(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,ierr)
                CALL MPI_RECV(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,Status,ierr)
                If (GlobalConvergence .neqv. LocalConvergence) then
                    GlobalConvergence = .False.
                     LocalConvergence = .False.
                end if
            End If
        end do
    End Do
    
Else

    Do o = 3, P
        If (R == o) then
            print *, "solution on ",R
            LocalConvergence = .False.
            GlobalConvergence = .False.
            
         do while (GlobalConvergence .eqv. .false.)
            CountParallel = 0
            
            A = Solver(A,l,m,W,iMax,eps,R)
            
            do x = 0, l
                FirstC(x+1) = A(2+(x*(m+2)))
                LastC(x+1)  = A((x+1)*(m+2)-1)
            end do
            
            MidC1 = LastC
            
            IF (o == 3) THEN
                CALL MPI_SEND(FirstC,l,MPI_real,1,6,MPI_COMM_WORLD,ierr)
                CALL MPI_RECV(LastC,l,MPI_real,1,6,MPI_COMM_WORLD,Status,ierr)
            ELSE
                CALL MPI_SEND(FirstC,l,MPI_real,R-1,6,MPI_COMM_WORLD,ierr)
                CALL MPI_RECV(LastC,l,MPI_real,R-1,6,MPI_COMM_WORLD,Status,ierr)
            End IF
            
            MidC2 = LastC
            LastC = MidC1
            
            IF (o == P) Then
                CALL MPI_SEND(LastC,l,MPI_real,2,6,MPI_COMM_WORLD,ierr)
                CALL MPI_RECV(FirstC,l,MPI_real,2,6,MPI_COMM_WORLD,Status,ierr)
            Else
                CALL MPI_SEND(LastC,l,MPI_real,R+1,6,MPI_COMM_WORLD,ierr)
                CALL MPI_RECV(FirstC,l,MPI_real,R+1,6,MPI_COMM_WORLD,Status,ierr)
            END IF
            
            LastC = FirstC
            FirstC = MidC2
            
            do x = 0, l
                If (abs(FirstC(x+1) - A((x*(m+2))+1))> eps) Then
                    A((x*(m+2))+1) = FirstC(x+1)
                Else 
                   CountParallel = CountParallel + 1
                End If
            end do
            
            do x = 1, l
                If (abs(LastC(x) - A(x*(m+2)))> 10.0**(-8)) Then
                    A(x*(m+2)) = LastC(x)
                Else 
                   CountParallel2 = CountParallel2 + 1
                End If
            end do
            
            If ((CountParallel > l) .and. (CountParallel2 > l)) then
                LocalConvergence = .True.
            else
                LocalConvergence = .False.
            end If
            
            GlobalConvergence = LocalConvergence
            
            Count = Count + 1
            
            do oo = 1, P
                If (oo /= o) then
                    CALL MPI_SEND(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,ierr)
                    CALL MPI_RECV(GlobalConvergence,1,MPI_logical,oo,7,MPI_COMM_WORLD,Status,ierr)
                    If (GlobalConvergence .neqv. LocalConvergence) then
                        GlobalConvergence = .False.
                        LocalConvergence = .False.
                    end if
                End If
            end do
          End Do
        End If
   End DO

End If

print *,"Parallel Convergence in:",R ,"is",LocalConvergence ,"Domain Decomposition iterations",Count
end subroutine

function Solver(A,l,m,W,iMax,eps,R) 
    integer                    :: W,l,m                !                Matrix Size, Matrix Rows, Matrix Columns
    real, dimension (W)        :: A,Solver,ANew, error !                Main Matrix, New Calculate Matrix, Difference Matrix
    !real, dimension (l,(W/l))  :: Test, TestNew                 !                (2D) printable Matrix
    integer                    :: R,i,j,k,n,nx,it,xt,ConElents
    real                       :: eps                  !                Convergance Criteria
    integer                    :: iMax
    logical                    :: Converged
    it = 1  ! -- The main loop iterations 
    ANew = A
    Converged = .False.
    k = 1
    
    do while ((it <= iMax) .and. (Converged .eqv. .false.))
       
        n = m
        nx = m + 1
        K = m + 3
        if (R > 2) then
           n = m + 1
           K = m + 4
           nx = m + 2
           xt = 1 + m*(R-2)
        end if
        
        if (R == 1)  xt = 2
        if (R == 2)  xt = 1 + m*(P-1)

        Do i = 2, l-1
            Do j = xt, n + xt - 2
                ANew(k) = (A(K+1) + A(K-1) + A(K+nx) + A(K-nx) + &
                          (h**2)*((pi**2) - 1)*exp(y(i))*sin(y(j)*pi))/4

                K = K + 1
            End Do
            K = K + 2
        End DO
        
        ConElents = 0
        do i = 1 , W
            error (i)  = abs(ANew(i)- A(i))
            if (error (i) <= eps) then
                ConElents = ConElents + 1
            End If
        end do
        

        If (ConElents == W) then
            Converged = .True.
        else
            Converged = .False.
        end If
    
        A = ANew
        it = it + 1
end do

    Solver = A
    
end function Solver

END Program Poisson_Equation

