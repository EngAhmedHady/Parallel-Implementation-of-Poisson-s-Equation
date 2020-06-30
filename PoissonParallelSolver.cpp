
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include "gnuplot.h"
#include "mpi.h"

using namespace std;

// ********** | Global Variables | ********
int TotalP;
static int P;
const float pi = 3.1416;
const int D = 3;
static float h;
MPI_Status status;

vector<double> initiate(int l,int m, vector<double> Y , int o, int W) 
{
    vector<double> LinearGrid;
    int n = 0;
    double pi = 3.1416; 
    if (o == 1)
    {
        for (int i = 0; i < l; ++i)
        {
            for (int j = 0; j < (m + 1); j++)
            {
                if ((j == 0) && (i < l-1)) LinearGrid.push_back(sin(pi * Y.at(i)));
                else LinearGrid.push_back(0);
            }
        }
    }
    else if (o == 2)
    {
        for (int i = 0; i < l; ++i)
        {
            for (int j = 0; j < (m + 1); j++)
            {
                if ((j == m) && (i < l-1)) LinearGrid.push_back(sin(pi * Y.at(i)) * exp(2));
                else LinearGrid.push_back(0);
                n++;
            }
        }
    }
    else
    {
        for (int i = 0; i < W; ++i) LinearGrid.push_back(0);
    }
    return (LinearGrid);
}

vector<double> solver(vector<double> A,vector<double> Y,int l,int m,int W,int R) 
{

    // Convergence Criteria
    double eps = pow(10,(-6)); // -- Maximum Error 
    int iMax = 100000; // ---------- Maximum Number of iterations
    double error;
    bool Converged = false;
    int it = 1; // ----------------- Iteration Counter

    vector<double> Anew;
    Anew.resize(W); 
    
    int k = 0; // ------------------ linear elements index
    int ConElents, n, nx, xt;
    double C,f;
    C = (pow(h, 2))*(pow(pi, 2) - 1);
    for (int i = 0; i < W; i++)	Anew[i] = A[i];
    while ((it <= iMax) && !Converged){
        n = m;
        nx = m + 1;
        k = m + 2;

        if (R > 2){
           n = m + 1;
           k = m + 3;
           nx = m + 2;
           xt = m*(R-2);
        }

        if (R == 1)  xt = 1;
        if (R == 2)  xt = m*(P-1);
        
        for (int i = 1; i < l-1; ++i){
            for (int j = xt; j < n + xt - 1; j++){
                f = exp(Y[i])*sin(Y[j]*pi);
                Anew[k] = (A[k+1] + A[k-1] + A[k+nx] + A[k-nx] + C * f)/4;
                k++;
            }
            k += 2;
        }
        
        ConElents = 0;
        for (size_t i = 0; i < W; i++)
        {
            error = abs(Anew[i]- A[i]);
            if (error <= eps) ConElents++;
        }
        
        if (ConElents == W) Converged = true;
        else Converged = false;

        for (int i = 0; i < W; i++)	A[i] = Anew[i];
        it++;
    }

    // if (R > 2) cout << R <<" Solved after \t" << it << endl;
    return (A);
}

vector<double> solution(int l, int R, int m, int W,vector<double> A,vector<double> Y)
{
    vector<double> FirstC, LastC, MidC1, MidC2;
    
    int CountParallel,CountParallel2,Count;
    bool LocalConvergence,GlobalConvergence = false;
    cout << "Solution on \t" << R << endl;
    float eps = pow(10,(-8));
    Count =0;

    if (R == 1)
    {
        FirstC.resize(l);
        while (!GlobalConvergence){
            CountParallel = 0;

            A = solver(A,Y,l,m,W,R);

            for (int i = 0; i < l; i++) LastC.push_back(A[(i+1)*(m+1)-2]);

            if (P > 2) 
            {
                MPI_Send(&LastC[0],l,MPI_DOUBLE,3,6,MPI_COMM_WORLD);
                MPI_Recv(&FirstC[0],l,MPI_DOUBLE,3,6,MPI_COMM_WORLD,&status);
            }
            else
            {
                MPI_Send(&LastC[0],l,MPI_DOUBLE,2,6,MPI_COMM_WORLD);
                MPI_Recv(&FirstC[0],l,MPI_DOUBLE,2,6,MPI_COMM_WORLD,&status);
            }

            for (int i = 0; i < l; i++)	LastC[i] = FirstC[i];

            for (int i = 1; i <= l; i++)
            {
                if (abs(LastC[i-1] - A[i*(m+1)-1]) > eps) A[i*(m+1)-1] = LastC[i-1];
                else CountParallel++;
            }	

            if (CountParallel >= l) LocalConvergence = true;
            else LocalConvergence = false;

            GlobalConvergence = LocalConvergence;
            Count++;
            
            for (size_t i = 2; i <= P; i++)
            {
                MPI_Send(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD);
                MPI_Recv(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD,&status);
                
                if (GlobalConvergence != LocalConvergence)
                {
                    GlobalConvergence = false;
                    LocalConvergence = false;
                }
            }
            LastC.clear();
        }
    } 
    else if (R == 2) 
    {
        LastC.resize(l);
        while (!GlobalConvergence){
            CountParallel = 0;
            A = solver(A,Y,l,m,W,R);

            for (int i = 0; i < l; i++) FirstC.push_back(A[1+(i*(m+1))]);

            if (P > 2) 
            {
                MPI_Send(&FirstC[0],l,MPI_DOUBLE,P,6,MPI_COMM_WORLD);
                MPI_Recv(&LastC[0],l,MPI_DOUBLE,P,6,MPI_COMM_WORLD,&status);
            }
            else
            {
                MPI_Send(&FirstC[0],l,MPI_DOUBLE,1,6,MPI_COMM_WORLD);
                MPI_Recv(&LastC[0],l,MPI_DOUBLE,1,6,MPI_COMM_WORLD,&status);
            }

            for (int i = 0; i < l; i++)	FirstC[i] = LastC[i];

            for (int i = 0; i < l; i++)
            {
                if (abs(FirstC[i] - A[i*(m+1)]) > eps) A[i*(m+1)] = LastC[i];
                else CountParallel++;
            }

            if (CountParallel >= l) LocalConvergence = true;
            else LocalConvergence = false;

            GlobalConvergence = LocalConvergence;
            Count++;

            for (size_t i = 1; i <= P; i++)
            {
                if ( i != 2)
                {
                    MPI_Send(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD);
                    MPI_Recv(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD,&status);

                    if (GlobalConvergence != LocalConvergence)
                    {
                        GlobalConvergence = false;
                        LocalConvergence = false;
                    }
                }
            }
            FirstC.clear();
        }
    }
    else
    {
        while (!GlobalConvergence){
            CountParallel = 0;
            CountParallel2 = 0;
            A = solver(A,Y,l,m,W,R);

            for (int i = 0; i < l; i++) FirstC.push_back(A[1+(i*(m+2))]);
            for (int i = 0; i < l; i++) 
            {
                LastC.push_back(A[(i+1)*(m+2)-2]);
                MidC1.push_back(A[(i+1)*(m+2)-2]);
            }

            if (R == 3) 
            {
                MPI_Send(&FirstC[0],l,MPI_DOUBLE,1,6,MPI_COMM_WORLD);
                MPI_Recv(&LastC[0],l,MPI_DOUBLE,1,6,MPI_COMM_WORLD,&status);
            }
            else
            {
                MPI_Send(&FirstC[0],l,MPI_DOUBLE,R-1,6,MPI_COMM_WORLD);
                MPI_Recv(&LastC[0],l,MPI_DOUBLE,R-1,6,MPI_COMM_WORLD,&status);
            }

            for (int i = 0; i < l; i++) 
            {
                MidC2.push_back(LastC[i]);
                LastC[i] = MidC1[i];
            }

            if (R == P) 
            {
                MPI_Send(&LastC[0],l,MPI_DOUBLE,2,6,MPI_COMM_WORLD);
                MPI_Recv(&FirstC[0],l,MPI_DOUBLE,2,6,MPI_COMM_WORLD,&status);
            }
            else
            {
                MPI_Send(&LastC[0],l,MPI_DOUBLE,R+1,6,MPI_COMM_WORLD);
                MPI_Recv(&FirstC[0],l,MPI_DOUBLE,R+1,6,MPI_COMM_WORLD,&status);
            }

            for (int i = 0; i < l; i++) 
            {
                LastC[i] = FirstC[i];
                FirstC[i] = MidC2[i];
            }
            

            for (int i = 0; i < l; i++)
            {
                if (abs(FirstC[i] - A[i*(m+2)]) > eps) A[i*(m+2)] = FirstC[i];
                else CountParallel++;

                if (abs(LastC[i] - A[(i+1)*(m+2)-1]) > eps) A[(i+1)*(m+2)-1] = LastC[i];
                else CountParallel2++;
            }

            if ((CountParallel >= l) && (CountParallel2 >= l))LocalConvergence = true;
            else LocalConvergence = false;

            GlobalConvergence = LocalConvergence;
            Count++;

            for (size_t i = 1; i <= P; i++)
            {
                if ( i != R)
                {
                    MPI_Send(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD);
                    MPI_Recv(&GlobalConvergence,1,MPI_BYTE,i,7,MPI_COMM_WORLD,&status);

                    if (GlobalConvergence != LocalConvergence)
                    {
                        GlobalConvergence = false;
                        LocalConvergence = false;
                    }
                }
            }
            FirstC.clear();
            LastC.clear();
            MidC1.clear();
            MidC2.clear();
        }       
    }
    cout <<"Parallel Convergence in: \t"<< R << "\t is \t"<< GlobalConvergence << "\t Domain Decomposition iterations \t"<< Count << endl;
    return(A);
}

int main(int argc, char **argv)
{ 
    double t; 
	t = clock();
    // =================| MPI INTIALIZATION |=================
    int R;
    MPI_Init(&argc, &argv); // ------------------------ Initialize the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &TotalP); // ----- Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &R); // ---------- Get the rank of the process
    P = TotalP - 1; // ------------------------------- Calculating Number of Slaves 
  
    //|Rule (1)|.......| Minimum number of processors is 3 |..
    if (TotalP < 3)
    {
        cout << "The minimum valid number of processors is 3 \n"
             << "Abort the program ... \n";
        MPI_Finalize();
        return 0;
    }

    // ==============| Main Prameters Defining |===============
    /* ------------------------------------------------\
	|l: Grid size (length and width for square domain) |
    |   Also the number of nodes for the real domain.  |
	|m: Size of the portion for each processor.        |
	|D: Domain Real Size (we considered the dimensions |
    |   of the selected case to be 3 X 3).             |
    |W: Linear Portion Size for Each Processor.        |
    |h: Step size in x,y directions (dx,dy).           |
	\-------------------------------------------------*/
    int l,m,W;   

    vector<double> Y; // ----------- Position value
    vector<double> A; // ----------- Portion Linear vector

    // ==============|| Master Processor Code || =============
    if (R == 0) 
    {
        // Ploting File INTIALIZATION
        std::ofstream MyGrid ("ParallelCGrid.dat");
        
        // Grid Size Request
	    cout << "Insert Grid Dimensions: ";
        cin >> l;

        // Calculating The main Parameters
        m = int(float(l)/float(P));
        h = D/(float(l) - 1);

        //--------------------------------
                
        //Sending the calculated values to the other processors
        for (size_t i = 1; i <= P; i++)
        {
            MPI_Send(&m,1,MPI_INT,i,3,MPI_COMM_WORLD);
            MPI_Send(&l,1,MPI_INT,i,2,MPI_COMM_WORLD);       
        }
        
        //|Rule (2)|.| Grid size Should be divisible on the slaves |..
        if ((abs((float(l)/float(P)) - m)) > 0)
            cout << "Unbalanced (m = (l / (np-1)) is not integer) " << abs((float(l)/float(P)) - m) << endl;

        if ((m > 1) && (abs((float(l)/float(P)) - m) == 0)) 
        {
            cout << "Solution initialized with m = \t" << m << endl;
            for (int i = 0; i < l; i++)  Y.push_back(double(i)*h);

            for (int o = 1; o <= P; o++)
            {
                W = l * (m + 1);
                if (o > 2) W = l * (m + 2);

                A = initiate(l,m,Y,o,W);

                MPI_Send(&W,1,MPI_INT,o,0,MPI_COMM_WORLD);
                MPI_Send(&A[0],W,MPI_DOUBLE,o,1,MPI_COMM_WORLD);
                MPI_Send(&Y[0],l,MPI_DOUBLE,o,4,MPI_COMM_WORLD);
            }
            A.clear();
            double** Grid = new double*[l];
	        for (int i = 0; i < l; ++i)	Grid[i] = new double[l];
            int k;
            for (int o = 1; o <= P; o++)
            {
                W = l * (m + 1);
                if (o > 2) W = l * (m + 2);
                A.resize(W);
                MPI_Recv(&A[0],W,MPI_DOUBLE,o,8,MPI_COMM_WORLD,&status);
                if (o == 1)
                {
                    k = 0;
                    for (size_t i = 0; i < l; i++)
                    {
                        for (size_t j = 0; j < m; j++)
                        {
                            Grid[i][j] = A[k];
                            k++;
                        }
                        k++;
                    }      
                }
                else if (o == 2)
                {
                    k = 1;
                    for (size_t i = 0; i < l; i++)
                    {
                        for (size_t j = l-m; j < l; j++)
                        {
                            Grid[i][j] = A[k];
                            k++;
                        }
                        k++;
                    }
                }
                else
                {
                    k = 1;
                    for (size_t i = 0; i < l; i++)
                    {
                        for (size_t j = m*(o-2); j < m*(o-1); j++)
                        {
                            Grid[i][j] = A[k];
                            k++;
                        }
                        k+=2;
                    }
                }
                A.clear();
            }
            // Print The Matrix
	        for (int i = 0; i < l; i++)
	        {
	    	    for (int j = 0; j < l; j++)	MyGrid << Grid[i][j] << "\t";
		        MyGrid << endl;
	        }
            gnuplot plot;
    	    plot("load 'ParallelCGrid.gnu'");      
        }
        else
        {
           cout << "Insufficient Grid \n Abort the program ... \n";
           MPI_Finalize();
           return 0;
        }
    }

    if (R != 0) 
    {
        MPI_Recv(&m,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
        MPI_Recv(&l,1,MPI_INT,0,2,MPI_COMM_WORLD,&status);

        if ((m > 1) && (abs((float(l)/float(P)) - m) == 0)) 
        {
            h = D/(float(l) - 1);
            MPI_Recv(&W,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
            A.resize(W); //Resizes the container (A) so that it contains ‘W’ elements.
            Y.resize(l); //Resizes the container (Y) so that it contains ‘l’ elements.
            MPI_Recv(&A[0],W,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);
            MPI_Recv(&Y[0],l,MPI_DOUBLE,0,4,MPI_COMM_WORLD,&status);
            A = solution(l,R,m,W,A,Y);
            MPI_Send(&A[0],W,MPI_DOUBLE,0,8,MPI_COMM_WORLD);
            A.clear();
        }
        else
        {
            MPI_Finalize();
            return 0;        
        }

    }

    t = double(clock()-t)/double(CLOCKS_PER_SEC);
    MPI_Finalize();  // ----------------------------- Finalize the MPI environment.
    cout << "Simulation Time: \t" << t << "\t Seconds in processor No. \t" << R <<endl;
	return 0;
}

