
#include <iostream>
#include <cmath>
#include <fstream>
#include "gnuplot.h"


using namespace std;

int main()
{
	int m, n;
	cout << "Insert Grid Dimensions: ";
	cin >> m >> n;

	// == grid size ===
	/* -----------------------------------\
	| we considered the dimensions of the |
	| selected case to be 3 X 3 and the   |
	| step size in x direction (dx) is: h,|
	| step size in y direction (dy) is: k |
	| and aspect ratio (dx/dy) is: B      |
	\------------------------------------*/
	float h; float k; float B; float C; float pi;
	h = 3 / ((float)m - 1);
	k = 3 / ((float)n - 1);
	B = pow(((float)h / (float)k), 2);
	C = 2 * (1 + B);
	pi = 3.1416;
	cout << h <<"\t" << k << "\t" << B << "\t" << pi <<"\n";
	
	//==========|Iteration Control|============
	int It = 1; // -------------- Iteration loops counter
	int iMax = 100000; // ---- Maximum Number of iterations
	bool Converged = false;
	int ConvPos = 0;
	float eps = pow(10,(-5));
	// ========= |Intialize Data File| =======
	std::ofstream MyGrid ("CGrid.dat");
	
	/*Allocate 2D Matrix with dynamic dimension*/
	// ========== Main Grid =================
	double** Grid = new double*[m];
	for (int i = 0; i < m; ++i)
		Grid[i] = new double[n];

	double** NewGrid = new double*[m];
	for (int i = 0; i < m; ++i)
		NewGrid[i] = new double[n];

	double** Error = new double*[m];
	for (int i = 0; i < m; ++i)
		Error[i] = new double[n];
	// =======================================
	double* X = new double[n];
	double* Y = new double[m];


	/*Intializing (Filling) The Matrix with a value*/
	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; j++)	
			Grid[i][j] = 0;

	// Intializing (Filling) X,Y Vectors
	for (int i = 0; i < n; i++)
		X[i] = double(i)*h;

	for (int i = 0; i < m; i++) 
		Y[i] = double(i)*k;


	// Boundary conditions
	for (int i = 0; i < m - 1; i++)
	{
		Grid[i][n - 1] = sin(pi * 2) * exp(2);
		Grid[i][0] = sin(pi * 2);
	}

	cout << endl << endl;

	// =====| Finit Difference Solution ====
	double t; 
	t = clock();

	while ((It <= iMax) && !Converged)
	{
		ConvPos = 0;
		for (int i = 1; i < m-1; i++)
		{
			for (int j = 1; j < n-1; j++)
			{
				NewGrid[i][j] = (Grid[i + 1][j] + Grid[i - 1][j] +
								(pow(pi, 2) - 1) * exp(X[i]) * sin(Y[j]*pi) +
				                B * (Grid[i][j + 1] + Grid[i][j - 1])) / C;
			}
		}
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; j++)	
			{
				Error[i][j] = abs(NewGrid[i][j] - Grid[i][j]);
				Grid[i][j] = NewGrid[i][j];
				if (Error[i][j] <= eps) ConvPos++;
			}		
		}
		cout<<(double(It)/double(iMax))*100<<"\%"<<"\r";
		if (ConvPos == (m*n)) Converged = true;
		else Converged = false;
		It++;
	}

	t = double(clock()-t)/double(CLOCKS_PER_SEC);

    // Print The Matrix
	for (int i = 0; i <= m - 1; i++)
	{
		for (int j = 0; j <= n - 1; j++)
			MyGrid << Grid[i][j] << "\t";
		MyGrid << endl;
	}
	cout << endl;

	if (Converged)
		cout <<"Simulation Time: "<< t << " Seconds" << "\t Convergence achieved after: "<< It << " iterations \n";
	else
		cout <<"Simulation Time: "<< t << " Seconds" << "\t No convergence after: "<< It << " iterations \n";
	
	gnuplot plot;
	plot("load 'CGrid.gnu'");

	return 0;
}

