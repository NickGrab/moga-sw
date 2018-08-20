#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void PartialPivoting(int N, double A[][N],double *B)
{
	int lambda[N], i, rho, k, j;
	double maxRhoVal = 0.0, temp1 = 0.0, temp2 = 0.0, l[N];

	for (i = 0; i < N; i++)
		lambda[i] = i;

	for (k = 0; k < N; k++)
	{
		maxRhoVal = A[k][k];

		// Pivoting step
		// 1. Find largest element in a column
		for (i = k; i < N; i++)
		{
			if (A[i][k] >= maxRhoVal) { maxRhoVal = A[i][k];  lambda[k] = i;}
		}

		// 2. Exchange row number "rho" and row number "k"
		if (lambda[k] != k)
		{
			for (j = 0; j < N; j++)
			{
				temp1 = A[k][j];
				A[k][j] = A[lambda[k]][j];
				A[lambda[k]][j] = temp1;
			}
			temp2 = B[k];
			B[k] = B[lambda[k]];
			B[lambda[k]] = temp2;
		}

		//printf("For row number %d \n", k);
		//printf("MaxValue = %.4f ", maxRhoVal);
		//printf("Exchanging row number %d with row number %d\n", k, lambda[k]);

		//printf("Row number %d is now %.4f %.4f %.4f %.4f %.4f\n", k, A[k][0], A[k][1], A[k][2], A[k][3], B[k]);
		//printf("Row number %d is now %.4f %.4f %.4f %.4f %.4f\n", lambda[k], A[lambda[k]][0], A[lambda[k]][1], A[lambda[k]][2], A[lambda[k]][3], B[lambda[k]]);

		// Elimination step
		if (A[k][k] == 0.0) continue;

		//printf("\n%.4f\n\n", A[k][k]);

		for (i = k+1; i < N; i++) 
		{
			if (A[i][k] == 0.0) continue;
	
			l[k] = A[i][k]/A[k][k];
			A[i][k] = A[i][k] - l[k]*A[k][k];

			B[i] = B[i] - l[k]*B[k];
			for (j = k+1; j < N; j++)
				A[i][j] = A[i][j] - l[k]*A[k][j];
		}

		//printf("Row number %d is now %.4f %.4f %.4f %.4f %.4f\n", k, A[k][0], A[k][1], A[k][2], A[k][3], B[k]);
                //printf("Row number %d is now %.4f %.4f %.4f %.4f %.4f\n", lambda[k], A[lambda[k]][0], A[lambda[k]][1], A[lambda[k]][2], A[lambda[k]][3], B[lambda[k]]);
		
	}

	return;
}

void BackSubstitution(int N, double A[][N], double *B, double *X)
{
	int i = 0, j = 0;
	double sum = 0.0;


	for (i = N-1; i >= 0; i--)
	{
		sum = 0.0;
		for (j = i+1; j < N; j++ ) sum += A[i][j]*X[j];
		
		X[i] = (B[i] - sum)/A[i][i];
	}
	return;
}

void TriDiagonalMatrixSol(int N, double *a, double *b, double *c, double *d)
{
	int i = 0;
	for (i = 0; i < N; i++)
	{
		if (i == 0)
		{
			c[0] = c[0]/b[0];
			d[0] = d[0]/b[0];
		}	
		else if (i == N-1)
		{
			d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
		}
		else
		{
			c[i] = c[i]/(b[i] - a[i]*c[i-1]);
			d[i] = (d[i] - a[i]*d[i-1])/(b[i] - a[i]*c[i-1]);
		}
	}
	
	return;
}

void SolveX(int N, double *a, double *b, double *c, double *d, double *X)
{
	int i=0;
	X[N] = d[N-1];
	for (i = N-2; i >= 0; i--)
		X[i+1] = d[i] - c[i]*X[i+2];
}

int main()
{
	int i = 0, j = 0;
	//double A[4][4], B[4], X[4];

	//A[0][0] = 0.02; A[0][1] = 0.01; A[0][2] = 0.00; A[0][3] = 0.00; B[0] = 0.02;
        //A[1][0] = 1.00; A[1][1] = 2.00; A[1][2] = 1.00; A[1][3] = 0.00; B[1] = 1.00;
	//A[2][0] = 0.00; A[2][1] = 1.00; A[2][2] = 2.00; A[2][3] = 1.00; B[2] = 4.00;
	//A[3][0] = 0.00; A[3][1] = 0.00; A[3][2] = 100.0;A[3][3] = 200.0;B[3] = 800.0;

	//A[0][0] = 3.00; A[0][1] = -1.0; A[0][2] = 0.00; A[0][3] = 0.00; B[0] = 2.00;
        //A[1][0] = -1.0; A[1][1] = 3.00; A[1][2] = -1.0; A[1][3] = 0.00; B[1] = 1.00;
        //A[2][0] = 0.00; A[2][1] = -1.0; A[2][2] = 3.00; A[2][3] = -1.0; B[2] = 1.00;
        //A[3][0] = 0.00; A[3][1] = 0.00; A[3][2] = -1.0; A[3][3] = 3.00; B[3] = 2.00;


	FILE *fp, *fw;
	fp = fopen("disp.txt" , "r");
	fw = fopen("interp.txt", "w");

	if (fp == NULL) { printf("Cannot open file.... exiting \n"); exit(1); }
	
	int N = 6;
	double x[N], f[N], e[N], h[N], d[N], r[N], z[N];
	double df, dx, xs, x1, y;
	char buf[100];

	while (1)
	{
		fgets(buf, 100, fp);
		if (feof(fp)) {printf("End of file reached \n"); break;}
		sscanf(buf, "%lf %lf\n", &x[i], &f[i]);
		i++;
	}
	fclose(fp);

	N = N-1;
	for (i = 0; i <= N-1; i++)
	{
		if (i < N-1)
		{
			e[i] = 2.0*(x[i+2] - x[i]);
			h[i] = x[i+2] - x[i+1];
			d[i] = x[i+1] - x[i];
			r[i] = 6*(f[i+2] - f[i+1])/h[i] - 6*(f[i+1] - f[i])/d[i];
			 printf("%10.4f %10.4f %10.4f %10.4f \n", e[i], d[i], h[i], r[i]);
		}
		else if (i == N-1)
		{
			d[i] = x[i+1] - x[i];
		}
	}	

	for (i = 0; i <= N-3; i++)
	{
		df = d[i]/e[i];
		e[i+1] = e[i+1] - df*h[i];
		r[i+1] = r[i+1] - df*r[i];
	}

	for (i = N-3; i >= 0; i--)
	{
		df = h[i]/e[i+1];
		r[i+1] = r[i] - r[i+1]*df;
	}

	for (i = 0; i <= N-2 ; i++)
		z[i+1] = r[i]/e[i];
	z[0] = 0.0; z[N] = 0.0;
	
	for (i = 0; i < N; i++)
	{
		dx = (x[i+1] - x[i])/100.0;
		xs = x[i] + 0.1;
		for (x1 = xs; x1 < x[i+1]; x1 = x1+dx)
		{
			y = -z[i]*pow( (x1-x[i+1]),3 )/(6.0*d[i]) +
			   z[i+1]*pow( (x1-x[i]),3 )/(6.0*d[i]) +
			 (f[i+1]/d[i] - z[i+1]*d[i]/6.0)*(x1 - x[i]) +
			 (-f[i]/d[i] + z[i]*d[i]/6.0)*(x1 - x[i+1]);
		
			fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
		}
	}

	fclose(fw);
	//printf("%10.4f %10.4f %10.4f %10.4f\n", X[0], X[1], X[2], X[3]);
	return 0;
}


