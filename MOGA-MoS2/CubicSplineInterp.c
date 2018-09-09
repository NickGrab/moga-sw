#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
	int i = 0;

	FILE *fp, *fw;
	fp = fopen("disp.txt" , "r");
	fw = fopen("interp.txt", "w");

	if (fp == NULL) { printf("Cannot open file to read.... exiting \n"); exit(1); }
	if (fw == NULL) { printf("Cannot open file to write.... exiting \n"); exit(1);}	

	int N = 6;
	double x[N], f[N], e[N], h[N], d[N], r[N], z[N];
	double df, dx, xs, x1, y;
	char buf[100];

	while (1)
	{
		fgets(buf, 100, fp);
		if (feof(fp)) {printf("End of file reached \n"); break;}
		sscanf(buf, "%lf %lf", &x[i], &f[i]);
		i++;
	}

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
		dx = (x[i+1] - x[i])/5.0;
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

	fclose(fp);
	fclose(fw);
	
	return 0;
}


