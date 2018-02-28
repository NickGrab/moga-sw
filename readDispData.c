#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Range
{
	double min_x, min_f, max_x, max_f;

};

void ReadData(int numBands, double x_vasp_data[][1000], double f_vasp_data[][1000], double x_gulp_data[][1000], double f_gulp_data[][1000])
{
	int regionsPerBand = 3, pointsPerRegion_vasp = 21; 

	int points_GM_gulp = 183, points_MK_gulp = 106, points_KG_gulp = 211;

	//double x_vasp_data[numBands][1000], f_vasp_data[numBands][1000];

	//double x_gulp_data[numBands][1000], f_gulp_data[numBands][1000];

	double temp_x, temp_y;

	double phonon_shift[3] = {0, 183, 289};
	double conv_x = 1000, conv_y = 33.35641;
	int band_offset = 26, band_segm_length = 21, gulp_offset = 3;

	int i = 0, j = 0, k = 0, lineNum_vasp = 0, lineNum_gulp = 0, bands_vasp = 0, bands_gulp = 0;

	const char *phonon_file[3] = {"GULP/G_M/phonon.disp", "GULP/M_K/phonon.disp", "GULP/K_G/phonon.disp"};
	const char *band_file = "DFT/band.dat";
	char lines_gulp[200], lines_band[200], lines_band_ph[200];

	FILE *vasp, *gulp[3];

	/* Reading Vasp File band.dat */
	vasp = fopen(band_file, "r");
	if (vasp == NULL) { printf("Unable to open file %s for reading dispersion data .... exiting \n", band_file); exit(1); }
	
	while (lineNum_vasp < band_offset) { fgets(lines_band, 200, vasp); lineNum_vasp++; }

	for (i = 0; i < numBands; i++)
	{
		for (j = 0; j < regionsPerBand; j++)
		{
			fgets(lines_band, 200, vasp);
			lineNum_vasp++;
			while(strcmp(lines_band,"\n") != 0)
			{
				sscanf(lines_band, "%lf %lf", &temp_x, &temp_y);
				x_vasp_data[i][k /*+ i*regionsPerBand*band_segm_length*/] = temp_x;
				f_vasp_data[i][k /*+ i*band_segm_length*/] = temp_y;
				//printf("i %d j %d = %10.6f %10.6f \n", i + 1, k + 1 + i*regionsPerBand*band_segm_length,
				//	x_vasp_data[i][k + i*regionsPerBand*band_segm_length], f_vasp_data[i][k + i*band_segm_length]);

				k++;
				fgets(lines_band, 200, vasp);
				lineNum_vasp++;
			}
			//printf("Line Number = %d \n", lineNum_vasp);
		}
		fgets(lines_band, 200, vasp);
		lineNum_vasp++;
		k=0;
	}

	//printf("Total lines read = %d \n", lineNum_vasp);
	/* Reading Vasp file complete*/

	/* Reading Gulp file */
	int curPoints = 0, prevPoint = 0;
	for (i = 0; i < regionsPerBand; i++)
        {
		lineNum_gulp = 0;
                gulp[i] = fopen(phonon_file[i], "r");
                if (gulp[i] == NULL) { printf("Unable to open file %s for reading dispersion data .... exiting \n", phonon_file[i]); exit(1); }
		
		while (lineNum_gulp < gulp_offset) {fgets(lines_gulp, 200, gulp[i]); lineNum_gulp++; }
		
		if (i == 0) curPoints = points_GM_gulp;
		else if (i == 1) curPoints = points_MK_gulp;
		else if (i == 2) curPoints = points_KG_gulp;
		
		for (j = 0; j < curPoints; j++)
		{
			for (k = 0; k < numBands; k++)
			{
				fgets(lines_gulp, 200, gulp[i]);
				sscanf(lines_gulp, "%lf %lf", &temp_x, &temp_y);
				x_gulp_data[k][prevPoint + j] = temp_x + prevPoint;
				f_gulp_data[k][prevPoint + j] = temp_y;		
				//printf("File = %s, i %d j %d, x = %4d y = %10.6f\n", phonon_file[i], k+1, prevPoint + j + 1, 
				//				x_gulp_data[k][prevPoint + j], f_gulp_data[k][prevPoint + j]);
				printf("File = %s, i = %d, j = %d, %10.4f %10.6f\n", phonon_file[i], k+1, prevPoint + j + 1, 
					x_gulp_data[k][prevPoint + j], f_gulp_data[k][prevPoint + j]);
			}	
		}

		prevPoint += curPoints;
        }	

	fclose(vasp); fclose(gulp[0]); fclose(gulp[1]); fclose(gulp[2]);

	/*FILE *fw1 = fopen("gulp.txt","w"), *fw2 = fopen("vasp.txt", "w");

	for(i = 0; i < 500; i++)
		fprintf(fw1, "%10.6f %10.6f \n", x_gulp_data[0][i], f_gulp_data[0][i]);

	for(i = 0; i < regionsPerBand*pointsPerRegion_vasp; i++)
		fprintf(fw2, "%10.6f %10.6f \n", x_vasp_data[0][i]*1000.0, f_vasp_data[0][i]*33.35641);


	fclose(fw1); fclose(fw2);	*/

	return;
}

double CubicSplineInterp(int N, double *x, double *f, struct Range *limit, int ref_N, double *ref_x, double *ref_f)
{


        int i = 0, j = 0, k = 0, rangeFound = 0;

	FILE *fw = fopen("interp.txt", "w");

        if (fw == NULL) { printf("Cannot open file to write.... exiting \n"); exit(1);}

        double e[N], h[N], d[N], r[N], z[N];
        double df, dx, xs, x1, y, sqErr = 0.0;
        char buf[100];

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

        for (j = 0; j < ref_N; j++)
        {
                x1 = ref_x[j];
		if (limit->min_x > x1 || limit->max_x < x1) continue;


		for (k = 0; k < N; k++)
		{
			if (x[k] < x1 && x[k+1] > x1)
			{
				i = k;
				rangeFound = 1; 
				break;
			}

			else if (x[k] == x1)
			{
				y = f[k]; break;
			}

			else if (x[k+1] == x1)
			{
				y = f[k+1]; break;
			}

		}
		if (rangeFound == 0) 
		{ 
			if (x[k] == x1 || x[k+1] == x1)
			{
				sqErr += pow(y-ref_f[j], 2);
				fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
			}
			continue; 
		}

		printf("j = %d, x1 = %10.6f , Range = %d, xleft = %10.6f, xright = %10.6f \n", j, x1, i+1, x[i], x[i+1]);

                y = -z[i]*pow( (x1-x[i+1]),3 )/(6.0*d[i]) +
                     z[i+1]*pow( (x1-x[i]),3 )/(6.0*d[i]) +
                    (f[i+1]/d[i] - z[i+1]*d[i]/6.0)*(x1 - x[i]) +
                    (-f[i]/d[i] + z[i]*d[i]/6.0)*(x1 - x[i+1]);

		printf("j = %d, x1 = %10.6f , Range = %d, xleft = %10.6f, xright = %10.6f, y = %10.6f, ref_f = %10.6f \n", j, x1, i+1, x[i], x[i+1], y, ref_f[j]);

		sqErr += pow(y-ref_f[j], 2);

                fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
        }
  
	fclose(fw);

	return sqrt(sqErr);
}
int main(int argc, char *argv[])
{
	int numBands = 9,i = 0;
	int regionsPerBand = 3, pointsPerRegion_vasp = 21;

	double x_vasp_data[numBands][1000], f_vasp_data[numBands][1000];
        double x_gulp_data[numBands][1000], f_gulp_data[numBands][1000];	

	double x[500], f[500];
	double ref_x[regionsPerBand*pointsPerRegion_vasp], ref_f[regionsPerBand*pointsPerRegion_vasp];

	double SquaredError = 0.0;

	double min_x = 10000000000.0, max_x = -10000000000.0, min_f = 10000000000.0, max_f = -10000000000.0;

	struct Range *lim = (struct Range*)malloc(sizeof(struct Range));

	ReadData(numBands, x_vasp_data, f_vasp_data, x_gulp_data, f_gulp_data);

	FILE *fw1 = fopen("gulp.txt","w"), *fw2 = fopen("vasp.txt", "w");

	int BandNum = 0;
	//sscanf(argv[1], "%s     \n", tempStr);
	sscanf(argv[1], "%d", &BandNum);

	printf("%d\n", BandNum);

        for(i = 0; i < 500; i++)
	{
                fprintf(fw1, "%10.6f %10.6f \n", x_gulp_data[BandNum][i], f_gulp_data[BandNum][i]);
		x[i] = x_gulp_data[BandNum][i];
		f[i] = f_gulp_data[BandNum][i];
		if (x[i] <= min_x) min_x = x[i];
		if (x[i] >= max_x) max_x = x[i];
		if (f[i] <= min_f) min_f = f[i];
		if (f[i] >= max_f) max_f = f[i];
	}

        for(i = 0; i < regionsPerBand*pointsPerRegion_vasp; i++)
	{
                fprintf(fw2, "%10.6f %10.6f \n", x_vasp_data[BandNum][i]*1000.0, f_vasp_data[BandNum][i]*33.35641);
		ref_x[i] = x_vasp_data[BandNum][i]*1000.0;
		ref_f[i] = f_vasp_data[BandNum][i]*33.35641;
	}


        fclose(fw1); fclose(fw2);

	lim->min_x = min_x; lim->min_f = min_f; lim->max_x = max_x; lim->max_f = max_f;

	printf("Min x = %10.6f, Min f = %10.6f, Max x = %10.6f, Max f = %10.6f \n", lim->min_x, lim->min_f, lim->max_x, lim->max_f);

	SquaredError = CubicSplineInterp(500, x, f, lim, regionsPerBand*pointsPerRegion_vasp, ref_x, ref_f);

	printf("%10.6f \n", SquaredError);

	//void ErrorEstimation(regionsPerBand*pointsPerRegion_vasp, lim, ref_x, ref_f, )

	return 0;
}







