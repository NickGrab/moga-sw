#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

struct Range
{
        double min_x, min_f, max_x, max_f;

};

void file_copy(char file[], char file_copy[], int bands) { // creates a copy of a file
    FILE *in_stream, *out_stream;
    char c;
    
    in_stream = fopen(file,"r");
    if (in_stream == NULL)
    {
        printf("Cannot open file %s \n", file); // null pointer handler
        exit(0);
    }
    
    out_stream = fopen(file_copy, "a");
    if (out_stream == NULL)
    {
        printf("Cannot open file %s \n", file_copy); // null pointer handler
        exit(0);
    }
    
    c = fgetc(in_stream);
    while (c != EOF)
    {
        fputc(c, out_stream);
        c = fgetc(in_stream);
    }
  
    if (bands == 1) { // I know this method of setup is dumb, but it was a lot easier than modifying my previous method.
	fprintf(out_stream,"dispersion 1 %d \n", 183);
        fprintf(out_stream,"0.0 0.0 0.0 to 0.5 0.0 0.0\n");
    } else if (bands == 2) {
	fprintf(out_stream,"dispersion 1 %d \n", 106);
        fprintf(out_stream,"0.5 0.0 0.0 to 0.33333 0.33333 0.0\n");
    } else if (bands == 3) {
	fprintf(out_stream,"dispersion 1 %d \n", 211);
        fprintf(out_stream,"0.33333 0.33333 0.0 to 0.0 0.0 0.0\n");
    }

    if (bands > 0) {
        fprintf(out_stream,"output phon phonon\n");
        fprintf(out_stream,"shrink 5 5 5\n\n");
    }  
   
   
    fclose(in_stream);
    fclose(out_stream);
    
}

void modify_input(double A[3], double rho[3], double B[3], double lambda[2], double gamma_3b[2]) {

    FILE *ff;
    char filename[] = "forcefield";
    ff = fopen(filename, "r+");

    double theta0 = 81.9966;

    char c, str[200];
    int line=0;
    while (line !=4) {
        if((c = fgetc(ff)) == '\n') line++;
    }

    int i, target;

    for (i=0;i<3;i++) { 

        fseek(ff,11,SEEK_CUR);
        sprintf(str,"%6.3lf   %6.3lf   %7.4lf",A[i],rho[i],B[i]); 

        fputs(str,ff); 

        target = line + 1;
        while (line != target) {
        if((c = fgetc(ff)) == '\n') line++; 
        }

    }

    while (line !=11) {
        if((c = fgetc(ff)) == '\n') line++; 
    }

    for (i=0;i<2;i++) { //modifies 3-body terms

        fseek(ff,11,SEEK_CUR);
        double gamma_dummy = gamma_3b[i];
        sprintf(str,"%7.4lf  %7.4lf %6.3lf %6.3lf",lambda[i],theta0,gamma_3b[i],gamma_dummy); 
        fputs(str,ff); 

        target = line + 1;
        while (line != target) {
            if((c = fgetc(ff)) == '\n') line++; // moves to next line
        }
    }

    fclose(ff); //close input file

}

void read_output(double* err_a,double* elast) {

    FILE *af;
   
    char output[] = "afterfit.out";
    af = fopen(output,"r"); // read output file
    
    if (af == NULL)
    {
        printf("Cannot open file %s \n", output); // null pointer handler
        exit(0);
    }
    
    
    double c_11, c, el_c;
    double lat_const = 12.29, elast_real = 238; // elastic constant in GPa
    double elast_calc;

    char tempStr1[100], tempStr2[100];
    double temp1, temp2, temp3;

    char lines[200], search_string[] = "  Comparison of initial and final structures : ", ch;
    while (fgets(lines,200,af) != NULL) { // search output file for ^ test string
        int i;
        if (strstr(lines,search_string)) {
	    printf(lines);
            break;
        }
    }
    
    int line_down=0;
    while (line_down !=5) {
        if((ch = fgetc(af)) == '\n') line_down++; // moves to constant "a"
    }
    fgets(lines,200,af); // reads line in
    int len = strlen(lines);
    double err_dum; // dummy error
    printf("Extracting err_a\n");
    sscanf(lines, "%s %lf %lf %lf %s %lf", tempStr1, &temp1, &temp2, &temp3, tempStr2, &err_dum); // assigns value of percent error in a
    *err_a = err_dum;
    printf("Extracting err_a done\n");    

    line_down = 0;
    while (line_down !=1) { // moves down to elastic constant "c" (for elastic fit)
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);

    printf("Extracting C\n");
    sscanf(lines,"%s %lf %lf %lf %s %lf", tempStr1, &temp1, &c, &temp2, tempStr2, &temp3); // assigns value of elastic constant "c" (c_final)
    printf("Extracting C done\n");

    line_down = 0;
    while (line_down !=21) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);

    printf("Extracting C11\n");
    sscanf(lines, "%*d %lf",&c_11); // assigns value of elastic constant C11
    printf("Extracting C11 done\n");    

    elast_calc = c_11*c*(2/lat_const); // elastic constant for comparison
    *elast = fabs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    
    return;
}

/*double phonon_disp() {
    
    FILE *phonon, *band;
    const char *phonon_file[3] = {"G_M/phonon.disp","M_K/phonon.disp","K_G/phonon.disp"}; // phonon dispersion segments created by gulp in folders created by processor_run()
    const char *band_file = "../band.dat"; // band.dat file location

    char lines_phonon[200], lines_band[200], lines_band_ph[200]; // used for line by line storage
    char c_phonon,c_band; // used for return counts for document navigation
    int j, j_comp, k, l_phonon = 0, l_band = 0; // loop counters
    
    int weights[9] = {10,10,10,3,3,3,1,1,1}; // weights for the chi^2 fit
    double conv_x = 1000; // reciprical position conversion rate
    double conv_y = 33.35641; // frequencies conversion rate (THz) -> 1/cm
    double phonon_shift[3] = {0,183,289};
    
    double x_band, y_band; // current position and frequency of empirical band data
    double chi_sq; // chi-squared accumulator for fitting
    double new_x_phonon, old_x_phonon, new_y_phonon, old_y_phonon, y_interp_phonon; // bracketing position and frequency values for eventually interpolation
    int target; // target line for current segment reading in band.dat
    
    int band_offset = 26; // number of lines in band.dat before first segment data
    int band_segm_lngth = 21; // number of data points in a given segment
    int phonon_offset = 3; // offset before first segment data in phonon.disp file
    
    band = fopen(band_file,"r"); // opens band.dat file for reading
    
    if (band == NULL)
    {
        printf("Cannot open file %s \n", band_file); // null pointer handler
        exit(0);
    }
    
    for (k=0;k<3;k++) { // loops over segments
        
        //printf("The value of k is : %d\n",k); TEST PRINT
        
        phonon = fopen(phonon_file[k],"r"); // reads from 'current segment' (G-M, M-K, K-G)

        if (phonon == NULL)
        {
            printf("Cannot open file %s \n", phonon_file[k]); // null pointer handler
            exit(0);
        }

        rewind(band); // resets band.dat data stream to beginning of file for loop reading consistency
        
        target = band_offset + k*(band_segm_lngth+1); // sets target for first line of current segment in band.dat file
        while (l_band != target) {
            if((c_band = fgetc(band)) == '\n') l_band++; // moves down to first line of 'current segment' in band file
        }
        
        for (j=0;j<9;j++) { // loops over each band
            
            fgets(lines_band,200,band); // gets first value for band.dat file in current segement
            sscanf(lines_band, "%lf %lf", &x_band, &y_band); // reads in first data point of band 1 for current segment
            
            rewind(phonon); // resets to top of phonon file for subsequent loops
            
            l_phonon = 0; // resets line counter
            
            while (l_phonon != phonon_offset + j) {
                if((c_phonon = fgetc(phonon)) == '\n') l_phonon++; // moves down to first line of 'current segment' in phonon file
            }
            
            fgets(lines_phonon,200,phonon); // gets first line in current segment of phonon.disp
            sscanf(lines_phonon, "%lf %lf", &old_x_phonon, &old_y_phonon); // reads in initial values
            
            old_x_phonon += phonon_shift[k];
            
            l_phonon = 0;
            while (l_phonon != 8) {
                if((c_phonon = fgetc(phonon)) == '\n') l_phonon++; // moves down nine lines (to the next line of current segment)
            }
            
            while (fgets(lines_phonon,200,phonon) != NULL) { // runs until phonon gives a null (end of document)
                
                sscanf(lines_phonon, "%lf %lf", &new_x_phonon, &new_y_phonon); // scans in new phonon values (from fgets above)
                new_x_phonon += phonon_shift[k];
                
                if (x_band*conv_x <= new_x_phonon) { // checks if current value of band (to be interpolated to) is beneath upper bound
                    y_interp_phonon = old_y_phonon + (x_band*conv_x-old_x_phonon)*(new_y_phonon-old_y_phonon)/(new_x_phonon-old_x_phonon+0.0000000001); // uses linear interpolation to find phonon.disp value for comparison, perturbation ensures no division by zero
                    
                    chi_sq += weights[j]*pow((y_interp_phonon-y_band*conv_y),2); // adds chi squared value
                    //printf("%f\n", chi_sq);
                    fgets(lines_band,200,band); // scans in new band.dat data for comparison
                    
                }
                
                
//                if (k>0) {
//                    printf("\n%f %f %f \n", old_x_phonon, x_band*conv_x, new_x_phonon);
//                    printf("%f %f %f \n\n", old_y_phonon, y_band*conv_y, new_y_phonon);
////                    printf("%f\n", y_interp_phonon); // TEST PRINTS (DELETE US)
//                }

                
                
                if (strncmp(lines_band,"\n",200) == 0) { // will end loop if band has reached last line of segment (empty line)
                    break;
                }
                
                sscanf(lines_band, "%lf %lf", &x_band, &y_band); // scans in values from band.dat
                
                old_x_phonon = new_x_phonon; // set lower bound of checked interval for interpolation
                old_y_phonon = new_y_phonon;
                
                l_phonon = 0;
                while (l_phonon != 8) {
                    if((c_phonon = fgetc(phonon)) == '\n') l_phonon++; // moves down nine lines (to the next line of current segment)
                }
                
            }
            
            rewind(band); // rewinds band.dat stream to beginning
            
            l_band = 0;
            while (l_band != band_offset) {
                if((c_band = fgetc(band)) == '\n') l_band++; // skips intro block (necessary to avoid weird catches/ensure consistency)
            }
            
            j_comp = 0;

            do {
                strcpy(lines_band_ph,lines_band);
                fgets(lines_band,200,band);
                if (strncmp(lines_band,lines_band_ph,200)==0) {
                    j_comp++;
                }
            } while (j_comp != j); // moves to beginning of next band (marked by double space)
            
            
            l_band = 0;
            while (l_band != k*(band_segm_lngth+1)) {
                if((c_band = fgetc(band)) == '\n') l_band++; // moves down to first line of current segment in current band of band file
            }
        }
        
        fclose(phonon); // closes current phonon file (in preparation of opening of next phonon file

    }
    //chi_sq = chi_sq/567;
    return chi_sq;

}*/

void ReadData(int numBands, double x_vasp_data[][10000], double f_vasp_data[][10000], double x_gulp_data[][10000], double f_gulp_data[][10000])
{
        int regionsPerBand = 3, pointsPerRegion_vasp = 21;
        int points_GM_gulp = 183, points_MK_gulp = 106, points_KG_gulp = 211;

	double temp_x, temp_y;

        double phonon_shift[3] = {0, 183, 289};
        double conv_x = 1000, conv_y = 33.35641;
        int band_offset = 26, band_segm_length = 21, gulp_offset = 3;

        int i = 0, j = 0, k = 0, lineNum_vasp = 0, lineNum_gulp = 0, bands_vasp = 0, bands_gulp = 0;

        const char *phonon_file[3] = {"G_M/phonon.disp", "M_K/phonon.disp", "K_G/phonon.disp"};
        const char *band_file = "../UTIL/band.dat";
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
				//      x_vasp_data[i][k + i*regionsPerBand*band_segm_length], f_vasp_data[i][k + i*band_segm_length]);
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
			}
		}

		prevPoint += curPoints;
        }

        fclose(vasp); fclose(gulp[0]); fclose(gulp[1]); fclose(gulp[2]);

	return;
}

double CubicSplineInterp(int N, double *x, double *f, struct Range *limit, int ref_N, double *ref_x, double *ref_f, int id)
{
        int i = 0, j = 0, k = 0, rangeFound = 0;
        char fileName[100];
        sprintf(fileName, "interp-%d.txt", id);

        FILE *fw = fopen(fileName, "w");

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
                y = -z[i]*pow( (x1-x[i+1]),3 )/(6.0*d[i]) +
                     z[i+1]*pow( (x1-x[i]),3 )/(6.0*d[i]) +
                    (f[i+1]/d[i] - z[i+1]*d[i]/6.0)*(x1 - x[i]) +
                    (-f[i]/d[i] + z[i]*d[i]/6.0)*(x1 - x[i+1]);

		sqErr += pow(y-ref_f[j], 2);

                fprintf(fw, "%10.4f \t %10.4f \n", x1, y);
	}

	fclose(fw);
        return sqrt(sqErr);

}

double ErrorPhononDispersion()
{
        int BandNum = 0, numBands = 9, i = 0;

        int regionsPerBand = 3, pointsPerRegion_vasp = 101;

        double x_vasp_data[numBands][10000], f_vasp_data[numBands][10000];
        double x_gulp_data[numBands][10000], f_gulp_data[numBands][10000];

        double x[5000], f[5000];
        double ref_x[regionsPerBand*pointsPerRegion_vasp], ref_f[regionsPerBand*pointsPerRegion_vasp];
        double Error[numBands], w[numBands];

	w[0] = 1.0; w[1] = 1.0; w[2] = 1.0; w[3] = 1.0; w[4] = 1.0; w[5] = 1.0; w[6] = 0.03; w[7] = 0.03; w[8] = 0.03;
        double SquaredError = 0.0;

        double min_x = 10000000000.0, max_x = -10000000000.0, min_f = 10000000000.0, max_f = -10000000000.0;

        struct Range *lim = (struct Range*)malloc(sizeof(struct Range));

        ReadData(numBands, x_vasp_data, f_vasp_data, x_gulp_data, f_gulp_data);
  
        for (BandNum=0; BandNum < numBands; BandNum++)
        {
        	for(i = 0; i < 500; i++)
        	{
        		x[i] = x_gulp_data[BandNum][i];
        		f[i] = f_gulp_data[BandNum][i];
        		if (x[i] <= min_x) min_x = x[i];
                        if (x[i] >= max_x) max_x = x[i];
			if (f[i] <= min_f) min_f = f[i];
                        if (f[i] >= max_f) max_f = f[i];
                }

                for(i = 0; i < regionsPerBand*pointsPerRegion_vasp; i++)
                {
                        ref_x[i] = x_vasp_data[BandNum][i]*1000.0;
                        ref_f[i] = f_vasp_data[BandNum][i]*33.35641;
                }

                lim->min_x = min_x; lim->min_f = min_f; lim->max_x = max_x; lim->max_f = max_f;

                Error[BandNum] = CubicSplineInterp(500, x, f, lim, regionsPerBand*pointsPerRegion_vasp, ref_x, ref_f, BandNum);

                SquaredError += w[BandNum]*Error[BandNum];
	}
        //for (BandNum=0; BandNum < numBands; BandNum++) printf("%10.6f \t ", Error[BandNum];
	return (SquaredError/(numBands*regionsPerBand*pointsPerRegion_vasp));
}

void initialize_population(int population_num) {
    
    FILE *population;
    
    population = fopen("ga.in", "w");
    
    double frac_perturb = 0.1;
    double rand_frac;
    
    fprintf(population, "%d\n", population_num);
    
  
    double variables[13] = {1,7,5,0.5,0.5,0.5,20,20,20,15,15,1,1};
    double rand_var[13]; // array for holding random perturbations of variables
 
    int i,j,k;
    
    srand(time(NULL));
    
    for (i=0;i<population_num;i++) {
        for (j=0;j<13;j++) {
            rand_frac = rand() / (double)(RAND_MAX)*(2*frac_perturb) - frac_perturb;
            rand_var[j] = (1-rand_frac)*variables[j];
            fprintf(population,"%lf ",rand_var[j]);
        }
        fprintf(population,"\n");
    }
    
    fclose(population);
}

void processor_run(int iter_num, char folder[], char inputs[]) {
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[3], rho[3], B[3], lambda[2], gamma_3b[2]; //allocates variables, note/recall equilvalences in gammas/rmax_3bs
    double err_a, elast, chi_sq; 
    char file[200], path[200], c,dc[20];
    int bands;

    //double weight_a = 5.0, weight_c11 = 1.0;

    mkdir(folder,S_IRWXU);
    // copy cell, forcefield, etc. into folder
    
    const char *gulp_in[2] = {"cell","forcefield"};
    int i;    
   
    for(i=0;i<2;i++) {      // copies "cell" and "forcefield" into new folder
        strcpy(file,gulp_in[i]);
        strcpy(path,folder);
        strcat(path,"/");
        strcat(path,file);
        
        file_copy(file,path,bands);
    }
 
    chdir(folder);   //enter folder


    // scans line for variables

    sscanf(inputs,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&A[0],&A[1],&A[2],&rho[0],&rho[1],&rho[2],&B[0],&B[1],&B[2],&lambda[0],&lambda[1],&gamma_3b[0],&gamma_3b[1]);
    modify_input(A,rho,B,lambda,gamma_3b); 
    
    FILE *gulp_input;

    gulp_input = fopen("in.gulp","w"); // create gulp input file
    
    fprintf(gulp_input, "optim relax conp comp phon nofreq\n");
    fclose(gulp_input);
    

    strcpy(file,"cell");
    strcpy(path,"in.gulp");
    bands = 0;
    file_copy(file,path,bands); // appends "cell" to in.gulp
    
    strcpy(file,"forcefield");
    file_copy(file,path,bands); // appends "forcefield" to in.gulp
   

    const char *segment[3] = {"G_M","M_K","K_G"};
    for(i=0;i<3;i++) {
        
        strcpy(dc,segment[i]);
        mkdir(dc,S_IRWXU); // creates directory for segment
        
        strcpy(file,"in.gulp");
        strcpy(path,dc);
        strcat(path,"/");
        strcat(path,file);  // copies in.gulp into new directory
        bands = i+1;
        
        file_copy(file,path,bands);
        
        chdir(dc); // enter segment directory
        
        system("/home/pv-02/hpc-23/kris658/SOFTWARE/GULP/Src/Linux/gulp < in.gulp > afterfit.out"); // runs gulp
        
        chdir(".."); // exits new directory
    }
    chdir(segment[0]); // arbitrary choice of segment for objectives 1 & 2 (same for all)
    printf("Extracting elastic constant \n");
    read_output(&err_a,&elast); // extracts objective 1 & 2: error in lattice constant a and error in elastic constant

    //err_a *= weight_a; elast *= weight_c11;

    printf("Extracting elastic constant --done\n");    

    chdir(".."); // moves back to 'folder' directory to read all segments for phonon dispersion chi squared calculation
    //chi_sq = phonon_disp(); 
    printf("Computing error in dispersion \n");
    chi_sq = ErrorPhononDispersion();

    printf("Computing error in dispersion --done\n");

    FILE *output;

    output = fopen("ga_line","w"); 
    fprintf(output,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
	            A[0],A[1],A[2],rho[0],rho[1],rho[2],B[0],B[1],B[2],lambda[0],lambda[1],gamma_3b[0],gamma_3b[1],err_a,elast,chi_sq);
    
    chdir("..");
    fclose(output);
}

int main(int argc, char *argv[]) {
    
    
    int myid; /* My rank */
    int nprocs; /* Number of processors */
    int iteration_num = 10; // number of training iterations
    int population_num = 200; // population size for ga.in training
    int initialized;

    initialize_population(population_num); // initializes population

    
    int i,j,line; // iterators
    char c,variables[500],folder[200];
    char test[200];
    
    for(j=0;j<iteration_num;j++) { // recursively optimizes for number of iterations specified


        MPI_Initialized(&initialized);
        if (!initialized) {
            MPI_Init(&argc, &argv);
        }
        MPI_Comm_rank(MPI_COMM_WORLD, &myid); // starts MPI
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      
        FILE *input;
        input = fopen("ga.in","r"); // opens ga.in file for reading
        
        if (input == NULL)
        {
            printf("Cannot open file \n"); // null handler
            exit(0);
        }


        for (i=myid;i<population_num;i+=nprocs) { // loops through population list, assigning members to processors

            sprintf(folder,"%d",i); // creates a folder to hold all analysis on current line

            // gets line specified by loop and reads it into memory
            
            line = 0;
            rewind(input);

            while (line!=(i+1)) {
                if((c = fgetc(input)) == '\n') line++;
            }
            fgets(variables,500,input);
            
            // end

            processor_run(j, folder,variables); // runs analysis on line input

        }
        fclose(input);

	MPI_Barrier(MPI_COMM_WORLD);

        if (myid == 0) {
            FILE *ga_input,*ga_line;
            //sprintf(test,"ga_%d.in", j);
            ga_input = fopen("ga.in","w");

            if (ga_input == NULL)
            {
                printf("Cannot open file \n");
                exit(0);
            }

            fprintf(ga_input,"%d\n",population_num);

            for (i=0;i<population_num;i++) {
                sprintf(folder,"%d",i);
                chdir(folder);

                ga_line = fopen("ga_line","r");
                fgets(variables,200,ga_line);
                fputs(variables,ga_input);

                fclose(ga_line);

                chdir("..");
                // code that removes directories files generated above
                char rm_dir[200];
                sprintf(rm_dir,"rm -rf %s", folder);
                system(rm_dir);
            }
            fclose(ga_input);

            char cmdexec1[200];
	    sprintf(cmdexec1,"./ga ga.in %d", j+1);
	    system(cmdexec1);
            char cmdexec[200];
	    sprintf(cmdexec, "cp value.d ga.in");
	    system(cmdexec);
        }
	printf("Iteration = %d done \n", j);
    MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();

    return 0;
}
