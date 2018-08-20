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
    
    if (bands == 1) { // X-G region
        fprintf(out_stream,"dispersion 1 %d \n", 144); 
        fprintf(out_stream,"0.5 0.0 0.0 to 0.0 0.0 0.0\n");
    } else if (bands == 2) { // G-Y region
        fprintf(out_stream,"dispersion 1 %d \n", 79);
        fprintf(out_stream,"0.0 0.0 0.0 to 0.0 0.5 0.0\n");
    } else if (bands == 3) { // Y-S region
        fprintf(out_stream,"dispersion 1 %d \n", 143); 
        fprintf(out_stream,"0.0 0.5 0.0 to 0.5 0.5 0.0\n");
    }
     else if (bands == 4) { // S-G region
        fprintf(out_stream,"dispersion 1 %d \n", 164); 
        fprintf(out_stream,"0.5 0.5 0.0 to 0.0 0.0 0.0\n");
    }
    
    if (bands > 0) {
        fprintf(out_stream,"output phon phonon\n");
        fprintf(out_stream,"shrink 5 5 5\n\n");
    }
    
    fclose(in_stream);
    fclose(out_stream);
    
}

void modify_input(double A[6], double rho[6], double B[6], double lambda[5], double gamma_MoTe1Te2[2], double gamma_rest[4]) {
    
    double theta0[5] = {82.50, 76.75, 115.75, 76.04, 64.22}; // geometry-dependent constant
    
    FILE *ff;
    char filename[] = "forcefield"; // Input file name
    ff = fopen(filename, "r+");
    
    char c, str[200]; //current character, string to hold number
    int line=0; // current line, current position in file
    
    // moves to 5th line for first modification

    while (line !=4) {
        if((c = fgetc(ff)) == '\n') line++;
    }
    int i, target;
    
    for (i=0;i<6;i++) { //modifies 2-body terms 6-pairs in total Te1-Te2, Te1-Te1, Te2-Te2, Mo-Te1, Mo-Te2, Mo-Mo
        
        fseek(ff,11,SEEK_CUR); //goes to first modified value
        sprintf(str,"%6.3lf   %6.3lf   %7.4lf",A[i],rho[i],B[i]); // Modifies A, rho and B two every pair corresponding to the tuple in ga.in file

        fputs(str,ff); // updates first line w/ correct values, current cursor position at end of B
    
        target = line+1;
        while (line != target) {
        if((c = fgetc(ff)) == '\n') line++; // Loops moves the counter to next line
        } 
    }
    
    while (line !=14) {
        if((c = fgetc(ff)) == '\n') line++; // moves to 15th line for 3-body terms
    }
    
    for (i=0;i<5;i++) { //modifies 3-body terms 5 in total Mo-Te1-Te2, Mo-Te1-Te1, Mo-Te2-Te2, Te1-Mo-Mo, Te2-Mo-Mo
        
        fseek(ff,11,SEEK_CUR); //goes to first modified value
        if (i==0) {
		// For Mo-Te1-Te2 gamma0 and gamma1 are different
        	sprintf(str,"%7.4lf  %7.4lf %6.3lf %6.3lf", lambda[i], theta0[i], gamma_MoTe1Te2[0], gamma_MoTe1Te2[1]);
	}
	else if (i>=1) {
		sprintf(str,"%7.4lf  %7.4lf %6.3lf %6.3lf", lambda[i], theta0[i], gamma_rest[i-1], gamma_rest[i-1]);
	}

        fputs(str,ff); // updates first line w/ correct values, current cursor position at end of rmax
        
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
    double lat_const = 13.97, elast_real = 128.7; // elastic constant in GPa
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
    while (line_down !=30) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);

    printf("Extracting C11\n");
    int temp;
    printf("%s\n", lines);
    sscanf(lines, "%*d %lf", &c_11); // assigns value of elastic constant C11
    printf("Extracting C11 done\n");    

    printf("c = %f, c_11 = %f\n", c, c_11);

    elast_calc = c_11*c*(2/lat_const); // elastic constant for comparison
    *elast = fabs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    
    return;
}

void ReadData(int numBands, double x_vasp_data[][10000], double f_vasp_data[][10000], double x_gulp_data[][10000], double f_gulp_data[][10000])
{
        int regionsPerBand = 4, pointsPerRegion_vasp = 101;
        int points_XG_gulp = 144, points_GY_gulp = 79, points_YS_gulp = 143, points_SG_gulp = 164;

	double temp_x, temp_y;

        double phonon_shift[4] = {0, 144, 223, 366};
        double conv_x = 1000, conv_y = 33.35641;
        int band_offset = 26, band_segm_length = 21, gulp_offset = 3;

        int i = 0, j = 0, k = 0, lineNum_vasp = 0, lineNum_gulp = 0, bands_vasp = 0, bands_gulp = 0;

        const char *phonon_file[4] = {"X_G/phonon.disp", "G_Y/phonon.disp", "Y_S/phonon.disp", "S_G/phonon.disp"};
        const char *band_file = "../UTIL/band.dat";
        char lines_gulp[200], lines_band[200], lines_band_ph[200];

        FILE *vasp, *gulp[4];

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

                if (i == 0) curPoints = points_XG_gulp;
                else if (i == 1) curPoints = points_GY_gulp;
                else if (i == 2) curPoints = points_YS_gulp;
		else if (i == 3) curPoints = points_SG_gulp;

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
        int BandNum = 0, numBands = 18, i = 0;

        int regionsPerBand = 4, pointsPerRegion_vasp = 101;

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

                Error[BandNum] = CubicSplineInterp(530, x, f, lim, regionsPerBand*pointsPerRegion_vasp, ref_x, ref_f, BandNum);

                SquaredError += w[BandNum]*Error[BandNum];
	}
        //for (BandNum=0; BandNum < numBands; BandNum++) printf("%10.6f \t ", Error[BandNum];
	return (SquaredError/(numBands*regionsPerBand*pointsPerRegion_vasp));
}

void initialize_population(int population_num) {
    
    FILE *population;
    
    population = fopen("ga.in", "w");
    
    double frac_perturb = 0.05;
    double rand_frac;
    
    fprintf(population, "%d\n", population_num);
    
    //initial guesses for the 6 As, 6rhos, 6 Bs, 5 lambdas, 2 gamma_MoTe1Te2 and 4 gamma values for other 3 body pairs
    double variables[29] = {3, 3, 3, 6.6, 6.6, 2.7, 2.0, 2.0, 2.0, 0.31, 0.31, 0.28, 50.0, 50.0, 50.0, 13.6, 13.6, 38.5, 32.9, 32.9, 32.9, 21.2, 21.2, 1.9, 1.9, 1.9, 1.9, 5.7, 5.7};
    double rand_var[29]; // array for holding random perturbations of variables
    
    int i,j,k;
    
    srand(time(NULL));
    
    for (i=0;i<population_num;i++) {
        for (j=0;j<29;j++) {
            rand_frac = rand() / (double)(RAND_MAX)*(2*frac_perturb) - frac_perturb;
            rand_var[j] = (1-rand_frac)*variables[j];
            fprintf(population,"%lf ",rand_var[j]);
        }
        fprintf(population,"\n");
    }
    
    fclose(population);
}

void processor_run(int id, int iter_num, char folder[], char inputs[]) {
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[6], rho[6], B[6], lambda[5], gamma_MoTe1Te2[2], gamma_rest[4]; //variables modified in the force field file
    double err_a, elast, chi_sq; 
    char file[200], path[200], c,dc[20];
    int bands;

    double weight_a = 1.0, weight_c11 = 1.0;

    mkdir(folder,S_IRWXU);
    // copy cell, forcefield, etc. into folder
    
    const char *gulp_in[2] = {"cell","forcefield"};
    int i;    
   
    for(i=0;i<2;i++) {      // copies "cell" and "forcefield" into new folder
        strcpy(file,gulp_in[i]);
        strcpy(path,folder);
        strcat(path,"/");
        strcat(path,file);
        
        file_copy(file,path,0);
    }
 
    chdir(folder);   //enter folder

    char cwd[100];
    getcwd(cwd, sizeof(cwd));
    printf("Current directory is %s\n", cwd);

    // scans line for variables
    
    sscanf(inputs,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &A[0],&A[1],&A[2],&A[3],&A[4],&A[5],&rho[0],&rho[1],&rho[2],&rho[3],&rho[4],&rho[5],&B[0],&B[1],&B[2],&B[3],&B[4],&B[5], &lambda[0],&lambda[1],&lambda[2],&lambda[3],&lambda[4],&gamma_MoTe1Te2[0],&gamma_MoTe1Te2[1], &gamma_rest[0], &gamma_rest[1], &gamma_rest[2], &gamma_rest[3]);

    modify_input(A,rho,B,lambda,gamma_MoTe1Te2, gamma_rest); // modifies "forcefield" file

    FILE *gulp_input;

    gulp_input = fopen("in.gulp","w"); // create gulp input file
    
    fprintf(gulp_input, "optim relax conp comp phon nofreq\n");
    fclose(gulp_input);
    
    char file1[100], file2[100], path1[200], path2[200];
   

    strcpy(file1,"cell");
    strcpy(path1,"in.gulp");
    bands = 0;
    file_copy(file1,path1,bands); // appends "cell" to in.gulp
    
    strcpy(file2,"forcefield");
    file_copy(file2,path1,bands); // appends "forcefield" to in.gulp
   

    const char *segment[4] = {"X_G","G_Y","Y_S","S_G"};
   
    char file0[4][100], path0[4][100];
 
    for(i=0;i<4;i++) {
        
        strcpy(dc,segment[i]);
        mkdir(dc,S_IRWXU); // creates directory for segment
        
        strcpy(file0[i],"in.gulp");
        strcpy(path0[i],dc);
        strcat(path0[i],"/");
        strcat(path0[i],file0[i]);  // copies in.gulp into new directory
        bands = i+1;
        
        file_copy(file0[i],path0[i],bands);
        
        chdir(dc); // enter segment directory
        
        system("/home/pv-02/hpc-23/kris658/SOFTWARE/GULP/Src/Linux/gulp < in.gulp > afterfit.out"); // runs gulp
        
        chdir(".."); // exits new directory
    }
    chdir(segment[0]); // arbitrary choice of segment for objectives 1 & 2 (same for all)
    printf("%d Extracting elastic constant \n", id);
    read_output(&err_a,&elast); // extracts objective 1 & 2: error in lattice constant a and error in elastic constant

    err_a *= weight_a; elast *= weight_c11;

    printf("Extracting elastic constant --done\n");    

    chdir(".."); // moves back to 'folder' directory to read all segments for phonon dispersion chi squared calculation
    //chi_sq = phonon_disp(); 
    printf("Computing error in dispersion \n");
    chi_sq = ErrorPhononDispersion();


    printf("Computing error in dispersion --done\n");

    FILE *output;

    output = fopen("ga_line","w"); 
    fprintf(output,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A[0],A[1],A[2],A[3],A[4],A[5],rho[0],rho[1],rho[2],rho[3],rho[4],rho[5],B[0],B[1],B[2],B[3],B[4],B[5],lambda[0], lambda[1], lambda[2],lambda[3],lambda[4],gamma_MoTe1Te2[0],gamma_MoTe1Te2[1], gamma_rest[0], gamma_rest[1], gamma_rest[2], gamma_rest[3],err_a,elast,chi_sq);
    
    chdir("..");
    fclose(output);
    
}

int main(int argc, char *argv[]) {
    
    
    int myid; /* My rank */
    int nprocs; /* Number of processors */
    int iteration_num = 100; // number of training iterations
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
            printf("%s\n", variables);
            processor_run(myid,j, folder,variables); // runs analysis on line input

        }
        fclose(input);

	MPI_Barrier(MPI_COMM_WORLD);

      /*  if (myid == 0) {
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
    MPI_Barrier(MPI_COMM_WORLD);*/
    }
    MPI_Finalize();

    return 0;
}
