#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"

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
        fprintf(out_stream,"dispersion 1 %d \n", 183); // what is 183 here? Better as variable?
        fprintf(out_stream,"0.0 0.0 0.0 to 0.5 0.0 0.0\n");
    } else if (bands == 2) {
        fprintf(out_stream,"dispersion 1 %d \n", 106); // what is 183 here? Better as variable?
        fprintf(out_stream,"0.5 0.0 0.0 to 0.33333 0.33333 0.0\n");
    } else if (bands == 3) {
        fprintf(out_stream,"dispersion 1 %d \n", 211); // what is 183 here? Better as variable?
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
    
    double theta0 = 81.9966; // geometry-dependent constant
    
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
    
    for (i=0;i<3;i++) { //modifies 2-body terms
        
        fseek(ff,11,SEEK_CUR); //goes to first modified value
        //sprintf(str,"%.3lf    %.3lf    %.4lf   0.00 %.6f",A[i],rho[i],B[i],rmax_2b[i]);
        sprintf(str,"%6.3lf   %6.3lf   %7.4lf",A[i],rho[i],B[i]); // assumes A, rho, B will not exceed 100

        fputs(str,ff); // updates first line w/ correct values, current cursor position at end of rmax
    
        target = line + 1;
        while (line != target) {
        if((c = fgetc(ff)) == '\n') line++; // moves to next line
        }
        
    }
    
    while (line !=11) {
        if((c = fgetc(ff)) == '\n') line++; // moves to 12th line for 3-body terms
    }
    
    for (i=0;i<2;i++) { //modifies 3-body terms
        
        fseek(ff,11,SEEK_CUR); //goes to first modified value
        //sprintf(str,"%.4lf  %.4lf  %.3lf  %.3lf   0.0 %.6lf   0.0 %.6lf   0.0 %.6lf",lambda[i],theta0,gamma_3b[i],gamma_3b[i],rmax_3b[0],rmax_3b[0],rmax_3b[1]);
        double gamma_dummy = gamma_3b[i];
        sprintf(str,"%7.4lf  %7.4lf %6.3lf %6.3lf",lambda[i],theta0,gamma_3b[i],gamma_dummy); // assumes lambda, theta0, gammas will not exceed 100
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
    double lat_const = 12.29, elast_real = 238; // elastic constant in GPa
    double elast_calc;

    char lines[200], search_string[] = "  Comparison of initial and final structures : ", ch;
    while (fgets(lines,200,af) != NULL) { // search output file for ^ test string
        int i;
        if (strstr(lines,search_string)) {
            break;
        }
    }
    
    int line_down;
    while (line_down !=5) {
        if((ch = fgetc(af)) == '\n') line_down++; // moves to constant "a"
    }
    fgets(lines,200,af); // reads line in
    int len = strlen(lines);
    double err_dum; // dummy error
    sscanf(&lines[len-8], "%lf",&err_dum); // assigns value of percent error in a
    *err_a = err_dum;
    
    line_down = 0;
    while (line_down !=1) { // moves down to elastic constant "c" (for elastic fit)
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    sscanf(lines,"%*s %*lf %lf", &c); // assigns value of elastic constant "c" (c_final)

    line_down = 0;
    while (line_down !=21) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    sscanf(lines, "%*d %lf",&c_11); // assigns value of elastic constant C11
    
    elast_calc = c_11/c*(2/lat_const); // elastic constant for comparison
    *elast = fabs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    

}

double phonon_disp() {
    
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

}

void initialize_population(int population_num) {
    
    FILE *population;
    
    population = fopen("ga.in", "w");
    
    double frac_perturb = 0.005;
    double rand_frac;
    
    fprintf(population, "%d\n", population_num);
    
    //initial guesses for the 3 As, 3 rhos, 3 Bs, 2 lambdas, and 2 independent gammas respectively modify manually as necessary
    double variables[13] = {5,5,5,0.5,0.5,0.5,20,20,20,15,15,1,1};
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

void processor_run(char folder[], char inputs[]) {
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[3], rho[3], B[3], lambda[2], gamma_3b[2]; //allocates variables, note/recall equilvalences in gammas/rmax_3bs
    double err_a, elast, chi_sq; // allocates objectives, namely error in a, elastic (multiply these to get single objective) and chi_sq
    char file[200], path[200], c,dc[20];
    int bands;

    // create folder for this instance of gulp run (now done in first copy cell step)

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

    modify_input(A,rho,B,lambda,gamma_3b); // modifies "forcefield" file

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
   
    gulp_input = fopen("in.gulp","a");
    // 183 105 211
    fclose(gulp_input);

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
    read_output(&err_a,&elast); // extracts objective 1 & 2: error in lattice constant a and error in elastic constant
    
    chdir(".."); // moves back to 'folder' directory to read all segments for phonon dispersion chi squared calculation
    chi_sq = phonon_disp();  // may need to modify paths (don't know where phonon_disp() is placed by gulp)

    FILE *output;

    output = fopen("ga_line","w"); // writes to a single line text document ga_line -> sloppy workaround to avoid race conditions. Find something more intellegent!

    fprintf(output,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A[0],A[1],A[2],rho[0],rho[1],rho[2],B[0],B[1],B[2],lambda[0],lambda[1],gamma_3b[0],gamma_3b[1],err_a,elast,chi_sq);
    
    chdir("..");
    fclose(output);
    
   
    
}

int main(int argc, char *argv[]) {
    
    
    int myid; /* My rank */
    int nprocs; /* Number of processors */
    int iteration_num = 1; // number of training iterations
    int population_num = 5; // population size for ga.in training
    int initialized;

    initialize_population(population_num); // initializes population

    
    int i,j,line; // iterators
    char c,variables[200],folder[200];
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
            fgets(variables,200,input);
            
            // end

            processor_run(folder,variables); // runs analysis on line input

        }
        fclose(input);

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
        }

        // ./ga.c
        // Work ga.c output
        //file_copy("ga.in","ga.out",)
    }
    MPI_Finalize();

    return 0;
}
