#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "mpi.h"



int population_num = 5;
double theta0 = 81.9966; // necessary dummy variable for unchanging variable. Worth knowing globally
int iteration_num = 1;



void file_copy(char file[], char file_copy[], int bands) {
    FILE *in_stream, *out_stream;
    char c;
    
    in_stream = fopen(file,"r");
    if (in_stream == NULL)
    {
        printf("Cannot open file %s \n", file);
        exit(0);
    }
    
    out_stream = fopen(file_copy, "a");
    if (out_stream == NULL)
    {
        printf("Cannot open file %s \n", file_copy);
        exit(0);
    }
    
    c = fgetc(in_stream);
    while (c != EOF)
    {
        fputc(c, out_stream);
        c = fgetc(in_stream);
    }
    
    if (bands == 1) { // I know this method of setup is dumb, but it was a lot easier than modifying my previous method.
        fprintf(out_stream,"0.0 0.0 0.0 to 0.5 0.0 0.0\n");
    } else if (bands == 2) {
        fprintf(out_stream,"0.5 0.0 0.0 to 0.33333 0.33333 0.0\n");
    } else if (bands == 3) {
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
    
    // NOTE: STANDARD SPACING IS CRUTIAL FOR PROPER CODE FUNCTION
    // we should come up with a standard, but for now the only thing that needs to be changed
        // is adding a space in the beginning of the three body term so all terms are in-line and other values at top

    
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

void read_output(double err_a,double elast) {

    FILE *af;
    char output[] ="afterfit.out";
    af = fopen(output,"r"); // read output file
    
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
    sscanf(&lines[len-8], "%lf",&err_a); // assigns value of percent error in a
    
    line_down = 0;
    while (line_down !=1) { // moves down to elastic constant "c" (for elastic fit)
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    char cw; // dummy throwaway
    double c_init; //dummy throwaway
    sscanf(lines,"%s %lf %lf", &cw, &c_init, &c); // assigns value of elastic constant "c" (c_final)

    line_down = 0;
    while (line_down !=21) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    int throw; //throwaway int
    sscanf(lines, "%d %lf",&throw,&c_11); // assigns value of elastic constant C11
    
    elast_calc = c_11/c*(2/lat_const); // elastic constant for comparison
    elast = fabs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    

}

void phonon_disp(double chi_sq) {
    
    // put each phonon.disp file in a different directory?
    FILE *pd, *dft;
    const char *phonon[3] = {"G_M/phonon.disp","M_K/phonon.disp","K_G/phonon.disp"};

    char lines_pd[200], lines_dft[200], lines_dft_ph[200], c_pd,c_dft;
    int j, j_comp, k, l_pd,l_dft;
    
    int weights[9] = {10,10,10,3,3,3,1,1,1}; // weights for the chi^2 fit
    double conv_x = 1000;
    double conv_y = 33.35641; // conversions b/t files (if needed)
    
    double x_dft, y_dft;
    double new_x_pd, old_x_pd, new_y_pd, old_y_pd, y_interp_pd;
    int breaking;
    
    int dft_offset = 26, dft_segm_lngth = 21;
    int pd_offset = 3;
    
    //double posit_dft[dft_segm_lngth], freq_dft[dft_segm_lngth];
    //double posit_pd[pd_segm_lngth], freq_pd[pd_segm_lngth];
    
    
    const char *band = "band.dat";
    dft = fopen(band,"r");
    
    for (k=0;k<3;k++) {

        pd = fopen(phonon[k],"r"); // reads from 'current segment'
        rewind(dft);

        while (l_dft != dft_offset + k*(dft_segm_lngth+1)) {
            if((c_dft = fgetc(dft)) == '\n') l_dft++; // moves down to first line of 'current segment' in dft file
        }

        fgets(lines_dft,200,dft);
        sscanf(lines_dft, "%lf %lf", &x_dft, &y_dft); // reads in initial values
        
        for (j=0;j<9;j++) { // loops over each band (parallelizable - worth it?)

            rewind(pd); // resets to top of pd file for subsequent loops
            
            l_pd = 0;

            while (l_pd != pd_offset + j) {
                if((c_pd = fgetc(pd)) == '\n') l_pd++; // moves down to first line of 'current segment' in pd file
            }
            fgets(lines_pd,200,pd);
            sscanf(lines_pd, "%lf %lf", &old_x_pd, &old_y_pd); // reads in initial values
            
            l_pd = 0;
            while (l_pd != 9) {
                if((c_pd = fgetc(pd)) == '\n') l_pd++; // moves down nine lines (to the next line of current segment)
            }
            
            while (fgets(lines_pd,200,pd) != NULL) { // runs until pd gives a null (end of document)
                sscanf(lines_pd, "%lf %lf", &new_x_pd, &new_y_pd); // scans in new pd values (from fgets above)
                //printf("%f\n", new_x_pd);
                do {
                    if (x_dft <= new_x_pd) { // checks if current value of dft (to be interpolated to) is below upper bound of pd interval
                        y_interp_pd = old_y_pd + (x_dft*conv_x-old_x_pd)*(new_y_pd-old_y_pd)/(new_x_pd-old_x_pd); // uses linear interpolation to find value to compare
                        //printf("%f\n", y_interp_pd);
                        chi_sq += weights[j]*pow((y_interp_pd-y_dft*conv_y),2); // adds chi squared value
                        //printf("%f\n", chi_sq);

                         //scan in new dft values, check if b/t again before getting new line
                        
                        fgets(lines_dft,200,dft);
                        if (strncmp(lines_dft,"\n",200) == 0) { // will end loop if dft has reached last line (prevents infinite loop)
                            breaking = 1;
                            break;
                        }
                        sscanf(lines_dft, "%lf %lf", &x_dft, &y_dft);
                    }
                    
                } while (x_dft <= new_x_pd);
                    
                if (breaking == 1) { // will end loop if dft has reached last line (prevents infinite loop)
                    breaking = 0;
                    break; // inelegant to repeat, I know, but it works, so sue me.
                }
                    
                old_x_pd = new_x_pd; // set lower bound of checked interval for interpolation
                old_y_pd = new_y_pd;
                l_pd = 0;
                while (l_pd != 9) {
                    if((c_pd = fgetc(pd)) == '\n') l_pd++; // moves down nine lines (to the next line of current segment)
                }
                
            }
            
            j_comp = 0;
            rewind(dft);
            
            l_dft = 0;
            while (l_dft != dft_offset) {
                if((c_dft = fgetc(dft)) == '\n') l_dft++; // skips intro block (necessary to avoid weird catches/ensure consistency)
            }
            
            do {
                strcpy(lines_dft_ph,lines_dft);
                fgets(lines_dft,200,dft);
                if (strncmp(lines_dft,lines_dft_ph,200) == 0) {
                    j_comp++;
                }
            } while (j_comp != j); // moves to beginning of next band (marked by double space)
            
            
            l_dft = 0;
            while (l_dft != k*(dft_segm_lngth+1)) {
                if((c_dft = fgetc(dft)) == '\n') l_dft++; // moves down to first line of current segment in current band of dft file
            }
        }
        
        fclose(pd); // closes current pd file (in preparation of opening of next pd file
    }

}

void initialize_population() {
    
    FILE *population;
    
    population = fopen("ga.in", "w");
    
    double frac_perturb = 0.05;
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
        }
        fprintf(population,"\n");
    }
}

void processor_run(char folder[], char inputs[]) {
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[3], rho[3], B[3], lambda[2], gamma_3b[2]; //allocates variables, note/recall equilvalences in gammas/rmax_3bs
    double err_a, elast, chi_sq; // allocates objectives, namely error in a, elastic (multiply these to get single objective) and chi_sq
    double theta0; //allocates constant (for convenience, -> read/specified once and set)
    char file[200], path[200], c,dc[20];
    int bands;
    //int test;
    
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
    fprintf(gulp_input,"dispersion 1 %d \n", 183); // what is 183 here? Better as variable?
    fclose(gulp_input);

    const char *segment[3] = {"G_M","M_K","K_G"};
    
    for(i=0;i<3;i++) {
        
        strcpy(dc,segment[i]);
        mkdir(dc,S_IRWXU); // creates directory for segment
        
        strcpy(file,"in.gulp");
        strcpy(path,dc);
        strcat(path,"/");
        strcat(path,file);  // copies in.gulp into new directory
        bands = 1;
        
        file_copy(file,path,bands);
        
        chdir(dc); // enter segment directory
        
        //gulp < in.gulp > afterfit.out; // runs gulp
        
        chdir(".."); // exits new directory
    }
    
    // file collation test
    err_a = 0.123;
    elast = 0.055;
    chi_sq = 1000;
    
//    chdir(segment[0]); // arbitrary choice of segment for objectives 1 & 2 (same for all)
//    read_output(err_a,elast); // extracts objective 1 & 2: error in lattice constant a and error in elastic constant
//    chdir("..");
//
//    phonon_disp(chi_sq);  // may need to modify paths (don't know where phonon_disp() is placed by gulp)

    FILE *output;

    output = fopen("ga_line","w");

    fprintf(output,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",A[0],A[1],A[2],rho[0],rho[1],rho[2],B[0],B[1],B[2],lambda[0],lambda[1],gamma_3b[0],gamma_3b[1],err_a,elast,chi_sq);
    
    chdir("..");
    fclose(output);
    
}



int main(int argc, char *argv[]) {
    
    // First, initialize population
//    initialize_population();
    
    
    int myid =1; /* My rank */
    int nprocs = 1; /* Number of processors */
    int i,j,line; // iterators
    char c,variables[200],folder[200];
    char test[200];

    for(j=0;j<iteration_num;j++) {

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


        FILE *input;
        input = fopen("ga.in","r");
        
        if (input == NULL)
        {
            printf("Cannot open file \n");
            exit(0);
        }


        for (i=myid;i<population_num;i+=nprocs) {

            sprintf(folder,"%d",i);

            line = 0;
            fgets(test,200,input);
            rewind(input);

            while (line!=(i+1)) {
                if((c = fgetc(input)) == '\n') line++;


            }
            fgets(variables,200,input);

            processor_run(folder,variables);

        }

        MPI_Finalize();

        FILE *ga_input,*ga_line;
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
        }



        // ./ga.c
        // Change ga.c output to ga.in -> modify if necessary
    }

    return 0;
}
