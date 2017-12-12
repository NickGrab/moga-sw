#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>




//float rho[3] = {0.345,0.3241,1.567};
//float B[3] = {3.4546,0.456,2.976};
//double rmax_2b[3];

//float lambda[2] = {1.34,7};
int population_num = 20;
double theta0 = 81.9966; // necessary dummy variable for unchanging variable. Worth knowing globally
//double gamma_3b[2] = {1.21,3.655};
//double rmax_3b[2];


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

double *read_output(double obj[]) {

    FILE *af;
    char output[] ="afterfit.out";
    af = fopen(output,"r"); // read output file
    
    double c_11, c, el_c;
    double lat_const = 12.29, elast_real = 238; // elastic constant in GPa

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
    sscanf(&lines[len-8], "%lf",&obj[0]); // assigns value of percent error in a
    
    line_down = 0;
    while (line_down !=1) { // moves down to elastic constant "c" (for elastic fit)
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    char cw; // dummy throwaway
    double c_init; //dummy throwaway
    sscanf(lines,"%s %lf %lf", cw, &c_init, &c); // assigns value of elastic constant "c" (c_final)

    line_down = 0;
    while (line_down !=21) { // moves down to elastic constant C11
        if((ch = fgetc(af)) == '\n') line_down++;
    }
    
    fgets(lines,200,af);
    int throw; //throwaway int
    sscanf(lines, "%d %lf",throw,&c_11); // assigns value of elastic constant C11
    
    elast_calc = c_11/c*(2/lat_const); // elastic constant for comparison
    obj[1] = abs(elast_calc-elast_real)/elast_real*100; // percent error in elastic constant
    
    fclose(af);
    
    return obj;

}

void phonon_disp() {
    
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
    
    population = fopen("population.txt", "w");
    
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

            fprintf(population, "%lf ",rand_var[j]);
        }
        fprintf(population,"\n");
    }
}

void processor_run(char folder[], char inputs[]) {
    
    //initialize variables (make global only if using private copies in this subroutine (OpenMPI))
    double A[3], rho[3], B[3], lambda[2], gamma_3b[2]; //allocates variables, note/recall equilvalences in gammas/rmax_3bs
    double *obj[3]; // allocates objectives, namely error in a, elastic (multiply these to get single objective) and chi_sq
    double obj_place[3];
    double err_a,c_final,elastic, chi_sq; // allocates objectives
    double theta0; //allocates constant (for convenience, -> read/specified once and set)
    char file[200], path[200], c,dc[20];
    int bands;
    int test;
    
    // create folder for this instance of gulp run (now done in first copy cell step)
    mkdir(folder,S_IRWXU);
    // copy cell, forcefield, etc. into folder
    
    strcpy(file,"cell");
    strcpy(path,folder);
    strcat(path,"/");
    strcat(path,file); // creates path for new file;
    
    file_copy(file,path,bands);
    
    strcpy(file,"forcefield");
    strcpy(path,folder);
    strcat(path,"/");
    strcat(path,file); // creates path for new file;
    
    file_copy(file,path,bands);

    //enter folder
    chdir(folder);


    // do any setup for gulp to run in folder
    //system("npm install gulp"); // for my computer only, different for hpc?, universal implementation possible?
    // read inputs line into variable (will need dummies or local)
    modify_input(A,rho,B,lambda,gamma_3b);

    FILE *gulp_input;

    gulp_input = fopen("in.gulp","w");
    fprintf(gulp_input, "optim relax conp comp phon nofreq\n");
    fclose(gulp_input);
    
    strcpy(file,"cell");
    strcpy(path,"in.gulp");
    bands = 0;

    file_copy(file,path,bands);
    
    strcpy(file,"forcefield");
    file_copy(file,path,bands);
   
    gulp_input = fopen("in.gulp","a");
    fprintf(gulp_input,"dispersion 1 %d \n", 183); // what is 183 here? Better as variable?
    fclose(gulp_input);


    strcpy(dc,"G_M");
    mkdir(dc,S_IRWXU);

    strcpy(file,path);
    strcpy(path,dc);
    strcat(path,"/");
    strcat(path,file); // creates path for new file;
    bands = 1;

    file_copy(file,path,bands);

    strcpy(dc,"M_K");
    mkdir(dc,S_IRWXU);

    strcpy(path,dc);
    strcat(path,"/");
    strcat(path,file); // creates path for new file;
    bands = 2;
    file_copy(file,path,bands);
    
    strcpy(dc,"K_G");
    mkdir(dc,S_IRWXU);

    strcpy(path,dc);
    strcat(path,"/");
    strcat(path,file); // creates path for new file;
    bands = 3;
    file_copy(file,path,bands);
    
    // run gulp in G_M, M_K, and K_G. Use OpenMP threads? // (LOOK AT ME PLEASE!)
    // inside loops/threads:
        // any necessary setup for in.gulp
            // gulp < in.gulp > afterfit.out
        // get objectives from afterfit.out -> choose any one afterfit.out from bands
            // obj = read_output(obj_place);
                // pretty sure this will work to get the array over. arrays in c are ugly :/ *Especially* when they have to work in parallel
        // get chi^2 value
            // phonon_disp();
                // may need to modify paths (don't know where phonon_disp() is placed by gulp)
                // 98% positive it will work otherwise
    // END
    
    
}



int main() {
// TESTING
    //read_output();
    //phonon_disp();
    //initialize_population();
    //processor_run("folder","3.4 5.768 3.45 0.35 1.67 2.334 4.345 0.345 0.3454 1.234 1.234 2.11 0.551");
// END
    
    // First, initialize population
    initialize_population();
    
    // Create ga.in file
    FILE *ga_in;
    ga_in = fopen("ga.in","w"); // creates first ga.in file (needs to be created outside parallelization loop)
    fprintf(ga_in,"%d\n", population_num); // first line: number of trials
    
    int i;
    
    for (i=0;i<(population_num),i++) {
        fprintf(ga_in,"\n"); // creates as many blank spaces in document as there will be lines (needed for parallelization loop)
    }
    
    
    // parallelize over processor_run
    // MPI assignment loop from 0 to population_num - 1
        // read line by file -> assign line to array
        // assign folder name
            // ex: run_1, run_2, etc.
        // call processor_run(run_x,array)
    // end parallelization
    // consoladate all lines into an input for ga.c (Ankit's file)
        // ^^ Need help figuring out how to do this off a parallel script. Can't figure out a solution that doesn't lead to race conditions
        // any attempt to write directly in the loop risks the overwriting each other. But we can't allocate one big 2D array/array of pointers because they aren't guaranteed to be one consistent size. We may need to use padding to change that to make this work, but I wanted to see if you guys had any ideas.
    // run ga.c
    // check if conditions are satisfied (I think ga.c does this?)
    // send output back to MPI for a second loop
    
    return 0;
}
