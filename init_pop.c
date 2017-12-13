//
// Identical to function of same name in proper .c file
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main initialize_population() {
    
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
    return 0;
}
