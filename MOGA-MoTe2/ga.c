#include <stdio.h>
#include<stdlib.h>
#include <math.h>

time_t t;
#define N 200
#define N_MAX 1
int nvar=32;
double sigmashare=0.2;
int alpha=2;
double crossover_P=0.9;
double mutation_P=0.1;

double fraction = 0.20;

int readinput(char fileName[], double x[][N],double y1[],double y2[], double y3[]);
int fitness(int p[],double a1[],double a2[], double a3[], int front[][N+1]);
int ranking(int p[],double r[],int nfront,int front[][N+1],double x[][N],double x0[][2]);
int selection(double r[],int ptemp[]);
int crossover(int matingpool[],double x[][N],double x0[][2],double xnew[][N]);
int mutation (int nstp,double xnew[][N],double x0[][2]);
int copy_population(double x[][N],double xnew[][N]);

int main(int argc, char *argv[]){

	FILE *fout;
	int i,j,nstp;
	int population[N],matingpool[N];
	double 	x[nvar-3][N],x0[nvar-3][2],xnew[nvar-3][N];
	double y1[N],y2[N],y3[N] ;  // objective functions =x1+x2 & y2=1/(x1+x2)
	double rank[N];
	int nfront;
	int front[N][N+1];
	double variables[29] = {3, 3, 3, 6.6, 6.6, 2.7, 2.0, 2.0, 2.0, 0.31, 0.31, 0.28, 50.0, 50.0, 50.0, 13.6, 13.6, 38.5, 32.9, 32.9, 32.9, 21.2, 21.2, 1.9, 1.9, 1.9, 1.9, 5.7, 5.7};
	srand((unsigned) time(&t));
	//Defining Limit
	for (i=0;i<nvar;i++)
	{
		if (i < 29)
		{
			x0[i][0]= variables[i]*(1 - fraction); 
			x0[i][1]= variables[i]*(1 + fraction); 
		}
	}
	fprintf(stdout, "Reading File \n"); 
	readinput(argv[1], x,y1,y2,y3);
	nstp=1;
	
	fprintf(stdout, "Performing non dominated sorting \n");
    	nfront=fitness(population,y1,y2,y3,front);

	fprintf(stdout, "Calculating ranks \n");
    	ranking(population,rank,nfront,front,x,x0);

	fprintf(stdout, "Performing binary tournament selection \n");
    	selection(rank,matingpool);

	fprintf(stdout, "Performing binary crossover \n");
    	crossover(matingpool,x,x0,xnew);

	fprintf(stdout, "Performing mutation \n");
    	mutation(nstp,xnew,x0);

	fprintf(stdout, "Copying new population to old \n");
    	copy_population(x,xnew);
    
    	printf("----GA STEP %d FINISHED\n",nstp);

	fprintf(stdout, "Updating file for further computation\n");

        char fileName[50];
	sprintf(fileName, "firstFront_%s.d", argv[2]);

        fout = fopen(fileName,"wo");
        for(i = 1; i <= front[0][0]; i++)
        {
            int n = front[0][i];
            fprintf(fout, "%d \t %f \t %f \t %f\n", n, y1[n], y2[n], y3[n]);
        }

	printf("Front printiting done\n");

	fout = fopen("value.d","wo");
	fprintf(fout,"%d\n",N);  
   	for (i=0;i<N;i++)
	{
		for (j = 0; j < nvar-3; j++)
       			fprintf(fout,"%f ",x[j][i]);
		fprintf(fout, "\n");
	}
  	fclose(fout);
	printf("Value.d file done\n");

	return 0;
}

int fitness(int p[],double a1[],double a2[],double a3[],int front[][N+1]){

	int  i,j,k,l,m,n,listN,var;
	int dominate[N],nfront;
	int partialP[N+1],Ptemp[N+1];
	int listofmember[N][N+1];

	/*----Initilazition-------*/

	for(i=0;i<N;i++)
	{
  		dominate[i]=0;
	}

	for(i=0;i<N;i++)
  		for(j=0;j<N+1;j++)
		{
    			listofmember[i][j]=0;
     			front[i][j]=0;
		}

	for(i=0;i<N+1;i++)
	{
  		Ptemp[i]=0;
  		partialP[i]=0;
	}

	/*--------------------------*/
	for(i=0;i<N;i++)
	{
   		for(j=i+1;j<N;j++)
		{
     			// Keeping third objective constant and comparing the first 2
     			if(a1[i]>=a1[j] && a2[i]==a2[j] && a3[i] == a3[j])
			{
       				dominate[i]+=1;
       				m=listofmember[j][0]+1;
       				listofmember[j][m]=i;
       				listofmember[j][0]+=1;
     			}
     			else if (a1[i]==a1[j] && a2[i]>=a2[j] && a3[i] == a3[j])
     			{
				dominate[i]+=1;
				m=listofmember[j][0]+1;
				listofmember[j][m]=i;
				listofmember[j][0]+=1;
     			}
     			else if (a1[i] <= a1[j] && a2[i] == a2[j] && a3[i] == a3[j])
     			{
				dominate[j]+=1;
				m=listofmember[i][0]+1;
				listofmember[i][m]=j;
				listofmember[i][0]+=1;
     			}
     			else if (a1[i]==a1[j] && a2[i] < a2[j] && a3[i] == a3[j])
     			{
				dominate[j]+=1;
				m=listofmember[i][0]+1;
				listofmember[i][m]=j;
				listofmember[i][0]+=1;
     			}
    			else if (a1[i] < a1[j] && a2[i] < a2[j] && a3[i] == a3[j])
    			{
       				dominate[j]+=1;
       				m=listofmember[i][0]+1;
       				listofmember[i][m]=j;
       				listofmember[i][0]+=1;
     			}
    			else if (a1[i] > a1[j] && a2[i] > a2[j] && a3[i] == a3[j])
    			{
       				dominate[i]+=1;
       				m=listofmember[j][0]+1;
       				listofmember[j][m]=i;
       				listofmember[j][0]+=1;
     			}
     			// Keeping the first objective constant and changing the other two
     			else if (a1[i]==a1[j] && a2[i]==a2[j] && a3[i] > a3[j])
     			{
        			dominate[i]+=1;
        			m=listofmember[j][0]+1;
        			listofmember[j][m]=i;
        			listofmember[j][0]+=1;
     			}
     			else if (a1[i]==a1[j] && a2[i] == a2[j] && a3[i] < a3[j])
     			{
        			dominate[j]+=1;
        			m=listofmember[i][0]+1;
        			listofmember[i][m]=j;
        			listofmember[i][0]+=1;
     			}
    			else if (a1[i] == a1[j] && a2[i] < a2[j] && a3[i] < a3[j])
    			{
       				dominate[j]+=1;
       				m=listofmember[i][0]+1;
       				listofmember[i][m]=j;
       				listofmember[i][0]+=1;
     			}
    			else if (a1[i] == a1[j] && a2[i] > a2[j] && a3[i] > a3[j])
			{
       				dominate[i]+=1;
       				m=listofmember[j][0]+1;
       				listofmember[j][m]=i;
       				listofmember[j][0]+=1;
     			}
    			// Keeping the second objective constant and changing the other two
    			else if (a1[i] < a1[j] && a2[i] == a2[j] && a3[i] < a3[j])
    			{
       				dominate[j]+=1;
       				m=listofmember[i][0]+1;
       				listofmember[i][m]=j;
       				listofmember[i][0]+=1;
     			}
    			else if (a1[i] > a1[j] && a2[i] == a2[j] && a3[i] > a3[j])
			{
       				dominate[i]+=1;
       				m=listofmember[j][0]+1;
       				listofmember[j][m]=i;
       				listofmember[j][0]+=1;
     			}
    			// Changing all the objectives
    			else if (a1[i] < a1[j] && a2[i] < a2[j] && a3[i] < a3[j])
    			{
    	 			dominate[j]+=1;
	 			m=listofmember[i][0]+1;
	 			listofmember[i][m]=j;
	 			listofmember[i][0]+=1;
    			}
    			else if (a1[i] > a1[j] && a2[i] > a2[j] && a3[i] > a3[j])
    			{
	 			dominate[i]+=1;
				m=listofmember[j][0]+1;
	 			listofmember[j][m]=i;
				listofmember[j][0]+=1;
    			}

   		}
	}

	for(i=0;i<N;i++)
	{
  		if(dominate[i]==0)
		{
       			m=Ptemp[0]+1;
       			Ptemp[m]=i;
       			Ptemp[0]+=1;
       		}
	}

	nfront=0;
	while(Ptemp[0]> 0)
	{
   		for(i=1;i<=Ptemp[0];i++)
		{
      			front[nfront][0]+=1;
      			m=front[nfront][0];
      			n=Ptemp[i];
      			front[nfront][m]=n;
      			listN=listofmember[n][0];
      			if(listN > 0)
			{
      				for(j=1;j<=listN;j++)
				{
         				var=listofmember[n][j];
         				dominate[var]-=1;
         				if(dominate[var]==0)
					{
            					k=partialP[0]+1;
            					partialP[k]=var;
            					partialP[0]+=1;
          				}
        			}
      			}
   		}
 		Ptemp[0]=partialP[0];
 		for(l=1;l<=partialP[0];l++)
		{
    			Ptemp[l]=partialP[l];
  		}
 		partialP[0]=0;
 		nfront+=1;
	}
	return(nfront);
}

int ranking(int p[],double r[],int nfront,int front[][N+1],double x[][N],double x0[][2])
{

	int i,j,k,l,m,n;
	double maxrank=N,min_rank;
	double dmn,niche[N];

	//printf("---Ranking starts---\n");
	for(i=0;i<N;i++)
   		niche[i]=0;

	for(i=0;i<nfront;i++)
	{
 		for(j=1;j<=front[i][0];j++)
		{
    			m=front[i][j];
    			r[m]=maxrank;
		//    printf("%d ",m,r[m]);
 		}
		//  printf("\n before adjustment\n");
 		for(j=1;j<=front[i][0];j++)
		{
     			m=front[i][j];
   			for(k=j;k<=front[i][0];k++)
			{
        			n=front[i][k];
        			dmn=0;
        			for(l=0;l<nvar - 3;l++)
				{
            				dmn+=pow(((x[l][m]-x[l][n])/(x0[l][1]-x0[l][0])),2);
      				}
      				dmn=sqrt(dmn);
				//      printf("%d %d %f\n",m,n,dmn);
      				if(dmn <= sigmashare)
				{
        				niche[m]+=(1-pow((dmn/sigmashare),alpha));
        				niche[n]+=(1-pow((dmn/sigmashare),alpha));
       				}
   			}
 		}
		//printf("Next front\n");
		min_rank=999999999999;
		for(j=1;j<=front[i][0];j++)
		{
   			m=front[i][j];
   			r[m]/=niche[m];
   			if(r[m]<min_rank) min_rank=r[m];
 		}
 		maxrank=min_rank;
	}
	
	return(1);
}
// Roulette Selection Operator
int selection(double r[],int ptemp[])
{

	int  i,j,n;
	int xlow,xhigh;
	double rsum=0,m,rprob[N],raprob[N];

	for (i=0;i<N;i++) rsum+=r[i];
	for(i=0;i<N;i++)  rprob[i]=r[i]/rsum;
	raprob[0]=rprob[0];
	for(i=1;i<N;i++) raprob[i]=raprob[i-1]+rprob[i];

	for (i=0;i<N;i++)
	{
  		m=(double)rand()/(double)RAND_MAX;
  		xlow=0;
  		xhigh=N;
  		n=N/2;
  		while(xlow < xhigh)
		{
  			if(m>raprob[n-1] && m<=raprob[n])
			{
    				break;
   			}
   			else if (m > raprob[n])
			{
    				xlow=n+1;
    				n=(int)(xlow+xhigh)/2;
   			}
   			else if(m < raprob[n])
			{
    				xhigh=n-1;
    				n=(int)(xlow+xhigh)/2;
   			}
  		}
 		ptemp[i]=n;
	}

	return(1);
}

int crossover(int matingpool[],double x[][N],double x0[][2],double xnew[][N])
{
	int i,j=0,count=0,k,l,m;
	double alpha=0.5,ui,pi,yi,tag;
	double temp,x1[nvar-3],x2[nvar-3];

 	for(m=0;m<N;m++)
	{
   		fprintf(stdout, "Population = %d\n", m+1);
   		i=matingpool[m];
   		pi=(double)rand()/(double)RAND_MAX;
   		if(pi<=crossover_P)
		{
     			count+=1;
     			if(count==1)
			{
     				for(l=0;l<nvar-3;l++) x1[l]=x[l][i];
     			}
     			else if(count==2)
			{
       				for(l=0;l<nvar-3;l++) x2[l]=x[l][i];
     
       				for(l=0;l<nvar-3;l++)
				{
         				if(x1[l]>x2[l])
					{
            					temp=x2[l];
            					x2[l]=x1[l];
            					x1[l]=temp;
         				}
       				}
      	
				k=0;
      				while(k<2)
				{ 
         				for(l=0;l<nvar-3;l++) 
					{
            					tag=0;
             					while(tag<1)
						{
                 					ui=(double)rand()/(double)RAND_MAX;
                 					yi=(1+2*alpha)*ui-alpha;
                 					xnew[l][j]=(1-yi)*x1[l]+yi*x2[l];
                 					if(xnew[l][j] >= x0[l][0] && xnew[l][j]<= x0[l][1]) {tag=1;}
           					}
         				}
         				j++;
         				k++;
      				}
     				count=0;
     			}
    		}
   		else 
		{
        		for(l=0;l<nvar-3;l++) xnew[l][j]=x[l][i];
        		j++;
     		}
  	}
	if(j==N-1)
	{
   		i=matingpool[0];
   		for(l=0;l<nvar-3;l++) xnew[l][j]=x[l][i];
   		j++;
	}
//printf("new size %d \n",j); 
	return 1;
}

int mutation (int nstp,double xnew[][N],double x0[][2])
{

	int i,l,tu;
	double pi,ui,ri,temp,b=2,n;

	n = pow((1-(nstp/N_MAX)),b);

	for(i =0;i < N;i++)
	{
   		pi=(double)rand()/(double)RAND_MAX;
   		if(pi < mutation_P)
		{
      			for(l=0;l<nvar-3;l++)
			{
          			ui=(double)rand()/(double)RAND_MAX;
          			ri=(double)rand()/(double)RAND_MAX;
          			if(ui <= 0.5) tu = -1;
          			if(ui  > 0.5) tu = 1;
          			temp=xnew[l][i]+tu*(x0[l][1]-x0[l][0])*(1-pow(ri,n));
         			if(temp>x0[l][0] && temp < x0[l][1])
				{
            				xnew[l][i]=temp;
         			}
        		}
   		}
	}
	return(1);
}

int copy_population(double x[][N],double xnew[][N])
{
	int i,l;

	for(i=0;i<N;i++)
  		for(l=0;l<nvar-3;l++)
      			x[l][i]=xnew[l][i];
	return(1);
}  

int readinput(char fileName[], double x[][N],double y1[],double y2[],double y3[]){

	FILE *fin;
	int i, j, num;
	float t[nvar];
	fin = fopen(fileName,"r");
	fscanf(fin, "%d",&num);
	for(i=0;i<N;i++){
		for (j = 0; j < nvar; j++)
			fscanf(fin, "%f", &t[j]);	

      		for (j = 0; j < nvar - 3; j++)
			x[j][i] = t[j];

      		y1[i]=t[nvar - 3];
      		y2[i]=t[nvar - 2];
		y3[i]=t[nvar - 1];
		for (j = 0; j < nvar-3; j++)
			printf("%f ", x[j][i]);
		printf("%f %f %f\n", y1[i], y2[i], y3[i]);
		//printf("%f %f %f %f %f %f %f \n", x[0][i], x[1][i], x[2][i],x[3][i], y1[i], y2[i], y3[i]);
  	}
	fclose(fin);

	return(0);
}    
