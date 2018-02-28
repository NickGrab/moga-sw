#include <stdio.h>
#include <stdlib.h>
#include<string.h>
int main()
{
	char a[200];
	FILE *bands = fopen("DFT/band.dat", "r");
	while(1)
	{
		fgets(a,200, bands);
		
		if (strcmp(a,"\n") < 5)
		{
			printf("%s", a);
			break;
		}
		printf("%s", a);
	}
	return 0;
}
