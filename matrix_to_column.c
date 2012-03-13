//Reads 2D square matrix data file, writes column (1D) data file
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main (int argc, char *argv[])
{
	int row, col;
	int side = atoi(argv[1]);
	float tmp_storage;
	char *file_twoD = argv[2];
	char *file_oneD = malloc(strlen("1D_") + strlen(file_twoD) + 1); 
	FILE* p_oned, *p_twod;
	strcpy(file_oneD, "1D_");
	strcat(file_oneD, file_twoD);
 
	p_twod = fopen(file_twoD, "r");
	p_oned = fopen(file_oneD, "w");
	
	for ( row = 0; row < side; row++){
		for ( col = 0; col < side; col++){
			fscanf(p_twod, "%f", &tmp_storage);
			fprintf(p_oned,"%d %d %7.2f\n ",row, col, tmp_storage);
		}
	}
	
	return(0);
}
