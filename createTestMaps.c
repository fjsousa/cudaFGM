///////////////////////////////////////////////////////////////////////////////
//	createTestMaps.c
// 
//  An executable is created to create aspect and slope maps for the 3 test 
//	cases:
//
//			Sc0, Sc1 - aspect/slope = 0/0
//			Sc1 - aspect/slope = 203/0.5
//
// Inputs: "RunSet.in" file - Rows, Cols.
//
// Outputs: "aspect.map" and "slope.map"
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

int main (void)
{
	int row, col;
	int N, M;     				//Number of Rows (N) and Columns (M)
	FILE *IN;             //input file for rows/cols, width/height
	FILE *aspect;					//aspect, slope output maps
	FILE *slope;					//
	char buffer[100];     //buffer to use when fgets skips lines
	int n;

	//////////////////
	//Reads setup data 
	IN = fopen("RunSet.in","r");
	//skips 22 text lines of RunSet.in
	for (n = 0; n < 26; n++)
		fgets(buffer, 100, IN);
		//Reads setup data
	fscanf(IN, "%d %d", &N, &M);
	
	/////////////////////////////////////////
	//Writes aspect slope maps for difference 
	//test cases
	
	//Sc0 and Sc1
	aspect = fopen("aspect_sc0_1.map", "w");
	slope = fopen("slope_sc0_1.map", "w");
	for (row = 0; row < N; row++)
	{
		for (col = 0; col < M; col++)
		{	
			fprintf(aspect, " %4.1f ", 0.0);
			fprintf(slope, " %4.1f ", 0);
		}
		fprintf(aspect,"\n");
		fprintf(slope,"\n");
	}
	fclose(aspect); fclose(slope);

	//Sc2
	aspect = fopen("aspect_sc2.map", "w");
	slope = fopen("slope_sc2.map", "w");
	for (row = 0; row < N; row++)
	{
		for (col = 0; col < M; col++)
		{	
			fprintf(aspect, " %4.1f ", 203.0);
			fprintf(slope, " %4.1f ", 0.5);
		}
		fprintf(aspect,"\n");
		fprintf(slope,"\n");
	}
	fclose(aspect); fclose(slope);

	return (0);
}
