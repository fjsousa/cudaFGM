#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fireCudaLib.h"

int Printmap ( float*, char*, int, int);

int errorStuff(int N_isoLines, float upperIsoLine, int Rows, float CellHt, int Scenario, float* ignExact, float* ignNumeric, float totalTime)
{
	int Cols = Rows;
	float* ignError;
	float* isoTime, *SEI;
	int* A, *B;
	int row, col, cell, n;
	float Sum, Num, Den, AvgError, isoLine;
	FILE* SEIdata;
	char tmpchar[10] ="", fileString[24] = "";

	ignError   = (float*)malloc(Rows*Cols*sizeof(float));
	A          = (int*)malloc(Rows*Cols*sizeof(int));
	B          = (int*)malloc(Rows*Cols*sizeof(int));
	isoTime    = (float*)malloc(N_isoLines*sizeof(float));
	SEI        = (float*)malloc(N_isoLines*sizeof(float));
	
	//Residue and Global error
	Sum = 0;
	for (row = 0; row < Rows; row++)
	{
		for(col = 0; col < Cols; col++)
		{
			ignError[col + Cols*row] = fabs( ignExact[col + Cols*row] - ignNumeric[col + Cols*row] );
			Sum += ignError[col + Cols*row];
		}
	}
	AvgError = Sum / Cols/Rows;

	//Iso time lines loop
	Num = Den  = 0;
	for( n = 0; n < N_isoLines; n++)
	{
		isoLine = isoTime[n] = upperIsoLine/(N_isoLines ) * (n + 1);
		for (cell = 0; cell < Rows*Cols; cell++)
		{
			 A[cell] = (ignExact[cell] <= isoLine) ? 1 : 0;
			 B[cell] = (ignNumeric[cell] <= isoLine) ? 1 : 0;
		}
		
		for (cell = 0; cell < Rows*Cols; cell++)
		{
			Num += A[cell]*B[cell];
			Den += A[cell]+B[cell];
		}

		SEI[n] = Num/Den*2;
	}

	//Print Error Map
	strcat(fileString, "ignError");
	strcat(fileString, TAG);
		
	if (Scenario == 0)	
		strcat(fileString, "Sc0_");
	if (Scenario == 1)	
		strcat(fileString, "Sc1_");
	if (Scenario == 2)	
		strcat(fileString, "Sc2_");
	if (Scenario == 3)	
		strcat(fileString, "Sc3_");
	
	if (Cols <= 99)
	{
		strcat(fileString, "00");
		sprintf(tmpchar, "%2d", Cols);
	}
	else if (Cols <= 999)
	{
		strcat(fileString, "0");
		sprintf(tmpchar, "%3d", Cols);
	}
	else if (Cols <= 9999)
		sprintf(tmpchar, "%4d", Cols);

	strcat(fileString, tmpchar);
	strcat(fileString, ".dat");
  Printmap(ignError, fileString, Rows, Cols);
	
	//Print SEI data
	strcpy(fileString, "SEI");
	strcat(fileString, TAG);
	if (Scenario == 0)	
		strcat(fileString, "Sc0");
	if (Scenario == 1)	
		strcat(fileString, "Sc1");
	if (Scenario == 2)	
		strcat(fileString, "Sc2");
	if (Scenario == 3)	
		strcat(fileString, "Sc3");
	strcat(fileString, ".dat");

	SEIdata = fopen(fileString,"a");
	fprintf(SEIdata, "Avg error = %7.3f | %7.2f [%5.3f]\n", AvgError, totalTime, AvgError/totalTime*100);
	fprintf(SEIdata, "Cells = %4d; %5.2f\n", Rows, FeetToMeters(CellHt));
	for (n = 0; n < N_isoLines; n++)
		fprintf(SEIdata, "%f %f\n",  isoTime[n], SEI[n]);
	fprintf(SEIdata, "\n");
	
	//Close stuff
	free(ignExact);
	free(ignNumeric);
	free(ignError);
	free(isoTime); 
	free(SEI);
	free(A); 
	free(B);
	fclose(SEIdata);
	
	return (0);

}

int Printmap ( float* map, char* fileName, int Rows, int Cols )
{
    FILE *fPtr;
    int cell, col, row;

    if ( (fPtr = fopen(fileName, "w")) == NULL )
    {
        printf("Unable to open output map \"%s\".\n", fileName);
        return (FIRE_STATUS_ERROR);
    }

    for ( row = 0; row < Rows; row++ )
    {
        for ( cell=row*Cols, col=0; col<Cols; col++, cell++ )
        {
            fprintf(fPtr, "  %5.2f ", (map[cell]==INF) ? 000.00 : map[cell]);
        }
        fprintf(fPtr, "\n");
    }
    fclose(fPtr);
    return (FIRE_STATUS_OK);
}
