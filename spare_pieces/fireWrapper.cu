#define WithCUBLAS 1

/*
 *******************************************************************************
 *
 *  fireWrapper.c
 *
 *  Description
 *      Cellular automata like algorithm using fireLib modified Cuda functions.
 *		Main cicle files, kernel wraper and test case of a fire spread model.
 *
 *******************************************************************************
 */

#include "fireCudaLib.h"
#include <time.h>
#if WithCUBLAS
#include <cublas.h>
#endif

////////////////////////////////////////////////////////////////////////////////
// Other Functions 
////////////////////////////////////////////////////////////////////////////////
float* ExactCircle( int, int, float, float);
float* ExactElipse( int, int, float, float, float, float);
int errorStuff(int, float, int, float, int, float*, float*, float);
int PrintMap ( float*, char* );
//float *getCPUign( void );

/*timing structures*/
typedef struct{
	float total, init, precond, fireSpreadNo, fireSpread, examiningNeighbor, cublas;
} timeStructure;

typedef struct{
	    clock_t total, init, precond, fireSpreadNo, fireSpread, examiningNeighbor, cublas;
} clockTic;

//Global variables
static size_t Model;      								/* NFFL 1 */
static float WindSpd;     								/* mph */
static float WindDir;   									/* degrees clockwise from north */
static float Slope;   										/* fraction rise / reach */
static float Aspect;  	 									/* degrees clockwise from north */
static float M1;
float mapW;
float mapH;
static float CellWd;				  				  	/* Cell width (E-W) in feet. */
static float CellHt;					    				/* Cell height (N-S) in feet. */
int Rows;   															/* Number of rows in each map. */
int Cols;    															/* Number of columns in each map. */
float* modelArray;					//Array for Models

int main ( int argc, char **argv )
{
	unsigned int device;
	unsigned int Scenario;
  FuelCatalogPtr catalog;     /* fuel catalog handle */
  float* modelArray_d;				//Array of continuous memory that holds the Model data in the device 
	static int n_itt = 0;      /* counter for number of time steps */
	int row, col, ignRow, ignCol; 
	int cell, cells;  					/* neighbor index, total number of map cells */
  int N_isoLines;
	float upperIsoLine;
	float SpreadMax_pre, Residue = INF, Residue_max, simTime;
	float DistHV;
	float* exactMap; /**cpuMap, *errorCPGP;*/
  size_t *fuelMap;            /* ptr to fuel model map */
  float *ignMap;             	/* ptr to ignition time map (minutes) */
  float *ignMap_new;             	/* ptr to ignition time map (minutes) */
  float *slpMap;             	/* ptr to slope map (rise/reach) */
  float *aspMap;             	/* ptr to aspect map (degrees from north) */
  float *wspdMap;            	/* ptr to wind speed map (ft/min) */
  float *wdirMap;            	/* ptr to wind direction map (deg from north) */
  float *m1Map;              	/* ptr to 1-hr dead fuel moisture map */
  size_t *fuelMap_d;         	// ptr to Cuda arrays
  float *ignMap_d;           	//...
  float *ignMap_new_d;           	//...
	float *slpMap_d;       	   	//...
  float *aspMap_d;           	//...
  float *wspdMap_d;          	//...
  float *wdirMap_d;          	//...
  float *m1Map_d;            	//...
	float *spread0Map_d;				//
	float *RxIntensityMap_d;		//
	float *spreadMaxMap_d;
	float *phiEffWindMap_d;
	float *eccentricityMap_d;
	float* diff_d;
	float *azimuthMaxMap_d;
	float Spread0, SpreadMax;   
	float eccentricity, azimuthMax;
	clockTic start, end;
	timeStructure times;
	char fileString[25] = "";
	char tmpchar[6]  = "";
	FILE* IN, *timeData;
	FILE *slope_file, *aspect_file;
  float* azimuthMaxMap;
  float* eccentricityMap;
  float* phiEffWindMap;
	float* RxIntensityMap;

	
	
	start.total = clock();
	
	IN = fopen("RunSet.in","r");
	fscanf(IN, "%d %d %d %d %f %d %f", &device, &Scenario, &Rows, &Cols, &Residue_max, &N_isoLines, &upperIsoLine);
	printf("\n\nStarting model "TAG" with Scn#%d, Node Size = %4dx%4d.\n",Scenario, Rows, Cols);
	printf("\nTurn messages off for optimal performance\n");
	
	cudaSetDevice(device);
	
	//Grid dimensions
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(Rows/BLOCK_SIZE, Cols/BLOCK_SIZE);
	
	
	start.cublas = clock();
	
	#if WithCUBLAS
	//Init Cublas
	if( cublasInit() != CUBLAS_STATUS_SUCCESS)
	{
		printf("\nERROR: CUBLAS initialization!\n");
		return (1);
	}
	else
		printf("\nCUBLAS initialized\n");
	#endif
	end.cublas = clock();
	times.cublas = ((float) (end.cublas - start.cublas))/CLOCKS_PER_SEC;

	/* NOTE 3: allocate all the maps. */
  cells = Rows * Cols;
	if ( (ignMap   					= (float *) calloc(cells, sizeof(float)))   == NULL
		|| (ignMap_new   			= (float *) calloc(cells, sizeof(float)))   == NULL
		|| (azimuthMaxMap   	= (float *) calloc(cells, sizeof(float)))   == NULL
		|| (eccentricityMap   = (float *) calloc(cells, sizeof(float)))   == NULL
		|| (phiEffWindMap   	= (float *) calloc(cells, sizeof(float)))   == NULL
		|| (RxIntensityMap   	= (float *) calloc(cells, sizeof(float)))   == NULL
  	|| (slpMap   = (float *) calloc(cells, sizeof(float)))   == NULL
    || (aspMap   = (float *) calloc(cells, sizeof(float)))   == NULL
    || (wspdMap  = (float *) calloc(cells, sizeof(float)))   == NULL
    || (wdirMap  = (float *) calloc(cells, sizeof(float)))   == NULL
    || (m1Map    = (float *) calloc(cells, sizeof(float)))   == NULL
    || (fuelMap  = (size_t *) calloc(cells, sizeof(size_t))) == NULL)
	{
        fprintf(stderr, "Unable to allocate maps with %d cols and %d rows"
		"in Host memory.\n", Cols, Rows);
        return (1);
	}

	if (Scenario == 0)
  {
		Model  = 1;      			
		WindSpd = 0;     			
		WindDir = 0;     			
		Slope   = 0;    				
		Aspect  = 0.0;    		
		M1      = 0.05;    			
		mapW = MetersToFeet(7000);
		mapH = MetersToFeet(7000);		
		CellWd = mapW/(Cols -1);				  				  	
		CellHt = mapH/(Rows -1);					    				
		
		/* NOTE 4: initialize all the maps -- modify them as you please. */
		for ( cell=0; cell<cells; cell++ )
  	{
			fuelMap[cell]  = Model;
			slpMap[cell]   = Slope;
			aspMap[cell]   = Aspect;
			wspdMap[cell]  = 88. * WindSpd;      //convert mph into ft/min 
			wdirMap[cell]  = WindDir;
			m1Map[cell]    = M1;
			ignMap[cell]   = 1000;
			ignMap_new[cell]   = 1000;
  	}

		/* NOTE 5: set an ignition time & pattern (this ignites the middle cell). */
		ignCol = Cols/2;
		ignRow = Rows/2;
		cell = Cols/2 + Cols*(Rows/2);
		ignMap[cell] = 0;	
		ignMap_new[cell] = 0;	
	}

	else if (Scenario == 1)
	{
		Model  = 1;      		
		WindSpd = 5;     		
		WindDir = 45;    	
		Slope   = 0;  		
		Aspect  = 0; 
		M1      = 0.05;   	
		mapW = MetersToFeet(7000);
		mapH = MetersToFeet(7000);
		CellWd = mapW/(Cols -1);				  				  	
		CellHt = mapH/(Rows -1);					    				

		for ( cell=0; cell<cells; cell++ )
  	{
			fuelMap[cell]  = Model;
			slpMap[cell]   = Slope;
			aspMap[cell]   = Aspect;
			wspdMap[cell]  = 88. * WindSpd;      //convert mph into ft/min 
			wdirMap[cell]  = WindDir;
			m1Map[cell]    = M1;
			ignMap[cell]   = 1000;
			ignMap_new[cell]   = 1000;
  	}

		ignCol = Cols/4;
		ignRow = Rows*3/4;
		cell = Cols*1/4 + Cols*(Rows*3/4);
		ignMap[cell] = 0;	
		ignMap_new[cell] = 0;	
	}
	else if (Scenario == 2)
	{
		Model  = 1;      		
		WindSpd = 5;     		
		WindDir = 45.;    	
		Slope   = 0.5;  		
		Aspect  = 23 + 180; 
		M1      = 0.05;   	
		mapW = MetersToFeet(7000);
		mapH = MetersToFeet(7000);
		CellWd = mapW/(Cols -1);				  				  	
		CellHt = mapH/(Rows -1);					    				

		for ( cell=0; cell<cells; cell++ )
  	{
			fuelMap[cell]  = Model;
			slpMap[cell]   = Slope;
			aspMap[cell]   = Aspect;
			wspdMap[cell]  = 88. * WindSpd;      //convert mph into ft/min 
			wdirMap[cell]  = WindDir;
			m1Map[cell]    = M1;
			ignMap[cell]   = 1000;
			ignMap_new[cell]   = 1000;
  	}

		ignCol = Cols/4;
		ignRow = Rows*3/4;
		cell = Cols*1/4 + Cols*(Rows*3/4);
		ignMap[cell] = 0;	
		ignMap_new[cell] = 0;	
	}

	//GOST SCENARIO only on Selial_Lean
	
	else if (Scenario == 3)
	{
		Model  = 1;      
		WindSpd = 0;//0;     
		WindDir = 0;   
		M1      = 0.1;
		mapW = MetersToFeet(7000);
		mapH = MetersToFeet(7000);	
		CellWd = mapW/(Cols -1);				  				  	
		CellHt = mapH/(Rows -1);					    				
	
		createMaps(Rows + 1, Cols + 1, (double) mapW, (double) mapH);
		slope_file = fopen("slope.map","r");
		 aspect_file = fopen("aspect.map","r");

		for ( row = 0; row < Rows; row++ )
		{
  		for ( col = 0; col < Cols; col ++)
			{
				cell = col + row*Cols;
		
        fscanf(slope_file, " %f ", &slpMap[cell]);
				fscanf(aspect_file, " %f ", &aspMap[cell]);

      	fuelMap[cell]  = Model;
				wspdMap[cell]  = 88. * WindSpd; //convert mph into ft/min
				wdirMap[cell]  = WindDir;
				m1Map[cell]    = M1;
      	ignMap[cell]   = INF;
  		}
		}
		ignCol = Cols/4;
		ignRow = Rows/4;
		cell = Cols/4 + Cols*Rows/4;
  	ignMap[cell] = 0.0;
  }
	
	/* NOTE 6: create a standard fuel model catalog and a flame length table. */
	catalog = Fire_FuelCatalogCreateStandard("Standard", 13);
	#if Clatering	
	printf("\nCatalog initialized on Host.\n");
	#endif
	///////////////////////////////////////////////////////////////////////////
	//Preprocesing Stage
 	Preprocessor(catalog);

	DistHV = CellWd ;
	///////////////////////////////////////////////////////////////////////////
	//Store catalog in Global Cuda memory
	modelArray_d = CudaModelData(catalog);
	#if Clatering	
	printf("\nCatalog passed to Device.\n");
	#endif
	///////////////////////////////////////////////////////////////////////////
	//Preconditioned ignamap
	start.precond = clock();
	if (Scenario > 1)
	{
		SpreadMax_pre = NoWindNoSlope(0, m1Map[0], fuelMap[0], modelArray, RxIntensityMap );
		SpreadMax_pre = WindAndSlope(	0,
																	modelArray,
																	SpreadMax_pre, 
																	wspdMap[0], 
																	wdirMap[0], 
																	0, 
																	0, 
																	fuelMap[0], 
																	RxIntensityMap,
																	phiEffWindMap,
																	eccentricityMap,
																	azimuthMaxMap);


		ignMap = ExactElipse( ignRow, ignCol, 0, SpreadMax_pre, eccentricityMap[0], azimuthMaxMap[0]);
	}
	end.precond = clock();
	times.precond = ((float) (end.precond - start.precond))/CLOCKS_PER_SEC;
	///////////////////////////////////////////////////////////////////////////
	//Spread No Wind No Slope
	start.fireSpreadNo = clock();
	
	
	cudaMalloc((void**)&m1Map_d,   cells*sizeof(float));
	cudaMalloc((void**)&fuelMap_d, cells*sizeof(size_t));
	
	cudaMalloc((void**)&spread0Map_d, 	  cells*sizeof(float));
	cudaMalloc((void**)&RxIntensityMap_d, cells*sizeof(float));

	cudaMemcpy( m1Map_d,    m1Map,   cells*sizeof(float),  cudaMemcpyHostToDevice);	
	cudaMemcpy( fuelMap_d,  fuelMap, cells*sizeof(size_t), cudaMemcpyHostToDevice);
	
	FireKernel_NoWindNoSlope<<<dimGrid, dimBlock>>>(m1Map_d, RxIntensityMap_d, spread0Map_d, fuelMap_d, modelArray_d);

	/*  DEBUG DEBUG DEBUG DEBUG DEBUG 
	cudaMemcpy( SpreadMap, spread0Map_d,   cells*sizeof(float),  cudaMemcpyDeviceToHost);	
	PrintMap(SpreadMap, "debugSpread0.txt");
	DEBUG DEBUG DEBUG DEBUG DEBUG    */
	
	cudaFree( m1Map_d );
	
	checkCUDAError("No Wind No Slope Kernel");

	end.fireSpreadNo = clock();
	times.fireSpreadNo = ((float) (end.fireSpreadNo - start.fireSpreadNo))/CLOCKS_PER_SEC;
	#if Clatering	
	printf("\nFire Spread No Wind No Slope Done.\n");
	#endif
	///////////////////////////////////////////////////////////////////////////
	//Spread Wind And Slope
	start.fireSpread = clock();
	cudaMalloc((void**)&slpMap_d,   cells*sizeof(float));
	cudaMalloc((void**)&aspMap_d,   cells*sizeof(float));
	cudaMalloc((void**)&wspdMap_d,  cells*sizeof(float));
	cudaMalloc((void**)&wdirMap_d,  cells*sizeof(float));
	
	cudaMalloc((void**)&spreadMaxMap_d,    cells*sizeof(float));
	cudaMalloc((void**)&phiEffWindMap_d,   cells*sizeof(float));
	cudaMalloc((void**)&eccentricityMap_d, cells*sizeof(float));
	cudaMalloc((void**)&azimuthMaxMap_d,   cells*sizeof(float));

	cudaMemcpy( slpMap_d,   slpMap,  cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( aspMap_d,   aspMap,  cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( wspdMap_d,  wspdMap, cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( wdirMap_d,  wdirMap, cells*sizeof(float), cudaMemcpyHostToDevice);

	FireKernel_WindAndSlope<<<dimGrid, dimBlock>>>(	wspdMap_d, 
																									wdirMap_d, 
																									slpMap_d, 
																									aspMap_d, 
																									spread0Map_d, 
																									spreadMaxMap_d,
																									RxIntensityMap_d,
																									phiEffWindMap_d,
																									eccentricityMap_d,
																									azimuthMaxMap_d,
																									fuelMap_d, 
																									modelArray_d);

	/*  DEBUG DEBUG DEBUG DEBUG DEBUG
	cudaMemcpy( azimuthMaxMap,   azimuthMaxMap_d,   cells*sizeof(float),  cudaMemcpyDeviceToHost);	
	cudaMemcpy( eccentricityMap, eccentricityMap_d, cells*sizeof(float),  cudaMemcpyDeviceToHost);	
	cudaMemcpy( phiEffWindMap,   phiEffWindMap_d,   cells*sizeof(float),  cudaMemcpyDeviceToHost);	
	cudaMemcpy( SpreadMap,       spreadMaxMap_d,   	cells*sizeof(float),  cudaMemcpyDeviceToHost);	
	PrintMap(SpreadMap, "debugSpreadMax.txt");
	PrintMap(azimuthMaxMap,   "azimuthmaxmap.txt");
	PrintMap(eccentricityMap, "eccentricitymap.txt");
	PrintMap(phiEffWindMap,   "phieffwindmap.txt");
	DEBUG DEBUG DEBUG DEBUG DEBUG    */
	
	cudaFree( wspdMap_d  );
	cudaFree( wdirMap_d  );
	cudaFree( slpMap_d  );
	cudaFree( aspMap_d  );
	cudaFree( RxIntensityMap_d  );
	cudaFree( fuelMap_d );

  checkCUDAError("Wind & Slope Kernel");

	end.fireSpread = clock();
	times.fireSpread = ((float) (end.fireSpread - start.fireSpread))/CLOCKS_PER_SEC;
	#if Clatering	
	printf("\nFire Spread With Wind and Slope Done.\n");
	#endif
	
	///////////////////////////////////////////////////////////////////////////
	//Spread at Neighbors
	start.examiningNeighbor = clock();
	cudaMalloc((void**)&ignMap_d,      cells*sizeof(float));
	cudaMalloc((void**)&ignMap_new_d,  cells*sizeof(float));
	cudaMalloc((void**)&diff_d,  			 cells*sizeof(float));
	
	cudaMemcpy( ignMap_d,     ignMap,     cells*sizeof(float),        cudaMemcpyHostToDevice);
	cudaMemcpy( ignMap_new_d, ignMap_new, cells*sizeof(float),        cudaMemcpyHostToDevice);
	

	//Main loop: iterates in time until infinity or fire reached the edge 
	#if Clatering	
	printf("\nStarting Fire growth...\n");
	#endif
	while ( Residue > Residue_max)
	{
		n_itt++;
		
		FireKernel_SpreadAtNeighbors<<<dimGrid, dimBlock>>>(	ignMap_d, 
																													ignMap_new_d,
    			  																							spread0Map_d,
    			  																							spreadMaxMap_d,
    			  																							azimuthMaxMap_d,
    			  																							eccentricityMap_d,
    			  																							phiEffWindMap_d,
    			  																							DistHV,
																													diff_d);//9

	
		/*cudaMemcpy(ignMap_new, ignMap_new_d, cells*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(ignMap, ignMap_d, cells*sizeof(float), cudaMemcpyDeviceToHost);
		
  	PrintMap(ignMap, "ignMap.dat");
 	 	PrintMap(ignMap_new, "ignMap_new.dat");
		
		char Continue;
		gets(&Continue);
*/
		#if WithCUBLAS
		Residue = cublasSasum(cells, diff_d, 1);
		cublasScopy(cells, ignMap_new_d, 1, ignMap_d, 1); 
		//printf("\nResidue = %f  @ itt %d\n", Residue,  n_itt);
		#endif
	}

	cudaMemcpy( ignMap,   ignMap_d,  cells*sizeof(float), cudaMemcpyDeviceToHost);
	
	end.examiningNeighbor = clock();
	times.examiningNeighbor = ((float) (end.examiningNeighbor - start.examiningNeighbor))/CLOCKS_PER_SEC;
	
	end.total = clock();
	times.total = ((float)(end.total - start.total ))/CLOCKS_PER_SEC;
	
	checkCUDAError("Spread At Neighbors Kernel");	
	/////////////////////////////////////////////////////////////////////////////
	//PRINT AND CLOSE
	#if Clatering
	printf("Done.\n");
	printf("\nPrinting And Closing.\n");
	#endif
 	
	#if WithCUBLAS
	//Close Cublas
	if(  cublasShutdown() != CUBLAS_STATUS_SUCCESS)
	{
		printf("ERROR: CUBLAS Shutdown!\n");
		return (1);
	}
	else
		printf("CUBLAS Shutdown.\n");
	#endif

	//fetch particluar values for exact functions
	cudaMemcpy(&Spread0,      &spread0Map_d[Cols/2 + Rows/2*Cols],      sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&SpreadMax,    &spreadMaxMap_d[Cols/2 + Rows/2*Cols],    sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&eccentricity, &eccentricityMap_d[Cols/2 + Rows/2*Cols], sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&azimuthMax,   &azimuthMaxMap_d[Cols/2 + Rows/2*Cols],   sizeof(float), cudaMemcpyDeviceToHost);
	
	//Maximum ignition time
	simTime = 0;
	for (cell = 0; cell < Rows*Cols; cell++)
		if (ignMap[cell] > simTime ) simTime = ignMap[cell];
 	
	//Error related stuff
	if (Scenario == 0)
	{
		exactMap = ExactCircle( ignRow, ignCol, 0, Spread0);
		errorStuff(N_isoLines, upperIsoLine, Rows, CellHt, Scenario, exactMap, ignMap, simTime);
  PrintMap(ignMap, "ignMap.dat");
	}
	else if (Scenario == 1)
	{
		exactMap = ExactElipse( ignRow, ignCol, 0, SpreadMax, eccentricity, azimuthMax);
		errorStuff(N_isoLines, upperIsoLine, Rows, CellHt, Scenario, exactMap, ignMap, simTime);
	}
	else if (Scenario == 2)
	{
		exactMap = ExactElipse( ignRow, ignCol, 0, SpreadMax, eccentricity, azimuthMax);
		errorStuff(N_isoLines, upperIsoLine, Rows, CellHt, Scenario, exactMap, ignMap, simTime);
	}

	//print ignMap for numerical solution
	strcat(fileString, "ign");
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
  PrintMap(ignMap, fileString);

	#if Clatering
	printf("\nPrinting for Ign map for Numerical solution... Done.\n");
	#endif
	//print time data for this model and this run
	strcpy(fileString, "time");
	strcat(fileString, TAG);
	strcat(fileString, "BS");
	if (BLOCK_SIZE < 9)	
	{	
		strcat(fileString, "0");
		sprintf(tmpchar, "%1d", BLOCK_SIZE);
	}
	else
		sprintf(tmpchar, "%2d", BLOCK_SIZE);
	
	strcat(fileString, tmpchar);
	
	if (Scenario == 0)	
		strcat(fileString, "_Sc0");
	if (Scenario == 1)	
		strcat(fileString, "_Sc1");
	if (Scenario == 2)	
		strcat(fileString, "_Sc2");
	if (Scenario == 3)	
		strcat(fileString, "_Sc3");	
		strcat(fileString, ".dat");
	
	timeData = fopen(fileString,"a");
	fprintf(timeData,"%5d %f %f %f %f %f %f %f %f %f %f %7d\n", Rows, 
	times.precond, times.fireSpreadNo, times.fireSpread, times.examiningNeighbor, 
	times.cublas,
	times.precond + times.fireSpreadNo + times.fireSpread + times.examiningNeighbor,
	times.fireSpreadNo + times.fireSpread + times.examiningNeighbor,
	times.total,
	(times.precond + times.fireSpreadNo + times.fireSpread + times.examiningNeighbor)/(times.total)*100 ,
	simTime, n_itt);

	#if Clatering
	printf("\nPrinting for timing data... Done.\n");
	#endif
	

	#if Clatering
	printf("%5d %f %f %f %f %f %f %f %f %f %f %7d\n", Rows, 
	times.precond, times.fireSpreadNo, times.fireSpread, times.examiningNeighbor, 
	times.cublas,
	times.precond + times.fireSpreadNo + times.fireSpread + times.examiningNeighbor,
	times.fireSpreadNo + times.fireSpread + times.examiningNeighbor,
	times.total,
	(times.precond + times.fireSpreadNo + times.fireSpread + times.examiningNeighbor)/(times.total)*100 ,
	simTime, n_itt);
	printf("\nPrinting for Error data. Done.\n");
	#endif
	/*
	//CPU GPU error
	strcpy(mapString, "e-CpuGpu_" );
	#if Scenario_0	
	strcat(mapString, "Sc0_");
	#endif
	#if Scenario_1	
	strcat(mapString, "Sc1_");
	#endif
	#if Scenario_2	
	strcat(mapString, "Sc2_");
	#endif
	sprintf(tmpchar, "%4d", Cols); 
	strcat(mapString, tmpchar);
	strcat(mapString, ".txt");
	//
	cpuMap    = getCPUign();
 	errorCPGP = fireError_b( ignMap, cpuMap);
  //
	PrintMap(errorCPGP, mapString);
  */

	fclose(IN);
	fclose(timeData);
	
	//se tirar os comentarios esta merda estoira...
	/*free(lockMap);
	free(exactMap);  
  free(fuelMap);            
  free(ignMap);             	
  free(slpMap);             	
  free(aspMap);             	
  free(wspdMap);            	
  free(wdirMap);            	
  free(m1Map);              	
	if (Scenario == 3)
	{
		fclose(slope_file);
		fclose(aspect_file);
	}
  cudaFree(modelArray_d);			 
	cudaFree(timeNext_d);					
	cudaFree(lockMap_d);
	cudaFree(lockTime_d);
  cudaFree(fuelMap_d);         
  cudaFree(ignMap_d);           	
	cudaFree(slpMap_d);       	   	
  cudaFree(aspMap_d);           	
 	cudaFree( wspdMap_d);          	
 	cudaFree( wdirMap_d);          	
  cudaFree(m1Map_d);            	
	cudaFree(spread0Map_d);				
	cudaFree(RxIntensityMap_d);		
	cudaFree(spreadMaxMap_d);
	cudaFree(phiEffWindMap_d);
	cudaFree(eccentricityMap_d);
	cudaFree(azimuthMaxMap_d);*/
	
	return (0);
}

///////////////////////////////////////////////////////////////////////////////
// Get Cpu matrix
///////////////////////////////////////////////////////////////////////////////
float *getCPUign( void)
{
	FILE* intxt;
	float *ignCPU;
	int row, col;
	char string[21] = "";
	char tmpchar[6] = "";

	ignCPU = (float*)malloc(Rows*Cols*sizeof(float));

	strcat(string, "ignLean_");
	#if Scenario_0	
	strcat(string, "Sc0_");
	#endif
	#if Scenario_1	
	strcat(string, "Sc1_");
	#endif
	#if Scenario_2	
	strcat(string, "Sc2_");
	#endif
	sprintf(tmpchar, "%4d", Cols);
	strcat(string, tmpchar);
	strcat(string, ".txt");

	intxt = fopen(string,"r");

	for (row = 0 ; row < Rows; row ++)
		for(col = 0 ; col < Cols; col++)
			fscanf(intxt,"%f", &ignCPU[col + row*Cols]);



	return(ignCPU);
}
///////////////////////////////////////////////////////////////////////////////
// Exact Solution Scenario 0
///////////////////////////////////////////////////////////////////////////////
float* ExactCircle( int ignRow, int ignCol, float ignTime, float spread0)
{
	float* ignExact;
	int row, col;

	ignExact = (float*) malloc(Rows*Cols*sizeof(float));

	for (row = 0 ; row < Rows; row++)
	{
		for (col = 0; col < Cols; col++)
		{
			ignExact[col + Cols*row] = ignTime + 
										sqrt( (col - ignCol)*(col - ignCol)*CellHt*CellHt + (row - ignRow)*(row - ignRow)*CellWd*CellWd)/
														spread0;
		}
	}
	
	return(ignExact);
}
///////////////////////////////////////////////////////////////////////////////
// Exact Solution Scenario 1
//
// Spread rate in a given direction is function of de ROS of the Wind and Slope 
// case , times a coefficient, function of the angle between the maximum spread
// rate (direction of the wind, etc) and the direction of the actual propagation
//
//R = RMax * F (azimuth, azimuthMax)
///////////////////////////////////////////////////////////////////////////////
float* ExactElipse( int ignRow, int ignCol, float ignTime, float Rmax, float Ecc, float azimuth_max)
{
	float* ignExact;
	int row, col;
	float Dist, DistH, DistW;
	float azimuth; 			//cell direction to north
	float dir; 					//angle between azimuth max and azimuth
	float F; 						//Factor to be apllied to Rmax


	ignExact = (float*)malloc(Rows*Cols*sizeof(float));

	for (row = 0 ; row < Rows; row++)
	{
		for (col = 0; col < Cols; col++)
		{
			
			DistW = (col - ignCol) * CellWd;
			DistH = (row - ignRow) * CellHt;

			//Special Cases
			if ((col - ignCol) == 0 && (row - ignRow) < 0)
				azimuth = 0.;
			else if ((col - ignCol) > 0 && (row - ignRow) == 0)
				azimuth = 90.;
			else if ((col - ignCol) == 0 && (row - ignRow) > 0)
				azimuth = 180.;
			else if ((col - ignCol) < 0 && (row - ignRow) == 0)
				azimuth = 270.;

			//1st Quadrant
			else if ( (col - ignCol) > 0 && (row - ignRow) < 0  )
			{	
				azimuth = fabs( atanf( DistW / DistH ) );
				azimuth = RadToDeg(azimuth);
			}
			//2nd Quadrant
			else if ( (col - ignCol) > 0 && (row - ignRow) > 0 )
			{	
				azimuth = atanf( DistH / DistW );
				azimuth = RadToDeg(azimuth) + 90.;
			}
			//3rd Quadrant
			else if ( (col - ignCol) < 0 && (row - ignRow) > 0 )
			{	
				azimuth = fabs(atanf( DistW / DistH ));
				azimuth = RadToDeg(azimuth) + 180.;
			}
			//4th Quadrant
			else if ( (col - ignCol) < 0 && (row - ignRow) < 0 )
			{	
				azimuth = atanf( DistH / DistW );
				azimuth = RadToDeg(azimuth) + 270.;
			}

			//printf("\nazimuth[%2d,%2d] = %5.2f\n", row, col, azimuth);

			if ((dir = fabs( azimuth_max - azimuth )) > 180. )
				dir = 360. - dir; // minimum distance between lines

			//printf("\ndir[%2d,%2d] = %5.2f\n", row, col, dir);
			dir = DegToRad(dir);
		
			F = (1- Ecc) / (1 - Ecc*cosf(dir));

			Dist = sqrt( DistH*DistH + DistW*DistW);

			ignExact[col + Cols*row] = ignTime + Dist / Rmax/F;
		}
	}
	
	return(ignExact);
}
///////////////////////////////////////////////////////////////////////////////
// Error
///////////////////////////////////////////////////////////////////////////////
float* fireError_a( float* One, float* Another)
{
	float* ignError;
	float AvgError, Sum = 0;
	int row, col;

	ignError = (float*) malloc(Cols*Rows*sizeof(float));

	for (row = 0 ; row < Rows; row++)
	{
		for (col = 0; col < Cols; col++)
		{
			ignError[col + Cols*row] = fabs( One[col + Cols*row] - Another[col + Cols*row] );
			Sum += ignError[col + Cols*row];
		}
	}
	
	AvgError = Sum / Cols/Rows;
	
	printf("  with Average Error of: %f", AvgError);

	return(ignError);
}


float* fireError_b( float* One, float* Another)
{
	float* ignError;
	float AvgError, Sum = 0;
	int row, col;

	ignError = (float*) malloc(Cols*Rows*sizeof(float));

	for (row = 0 ; row < Rows; row++)
	{
		for (col = 0; col < Cols; col++)
		{
			ignError[col + Cols*row] = fabs( One[col + Cols*row] - Another[col + Cols*row] );
			Sum += ignError[col + Cols*row];
		}
	}
	
	AvgError = Sum / Cols/Rows;
	
	printf(" / %f", AvgError);

	return(ignError);
}

///////////////////////////////////////////////////////////////////////////////
//Preprocessor Stage
//
//Computes specific Model variables outside the kernel, for all models 
///////////////////////////////////////////////////////////////////////////////
void Preprocessor(FuelCatalogPtr catalog)
{
	for (int m = 1; m <= FuelCat_MaxModels(catalog); m ++)
		Fire_FuelCombustion(catalog, (size_t) m);
}
///////////////////////////////////////////////////////////////////////////////
//Cuda Catalog function. 
//
//Copies catalog data to device. Recives the catalog handle, creates an array
//for the Models and another for the flames. Sends this arrays to the device.
//
//The Models array begins with the stride of each Model (Model + particles),
//then the model data, followed by the model's particle.
//
//Catalog structure is regarded as this				
//																	-> PartPtr[0]
//																-- 	
//																	-> PartPtr[1]
//						 	 	->ModelPtr[1]   -> PartArray[] 	--   				 
//						 	--										-> PartPtr[2]
//							 	->ModelPtr[2]   -> etc			--	
//							--										-> PartPtr[3]
//			-> ModelArray[]     ->ModelPtr[3]   -> etc
//						 	--
//			 				 	->ModelPtr[4]   -> etc
//catalog[]						 	--
//							 	->ModelPtr[etc] -> etc
//
//			
//			->flamePtr[] (size of flame classes)
//
///////////////////////////////////////////////////////////////////////////////
float* CudaModelData(FuelCatalogPtr catalog)
{
	float* modelArray_d;
	int sizeModels;				//Size of Models Array
	int sizeUPLastModel;	//offset that saves position ocupied by last Model + particle of model

	//Get Models size
	sizeModels = FuelCat_MaxModels(catalog);   //Begins with size of Models
	for (int now_model = 0; now_model < FuelCat_MaxModels(catalog); now_model++)
	{
		sizeModels += sizeModel;
	
		if (now_model != 0)
		{
			for (int now_particle = 0; now_particle < Fuel_MaxParticles(catalog,now_model); now_particle++)
			{
				//Model 0 has no particles
				sizeModels += sizePart;
			}
		}
	}

	#if Clatering
	printf("\nSize of Model Array: %d\n", sizeModels);
	#endif

	//malloc model array
	modelArray = (float*)malloc(sizeModels*sizeof(float));

	//put catalog data in Models Array
	sizeUPLastModel = FuelCat_MaxModels(catalog);
	for (int now_model = 0; now_model < FuelCat_MaxModels(catalog); now_model++)
	{
		modelArray[sizeUPLastModel + 0]   = Fuel_Mext (catalog,now_model);   			
 		modelArray[sizeUPLastModel + 1]   = Fuel_LifeAreaWtg(catalog,now_model ,0);  
 		modelArray[sizeUPLastModel + 2]   = Fuel_LifeRxFactor(catalog,now_model ,0);  
 		modelArray[sizeUPLastModel + 3]   = Fuel_BulkDensity(catalog,now_model );  
 		modelArray[sizeUPLastModel + 4]  = Fuel_PropFlux(catalog,now_model );  
 		modelArray[sizeUPLastModel + 5]  = Fuel_SlopeK(catalog,now_model );  
 		modelArray[sizeUPLastModel + 6]  = Fuel_WindB(catalog,now_model );  
 		modelArray[sizeUPLastModel + 7]  = Fuel_WindE(catalog,now_model );  
 		modelArray[sizeUPLastModel + 8]  = Fuel_WindK(catalog,now_model );  
		sizeUPLastModel += sizeModel;
			                         
		//Model 0 has no particles   
		if (now_model != 0 )
		{
			for (int now_particle = 0; now_particle < Fuel_MaxParticles(catalog,now_model); now_particle++)
			{

        modelArray[sizeUPLastModel + 0]  = Fuel_AreaWtg(catalog,now_model,now_particle); 
        modelArray[sizeUPLastModel + 1]  = Fuel_SigmaFactor(catalog,now_model,now_particle); 
				sizeUPLastModel += sizePart ;                           
			}                                                           
		}
		else
			modelArray[0] = FuelCat_MaxModels(catalog);  //Stride of 1st model in modelArray

		//Write each model (model + particle) size
		if(now_model != FuelCat_MaxModels(catalog) - 1)
			modelArray[now_model+1] = sizeUPLastModel;
	}

	//Cudamalloc and copy to device
	cudaMalloc((void**)&modelArray_d,    sizeModels*sizeof(float));
	cudaMemcpy(modelArray_d, modelArray, sizeModels*sizeof(float), cudaMemcpyHostToDevice);

	checkCUDAError("Spread At Neighbors Kernel");
	
	return (modelArray_d);

}	

/*
*******************************************************************************
Testing catalog structure. Prints info about particles and Models.

Take only catalog pointer
*******************************************************************************
 */

int Print_CatalogStruct(FuelCatalogPtr catalog)
{
	int now_model, now_particle;
	long partdescription;

	printf("\n->>There is a total of %2u Models in the catalog %s:\n\n", catalog->maxModels, catalog->name);

	printf("      Model ID | name   |               description      | mext  | MaxPartc.\n");
   	printf("      ----------------------------------------------------------------------\n");
	for (now_model = 0; now_model < catalog->maxModels; now_model++)
	{
		printf("      %2u       | %6s | %30s | %5.3f |     %1u\n",catalog->modelPtr[now_model]->modelId,
				  		       							   	      catalog->modelPtr[now_model]->name,
				  								  		 		  catalog->modelPtr[now_model]->desc,
				  								  		       	  catalog->modelPtr[now_model]->mext,
				  								  				  catalog->modelPtr[now_model]->maxParticles);
	}


	printf("\n\nThere is a total of 39 number of particles in the catalog:\n\n");

	for (now_model = 1; now_model < catalog->maxModels; now_model++)
	{
		printf(">>For model %2u: %s\n ", catalog->modelPtr[now_model]->modelId, 
										 catalog->modelPtr[now_model]->desc);

		printf("   Prt. ID |     Type      | load  |   S/V  | dens  | Silica | Si_eff  |  Area  | Sigma  \n");
		printf("    ------------------------------------------------------------------------------------  \n");
		for (now_particle = 0; now_particle < catalog->modelPtr[now_model]->maxParticles ; now_particle++)
		{
			if (catalog->modelPtr[now_model]->partPtr[now_particle]->type == 1)
			{
				partdescription = (long)"Dead particle";
			}
			else if (catalog->modelPtr[now_model]->partPtr[now_particle]->type == 2)
			{
				partdescription = (long)"Live Herb    ";
		 	}
			else
			{
				partdescription = (long)"Live Wood    ";
			}

			printf("    %2d      | %s | %5.3f | %6.1f | %3.1f  | %5.4f | %5.3f   | %6.3f | %5.3f\n",
				   now_particle,
				   partdescription,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->load,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->savr,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->dens,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->stot,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->seff,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->area,
				   catalog->modelPtr[now_model]->partPtr[now_particle]->sigma);
		
		}
		
		putchar('\n');
	}

	return(0);
}

int PrintMap ( float* map, char* fileName )
{
    FILE *fPtr;
    int cell, col, row;

    if ( (fPtr = fopen(fileName, "w")) == NULL )
    {
        printf("Unable to open output map \"%s\".\n", fileName);
        return (FIRE_STATUS_ERROR);
    }

    /*fprintf(fPtr, "north: %1.0f\nsouth: %1.0f\neast: %1.0f\nwest: %1.0f\nrows: %d\ncols: %d\n",
        (Rows*CellHt), 0., (Cols*CellWd), 0., Rows, Cols);*/
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

int fillMap(float* mapa/*, int Rows, int Cols*/)
{
	int rows, cols;
	for (rows = 0; rows < Rows; rows++)
	{
		for (cols = 0; cols < Cols; cols++)
		{
			mapa[cols + Cols*rows] = rows/Rows + cols/Cols;
		}
	}

	return(FIRE_STATUS_OK);
}


//*****************************************************************************
//                                Error check routine
//
//*****************************************************************************

void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                             cudaGetErrorString( err) );
        //exit(EXIT_FAILURE);
    }                         
}


