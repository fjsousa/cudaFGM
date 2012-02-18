////////////////////////////////////////////////////////////////////////
// *cudaFGM*
//    A fire growth model (FGM) that runs on NVIDIA GPUs. 
//
//firelib is used to create most of the fire properties. Main firelib 
//funcions used are "Fire_FuelCatalogCreateStandard" to create the fire 
//catalog, "Fire_SpreadNoWindNoSlope" and "Fire_SpreadWindSlopeMax" are 
//used to compute fire ellipse properties. 
//
//The kernel "FGM" is launched in each iteration of the fire growth model. 
//
//Main output is the ignition map "ignMap_###.dat", where ### is a user 
//defined tag
//
//Grass aspect and slope files can be read with "./cudaFGM 1". Grass slope files must
//be written in "percentage".
//
//Inputs provided in "RunSet.in": Map width (meters)
//																Map height
//																Number of rows
//																Number of cols	
//																Fuel Model (NFFL models) or 14-custom	
//																Wind speed
//																Wind Direction
//																Moisture (M1, M10, M100, Mherb, Mwood)
//																Custom Particle load (Not relevant if using 
//																	one of theNFFL 13) 
//																ignition point (X) as a %(0-1) of map width
//																ignition point (Y) as a %(0-1) of map height
//																GPU device
//                                ignition map file name (no spaces)
//																Verbosity (1 - more 0 - less)
//                                init ignMap to BEHAVE elipse for faster solution (1 - Yes, 0 - No)
//
//Change log: 17/02/2012          cudaFGM reads Grass aspect and slope file formats.  
//

#include "fireLib_float.h"
#include "header.h"
#include <cublas.h>
#include <time.h>


//////////////////
//Function headers
//
int PrintMap ( float*, char* );
//
float* BEHAVEelipse( float, int, int, float, float, float);
//
int Print_CatalogStruct(FuelCatalogPtr);
//
int FGM_cycle( 	float*,	float*,
								float*,	float*,
								float*,	float*,
								float*,	float*,
								float*,	float*,
								float*, float*,
								float*,	float*,
								float*,
								float,
								dim3,
								dim3,
								float,
								float,
								int);

//////////////////
//Global variables
int    Rows;                    //Map dimensions are global variables
int    Cols;										//
float mapW_m, mapH_m;						//map width height (meters) 
float mapW, mapH;    						//map width height (feet)

size_t Model;   
float WindSpd; 
float WindDir;  
float M1;      
float M10;     
float M100;    
float Mherb;   
float Mwood;   

//////
//Main
int main ( int argc, char *argv[] )
{
	int grass = atoi(argv[1]);  //Grass files are used if grass==1
	float slp_tmp, asp_tmp;     //slope and aspect temporary values
	unsigned int device;        //Assigned GPU device 
	//
	float Residue = INFINITY;		//Residue value
	float Residue_max = Smidgen;//Maximum residue			
  //
  int    row, col, cell;      /* row, col, and index of current cell */
  int    Cells;            		/* total number of map cells */
	float  ignX, ignY;					// ignition points
	float CellWd;
	//float CellHt;							//!Not used yet!
	//
	FuelCatalogPtr catalog;     /* fuel catalog handle */
 	float moisture[6];         /* fuel moisture content at current cell */
  //
	float particle_load;				//*CUSTOM FUEL MODEL* - particle load
	//
	float *initialMap;         //BEHAVE eliptical ignition map 
	size_t *fuelMap;           /* ptr to fuel model map */
  float *ignMap;             /* ptr to ignition time map (minutes) */
  float  *ignMap_new;        /* ptr to ignition time map (minutes) */
  float *slpMap;             /* ptr to slope map (rise/reach) */
  float *aspMap;             /* ptr to aspect map (degrees from north) */
  float *wspdMap;            /* ptr to wind speed map (ft/min) */
  float *wdirMap;            /* ptr to wind direction map (deg from north) */
  float *m1Map;              /* ptr to 1-hr dead fuel moisture map */
  float *m10Map;             /* ptr to 10-hr dead fuel moisture map */
  float *m100Map;            /* ptr to 100-hr dead fuel moisture map */
  float *mherbMap;           /* ptr to live herbaceous fuel moisture map */
  float *mwoodMap;           /* ptr to live stem fuel moisture map */
	float *spread0Map;				
	float *spreadMaxMap;
	float *azimuthMaxMap;
 	float *eccentricityMap;
	float *phiEffWindMap;
	//
	float *ignMap_d;
	float *ignMap_new_d;
	float *spread0Map_d;				
	float *spreadMaxMap_d;
	float *phiEffWindMap_d;
	float *eccentricityMap_d;
	float *azimuthMaxMap_d;
	float	*diff_d;
	//
	FILE *IN, *slope_file, *aspect_file;
	char buffer[100];     			//buffer to use when fgets skips lines
	char ignFileName[40];
	int n;
	int verbosity;              //level of verbosity of shell output
	int init_elipse;						//init ignMap to BEHAVE elipse for faster solution
	//
	clock_t start, end;
	double time;
	
	////////////////
	//Read RunSet.in
	IN = fopen("RunSet.in", "r");
	//skips 22 text lines of RunSet.in
	for (n = 0; n < 24; n++)
		fgets(buffer, 100, IN);
	fscanf(IN, "%f %f %d %d %d %f %f %f %f %f %f %f %f %f %f %d %s %d %d",  
	  	&mapW_m, &mapH_m, &Rows, &Cols, &Model, &WindSpd, &WindDir, 
			&M1, &M10, &M100, &Mherb, &Mwood, &particle_load, &ignX, &ignY, &device, &ignFileName, &verbosity, &init_elipse);
	//Input Checks
	if (  ignX > 1 || ignX < 0 || ignY > 1 || ignY < 0 )
	{
		printf("\nERROR: Runset.in - ignition point must be 0 <= ign(X,Y) <=1!\n");
		return (0);
	}
	if (  Cols != Rows )
	{
		printf("\nERROR: Runset.in - Rows must be equal to Cols!\n");
		return (0);
	}
	if (  mapW_m != mapH_m )
	{
		printf("\nERROR: Runset.in - Width must be equal to height!\n");
		return (0);
	}
	if (  Cols%BLOCK_SIZE != 0 ||  Rows%BLOCK_SIZE != 0)
	{
		printf("\nERROR: Cols and Rows must be multiples of BLOCK_SIZE! (cuda related restriction)\n");
		return (0);
	}

	////////////////
	//set GPU device
	cudaSetDevice(device);
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(Rows/BLOCK_SIZE, Cols/BLOCK_SIZE);

	//////////////
	//clublas init
	printf("\n>>Initializing CUBLAS...");
	if( cublasInit() != CUBLAS_STATUS_SUCCESS)
	{
		printf("\nERROR: CUBLAS initialization!\n");
		return (1);
	}
	else
		printf("Done.\n");
	
	///////////////////////
  //Allocate all the maps
  Cells = Rows * Cols;
  if ( (ignMap   					= (float *) calloc(Cells, sizeof(float))) == NULL
		|| (ignMap_new   			= (float *) calloc(Cells, sizeof(float)))   == NULL
		|| (azimuthMaxMap   	= (float *) calloc(Cells, sizeof(float)))   == NULL
		|| (spread0Map   			= (float *) calloc(Cells, sizeof(float)))   == NULL
		|| (spreadMaxMap     	= (float *) calloc(Cells, sizeof(float)))   == NULL
		|| (eccentricityMap   = (float *) calloc(Cells, sizeof(float)))   == NULL
		|| (phiEffWindMap   	= (float *) calloc(Cells, sizeof(float)))   == NULL
    || (slpMap   = (float *) calloc(Cells, sizeof(float))) == NULL
    || (aspMap   = (float *) calloc(Cells, sizeof(float))) == NULL
    || (wspdMap  = (float *) calloc(Cells, sizeof(float))) == NULL
    || (wdirMap  = (float *) calloc(Cells, sizeof(float))) == NULL
    || (m1Map    = (float *) calloc(Cells, sizeof(float))) == NULL
    || (m10Map   = (float *) calloc(Cells, sizeof(float))) == NULL
    || (m100Map  = (float *) calloc(Cells, sizeof(float))) == NULL
    || (mherbMap = (float *) calloc(Cells, sizeof(float))) == NULL
    || (mwoodMap = (float *) calloc(Cells, sizeof(float))) == NULL
    || (fuelMap  = (size_t *) calloc(Cells, sizeof(size_t))) == NULL )
  {
      fprintf(stderr, "Unable to allocate maps with %d cols and %d rows.\n",
          Cols, Rows);
      return (1);
  }
	//Cuda maps
	cudaMalloc((void**)&spread0Map_d, 	  Cells*sizeof(float));
	cudaMalloc((void**)&spreadMaxMap_d,    Cells*sizeof(float));
	cudaMalloc((void**)&phiEffWindMap_d,   Cells*sizeof(float));
	cudaMalloc((void**)&eccentricityMap_d, Cells*sizeof(float));
	cudaMalloc((void**)&azimuthMaxMap_d,   Cells*sizeof(float));
	cudaMalloc((void**)&ignMap_d,      Cells*sizeof(float));
	cudaMalloc((void**)&ignMap_new_d,  Cells*sizeof(float));
	cudaMalloc((void**)&diff_d,  			 Cells*sizeof(float));
	
	////////////
  //Map set up
	mapW = MetersToFeet(mapW_m);
	mapH = MetersToFeet(mapH_m);
	CellWd = mapW/(Cols -1);				  				  	
	//CellHt = mapH/(Rows -1);					    				
	//slope and aspect file read
	slope_file = fopen("slope.map","r");
	aspect_file = fopen("aspect.map","r");
	//If using grass file format for aspect and slope
	if (grass == 1)  
	{
		//fgets to skip header in aspect and slope files
		for (n = 0; n < 6; n++)
		{
			fgets(buffer, 100, slope_file);
			fgets(buffer, 100, aspect_file);
		}
		//data assignment to maps
		for ( row = 0; row < Rows; row++ )
		{
 		 	for ( col = 0; col < Cols; col ++)
			{
				cell = col + row*Cols;
      
				fscanf(aspect_file, "%f", &asp_tmp);
				fscanf(slope_file, "%f", &slp_tmp);
				slpMap[cell]     = slp_tmp/100; 										//Slope in firelib is a fraction
				asp_tmp = (asp_tmp - 90 < 0) ?                      //while in Grass is percentage rise/reach.
											asp_tmp - 90 + 360	: asp_tmp - 90 ;  //Aspect in firelib is N=0 and clockwise 
				aspMap[cell]	   = 360 - asp_tmp; 						      //while aspect in Grass is E=0 counter-clockwise
				fuelMap[cell]    = Model;
   		  wspdMap[cell]    = 88. * WindSpd;     							/* convert mph into ft/min */
   	 	  wdirMap[cell]    = WindDir;
   	 	  m1Map[cell]      = M1;
   	 	  m10Map[cell]     = M10;
     	  m100Map[cell]    = M100;
   		  mherbMap[cell]   = Mherb;
     		mwoodMap[cell]   = Mwood;
				ignMap[cell] 		 = 500;
  			ignMap_new[cell] = 500;
  		}
		}
		PrintMap(aspMap,"aspectTest.map");
	}
	else 
	{
		//data assignment to maps
		for ( row = 0; row < Rows; row++ )
		{
 		 	for ( col = 0; col < Cols; col ++)
			{
				cell = col + row*Cols;
      
				fscanf(aspect_file, "%f", &aspMap[cell] );
				fscanf(slope_file, "%f", &slpMap[cell] );
				fuelMap[cell]    = Model;
   		  wspdMap[cell]    = 88. * WindSpd;     /* convert mph into ft/min */
   	 	  wdirMap[cell]    = WindDir;
   	 	  m1Map[cell]      = M1;
   	 	  m10Map[cell]     = M10;
     	  m100Map[cell]    = M100;
   		  mherbMap[cell]   = Mherb;
     		mwoodMap[cell]   = Mwood;
				ignMap[cell] 		 = 500;
  			ignMap_new[cell] = 500;
  		}
		}
	}

	//ignition point - ignX and ignY is a percentage of the map height and width
	cell = Cols*ignX + Cols*Rows*ignY;
	ignMap[cell] 		 = 0;
  ignMap_new[cell] = 0;
	

	////////////////////////////////
  //Create fuel catalog
  
	//Create 13 + 0 (no fuel model) standard NFFL models and creates space for 
	//aditional custom model
	printf ("\n>>Creating standard fire models...");
	catalog = Fire_FuelCatalogCreateStandard("Standard", 14);
	printf ("Done.\n");
	//Create aditional custom model based on NFFL1
	//Only the PARTICLE LOAD is customized at the moment
	if ( Fire_FuelModelCreate (
  	catalog,     														//FuelCatalogData instance
    14,              												//fuel model number
    "CUSTOM",																//Name
    "Custom Fuel model", //longer description
    Fuel_Depth(catalog, 1),              		//bed depth (ft)
    Fuel_Mext(catalog, 1),  	            	//moisture of extinction (dl)
    Fuel_SpreadAdjustment(catalog, 1),      //spread adjustment factor (dl)
    1) != FIRE_STATUS_OK )									//maximum number of particles
	{
  	fprintf(stderr, "%s\n", FuelCat_Error(catalog));
  	Fire_FuelCatalogDestroy(catalog);
  	return (NULL);
	}
	//Add a particle to the custom model nº 14
	printf ("\n>>Creating custom fire model...");
  start = clock();
	if ( Fire_FuelParticleAdd (
  	catalog,     									// FuelCatalogData instance pointer
    14,              							//Custom fuel model id
    Fuel_Type(catalog,1,0),   
    particle_load,            		// Custom particle load              (lbs/ft2)
    Fuel_Savr(catalog,1,0),   		// surface-area-to-volume ratio     (ft2/ft3)
    Fuel_Density(catalog,1,0), 		//density                          (lbs/ft3)
    Fuel_Heat(catalog,1,0),  			//heat of combustion               (btus/lb)
    Fuel_SiTotal(catalog,1,0),    //total silica content               (lb/lb)
    Fuel_SiEffective(catalog,1,0))//effective silica content           (lb/lb)
				!= FIRE_STATUS_OK )
  {
    fprintf(stderr, "%s\n", FuelCat_Error(catalog));
    Fire_FuelCatalogDestroy(catalog);
    return (NULL);
  }
	else
	{
		end = clock();
		time = ((double) (end - start))/CLOCKS_PER_SEC;
		printf("Done with %lf seconds.\n",time);
	}
	
	/////////////////////////
	//Print catalog structure
	if (verbosity == 1) Print_CatalogStruct(catalog);

	//////////////////////////////////////////////////////////////
	//Preprocessing stage:Create sprea0 and spreadMax maps
	//Initialize ignMap and ignMap_new to BEHAVE eliptical ign map
	printf("\n>>Running preprocessor...");
	start = clock();
	for (cell = 0; cell < Cells; cell++)
	{
  	Model = fuelMap[cell];
  	moisture[0] = m1Map[cell];
  	moisture[1] = m10Map[cell];
  	moisture[2] = m100Map[cell];
  	moisture[3] = m100Map[cell];
  	moisture[4] = mherbMap[cell];
  	moisture[5] = mwoodMap[cell];
  	Fire_SpreadNoWindNoSlope(catalog, Model, moisture);
  	Fire_SpreadWindSlopeMax(catalog, Model, wspdMap[cell],
     wdirMap[cell], slpMap[cell], aspMap[cell]);
		
		spread0Map[cell]      = Fuel_Spread0(catalog,Model); 				
		spreadMaxMap[cell]    = Fuel_SpreadMax(catalog,Model);
		azimuthMaxMap[cell]   = Fuel_AzimuthMax(catalog,Model); 
		eccentricityMap[cell] = Fuel_Eccentricity(catalog,Model);
		phiEffWindMap[cell]  	= Fuel_PhiEffWind(catalog,Model);
	}
	//Initialize BEHAVE Elipse and update ignMap and ignMap_new 
	if (init_elipse == 1)
	{
		cell = Cols*ignX + Cols*Rows*ignY; //elipse is created with ignition point values 
		initialMap = BEHAVEelipse( CellWd, Rows*ignY, Cols*ignX, spreadMaxMap[cell], eccentricityMap[cell], azimuthMaxMap[cell]);
		for (cell = 0; cell < Cells; cell++)
			ignMap[cell] = ignMap_new[cell] = initialMap[cell];
	}
	end = clock();
	time = ((double) (end - start))/CLOCKS_PER_SEC;
	printf("Done with %lf seconds.\n",time);

	///////////////////////////////////
	//Fire growth model iterative cycle
	printf("\n>>Running FGM cycle...");
	start = clock();
	if ( FGM_cycle( ignMap, 				 ignMap_d,
								  ignMap_new, 		 ignMap_new_d,
								  spread0Map, 		 spread0Map_d,
								  spreadMaxMap, 	 spreadMaxMap_d,
								  phiEffWindMap, 	 phiEffWindMap_d,
								  eccentricityMap, eccentricityMap_d,
								  azimuthMaxMap, 	 azimuthMaxMap_d,
								 	diff_d,
								  CellWd,
								 	dimGrid,
								 	dimBlock,
								  Residue,
								  Residue_max,
								  Cells) !=1)
	{
		printf("\nERROR: FGM cycle.\n");
	}
	else
	{
		end = clock();
		time = ((double) (end - start))/CLOCKS_PER_SEC;
		printf("Done with %lf seconds.\n",time);
	}

	//////////////////////////
	//Mem copy of ignition map
	cudaMemcpy( ignMap,   ignMap_d,  Cells*sizeof(float), cudaMemcpyDeviceToHost);
	
	////////////////////
	//Print Ignition Map
	PrintMap(ignMap, ignFileName);

	//////////////
	//Close Cublas
	printf ("\n>>Closing CUBLAS...");
	if(  cublasShutdown() != CUBLAS_STATUS_SUCCESS)
	{
		printf("\nERROR: CUBLAS Shutdown!\n");
		return (1);
	}
	else
		printf("Done.\n");

}
/////////
//The END








////////////////
//More Functions
////////////////


///////////////////////////////////////////////////////////////////////////////
// Exact Solution Scenario 1
//
// Spread rate in a given direction is function of de ROS of the Wind and Slope 
// case , times a coefficient, function of the angle between the maximum spread
// rate (direction of the wind, etc) and the direction of the actual propagation
//
//R = RMax * F (azimuth, azimuthMax)
///////////////////////////////////////////////////////////////////////////////
float* BEHAVEelipse( float CellWd, int ignRow, int ignCol, float Rmax, 
																									float Ecc, float azimuth_max)
{
	float* elipseMap;
	int row, col;
	float Dist, DistH, DistW;
	float azimuth; 							//cell direction to north
	float dir; 									//angle between azimuth max and azimuth
	float F; 										//Factor to be apllied to Rmax
	float ignTime = 0;
	float CellHt  = CellWd;     //CellHt is equal to cell width
															//Code dos not handle yet rectangular domains

	elipseMap = (float*)malloc(Rows*Cols*sizeof(float));

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

			if ((dir = fabs( azimuth_max - azimuth )) > 180. )
				dir = 360. - dir; // minimum distance between lines

			dir = DegToRad(dir);
		
			F = (1- Ecc) / (1 - Ecc*cosf(dir));

			Dist = sqrt( DistH*DistH + DistW*DistW);

			elipseMap[col + Cols*row] = ignTime + Dist / Rmax/F;
		}
	}
	
	return(elipseMap);
}


////////////
//Print Maps
int PrintMap ( float* map, char* fileName )
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
            fprintf(fPtr, "  %5.2f ", (map[cell]==INFINITY) ? 000.00 : map[cell]);
        }
        fprintf(fPtr, "\n");
    }
    fclose(fPtr);
    return (FIRE_STATUS_OK);
}

/////////////////////////
//Print catalog structure
int Print_CatalogStruct(FuelCatalogPtr catalog)
{
	int now_model, now_particle;
	long partdescription;

	printf("\n>>There is a total of %2u Models in the catalog %s:\n\n", catalog->maxModels, catalog->name);

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
