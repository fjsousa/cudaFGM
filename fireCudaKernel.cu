///////////////////////////////////////////////////////////////////////////////
//fireCudaKernel.cu
//
//	Iterative cycle of the fire growth model, kernel caller and kernel.
////////////////////////////////////////////////////////////////////////////////
#include "header.h"
#include <cublas.h>
#include <stdio.h>

#define Cols (gridDim.x*blockDim.x)
#define Rows (gridDim.y*blockDim.y)

#define thx (threadIdx.x)
#define thy (threadIdx.y)

#define DistD (1.414213562*DistHV)

void checkCUDAError(const char*);
//
__global__ void FireKernel_SpreadAtNeighbors( float*, 
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float,
																							float*);
//////////////////////////////////////////////
//Kernel warper and iterative cycle of the FGM
int FGM_cycle( 	float *ignMap, 					float *ignMap_d,
								float *ignMap_new, 			float *ignMap_new_d,
								float *spread0Map, 			float *spread0Map_d,
								float *spreadMaxMap, 		float *spreadMaxMap_d,
								float *phiEffWindMap, 	float *phiEffWindMap_d,
								float *eccentricityMap, float *eccentricityMap_d,
								float *azimuthMaxMap, 	float *azimuthMaxMap_d,
								float* diff_d,
								float CellWd,
								dim3 dimGrid,
								dim3 dimBlock,
								float Residue,
								float Residue_max,
								int Cells)
								
								
{

int n_itt = 0; 							// counter for number of time steps

	cudaMemcpy( ignMap_d, ignMap, Cells*sizeof(float),cudaMemcpyHostToDevice);
	cudaMemcpy( ignMap_new_d, ignMap_new, Cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( spread0Map_d, spread0Map, Cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( spreadMaxMap_d, spreadMaxMap, Cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( phiEffWindMap_d, phiEffWindMap, Cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( eccentricityMap_d, eccentricityMap, Cells*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( azimuthMaxMap_d, azimuthMaxMap, Cells*sizeof(float), cudaMemcpyHostToDevice);
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
    			  																							CellWd,
																													diff_d);



		Residue = cublasSasum(Cells, diff_d, 1);
		cublasScopy(Cells, ignMap_new_d, 1, ignMap_d, 1); 
	}


	//////////////////
	//Cuda Error check	
	checkCUDAError(" FGM Kernel");	

	return (1);

}

/////////////
//Cuda Kernel 
__global__ void FireKernel_SpreadAtNeighbors( float* ignMap, 
											  											float* ignMap_new,
											  											float* spread0Map,
											  											float* spreadMaxMap,
											  											float* azimuthMaxMap,
											  											float* eccentricityMap,
											  											float* phiEffWindMap,
											  											float  DistHV,
																							float* diff)
{

	int nrow, ncol, ncell;
	float ignCell, ignTime_min, eccentricity;
	__shared__ float spreadAny_sh[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float ignTime_sh[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float dir_sh[BLOCK_SIZE][BLOCK_SIZE];
	
	/*int row  =threadIdx.y + blockIdx.y*blockDim.y;
	int col  =threadIdx.x + blockIdx.x*blockDim.x;
	int cell =col + Cols*row;*/
	#define row  ( threadIdx.y + blockIdx.y*blockDim.y)
	#define col  (threadIdx.x + blockIdx.x*blockDim.x)
	#define cell (col + Cols*row)
	
	//Updates ign Map from last iteration
	//ignCell = ignMap[cell];
	//ignMap_new[cell] = 101;
	//ignMap[cell] = 1;

	ignCell = ignMap[cell];
	if (ignCell > 0)
	{	
		ignTime_min = INFINITY;

		///////////////////////////////////////////////////////////////////////////
		//North Neighbor 
		#define RowN    (-1)
		#define ColN    (0)
		#define azimuth (180)

		nrow = row +  RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistHV / spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth

		///////////////////////////////////////////////////////////////////////////
		//North East Neighbor 
		#define RowN    (-1)
		#define ColN    (1)
		#define azimuth (225.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)	
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistD / spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		
		///////////////////////////////////////////////////////////////////////////
		//East Neighbor 
		#define RowN    (0)
		#define ColN    (1)
		#define azimuth (270.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
			
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistHV / spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min*      ( ignTime_sh[thx][thy] >= ignTime_min);
		}	
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South East Neighbor 
		#define RowN    (1)
		#define ColN    (1)
		#define azimuth (315)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistD/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South Neighbor 
		#define RowN    (1)
		#define ColN    (0)
		#define azimuth (0.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
			
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistHV / spreadAny_sh[thx][thy];

			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
			
		}	
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South West Neighbor 
		#define RowN    (1)
		#define ColN    (-1)
		#define azimuth (45)

		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		   
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistD / spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}	 
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//West Neighbor 
		#define RowN    (0)
		#define ColN    (-1)
		#define azimuth (90.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistHV / spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//North West Neighbor 
		#define RowN    (-1)
		#define ColN    (-1)
		#define azimuth (135)
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + DistD/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//a Neighbor 
		#define RowN    (-2)
		#define ColN    (-1)
		#define azimuth (153.434948823)//(333.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		
		
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		///////////////////////////////////////////////////////////////////////////
		//b Neighbor 
		#define RowN    (-2)
		#define ColN    (1)
		#define azimuth  (206.565051177)//  (26.56505117707799)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}

		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist

		///////////////////////////////////////////////////////////////////////////
		//c Neighbor 
		#define RowN    (-1)
		#define ColN    (-2)
		#define azimuth (116.565051177)//  (296.56505117707798)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//d Neighbor 
		#define RowN    (-1)
		#define ColN    (2)
		#define azimuth  (243.434948823)// (63.43494882292201)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//e Neighbor 
		#define RowN    (1)
		#define ColN    (-2)
		#define azimuth (63.434948823) //(243.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//f Neighbor 
		#define RowN    (1)
		#define ColN    (2)
		#define azimuth (296.565051177)// (116.56505117707799)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//g Neighbor 
		#define RowN    (2)
		#define ColN    (-1)
		#define azimuth (26.565051177)//(206.56505117707798)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;


		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//h Neighbor 
		#define RowN    (2)
		#define ColN    (1)
		#define azimuth  (333.434948823) //(153.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;

		
		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols)
		
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[ncell] < Smidgen && azimuthMaxMap[ncell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[ncell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[ncell];
				spreadAny_sh[thx][thy] = spreadMaxMap[ncell]*(1-eccentricity)/
																	(1- eccentricity*__cosf(dir_sh[thx][thy]));

				if (spreadAny_sh[thx][thy] > INFINITY)
					spreadAny_sh[thx][thy] = spread0Map[ncell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = ignMap[ncell] + Dist/ spreadAny_sh[thx][thy];

			
			ignTime_min = ignTime_sh[thx][thy]*( ignTime_sh[thx][thy] < ignTime_min)
										+ ignTime_min *      ( ignTime_sh[thx][thy] >= ignTime_min);
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist

	
	}	

	ignMap_new[cell] = ignTime_min;

	diff[cell] = ignTime_min - ignCell;

}

//////////////////////////
//Cuda Error check routine
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

