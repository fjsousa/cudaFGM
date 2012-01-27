///////////////////////////////////////////////////////////////////////////////
//fireCudaKernel.cu
//
//Kernel and other functions called from the kernel
////////////////////////////////////////////////////////////////////////////////
#include "header.h"

#define Cols (gridDim.x*blockDim.x)
#define Rows (gridDim.y*blockDim.y)

#define thx (threadIdx.x)
#define thy (threadIdx.y)

#define DistD (1.414213562*DistHV)
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

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
