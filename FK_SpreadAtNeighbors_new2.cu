///////////////////////////////////////////////////////////////////////////////
//fireCudaKernel.cu
//
//Kernel and other functions called from the kernel
////////////////////////////////////////////////////////////////////////////////
#include "fireCudaLib.h"

#define Cols (gridDim.x*blockDim.x)
#define Rows (gridDim.y*blockDim.y)

#define DistD (sqrtf(2.)*DistHV)
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

__global__ void FireKernel_SpreadAtNeighbors( float timeNow,
											  											float* ignMap, 
											  											float* timeNext,
											  											float* spread0Map,
											  											float* spreadMaxMap,
											  											float* azimuthMaxMap,
											  											float* eccentricityMap,
											  											float* phiEffWindMap,
											  											float  DistHV,
											  											unsigned int*  lockMap,
																							unsigned int*  lockTime)//11
{

	int nrow, ncol, ncell;
	float ignCell, ignNcell, eccentricity;
	bool looping;
	__shared__ float spreadAny_sh[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float ignTime_sh[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float dir_sh[BLOCK_SIZE][BLOCK_SIZE];

	#define row  ( threadIdx.y + blockIdx.y*blockDim.y)
	#define col  (threadIdx.x + blockIdx.x*blockDim.x)
	#define cell (col + Cols*row)
	
	//Updates ign Map from last iteration
	ignCell = ignMap[cell];

	if( ignCell > timeNow  )
	{
		if (*timeNext > ignCell)
		{	
			looping = true;
			while(looping)
			{
				if (atomicExch(lockTime, 1u) == 0u)
				{
					looping = false;
					if(*timeNext > ignCell)
					{
						*timeNext = ignCell;
					}
					atomicExch(lockTime, 0u);
				}
			}
		}
	}
	else if(ignCell == timeNow )
	{
		//Edge Flag update - N esta a correr com edge flags 
		/*if (row == 0 || row == Rows -1 || col == 0 || col == Cols - 1)
		{
			looping = true;
			while(looping)
			{
				if (atomicExch(lock, 1u) == 0u)
				{
					looping = false;
					*atEdge = 1;
					atomicExch(lock, 0u);
				}
			}
		}*/

		///////////////////////////////////////////////////////////////////////////
		//North Neighbor 
		#define RowN    (-1)
		#define ColN    (0)
		#define azimuth (0)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] =   ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth

		///////////////////////////////////////////////////////////////////////////
		//North East Neighbor 
		#define RowN    (-1)
		#define ColN    (1)
		#define azimuth (45.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth
		
		///////////////////////////////////////////////////////////////////////////
		//East Neighbor 
		#define RowN    (0)
		#define ColN    (1)
		#define azimuth (90.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South East Neighbor 
		#define RowN    (1)
		#define ColN    (1)
		#define azimuth (135)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South Neighbor 
		#define RowN    (1)
		#define ColN    (0)
		#define azimuth (180.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//South West Neighbor 
		#define RowN    (1)
		#define ColN    (-1)
		#define azimuth (225)

		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//West Neighbor 
		#define RowN    (0)
		#define ColN    (-1)
		#define azimuth (270.)

		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		///////////////////////////////////////////////////////////////////////////
		//North West Neighbor 
		#define RowN    (-1)
		#define ColN    (-1)
		#define azimuth (315)
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
	
		#if Stencil16
		///////////////////////////////////////////////////////////////////////////
		//a Neighbor 
		#define RowN    (-2)
		#define ColN    (-1)
		#define azimuth (333.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		///////////////////////////////////////////////////////////////////////////
		//b Neighbor 
		#define RowN    (-2)
		#define ColN    (1)
		#define azimuth (26.56505117707799)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist

		///////////////////////////////////////////////////////////////////////////
		//c Neighbor 
		#define RowN    (-1)
		#define ColN    (-2)
		#define azimuth (296.56505117707798)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//d Neighbor 
		#define RowN    (-1)
		#define ColN    (2)
		#define azimuth (63.43494882292201)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//e Neighbor 
		#define RowN    (1)
		#define ColN    (-2)
		#define azimuth (243.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//f Neighbor 
		#define RowN    (1)
		#define ColN    (2)
		#define azimuth (116.56505117707799)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//g Neighbor 
		#define RowN    (2)
		#define ColN    (-1)
		#define azimuth (206.56505117707798)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist
		
		///////////////////////////////////////////////////////////////////////////
		//h Neighbor 
		#define RowN    (2)
		#define ColN    (1)
		#define azimuth (153.43494882292202)
		#define Dist (sqrtf((RowN*RowN + ColN*ColN)*DistHV*DistHV))
		
		nrow = row + RowN;
		ncol = col + ColN;
		ncell = ncol + nrow*Cols;
		ignNcell = ignMap[ncell];

		if(nrow >= 0 && nrow < Rows && ncol >= 0 && ncol < Cols &&
		   ignNcell > timeNow && spreadMaxMap[cell] >= Smidgen)
		{
			//"SpreadAtAzimuth"
			if (phiEffWindMap[cell] < Smidgen && azimuthMaxMap[cell] == azimuth)
				spreadAny_sh[thx][thy] = spreadMaxMap[cell];
			else
			{
				if ((dir_sh[thx][thy] = fabsf(azimuthMaxMap[cell] - azimuth)) > 180)
					dir_sh[thx][thy] = 360. - dir_sh[thx][thy];
			
				dir_sh[thx][thy] = DegToRad(dir_sh[thx][thy]);

				eccentricity = eccentricityMap[cell];
				spreadAny_sh[thx][thy] = __fdividef(spreadMaxMap[cell]*(1-eccentricity),(1- eccentricity*__cosf(dir_sh[thx][thy])));

				if (spreadAny_sh[thx][thy] > INF)
					spreadAny_sh[thx][thy] = spread0Map[cell];
			}//"SpreadAtAzimuth"

			ignTime_sh[thx][thy] = timeNow + __fdividef(DistHV,spreadAny_sh[thx][thy]);

			if(ignTime_sh[thx][thy] < ignNcell)
			{
				looping = true;
				while(looping)
				{
					if (atomicExch(&(lockMap[ncell]), 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < ignMap[ncell])
						{
							ignMap[ncell] = ignTime_sh[thx][thy];
						}	
						atomicExch(&(lockMap[ncell]), 0u);
					}
				}
			}
			//Update timeNext
			if( ignTime_sh[thx][thy] < *timeNext )
			{	
				looping = true;
				while(looping)
				{
					if (atomicExch(lockTime, 1u) == 0u)
					{
						looping = false;
						if(ignTime_sh[thx][thy] < *timeNext)
						{
							*timeNext = ignTime_sh[thx][thy];
						}
						atomicExch(lockTime, 0u);
					}
				}
			}	
		}
		#undef RowN
		#undef ColN
		#undef azimuth
		#undef Dist

		#endif //For Stencil 16
	}	

}
