///////////////////////////////////////////////////////////////////////////////
//fireCudaKernel.cu
//
//Kernel and other functions called from the kernel
////////////////////////////////////////////////////////////////////////////////
#include "fireCudaLib.h"
								
#define Cols (gridDim.x*blockDim.x)
#define Rows (gridDim.y*blockDim.y)
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////

__global__ void FireKernel_WindAndSlope(	float* windUMap, 
											  									float* windDirMap,
											  									float* slopeMap,
											  									float* aspectMap,
											  									float* spread0Map,
											  									float* spreadMaxMap,
											  									float* RxIntensityMap,
											  									float* phiEffWindMap,
																					float* eccentricityMap,
																					float* azimuthMaxMap,
																					size_t* fuelMap,
																					float* modelArray) //14 
{

	int row, col, cell, Stride;
	float windB, windK,  phiSlope, phiWind, phiEw, upSlope, spreadMax;
	__shared__ float spread0_sh      		[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float slope_sh       		[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float effectiveWind_sh		[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float maxWind_sh			 		[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float lwRatio_sh					[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float split_sh						[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float x_sh								[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float y_sh								[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float Rv_sh							[BLOCK_SIZE][BLOCK_SIZE];
	__shared__ float a_sh								[BLOCK_SIZE][BLOCK_SIZE];

	row  = threadIdx.y + blockIdx.y*blockDim.y;
	col  = threadIdx.x + blockIdx.x*blockDim.x;
	cell = col + Rows*row; 
	
	Stride = modelArray[fuelMap[cell]];

	slope_sh[thx][thy] = slopeMap[cell];

	windB = modelArray[Stride + FUEL_WINDB];
	windK = modelArray[Stride + FUEL_WINDK];
	
	phiSlope = modelArray[Stride + FUEL_SLOPEK]*slope_sh[thx][thy]*slope_sh[thx][thy];
	phiWind  = modelArray[Stride + FUEL_WINDK] *powf(windUMap[cell],windB);
	
	//PhiWind tem um teste < smidgen em relacao a velocidade do vento WindUMap... 
	phiEw = phiSlope + phiWind;

	if((upSlope = aspectMap[cell]) >= 180.)
		upSlope = upSlope - 180;
	else
		upSlope = upSlope + 180;

	spread0_sh[thx][thy] = spread0Map[cell];

	//Situation 1 No fire Spread or reaction Intensity
	if(spread0Map[cell] < Smidgen)
	{	
		spreadMaxMap[cell]    = 0;
		eccentricityMap[cell] = 0;
		azimuthMaxMap[cell]   = 0;
		phiEffWindMap[cell]   = phiEw;
	}

	//Situation 2 No Wind and No Slope
	else if (phiEw < Smidgen)
	{
		phiEffWindMap[cell]   = 0;
		spreadMaxMap[cell]    = spread0_sh[thx][thy];
		eccentricityMap[cell] = 0;
		azimuthMaxMap[cell]   = 0;
	}

	//Situation 3 Wind with No Slope
	else if (slope_sh[thx][thy] < Smidgen)
	{
		effectiveWind_sh[thx][thy] = windUMap[cell];
		azimuthMaxMap[cell] = windDirMap[cell];	
		
		maxWind_sh[threadIdx.x][threadIdx.y] = 0.9*RxIntensityMap[cell];
		if(effectiveWind_sh[thx][thy] >  maxWind_sh[thx][thy])
		{
			phiEw = windK*__powf(maxWind_sh[thx][thy], windB);
			effectiveWind_sh[thx][thy] = maxWind_sh[thx][thy];
		}
		spreadMaxMap[cell] = spread0_sh[thx][thy]*(1 + phiEw);
		
		if(effectiveWind_sh[thx][thy] >  Smidgen)
		{
			lwRatio_sh[thx][thy] = 1. + 0.002840909 * effectiveWind_sh[thx][thy];
			if (lwRatio_sh[thx][thy] > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio_sh[thx][thy]*lwRatio_sh[thx][thy] - 1)/lwRatio_sh[thx][thy];
		}

		phiEffWindMap[cell] = phiEw;
	}

	//Situation 4 and 5 - slope with no wind and wind blows upSlope
	else if(windUMap[cell] < Smidgen || _Equal(upSlope, windDirMap[cell]))
	{
		azimuthMaxMap[cell] = upSlope;
		effectiveWind_sh[thx][thy] = powf(phiEw*modelArray[Stride + FUEL_WINDE], 1/windB);
		
		maxWind_sh[thx][thy] = 0.9*RxIntensityMap[cell];
		if(effectiveWind_sh[thx][thy] >  maxWind_sh[thx][thy])
		{
			phiEw = windK*__powf(maxWind_sh[thx][thy], windB);
			effectiveWind_sh[thx][thy] = maxWind_sh[thx][thy];
		}

		if(effectiveWind_sh[thx][thy] >  Smidgen)
		{
			lwRatio_sh[thx][thy] = 1. + 0.002840909 * effectiveWind_sh[thx][thy];
			if (lwRatio_sh[thx][thy] > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio_sh[thx][thy]*lwRatio_sh[thx][thy] - 1)/lwRatio_sh[thx][thy];
		}

		spreadMaxMap[cell] = spread0_sh[thx][thy]*(1 + phiEw);
		phiEffWindMap[cell] = phiEw;
	}
	//Situation 6 - Wind Blows cross Slope
	else
	{
		split_sh[thx][thy] = windDirMap[cell]; //split_sh here works as a buffer to winDirMap
		if (upSlope <= split_sh[thx][thy])
			split_sh[thx][thy] = split_sh[thx][thy] - upSlope;
		else
			split_sh[thx][thy] = 360. - upSlope + split_sh[thx][thy];

		split_sh[thx][thy] = DegToRad(split_sh[thx][thy]);
		x_sh[thx][thy]  = spread0_sh[thx][thy]*(phiSlope + phiWind*__cosf(split_sh[thx][thy]));
		y_sh[thx][thy]  = spread0_sh[thx][thy]*(phiWind*__sinf(split_sh[thx][thy]));
		Rv_sh[thx][thy] = sqrtf(x_sh[thx][thy]*x_sh[thx][thy] + y_sh[thx][thy]*y_sh[thx][thy]);

		spreadMax = spread0_sh[thx][thy] + Rv_sh[thx][thy];
		phiEw = spreadMax / spread0_sh[thx][thy] - 1;
		a_sh[thx][thy] = asinf(fabsf(y_sh[thx][thy]) / Rv_sh[thx][thy]);
		if(x_sh[thx][thy] >= 0.)
			a_sh[thx][thy] = (y_sh[thx][thy] >= 0.) ? a_sh[thx][thy]          : M_PI + M_PI - a_sh[thx][thy];
		else
			a_sh[thx][thy] = (y_sh[thx][thy] >= 0.) ? (M_PI - a_sh[thx][thy]) : M_PI + a_sh[thx][thy];
		
		split_sh[thx][thy] = RadToDeg(a_sh[thx][thy]);
		if (upSlope + split_sh[thx][thy] > 306.)
			azimuthMaxMap[cell] = upSlope + split_sh[thx][thy] - 360.;
		else
			azimuthMaxMap[cell] = upSlope + split_sh[thx][thy];

		effectiveWind_sh[thx][thy] = powf(phiEw*modelArray[Stride + FUEL_WINDE], 1/windB);
		//Do effective wind only if phiEw > Smidgen
		if(phiEw > Smidgen)
		{		
			maxWind_sh[thx][thy] = 0.9*RxIntensityMap[cell];
			if(effectiveWind_sh[thx][thy] >  maxWind_sh[thx][thy])
			{
				phiEw = windK*__powf(maxWind_sh[thx][thy], windB);
				effectiveWind_sh[thx][thy] = maxWind_sh[thx][thy];
				spreadMax = spread0_sh[thx][thy]*(1 + phiEw);
			}
		}

		if(effectiveWind_sh[thx][thy] >  Smidgen)
		{
			lwRatio_sh[thx][thy] = 1. + 0.002840909 * effectiveWind_sh[thx][thy];
			if (lwRatio_sh[thx][thy] > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio_sh[thx][thy]*lwRatio_sh[thx][thy] - 1)/lwRatio_sh[thx][thy];
		}

		spreadMaxMap[cell] = spreadMax;
		phiEffWindMap[cell] = phiEw;
	}
}


