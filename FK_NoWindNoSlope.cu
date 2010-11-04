///////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
#include "fireCudaLib.h"

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
__global__ void FireKernel_NoWindNoSlope(			float* M1Map,
																							float* RxIntensityMap,
																							float* spread0Map,
																							size_t* fuelMap,
											  											float* modelArray)
{

	int Stride, row, col, cell;
	float moisturePart, AreaWtg, ratio, RxIntensity;

	row  = threadIdx.y + blockIdx.y*blockDim.y;
	col  = threadIdx.x + blockIdx.x*blockDim.x;
	cell = col + (blockDim.x*gridDim.x)*row; 

	Stride = modelArray[fuelMap[cell]];
	moisturePart = M1Map[cell];

	AreaWtg = modelArray[Stride + FUEL_AREAWTG];

	ratio = AreaWtg*moisturePart/(modelArray[Stride + FUEL_MEXT]);

	RxIntensity =  modelArray[Stride + FUEL_LIFERXFACTOR]												

								*(1-2.59*ratio + 5.11*ratio*ratio - 3.52*ratio*ratio*ratio); //EtaM
	RxIntensityMap[cell] = RxIntensity;

	spread0Map[cell] = modelArray[Stride + FUEL_PROPFLUX]*RxIntensity /
											((250. + 1116.*moisturePart)*AreaWtg*    //Qig - Heat of pre Ignition
											 modelArray[Stride + FUEL_LIFEAREAWTG]*
											 modelArray[Stride + FUEL_SIGMAFACTOR]*
											 modelArray[Stride + FUEL_BULKDENSITY]);

}
