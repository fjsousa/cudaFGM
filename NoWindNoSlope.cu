///////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
#include "fireCudaLib.h"

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
float NoWindNoSlope(int cell, float M1Cell, size_t fuelCell, float* modelArray, float* RxIntensityMap )
{

	int Stride;
	float moisturePart, AreaWtg, ratio, RxIntensity, Spread0Cell;

	Stride = modelArray[fuelCell];
	moisturePart = M1Cell;

	AreaWtg = modelArray[Stride + FUEL_AREAWTG];

	ratio = AreaWtg*moisturePart/(modelArray[Stride + FUEL_MEXT]);

	RxIntensity =  modelArray[Stride + FUEL_LIFERXFACTOR]												

								*(1-2.59*ratio + 5.11*ratio*ratio - 3.52*ratio*ratio*ratio); //EtaM
	
	RxIntensityMap[cell] = RxIntensity;

	Spread0Cell = modelArray[Stride + FUEL_PROPFLUX]*RxIntensity /
											((250. + 1116.*moisturePart)*AreaWtg*    //Qig - Heat of pre Ignition
											 modelArray[Stride + FUEL_LIFEAREAWTG]*
											 modelArray[Stride + FUEL_SIGMAFACTOR]*
											 modelArray[Stride + FUEL_BULKDENSITY]);

	return (Spread0Cell);
}
