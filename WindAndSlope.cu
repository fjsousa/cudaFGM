///////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
#include "fireCudaLib.h"

float WindAndSlope(	int cell,
										float* modelArray,
										float spread0Cell,
										float windUCell, 
										float windDirCell,
										float slopeCell,
										float aspectCell,
										size_t fuelCell,
										float* RxIntensityMap,
										float* phiEffWindMap,
										float* eccentricityMap,
										float* azimuthMaxMap) 
{

	int Stride;
	float windB, windK,  phiSlope, phiWind, phiEw, upSlope, spreadMax, spreadMaxCell;
	float slope;       	
	float effectiveWind;
	float maxWind;			
	float lwRatio;			
	float split;				
	float x;					
	float y;					
	float Rv;				
	float a;					

	Stride = modelArray[fuelCell];

	slope  = slopeCell;

	windB = modelArray[Stride + FUEL_WINDB];
	windK = modelArray[Stride + FUEL_WINDK];
	
	phiSlope = modelArray[Stride + FUEL_SLOPEK]*slope *slope ;
	phiWind  = modelArray[Stride + FUEL_WINDK] *powf(windUCell,windB);
	
	//PhiWind tem um teste < smidgen em relacao a velocidade do vento WindUMap... 
	phiEw = phiSlope + phiWind;

	if((upSlope = aspectCell) >= 180.)
		upSlope = upSlope - 180;
	else
		upSlope = upSlope + 180;


	//Situation 1 No fire Spread or reaction Intensity
	if(spread0Cell < Smidgen)
	{	
		spreadMaxCell    = 0;
		eccentricityMap[cell] = 0;
		azimuthMaxMap[cell]   = 0;
		phiEffWindMap[cell]   = phiEw;
	}

	//Situation 2 No Wind and No Slope
	else if (phiEw < Smidgen)
	{
		phiEffWindMap[cell]   = 0;
		spreadMaxCell    = spread0Cell ;
		eccentricityMap[cell] = 0;
		azimuthMaxMap[cell]   = 0;
	}

	//Situation 3 Wind with No Slope
	else if (slope  < Smidgen)
	{
		effectiveWind  = windUCell;
		azimuthMaxMap[cell] = windDirCell;	
		
		maxWind = 0.9*RxIntensityMap[cell];
		if(effectiveWind  >  maxWind )
		{
			phiEw = windK*powf(maxWind , windB);
			effectiveWind  = maxWind ;
		}
		spreadMaxCell = spread0Cell *(1 + phiEw);
		
		if(effectiveWind  >  Smidgen)
		{
			lwRatio  = 1. + 0.002840909 * effectiveWind ;
			if (lwRatio  > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio *lwRatio  - 1)/lwRatio ;
		}

		phiEffWindMap[cell]  = phiEw;
	}

	//Situation 4 and 5 - slope with no wind and wind blows upSlope
	else if(windUCell < Smidgen || _Equal(upSlope, windDirCell))
	{
		azimuthMaxMap[cell] = upSlope;
		effectiveWind  = powf(phiEw*modelArray[Stride + FUEL_WINDE], 1/windB);
		
		maxWind  = 0.9*RxIntensityMap[cell];
		if(effectiveWind  >  maxWind )
		{
			phiEw = windK*powf(maxWind , windB);
			effectiveWind  = maxWind ;
		}

		if(effectiveWind  >  Smidgen)
		{
			lwRatio  = 1. + 0.002840909 * effectiveWind ;
			if (lwRatio  > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio *lwRatio  - 1)/lwRatio ;
		}

		spreadMaxCell = spread0Cell *(1 + phiEw);
		phiEffWindMap[cell] = phiEw;
	}
	//Situation 6 - Wind Blows cross Slope
	else
	{
		split  = windDirCell; //split_sh here works as a buffer to winDirMap
		if (upSlope <= split )
			split  = split  - upSlope;
		else
			split  = 360. - upSlope + split ;

		split  = DegToRad(split );
		x   = spread0Cell *(phiSlope + phiWind*cosf(split ));
		y   = spread0Cell *(phiWind*sinf(split ));
		Rv  = sqrtf(x *x  + y *y );

		spreadMax = spread0Cell  + Rv ;
		phiEw = spreadMax / spread0Cell  - 1;
		a  = asinf(fabsf(y ) / Rv );
		if(x  >= 0.)
			a  = (y  >= 0.) ? a           : M_PI + M_PI - a ;
		else
			a  = (y  >= 0.) ? (M_PI - a ) : M_PI + a ;
		
		split  = RadToDeg(a );
		if (upSlope + split  > 306.)
			azimuthMaxMap[cell] = upSlope + split  - 360.;
		else
			azimuthMaxMap[cell] = upSlope + split ;

		effectiveWind  = powf(phiEw*modelArray[Stride + FUEL_WINDE], 1/windB);
		//Do effective wind only if phiEw > Smidgen
		if(phiEw > Smidgen)
		{		
			maxWind  = 0.9*RxIntensityMap[cell];
			if(effectiveWind  >  maxWind )
			{
				phiEw = windK*powf(maxWind , windB);
				effectiveWind  = maxWind ;
				spreadMax = spread0Cell *(1 + phiEw);
			}
		}

		if(effectiveWind  >  Smidgen)
		{
			lwRatio  = 1. + 0.002840909 * effectiveWind ;
			if (lwRatio  > 1.00001)
				eccentricityMap[cell] = sqrtf(lwRatio *lwRatio  - 1)/lwRatio ;
		}

		spreadMaxCell = spreadMax;
		phiEffWindMap[cell] = phiEw;
	}

	return ( spreadMaxCell);
}


