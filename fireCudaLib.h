/*
 *******************************************************************************
 *
 *  fireCudaLib.h
 *
 *  Description
 *      Library of BEHAVE (Andrews 1986) fire behavior algorithms
 *      encapsulated and optimized for fire behavior simulation.
 *
 *  Legalities
 *      Copyright (C) 1996-2004 by Collin D. Bevins.  All rights reserved.
 *
 *  Description
 *      This header file describes the externally-visible facilities of
 *      the Fire Behavior Library C API.
 *
 *      This file really needs to be split into public and private portions.
 *
 *  History
 *      1996/09/04  Release 1.0.0
 *
 *      1999/03/05  Fixed NFFL07 live SAVR from 1500 to 1550.
 *      1999/03/05  Release 1.0.1
 *
 *      1999/04/22  Added "backing" to FuelModelData structure.
 *                  Added Fuel_BackingRate(catalog, model) access macro.
 *                  Added calculation of backing spread rate to function
 *                      Fuel_SpreadWindSlopeMax().
 *                  Removed ANSI_ARG and VOID constructs.
 *                  Added conditional M_PI definition.
 *                  Removed license.txt from distribution.
 *      1999/04/22  Release 1.0.2
 *
 *      1999/12/01  Fixed eccentricity sqrt(0) bug around line 796.
 *                  Fixed possible division by zero when determining SpreadAny
 *                  around lines 895-905.
 *      1999/12/01  Release 1.0.3
 *
 *      2002/09/30  Added 1 to number of max fuel models allocated
 *                  around line 1280
 *                  calloc(FuelCat_MaxModels(catalog)+1, sizeof(FuelModelPtr)))
 *                  otherwise it was freeing to many from the array
 *		2002/09/30	Release 1.0.4
 *******************************************************************************
 */

#ifndef FIRE_LIB_
#define FIRE_LIB_ 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
//#include <math_functions.h>
#define FIRELIB_VERSION "1.0.4"
#define FIRELIB_MAJOR_VERSION 1
#define FIRELIB_MINOR_VERSION 0
#define FIRELIB_PATCH_LEVEL   4

/*
 *------------------------------------------------------------------------------
 * Disable C++ name mangling.
 *------------------------------------------------------------------------------
 */

/*
#ifdef __cplusplus
#   define EXTERN extern "C"
#else
#   define EXTERN extern
#endif
*/

////////////////////////////////////////////////////////////////////////////////
//Problem parameters
//
////////////////////////////////////////////////////////////////////////////////

#define TAG "CUDA8_"
#define Stencil16 0
#define Clatering  1

#define INF 			9999999. /* or close enough */
#define BLOCK_SIZE 16

#define sizeModel (9)
#define sizePart  (2)


#define FUEL_MEXT         (0)
#define FUEL_LIFEAREAWTG  (1)
#define FUEL_LIFERXFACTOR (2)
#define FUEL_BULKDENSITY  (3)
#define FUEL_PROPFLUX     (4)
#define FUEL_SLOPEK       (5)
#define FUEL_WINDB        (6)
#define FUEL_WINDE        (7)
#define FUEL_WINDK        (8)

#define FUEL_AREAWTG      (9)
#define FUEL_SIGMAFACTOR  (10)

#define thx (threadIdx.x)
#define thy (threadIdx.y)
/*
 *------------------------------------------------------------------------------
 *  Macro pseudo functions.
 *------------------------------------------------------------------------------
 */

#define Smidgen                 (0.000001)
#define MetersToFeet(x)					((x)/0.3048 )
#define FeetToMeters(x) 				((x)*0.3048 )
#define DegToRad(x)				      ((x)*0.017453293)
#define RadToDeg(x)     				((x)*57.29577951)
#define IsZero(x)               (fabs(x)<Smidgen)
#define Equal(x,y)              (fabs((x)-(y))<Smidgen)
#define _Equal(x,y)             (fabsf((x)-(y))<Smidgen)

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

/*
 *------------------------------------------------------------------------------
 * Firelib return status codes.
 *------------------------------------------------------------------------------
 */

#define  FIRE_STATUS_OK         (0)
#define  FIRE_STATUS_ERROR      (-1)
#define  FIRE_STATUS_EOF        (1)

/*
 *------------------------------------------------------------------------------
 *  Fuel moisture and mass weighting classes.
 *------------------------------------------------------------------------------
 */

#define  FIRE_LIFE_CATS     (2) /* Number of fuel particle life categories */
#define  FIRE_LIFE_DEAD     (0)
#define  FIRE_LIFE_LIVE     (1)

#define  FIRE_SIZE_CLASSES  (6) /* Number of mass weighting classes. */

#define  FIRE_MCLASSES      (6) /* Number of fuel moisture classes. */
#define  FIRE_MCLASS_1HR    (0)
#define  FIRE_MCLASS_10HR   (1)
#define  FIRE_MCLASS_100HR  (2)
#define  FIRE_MCLASS_1000HR (3)
#define  FIRE_MCLASS_HERB   (4)
#define  FIRE_MCLASS_WOOD   (5)

/*
 *------------------------------------------------------------------------------
 *  FuelParticleData structure: fuel particle input and intermediate attributes.
 *------------------------------------------------------------------------------
 */

typedef struct fuelParticleDataStruct
{
    /* INPUT */
    float load;                	/* fuel loading                     (lb/sqft) */
    float savr;                	/* surface area-to-volume ratio        (1/ft) */
    float dens;                	/* particle density                 (lb/cuft) */
    float heat;                	/* heat of combustion                (BTU/lb) */
    float stot;                	/* total silica content        (fraction odw) */
    float seff;                	/* effective silica content    (fraction odw) */
    /* PARTICLE_DEPENDENT */
    float area;                	/* surface area */
    float sigma;               	/* exp(-138./sigma)                      (dl) */
    /* MODEL-DEPENDENT */
    float awtg;                	/* surface area derived weighting factor (dl) */
    float gwtg;                	/* size class area weighting factor */
    /* ENVIRONMENT-DEPENDENT */
    float mois;       	        /* particle moisture content       (fraction) */
    size_t live;                /* life category 0=dead, 1=live               */
    size_t type;              	/* type category 0=dead, 1=herb, 2=live woody */
    size_t sizeClass;           /* fuel moisture size class                   */
} FuelParticleData, *FuelParticlePtr, *PartPtr;

#define FIRE_TYPE_DEAD   (1)
#define FIRE_TYPE_HERB   (2)
#define FIRE_TYPE_WOOD   (3)

/* FuelParticleData structure access macros. */

#define Fuel_Live(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->live)

#define Fuel_Type(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->type)

#define Fuel_SizeClass(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->sizeClass)

#define Fuel_Load(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->load)

#define Fuel_Savr(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->savr)

#define Fuel_Heat(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->heat)

#define Fuel_Density(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->dens)

#define Fuel_SiTotal(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->stot)

#define Fuel_SiEffective(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->seff)

#define Fuel_SurfaceArea(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->area)

#define Fuel_AreaWtg(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->awtg)

#define Fuel_SizeAreaWtg(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->gwtg)

#define Fuel_SigmaFactor(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->sigma)

#define Fuel_Moisture(catalog,model,particle) \
                ((catalog)->modelPtr[(model)]->partPtr[(particle)]->mois)

/*
 *------------------------------------------------------------------------------
 *  FuelModelData structure: fuel model bed input attributes.
 *------------------------------------------------------------------------------
 */

typedef struct fuelModelDataStruct
{
    /* Input variables. */
    size_t modelId;             /* ! fuel model number                          */
    size_t combustion;          /* 0 if combustion not yet calculated         */
    size_t maxParticles;        /* maximum number of FuelParticles            */
    size_t particles;           /* current number of FuelParticles            */
    PartPtr *partPtr;           /* array of pointers to Fuel Particles        */
    char  *name;                /* ! fuel model short name                      */
    char  *desc;                /* ! fuel model description text                */
    char  *reserved1;           /* ! used for alignment                         */
    float depth;        	    /* fuel bed depth                        (ft) */
    float mext;                 /* dead fuel extinction moisture   (fraction) */
    float adjust;               /* spread rate adjustment factor         (dl) */
    /* Combustion intermediates. */
    float awtg[2];             	/* dead & live fuel area weighting factors    */
    float rxFactor[2];         	/* dead and live fuel rx factors              */
    float fineDead;            	/* fine dead fuel ratio                       */
    float liveFactor;          	/* live fuel moisture extinction factor       */
    float rhob;                	/* fuel bed bulk density                      */
    float taur;                	/* residence time                       (min) */
    float propFlux;            	/* propagating flux ratio                     */
    float slopeK;              	/* slope parameter 'k'                        */
    float windB;               	/* wind parameter 'b'                         */
    float windE;               	/* wind parameter (ratio**e/c)                */
    float windK;               	/* wind parameter (c * ratio**-e)             */
    /* Current environment. */
    float moisture[FIRE_MCLASSES]; /* array of fuel moistures (fraction odw) */
    float windFpm;             	/* wind speed                        (ft/min) */
    float windDeg;             	/* wind vector         (degrees from upslope) */
    float slope;               	/* slope                         (rise/reach) */
    float aspect;              	/* aspect (downslope) azimuth  (compass degs) */
    /* Updated by Fire_SpreadNoWindNoSlope() */
    float rxInt;               	/* reaction intensity          (BTU/sqft/min) */
    float spread0;             	/* no-wind, no-slope spread rate     (ft/min) */
    float hpua;                	/* heat per unit area              (BTU/sqft) */
    /* Updated by Fire_SpreadWindSlopeMax() */
    float spreadMax;           	/* spread in direction of max spread (ft/min) */
    float azimuthMax;          	/* direction of maximum spread      (degrees) */
    float effWind;             	/* effective windspeed                        */
    float lwRatio;             	/* length-to-width ratio for eff windspeed    */
    float eccentricity;        	/* eccentricity of ellipse for eff windspeed  */
    float backing;             	/* backing spread rate               (ft/min) */
    float phiW;                	/* wind factor                                */
    float phiS;                	/* slope factor                               */
    float phiEw;               	/* combined wind-slope factor                 */
    size_t wLimit;              /* wind limit 0=not reached, 1=reached        */
    size_t reserved2;           /* ! used for alignment                         */
    /* Updated by Fire_SpreadAtAzimuth() */
    float spreadAny;           	/* spread rate at arbitrary azimuth  (ft/min) */
    float azimuthAny;          	/* direction of arbitrary spread    (degrees) */
    float byrams;              	/* fireline intensity              (BTU/ft/s) */
    float flame;               	/* flame length                          (ft) */
    float scorch;              	/* scorch height                         (ft) */
} FuelModelData, *FuelModelPtr;

/* Fuel model input variable macros. */
#define Fuel_Model(catalog,model) \
                    ((catalog)->modelPtr[(model)]->modelId)

#define Fuel_Name(catalog,model) \
                    ((catalog)->modelPtr[(model)]->name)

#define Fuel_Desc(catalog,model) \
                    ((catalog)->modelPtr[(model)]->desc)

#define Fuel_Depth(catalog,model) \
                    ((catalog)->modelPtr[(model)]->depth)

#define Fuel_Mext(catalog,model) \
                    ((catalog)->modelPtr[(model)]->mext)

#define Fuel_SpreadAdjustment(catalog,model) \
                    ((catalog)->modelPtr[(model)]->adjust)

#define Fuel_CombustionFlag(catalog,model) \
                    ((catalog)->modelPtr[(model)]->combustion)

#define Fuel_MaxParticles(catalog,model) \
                    ((catalog)->modelPtr[(model)]->maxParticles)

#define Fuel_Particles(catalog,model) \
                    ((catalog)->modelPtr[(model)]->particles)

#define Fuel_ParticleArray(catalog,model) \
                    ((catalog)->modelPtr[(model)]->partPtr)

#define Fuel_ParticlePtr(catalog,model,particle) \
                    ((catalog)->modelPtr[(model)]->partPtr[(particle)])

/* Fuel model combustion intermediates macros. */
#define Fuel_LifeAreaWtg(catalog,model,life) \
                    ((catalog)->modelPtr[(model)]->awtg[(life)])

#define Fuel_LifeRxFactor(catalog,model,life) \
                    ((catalog)->modelPtr[(model)]->rxFactor[(life)])

#define Fuel_FineDead(catalog,model) \
                    ((catalog)->modelPtr[(model)]->fineDead)

#define Fuel_LiveMextFactor(catalog,model) \
                    ((catalog)->modelPtr[(model)]->liveFactor)

#define Fuel_BulkDensity(catalog,model) \
                    ((catalog)->modelPtr[(model)]->rhob)

#define Fuel_ResidenceTime(catalog,model) \
                    ((catalog)->modelPtr[(model)]->taur)

#define Fuel_PropFlux(catalog,model) \
                    ((catalog)->modelPtr[(model)]->propFlux)

#define Fuel_SlopeK(catalog,model) \
                    ((catalog)->modelPtr[(model)]->slopeK)

#define Fuel_WindB(catalog,model) \
                    ((catalog)->modelPtr[(model)]->windB)

#define Fuel_WindE(catalog,model) \
                    ((catalog)->modelPtr[(model)]->windE)

#define Fuel_WindK(catalog,model) \
                    ((catalog)->modelPtr[(model)]->windK)

/* Fuel model fire behavior variable macros. */
#define Fuel_RxIntensity(catalog,model) \
                    ((catalog)->modelPtr[(model)]->rxInt)

#define Fuel_Spread0(catalog,model) \
                    ((catalog)->modelPtr[(model)]->spread0)

#define Fuel_HeatPerUnitArea(catalog,model) \
                    ((catalog)->modelPtr[(model)]->hpua)

#define Fuel_SpreadMax(catalog,model) \
                    ((catalog)->modelPtr[(model)]->spreadMax)

#define Fuel_AzimuthMax(catalog,model) \
                    ((catalog)->modelPtr[(model)]->azimuthMax)

#define Fuel_SpreadAny(catalog,model) \
                    ((catalog)->modelPtr[(model)]->spreadAny)

#define Fuel_AzimuthAny(catalog,model) \
                    ((catalog)->modelPtr[(model)]->azimuthAny)

#define Fuel_EffectiveWind(catalog,model) \
                    ((catalog)->modelPtr[(model)]->effWind)

#define Fuel_LwRatio(catalog,model) \
                    ((catalog)->modelPtr[(model)]->lwRatio)

#define Fuel_Eccentricity(catalog,model) \
                    ((catalog)->modelPtr[(model)]->eccentricity)

#define Fuel_BackingRate(catalog,model) \
                    ((catalog)->modelPtr[(model)]->backing)

#define Fuel_PhiWind(catalog,model) \
                    ((catalog)->modelPtr[(model)]->phiW)

#define Fuel_PhiSlope(catalog,model) \
                    ((catalog)->modelPtr[(model)]->phiS)

#define Fuel_PhiEffWind(catalog,model) \
                    ((catalog)->modelPtr[(model)]->phiEw)

#define Fuel_WindLimit(catalog,model) \
                    ((catalog)->modelPtr[(model)]->wLimit)

#define Fuel_ByramsIntensity(catalog,model) \
                    ((catalog)->modelPtr[(model)]->byrams)

#define Fuel_FlameLength(catalog,model) \
                    ((catalog)->modelPtr[(model)]->flame)

#define Fuel_ScorchHeight(catalog,model) \
                    ((catalog)->modelPtr[(model)]->scorch)

/* Fuel model environment variable macros. */
#define Fuel_EnvMoisture(catalog,model,mclass) \
                    ((catalog)->modelPtr[(model)]->moisture[(mclass)])

#define Fuel_WindSpeed(catalog,model) \
                    ((catalog)->modelPtr[(model)]->windFpm)

#define Fuel_WindDir(catalog,model) \
                    ((catalog)->modelPtr[(model)]->windDeg)

#define Fuel_Slope(catalog,model) \
                    ((catalog)->modelPtr[(model)]->slope)

#define Fuel_Aspect(catalog,model) \
                    ((catalog)->modelPtr[(model)]->aspect)

/*
 *------------------------------------------------------------------------------
 *  FuelCatData structure; provides a complete fuel catalog.
 *------------------------------------------------------------------------------
 */

#define FIRE_CATALOG_MAGIC      (19520904L)
#define FIRE_ERROR_BUFFER_SIZE  (1024)

typedef struct fuelCatalogStruct
{
    long      magicCookie;      /* ! magic cookie for sanity checking           */
    int       status;           /* ! return status of most recent call          */
    size_t    maxModels;        /* maximum number of models in this catalog   */
    size_t    flameClasses;     /* size of the flame length array             */
    char         *name;         /* ! name for this catalog instance             */
    char         *error;        /* ! error message buffer                       */
    FuelModelPtr *modelPtr;     /* array of ModelPtr[maxModels+1]             */
    float       *flamePtr;   	/* flame length lookup array                  */
    float        flameStep;    	/* size of each flame length table class (ft) */
} FuelCatalogData, *FuelCatalogPtr;

#define FuelCat_MagicCookie(catalog)    (catalog->magicCookie)
#define FuelCat_MaxModels(catalog)      (catalog->maxModels)
#define FuelCat_Status(catalog)         (catalog->status)
#define FuelCat_FlameClasses(catalog)   (catalog->flameClasses)
#define FuelCat_FlameStep(catalog)      (catalog->flameStep)
#define FuelCat_FlameArray(catalog)     (catalog->flamePtr)
#define FuelCat_Name(catalog)           (catalog->name)
#define FuelCat_Error(catalog)          (catalog->error)
#define FuelCat_ModelArray(catalog)     (catalog->modelPtr)
#define FuelCat_ModelPtr(catalog,model) (catalog->modelPtr[model])


//------------------------------------------------------------------------------
//Header of new functions. Compiler complaning - don't know why declaration
//in file were the function is used isn't enough.
//------------------------------------------------------------------------------

int	createMaps(int, int, double, double );

float* ExactCircle( int, int, float, float);

float* ExactElipse( int, int, float, float, float, float);

int PrintMap ( float*, char* );

int errorStuff(int, float, int, float, int, float*, float*, float);

void Preprocessor(FuelCatalogPtr);

int Print_CatalogStruct(FuelCatalogPtr);

void checkCUDAError(const char*);

float* CudaModelData(FuelCatalogPtr);

float NoWindNoSlope(	int, float, size_t, float*, float* );

float WindAndSlope(	int,
										float*, 
										float, 
										float, 
										float, 
										float, 
										float, 
										size_t, 
										float*, 
										float*, 
										float*, 
										float*); 


__global__ void FireKernel_NoWindNoSlope(			float*,
																							float*,
																							float*,
																							size_t*,
											  											float* );

__global__ void FireKernel_WindAndSlope(	float*, 
											  									float*, 
											  									float*, 
											  									float*, 
											  									float*, 
											  									float*, 
											  									float*, 
											  									float*, 
																					float*, 
																					float*, 
																					size_t*,
																					float*); 


__global__ void FireKernel_SpreadAtNeighbors( float*, 
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float*,
											  											float,
																							float*);
/*
 *------------------------------------------------------------------------------
 *  Function prototypes for fire behavior computations.
 *------------------------------------------------------------------------------
 */

#define FIRE_NONE       (0)
#define FIRE_BYRAMS     (1)
#define FIRE_FLAME      (2)
#define FIRE_SCORCH     (4)

int Fire_FuelCombustion ( 
	FuelCatalogPtr,
	size_t
) ;


/*
 *------------------------------------------------------------------------------
 *  Function prototypes for creating and destroying fuel catalogs, fuel models,
 *  fuel particles, and flame length tables.
 *------------------------------------------------------------------------------
 */


int Fire_FlameLengthTable (
    FuelCatalogPtr catalog,     /* FuelCatalogData instance pointer           */
    size_t  flameClasses,       /* number of flame length classes             */
    float  flameStep	        /* flame length step value per class          */
) ;

FuelCatalogPtr Fire_FuelCatalogCreate (
    char  *name,                /* FuelCatalogData instance name              */
    size_t maxModels            /* maximum modelId allowed in this catalog    */
) ;

FuelCatalogPtr Fire_FuelCatalogCreateStandard (
    char  *name,                /* FuelCatalogData instance name              */
    size_t maxModels            /* maximum modelId allowed in this catalog    */
) ;

int Fire_FuelCatalogDestroy (
    FuelCatalogPtr catalog      /* FuelCatalogData instance pointer           */
) ;

int Fire_FuelModelCreate (
    FuelCatalogPtr catalog,     /* FuelCatalogData instance                   */
    size_t  model,              /* fuel model number            [0-maxModels] */
    char   *name,               /* short name                                 */
    char   *desc,               /* longer description                         */
    float  depth,	            /* bed depth                             (ft) */
    float  mext,    			/* moisture of extinction                (dl) */
    float  adjust,             	/* spread adjustment factor              (dl) */
    size_t  maxParticles        /* maximum number of fuel model particles     */
) ;

int Fire_FuelModelDestroy (
    FuelCatalogPtr catalog,     /* FuelCatalogData instance pointer           */
    size_t         model        /* fuel model id number         [0-maxModels] */
) ;

int Fire_FuelModelExists (
    FuelCatalogPtr catalog,     /* FuelCatalogData instance pointer           */
    size_t         model        /* fuel model id number         [0-maxModels] */
) ;

int Fire_FuelParticleAdd (
    FuelCatalogPtr catalog,     /* FuelCatalogData instance pointer           */
    size_t  model,              /* fuel model id number         [0-maxModels] */
    size_t  type,               /* FIRE_TYPE_DEAD, _TYPE_HERB, or _TYPE_WOOD  */
    float  load,      	        /* fuel load                        (lbs/ft2) */
    float  savr,         	    /* surface-area-to-volume ratio     (ft2/ft3) */
    float  dens,            	/* density                          (lbs/ft3) */
    float  heat,               	/* heat of combustion               (btus/lb) */
    float  stot,               	/* total silica content               (lb/lb) */
    float  seff                	/* effective silica content           (lb/lb) */
) ;

#ifdef NEED_STRDUP
char *strdup ( const char *str ) ;
#endif

#endif

/*
 *******************************************************************************
 * End of fireLib.h
 *******************************************************************************
 */
