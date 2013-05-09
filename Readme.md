#cudaFGM

A fire growth model (FGM) that runs on NVIDIA GPUs. 

 
## fireLib 

firelib is used to create most of the fire properties. Main firelib 
funcions used are "Fire_FuelCatalogCreateStandard" to create the fire 
catalog, "Fire_SpreadNoWindNoSlope" and "Fire_SpreadWindSlopeMax" are 
used to compute fire ellipse properties. 



##Main Inputs

Grass aspect and slope files can be read with "./cudaFGM 1". Grass slope files must be written in "percentage" of rise/reach.

Inputs provided in "RunSet.in": 

* Map width (meters)
* Map height
* Number of rows
* Number of cols
* Fuel Model (NFFL models) or 14-custom	
* Wind speed
* Wind Direction
* Moisture (M1, M10, M100, Mherb, Mwood)
* Custom Particle load (Not relevant if using 
		one of theNFFL 13) 
* ignition point (X) as a %(0-1) of map width
* ignition point (Y) as a %(0-1) of map height
* GPU device
* ignition map file name (no spaces)
* Verbosity (1 - more 0 - less)
init ignMap to BEHAVE elipse for faster solution (1 - Yes, 0 - No)

##Outputs

Main output is the ignition map "ignMap_###.dat", where ### is a user 
defined tag

##Stuff to do

Organise and trash unnecessary test cases cluttering the project
