#!/bin/bash

#This script creates the deterministic runs for the stochastic solution.
#Argument $1 povides the file with the pdf points.
declare -a ARRAY

#Associate $1 to stdin
exec 10<&0
exec < $1

#Copy real maps from maps folder to local folder
cp malcata_maps/malcata_aspect.grass aspect.map
cp malcata_maps/malcata_slope.grass slope.map

#Read stochastic data and built ARRAY
let count=0
while read line
do
	ARRAY[$count]=$line
	((count++))
done

#Lauch Deterministic FGM 
for ((i=1;i <= ${ARRAY[0]};i++)); do
	#Compute varibles to write into RunSet.in
	windspd=${ARRAY[$((4*($i-1)+1))]}
	winddir=${ARRAY[$((4*($i-1)+2))]}
	fueload=${ARRAY[$((4*($i-1)+3))]}
	fuelMoi=${ARRAY[$((4*($i-1)+4))]}
	if [ $i -lt 10 ]; then 
		tag=00$i
	elif [ $i -lt 100 ]; then  
			tag=0$i;
	else
		tag=$i
	fi
	
	string=$1
	ignmap=ignMap_${string:17:1}_Real_$tag.map
	echo '___________________________________ '
	echo '|Running Simulation #' $tag 
	echo '|Wind Speed=' $windspd 
	echo '|Wind Direction=' $winddir
	echo '|Fuel Load=' $fueload
	echo '|Fuel Moisture=' $fuelMoi

	
	#Rewrite RunSet.in wind new wind speed, direction, and ignMap file name
	sed '30 s/.*/'$windspd'/' RunSet.in >RunSet_tmp.in 
	sed '31 s/.*/'$winddir'/' RunSet_tmp.in >RunSet_tmp2.in
	sed '32 s/.*/'$fuelMoi'/' RunSet_tmp2.in >RunSet_tmp3.in
	sed '37 s/.*/'$fueload'/' RunSet_tmp3.in >RunSet_tmp4.in
	sed '41 s/.*/'$ignmap'/' RunSet_tmp4.in >RunSet.in

	#Run fgm 
	./cudaFGM 1

	./matrix_to_column 400 $ignmap
done

#Remove temporary files
rm RunSet_tmp*

#move ignition maps do folder ignMap_data
mv ignMap_* ignMapData/
mv 1D_ignMap* ignMapData/

scp ignMapData/1D_ignMap_${string:17:1}_Real* fsousa@zappa:~/Dropbox/shared_Rita/

#restore filedescriptor 10
exec 0<&10 10<&-
