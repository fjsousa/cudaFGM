#!/bin/bash
#Creates the maps used in the test cases for my master thesis

./createTestMaps

if [ $1 -lt 2 ]; then
	mv aspect_sc0_1.map aspect.map
	mv slope_sc0_1.map slope.map
	rm slope_* aspect_*
	time ./cudaFGM 0
fi


if [ $1 -eq 2 ]; then
	mv aspect_sc2.map aspect.map
	mv slope_sc2.map slope.map
	rm slope_* aspect_*
	time ./cudaFGM 0
fi

./matrix2col_plus_rowCol 600 ignMap_V.map
