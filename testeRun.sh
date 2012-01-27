#!/bin/bash

./createTestMaps

if [ $1 -lt 2 ]; then
	mv aspect_sc0_1.map aspect.map
	mv slope_sc0_1.map slope.map
	rm slope_* aspect_*
	time ./cudaFGM
fi


if [ $1 -eq 2 ]; then
	mv aspect_sc2.map aspect.map
	mv slope_sc2.map slope.map
	rm slope_* aspect_*
	time ./cudaFGM
fi

