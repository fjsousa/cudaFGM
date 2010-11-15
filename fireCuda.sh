#!/bin/bash

#Rows Cols
RUNS=11
NODES[0]=128
NODES[1]=256
NODES[2]=400
NODES[3]=512
NODES[4]=800
NODES[5]=912
NODES[6]=1024
NODES[7]=1248
NODES[8]=1520
NODES[9]=1696
NODES[10]=2048

Residue=0.00000001
nScn=4

N_isoLines=5
upperIsoLines[0]=3000
upperIsoLines[1]=1000
upperIsoLines[2]=1000
upperIsoLines[3]=3000

for ((i=3; i<$nScn; i++)); do	
	for ((j=0; j<$RUNS; j++)); do
		echo 0 $i ${NODES[$j]} ${NODES[$j]}  $Residue $N_isoLines ${upperIsoLines[$i]} > RunSet.in
		./fireCuda
		done
	done

mv *.dat ./data
