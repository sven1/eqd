#!/bin/bash

DATE=$(date +"%m-%d-%y")
DIR=results_${DATE}

#mkdir ${DIR}

if [ $# -eq 2 ]
then
	N_RANDOM_GRAPHS=200
	N=$1
	P=$2

	for ((i=1; i<=200; i++))
	do
		./eqDsatur_FIC R C 1 300 3 $N_RANDOM_GRAPHS $N $P $i >> TEST_65_07_FIC
	done
fi
