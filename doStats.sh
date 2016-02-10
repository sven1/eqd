#!/bin/bash

printHeader () {
	echo "Variant;Random Graph;N;P;Time_Total;Start_LB;Start_UB;LB;UB;Visited Nodes;Number_Of_Backtracks;Number_Of_FF;Number_Of_New_Cliques;Time_FF;Time_Find_Indep_Cliques;Time_Update_Indep_Cliques;Timeout;" >> ${1}
}

checkExistDir () {
	if [ ! -d ${DIR}/${TIMELIMIT}/${N} ]
	then
		mkdir -p ${DIR}/${TIMELIMIT}/${N}
	fi
}

startStudy () {
	if [ ${1} == N ]
	then
		FILE_EQD=${DIR}/${5}/${2}/EQD_${2}_${3}_${5}

		printHeader ${FILE_EQD}

		for ((i=1; i<=${4}; i++))
		do
			./eqDsatur R N 1 ${5} 3 ${4} ${2} ${3} ${i} 0 1 >> ${FILE_EQD}
		done
	elif [ ${1} == C ]
	then
		FILE_EQDC=${DIR}/${5}/${2}/EQDC_${N}_${3}_${5}

		printHeader ${FILE_EQDC}

		for ((i=1; i<=${4}; i++))
		do
			./eqDsatur R C 1 ${5} 3 ${4} ${2} ${3} ${i} 0 1 >> ${FILE_EQDC}
		done
	else
		FILE_EQD=${DIR}/${5}/${N}/EQD_${2}_${3}_${5}
		FILE_EQDC=${DIR}/${5}/${N}/EQDC_${2}_${3}_${5}

		printHeader ${FILE_EQD}
		printHeader ${FILE_EQDC}
		
		for ((i=1; i<=${4}; i++))
		do
			./eqDsatur R N 1 ${5} 3 ${4} ${2} ${3} ${i} 0 1 >> ${FILE_EQD}
			./eqDsatur R C 1 ${5} 3 ${4} ${2} ${3} ${i} 0 1 >> ${FILE_EQDC}
		done
	fi
}

DATE=$(date +"%m-%d-%y")
DIR=results_${DATE}

if [ $# -eq 0 ]
then
	echo "Usage: ./doStats [ N | C | B ] N [ P | P_MIN P_MAX P_STEP ] N_RANDOM_GRAPHS TIMELIMIT"
elif [ $# -eq 5 ]
then
	VARIANT=$1
	N=$2
	P=$3
	N_RANDOM_GRAPHS=$4
	TIMELIMIT=$5

	checkExistDir

	startStudy ${VARIANT} ${N} ${P} ${N_RANDOM_GRAPHS} ${TIMELIMIT} 
elif [ $# -eq 7 ]
then
	VARIANT=$1
	N=$2
	P_MIN=$3
	P_MAX=$4
	P_STEP=$5
	N_RANDOM_GRAPHS=$6
	TIMELIMIT=$7

	checkExistDir

	for j in `seq ${P_MIN} ${P_STEP} ${P_MAX}`
	do
		startStudy ${VARIANT} ${N} ${j} ${N_RANDOM_GRAPHS} ${TIMELIMIT} 
	done
else
	echo "error: invalud number of arguments"
fi
