#!/bin/bash

for ((i=40; i<=65; i+=5))
do
	./doStats.sh B ${i} 0.1 0.9 0.1 200 3600
done
