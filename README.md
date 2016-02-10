# EqDsatur

## Short Introduction

EqDsatur is an algorithm to color equitable a graph. Equitable means that the
color classes differs at most by one.

## Usage

To start EqDsatur simply type
	./eqDsatur
in the terminal. Without additional parameters you get some help which
parameters are missing (and necessary) and which you can add if its desired.

## Parameters

1. You can choose between N for normal mode, to select a specific graph as
   input, or R for random graph mode, to generate some random graphs as input.

2. You can choose between N for normal mode, to select the normal dsatur
   algorithm, or C for clique mode, to select **our** implementation with a
   pruning rule which uses cliques.

3. Doesn't matter anymore.

4. It's the time limit. Mostly 3600 or 300 seconds.

5. It's a treshold which is used by the vertex selection strategy method. It's
   recommanded to use *3* as value.

6. If you are in normal mode to select a graph, then the 6th parameter is the
   ressource file of the graph. Otherwise, if you are in random mode, you should
   add 4 other parameters.

   The first one doesn't matter anymore, you could put any number in there. The
   second one describes how much vertices you want. The 3th describes the density
   (global) for each edge. The 4th describes the specific random graph, it's a
   seed value. 

   If you want for instance 200 different random graphs you should put in
   different values each time, for example 1, 2, ..., 200.
   If you want always the same random graph, then put the same value in.
   With this value you can refer to a specific random graph.

7.  It's a percentage value, between 0 and 1. If our clique is by this value
	covered, then we are looking for a better one, if it's exist, otherwise we
	keep the old cliques.

8. It's a boolean value which describes if we use a start clique or not. 1 or 0
   are allowed values.

## References
