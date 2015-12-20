C = gcc
CXX = g++
CXXFLAGS = -std=c++11 -Wall $(OPT)
DEBUG = -g
OPT = -O3
NOOPT = -O0
DATE = $(shell date +"%m-%d-%y")

all : eqDsatur

eqDsatur : EqColoring.o Coloring.o Input.o Heuristic.o main.cpp
	$(CXX) $(CXXFLAGS) -o eqDsatur main.cpp EqColoring.o Coloring.o Input.o Heuristic.o

EqColoring.o : EqColoring.cpp EqColoring.hpp Coloring.hpp
	$(CXX) $(CXXFLAGS) -c EqColoring.cpp -o EqColoring.o

Coloring.o : Coloring.cpp Coloring.hpp
	$(CXX) $(CXXFLAGS) -c Coloring.cpp -o Coloring.o

Input.o : Input.cpp Input.hpp
	$(CXX) $(CXXFLAGS) -c Input.cpp -o Input.o

Heuristic.o : Heuristic.cpp Heuristic.hpp
	$(CXX) $(CXXFLAGS) -c Heuristic.cpp -o Heuristic.o

clean :
	rm eqDsatur *.o

ctar :
	rm *.tar

tar :
	tar -cvf eqDsatur_$(DATE).tar *.hpp *.cpp makefile Study.sh StudyOld.sh
