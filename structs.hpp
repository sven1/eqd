/**
 * \file structs.hpp
 * \author Sven Förster
 * \date 29.12.2015
 *
 * \version 1.0.0
 * cleaning code
 */

#ifndef STRUCTS_HPP
#define STRUCTS_HPP

#include "boost.hpp"
#include <string>
#include <ctime>

struct Current {
  long color;
  long nColors;
  long rank;
  Vertex node;
  long uncoloredVertices;
  int T;
  int M;
  bool createNewGraphs;
  int rankNC;
  std::vector<long> minBackupGraph;
  long minBP;
};

struct prevGraphs {
  GraphFord g;
  Vertex dlNode;
  long dlColor;
  int uncoloredVertices;
  std::vector<VertexFord> vert;
  int n;
  int nNew;
  int sumLB;
};

struct Bounds {
  long UB;
  long LB;
};

struct Backtracking {
  bool status;
  Vertex toNode;
  int toRank;
  std::vector<int> fRC;
  bool fastEnd;
};

struct Count {
  long visitedNodes;
  long newCliques;
  long nFF;
  long backtracks;
  long noNewGraphs;
};

struct Parameters {
  long nPruningRule;
  long timeLimit;
  long threshold;
  std::string ressource;
  long nRandomGraphs;
  long n;
  double p;
  char variant;
  char eqDsatur;
  int INIT_SRAND;
  double pNewCl;
  bool wStartCl;
};

struct Time {
  double timeFF;
  double timeUIC;
  double timeFindIndCliques;
  double timeTotal;
  std::clock_t start;
  std::clock_t startFF;
  std::clock_t startFIC;
  std::clock_t startUIC;
  bool timeout;
};

struct PropertyMap {
  NachbarMap n;
  CliqueMap cl;
  ColorMap c;
  IndexMap i;
  FBCMap fbc;
  GradMap g;
  RankMap r;
};

struct PropertyMapFF {
  CapacityMap c;
  ReverseEdgeMap re;
  ResidualCapacityMap rc;
  RefVertexMap rf;
};

struct Cliques {
  long nodesInClique;
  long nCliques;
  bool newClique;
  long nNodesSameCl;
  long nInCliqueBp;
};

struct Colors {
  std::vector<int> n;
};

#endif
