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

/**
 * \brief struct to represent some useful informations for the current step (color a specific node)
 */

struct Current {
  long color; /**< actual used color */
  long nColors; /**< amount of used colors */
  long rank; /**< actual rank */
  Vertex node; /**< actual vertex */
  long uncoloredVertices; /**< number of uncolored vertices */
  int T; /**< number of color classes with M nodes in it */
  int M; /**< number of nodes in the color class with the highest amount of nodes */
  bool createNewGraphs; /**< TRUE, if new backup graphs should be created */
  int rankNC; /**< rank of the last critical step (where new indep. cliques are created or backtracked) */
  std::vector<long> minBackupGraph; 
  long minBP;
};

/**
 * \brief structure for the backup graphs (for FF's) and some useful variables
 */

struct prevGraphs {
  GraphFord g; /**< backup graph */
  Vertex dlNode; /**< last colored vertices */
  long dlColor; /**< last used color */
  int uncoloredVertices; /**< uncolored vertices */
  std::vector<std::vector<VertexFord> > vertCliques; /**< vector with all vertices in cliques for FF */
  std::vector<std::vector<VertexFord> > vertCliquesColors; /**< vector with all nodes between cliques and colors for FF */
  std::vector<VertexFord> vertColors; /**< vector with all colors for FF */
  std::vector<VertexFord> vertRemaining; /**< vector with all reaming nodes i.e. all nodes that are not in a clique */
  VertexFord source1; /**< vertex for the source (2 FF) */
  VertexFord source2; /**< vertex for the source (1 FF) */
  VertexFord target1; /**< vertex for the target (2 FF) */
  VertexFord target2; /**< vertex for the target (1 FF) */
  
  int n; /**< number of vertices for the first FF */
  int nNew; /**< number of vertices for the second FF */
  int sumLB; /**< sum of the lower bounds */
};

/**
 * \brief structure for the upper and lower bound
 */

struct Bounds {
  long UB; /**< upper bound */
  long LB; /**< lower bound */
};

/**
 * \brief structure to store all relevant backtracking information
 */

struct Backtracking {
  bool status; /**< TRUE, if the algorithm backtracks at the moment */
  Vertex toNode; /**< node where the algorithm heads to */
  int toRank; /**< rank of the node where the algorithm heads to */
  std::vector<int> fRC; 
  bool fastEnd; 
};

/**
 * \brief some useful counts to create studies
 */

struct Count {
  long visitedNodes; /**< amount of visited nodes */
  long newCliques; /**< new created indep. cliques */
  long nFF; /**< how many FF (Ford Fulkerson) are executed */
  long backtracks; /**< how many backtracks the algorithm executed */
  long noNewGraphs; /**< how much new backup graphs are created */
};

/**
 * \brief structure that stores some information to begin the algorithm, or to initialize other components
 */

struct Parameters {
  long nPruningRule; /**< frequency to check for the pruning rule */
  long timeLimit; /**< time limit for the algorithm */
  long threshold; /**< threshold, variable for the PASS-VSS vertex selection strategy */
  std::string ressource; /**< filename where the algorithm should read the graph from */
  long nRandomGraphs; /**< for variant R (random graphs) amount of used random graph instances */
  long n; /**< number of vertices */
  double p; /**< density for an edge between two vertices */
  char variant; /**< R, to use random graphs or N to use a specific graph from the ressource variable */
  char eqDsatur; /**< N, to use the normal eqdsatur without clique, C to use the eqdsatur with clique. */
  int INIT_SRAND; /**< initializer to get a specific random graph instance */
  double pNewCl; /**< frequency to check for a new independent cliques */
  bool wStartCl; /**< TRUE (or 1) to use a start clique, FALSE (or 0) to use NOT a start clique */
};

/**
 * \brief structure that stores all time relevant variables and the timers for executing FF, the whole program, to find indep. cliques
 */

struct Time {
  double timeFF; /**< time for executing FF (the first and second one) */
  double timeUIC; /**< time time to update the independent cliques */
  double timeFindIndCliques; /**< time to find indep. cliques */
  double timeTotal; /**< total time used by the program */
  double timeBuildNetwork; /**< time to build the network */
  std::clock_t start; /**< timer for the whole program */
  std::clock_t startFF; /**< timer for FF */
  std::clock_t startBuildNetwork; /**< timer for building network */
  std::clock_t startFIC; /**< timer for finding indep. cliques */
  std::clock_t startUIC; /**< timer for updating indep. cliques */
  bool timeout; /**< TRUE, if the program timeouted (overstep the time limit) */
};

/**
 * \brief structure that stores the property maps for the given graph, i.e. property maps are needed to access the different properties for a vertex or an edge 
 */

struct PropertyMap {
  NachbarMap n; /**< property map to acess the neighbours of a vertex */
  CliqueMap cl; /**< property map to acess the clique of a vertex */
  ColorMap c; /**< property map to acess the color of a vertex */
  IndexMap i; /**< property map to acess the index of a vertex */
  FBCMap fbc; /**< property map to acess the forbidden colors for a vertex */
  GradMap g; /**< property map to acess the saturation degree for a vertex */
  RankMap r; /**< property map to acess the rank for a vertex */
};

/**
 * \brief structure that stores the property maps for the backup graph, i.e. property maps are needed to acess the different properties for a vertex or an edge
 */

struct PropertyMapFF {
  CapacityMap c; /**< property map to acess the capacity of an edge */
  ReverseEdgeMap re; /**< property map to acess the reverse edge of an edge */
  ResidualCapacityMap rc; /**< property map to acess the capacity of an edge in the corresponing residual network */
  RefVertexMap rf; /**< property map to acess the corresponding vertex in the actual graph (if existant for this vertex in the backup graph) */
};

/**
 * \brief structure that stores informations for the start und indep. cliques
 */

struct Cliques {
  long nodesInClique; /**< number of nodes in the independent cliques, plus the nodes in the start clique */
  long nCliques; /**< number of cliques (inc. start clique) */
  bool newClique; /**< TRUE, if we use a new clique in the actual step */
  long nNodesSameCl; /**< number of nodes that are colored from the indep. cliques */
  long nInCliqueBp; /**< number of nodes in the independent cliques */
  long vertInMinCl;
};

/**
 * \brief structure that represent the color classes
 */

struct Colors {
  std::vector<int> n; /**< vector that represent the color classes */
};

#endif
