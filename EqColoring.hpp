/**
 * @file EqColoring.hpp
 * @brief equitable coloring for a given graph
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#ifndef EQCOLORING_HPP
#define EQCOLORING_HPP

#include "structs.hpp"
#include "Coloring.hpp"

/**
 * \class EqColoring
 * \brief class for coloring a graph instance in equitable color classes (differ at most by one)
 */
class EqColoring : Coloring{
  private:
	  /**
	   * @brief backup graph
	   */
    prevGraphs gf;
	/**
	 * @brief property map for the backup graph
	 */
    PropertyMapFF pmf;

	/**
	 * @brief connect property map to the backup graph
	 * @return 
	 */
    bool initPrevGraphsFF();
	/**
	 * @brief create the backup graph for a specific color
	 * @param color color
	 * @return 
	 */
    bool initBackupGraphs(int color);

	/**
	 * @brief init. the network for the first FF 
	 * @param colors colors
	 * @return 
	 */
	bool initNetwork(int colors);
	/**
	 * @brief init. the network for the second FF
	 * @param colors colors
	 * @return 
	 */
	bool initNetwork2(int colors);
    
	/**
	 * @brief perform FF (i.e. edmunds karp max flow)
	 * @param vs source
	 * @param vt sink
	 * @return max flow
	 */
    long performEKMF(VertexFord &vs, VertexFord &vt);

	/**
	 * @brief update independent cliques if better indep. cliques are found after coloring vertex v
	 * @param v vertex
	 * @return TRUE, if better cliques are found, otherwise FALSE
	 */
    bool updateIndepCliques(Vertex &v);

	/**
	 * @brief reset capacity network
	 */
    void resetCap();

  public:
    EqColoring();
    EqColoring(const Parameters &parm);
    EqColoring(const Graph &g);
    EqColoring(const Graph &g, const Parameters &parm);

    ~EqColoring();

	/**
	 * @brief recursive (main) function for the normal dsatur algorithm
	 * @return 
	 */
    bool node();
	/**
	 * @brief recursive (main) function for the clique dsatur algorithm
	 * @return 
	 */
    bool nodeClique();

	/**
	 * @brief starts normal dsatur algorithm
	 * @return 
	 */
    bool dsatur();
	/**
	 * @brief starts clique dsatur algorithm
	 * @return 
	 */
    bool dsaturClique();

	/**
	 * @brief pruning rule from the paper
	 * @return 
	 */
    bool pruningRulePaper();
	/**
	 * @brief pruning rule with cliques (FF)
	 * @return 
	 */
    bool pruningRuleFF();

	/**
	 * @brief check equitability for this graph
	 * @return TRUE, if the color classes differs at most by one, otherwise FALSE
	 */
    bool checkEquitability() const;
	/**
	 * @brief check if it is an equitable coloring
	 * @return TRUE, if its an eq. coloring, otherwise FALSE
	 */
    bool checkEqColoring() const;

	/**
	 * @brief search for better indep. cliques
	 * @param sBetterClique TRUE, search for BETTER cliques, otherwise search for ANY indep. clique
	 * @return 
	 */
    bool useNewIndepCliques(bool sBetterClique);

	/**
	 * @brief check if the algorithm could prune here
	 * @return TRUE, then prune, otherwise continue
	 */
    bool pruneFF();
	/**
	 * @brief check if the algorithm could prune for this color here
	 * @param color color
	 * @return TRUE, then check for the next color, otherwise continue 
	 */
    bool pruneFF(int color);
    
	/**
	 * @brief calculate an upper bound
	 * @return upper bound
	 */
    int calcUB();
	/**
	 * @brief calculate an upper bound by the naive approach
	 * @return upper bound
	 */
    int naiveUB();

	/**
	 * @brief give a node a specific color
	 * @param v vertex
	 * @param color color
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool swapNodeToColor(Vertex v, int color);
	/**
	 * @brief find a color class with the min. amount of vertices and one with the max. amount of vertices in it.
	 * @param cMin minimum
	 * @param cMax maximum
	 * @return pair of indices of the color class
	 */
    std::pair<int, int> findMinMaxColorClass(int cMin, int cMax);
};

#endif
