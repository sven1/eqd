/**
 * @file Heuristic.hpp
 * @brief provides some useful methods to integrate with graphs
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#include "boost.hpp"
#include "structs.hpp"

/**
 * @class Heuristic
 * @brief provides some useful methods to integrate with graphs
 */

class Heuristic{
  private:
	  /**
	   * @brief create a random number between LB and UB
	   * @param LB lower bound
	   * @param UB upper bound
	   * @return random number between UB and LB
	   */
    static double rn(const double LB, const double UB);
  
  public:
	/**
	 * @brief find independent cliques in the graph g and store it 
	 * @param g graph 
	 * @param pm property map to access the graph
	 * @return TRUE, if it works, otherwise FALSE
	 */
    static bool findIndepCliques(Graph &g, PropertyMap &pm);

	/**
	 * @brief construct a random graph
	 * @param g graph 
	 * @param n number of vertices
	 * @param p density to construct an edge (between 0 and 1)
	 * @param INIT_SRAND seed value for the random graph, or a number to identify a graph
	 * @return TRUE, if it works, otherwise FALSE
	 */
    static bool constructRandomGraph(Graph &g, int n, double p, int INIT_SRAND);
};
