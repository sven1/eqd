/**
 * @file Input.hpp
 * @brief read a graph
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#include "boost.hpp"
#include "structs.hpp"
#include <string>

/**
 * \class Input
 * \brief class to read a graph or input args
 */

class Input{
  private:
	  /**
	   * @brief read the graph from ressource file
	   * @param filename ressource file
	   * @param sType type (for instance .col)
	   * @param g graph
	   * @return 
	   */
    static bool readGraph(const std::string &filename, const std::string &sType , Graph &g);
	/**
	 * @brief read graph in .col format
	 * @param filename ressource file
	 * @param g graph
	 * @return 
	 */
    static bool readGraphCol(const std::string &filename, Graph &g);

  public:
	/**
	 * @brief read the graph
	 * @param filename ressource file
	 * @param g graph
	 * @return 
	 */
    static bool readGraph(const std::string &filename, Graph &g);
	/**
	 * @brief read clique from graph
	 * @param g grapj
	 * @return 
	 */
    static bool readClique(Graph &g);
	/**
	 * @brief read and apply the input arguments
	 * @param argc number of arguments
	 * @param argv string of the arguments
	 * @param p struct where to store the input args
	 * @return 
	 */
    static bool readInputArgs(int argc, char** argv, Parameters &p);
};
