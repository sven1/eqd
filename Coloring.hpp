/**
 * \file Coloring.hpp
 * \author Sven Förster
 * \date 29.12.2015
 *
 * \version 1.0.0
 * cleaning code
 */

#ifndef COLORING_HPP
#define COLORING_HPP

#include "structs.hpp"
#include "Input.hpp"
#include "Heuristic.hpp"

/**
 * \class Coloring
 * \brief class for coloring a graph instance
 */

class Coloring {
  protected:
    Graph g; /**< the given graph */
    Bounds b; /**< struct for the upper and lower bound */
    PropertyMap pm; /**< struct with the property maps for accessing the graph propertys */
    Time t; /**< struct with all time informations and timers */
    Count c; /**< struct with counting informations, amount of uncolored vertices, used colors, etc. */
    std::vector<Vertex> startClique; /**< start clique stored in a vector */
    std::vector<std::vector<Vertex> > indClq; /**< independent cliques stored in a 2D vector */
    Current curr; /**< relevant informations for the current step in the search tree */
    Parameters parm; /**< relevant informations for the graph. number of nodes, density, etc. */
    Backtracking bt; /**< relevant informations for backtracking. status, rank of the node where to backtrack */
    Colors cc; /**< color classes are stored in this variable */
    Cliques cl; /**< useful informations for the cliques. Number of cliques, nodes in cliques, etc. */

	/**
	 * \brief initialize the property maps for the given graph
	 * \return TRUE, if no error occurs
	 */

    bool initPropMap();

	/**
	 * \brief initialize the neighbours in a vector for every node in the given graph as a property
	 * \return TRUE, if no error occurs
	 */

    bool initNeighbours();

	/**
	 * \brief resize needed vectors to an appropriate size (upper bound)
	 * \return TRUE, if no error occurs
	 */

    bool initResize();

	/**
	 * \brief initialize independent cliques and one start clique, if its wanted
	 * \return TRUE, if no error occurs
	 */

    bool initCliques();

	/**
	 * \brief initialize a struct with counting variables with default values
	 * \return TRUE, if no error occurs
	 */

    bool initCounts();
	/**
	 * \brief starts all other initialize methods
	 * \return TRUE, if no error occurs
	 */

    bool initVar();

	/**
	 * \brief initialize the timers with default value
	 * \return TRUE, if no error occurs
	 */

    bool setTime();

	/**
	 * \brief set upper and lower bounds for the equitable chromatic number
	 * \param[in] LB lower bound
	 * \param[in] UB upper bound
	 * \return TRUE, if no error occurs
	 */

    bool setBounds(int LB, int UB);

	/**
	 * \brief initialize the variables (when do i have to choose for new indep. cliques)
	 * \param[in] nNodesSameCl how many nodes are colored in the actual indep. cliques
	 * \param[in] nInCliqueBp how many nodes are in the indep. cliques
	 * \return TRUE, if no error occurs
	 */

	bool setBackupVar(int nNodesSameCl, int nInCliqueBp);

	/**
	 * \brief set informations for the cliques
	 * \param[in] nodesInClique how many nodes are in the indep cliques + start clique
	 * \param[in] nCliques how many cliques do we have
	 * \param[in] newClique TRUE, if we found a better clique
	 * \return TRUE, if no error occurs
	 */

    bool setClique(long nodesInClique, long nCliques, bool newClique);

	/**
	 * \brief set informations for the current step in the coloring branch
	 * \param[in] c actual color
	 * \param[in] r actual rank
	 * \param[in] node actual vertex
	 * \param[in] uncoloredVertices amount of uncolored vertices
	 * \param[in] T how many color classes do we have with the highest amount of nodes in it
	 * \param[in] M size of the color class with the highest amount of nodes in it
	 * \param[in] nColors how much colors are actual in use for the coloring
	 * \param[in] createNewGraphs true, if we need new backup graphs
	 * \param[in] rankNC last rank (position) where a critical change was made (new indep. cliques or backtrack)
	 * \return TRUE, if no error occurs
	 */

    bool setCurr(int c, int r, Vertex node, int uncoloredVertices, int T, int M, long nColors, bool createNewGraphs);
    
	/**
	 * \brief set the informations for backtracking
	 * \param[in] status TRUE, if the algorithm backtracks to a specific node
	 * \param[in] node the vertex where the algorithm backtracks to
	 * \param[in] toRank its the rank of the vertex where the algorithm backtracks to
	 * \return TRUE, if no error occurs
	 */

	bool setBacktracking(bool status, Vertex node, int toRank);
    
	/**
	 * \brief set the start informations for the algorithm
	 * \param[in] n number of nodes in the graphs
	 * \param[in] p density for an edge between two nodes in a random graph
	 * \param[in] npr frequency for checking the pruning rule
	 * \param[in] tl time limit for the algorithm to run
	 * \param[in] th treshold for the VSS-PSS strategy
	 * \param[in] res ressource graph (for instance a DIMACS graph)
	 * \param[in] nrg number of random graphs
	 * \param[in] variant R stands for random graphs, N stands for using a given, specific graph
	 * \param[in] pNewCl frequency for checking for new indep cliques
	 * \return TRUE, if no error occurs
	 */
	
	bool setParm(long n = 40, double p = 0.5, long npr = 1, long tl = 3600, long th = 2, std::string res = "res/queen7_7.col", long nrg = 200, char variant = 'R', double pNewCl = 0.5);
    
	/**
	 * \brief set T and M
	 * \return TRUE, if no error occurs
	 */
	
	bool setTandM(int T, int M);

	/**
	 * \brief search for a max clique containing vertex v
	 * \param[out] clq store the max clique in this variable
	 * \param[in] v this vertex must be in the clique
	 * \param[in] uncolored if its TRUE vertices could be choosen even if they are already colored, otherwise only uncolored vertices could be choosen
	 * \param[in] inNoOtherClique if its TRUE vertices could only be choosen if they are not in a other clique, otherwise it doesnt matter
	 * \return TRUE, if no error occurs
	 */

    bool findMaxClique(std::vector<Vertex> &clq, Vertex v, bool uncolored, bool inNoOtherClique);

	/**
	 * \brief calculate the maximal saturation degree of the graph (i.e. for everey node and then take the max. value)
	 * \return the maximal saturation degree of the graph
	 */

    int findMaxSatDeg();

	/**
	 * \brief search all vertices with the given saturation degree
	 * \param[in] satDeg the satuation degree
	 * \return a vector with the vertices which have the given saturation degree
	 */

    std::vector<Vertex> findVertexSatDeg(int satDeg);

	/**
	 * \brief calculate the running time
	 */

    void calcTime();
    
	/**
	 * \brief compare the degree of two vectors
	 * \param[in] v vertex one
	 * \param[in] w vertex two
	 * \return TRUE, if the degree from v is lower than the degree from w
	 */

	bool compareDegree(Vertex v, Vertex w);

	/**
	 * \brief put node v in vector tmp w.r.t. the variables uncolor and inOtherClique
	 * \param[in] v the vertex
	 * \param[out] tmp the vector where v should probably to put in
	 */

    bool putNodeWithParm(Vertex v, std::vector<Vertex> &tmp, bool uncolored, bool inOtherClique);
    
	/**
	 * \brief apply the clique to the graph
	 * \param[in] clq vector with nodes
	 */
	
	bool putInClique(std::vector<Vertex> &clq);

	/**
	 * \brief print neighbours from vertex v
	 * \param[in] v vertex v
	 */

    void printNeighbours(const Vertex &v) const;
    
	/**
	 * \brief show some useful information for vertex v
	 * \param[in] v vertex v
	 */
	
	void printVertexInfo(const Vertex &v) const;
    
	/**
	 * \brief print the forbidden colors for the vertex v
	 * \param[in] vertex v
	 */
	
	void printFBC(const Vertex &v) const;
    
	/**
	 * \brief print the counting (amount of new cliques, used colors, nodes, etc.)
	 */
	
	void printCounts() const;

	/**
	 * \brief check if the color for vertex v is allowed
	 * \param[in] v vertex v
	 * \return TRUE, if its a coloring for vertex v
	 */

    bool checkColoring(const Vertex &v) const;
    
	/**
	 * \brief check if clq is a clique
	 * \param[in] clq clique
	 * \return TRUE, if its a clique
	 */
	
	bool checkClique(const std::vector<Vertex> &clq) const;
    
	/**
	 * \brief calculate the rank where v have to batrack to
	 * \param[in] v vertex v
	 * \return TRUE, if the algorithm should backtrack
	 */
	
	bool checkForBacktracking(Vertex v);
    
	/**
	 * \brief check for vertex v the variables uncolord and inNoOtherClique
	 * \param[in] v vertex v
	 * \param[in] uncolored TRUE, if the vertex should be uncolored
	 * \param[in] inNoOtherClique TRUE, if the vertex should be in no clique
	 * \return TRUE, if v is acceptable w.r.t. the variables uncolored and inNoOtherClique
	 */
	
	bool checkNodeParm(Vertex v, bool uncolored, bool inNoOtherClique);
    
	/**
	 * \brief check if the run time is not over the time limit
	 */
	
	bool checkOvertime();

	/**
	 * \brief update the variables T and M (informations for color classes)
	 * \return TRUE, if no error occurs
	 */

    bool updateTandM(int lastColor, bool inc);
    
	/**
	 * \brief its a helper function for updateTandM
	 */
	
	int helpUpdateTandM(int M);

	/**
	 * \brief (maybe) increment the saturation degree for vertex v
	 * \return TRUE, if no error occurs
	 */

    bool incSatDeg(Vertex v, int color);

	/**
	 * \brief (maybe) decrease the saturation degree for vertex v
	 * \return TRUE, if no error occurs
	 */
    
	bool decSatDeg(Vertex v, int color);

	/**
	 * \brief backtrack to a specific rank for the given vertex v
	 * \return the rank where the algorithm have to backtrack
	 */

    int backtrackToRank(Vertex v);
    
	/**
	 * \brief backtracking for the leafs of the search tree
	 */

	void newUBBacktracking();

    Colors cls;

  public:
    Coloring();
    Coloring(const Parameters &parm);
    Coloring(const Graph &g);
    Coloring(const Graph &g, const Parameters &parm);
  
    ~Coloring();

	/**
	 * \brief algorithm to search for an upper bound for the equitable chromatic number
	 * \param[in] g graph
	 * \return an upper bound for the given graph
	 */

    long naiveUB(const Graph &g);
    
	/**
	 * \brief algorithm to search for a lower bound for the equitable chromatic number
	 * \param[in] g graph
	 * \return a lower bound for the given graph
	 */
	
	long eqlLB(const Graph &g);

	/**
	 * \brief algorithm to select the next vertex to colorize
	 * \return the next vertex
	 */

    Vertex passVSS();
    
	/**
	 * \brief helper function for the passVSS method
	 */
	
	int helpPassVSS(Vertex v, int maxSatDeg);

	/**
	 * \brief main function to traverse the search tree to find an equitable chromatic number with no cliques
	 * \return TRUE, if its terminated
	 */

    bool node();

	/**
	 * \brief print some information for all vertices
	 */

    void printVertexInfo() const;
    
	/**
	 * \brief print the adjacency matrix for the given graph
	 */
	
	void printAdjMatrix() const;
    
	/**
	 * \brief print some information for the cliques
	 */
	
	void printCliqueInfo() const;
    
	/**
	 * \brief print the forbidden colors, i.e. the colors that are forbidden to use for a vertex
	 */

	void printFBC() const;
    
	/**
	 * \brief print all nodes for the clique
	 * \param[in] clq clique
	 */
	
	void printClique(const std::vector<Vertex> &clq) const;
    
	/**
	 * \brief print the independent cliques
	 * \param[in] indClq the independent cliques
	 */
	
	void printIndepCliques(const std::vector<std::vector<Vertex> > &indClq) const;
    
	/**
	 * \brief print some relevant information for the graph
	 */
	
	void printGraphHeaders() const;
    
	/**
	 * \brief print ALL informations for the graph, countings, timers, etc.
	 */
	
	void printAll() const;
    
	/**
	 * \brief print the upper and lower bound
	 */
	
	void printBounds() const;
    
	/**
	 * \brief print some useful informations for the current step in the search tree
	 */
	
	void printCurrent() const;
    
	/**
	 * \brief print the color class for color (i+1)
	 */

	void printColorClass(int i) const;
    
	/**
	 * \brief print ALL color classes
	 */
	
	void printColorClasses() const;
    
	/**
	 * \brief format the results and afterwards they are printed on the screen
	 */
	
	void printStudyRandomGraphs() const;
    
	/**
	 * \brief find a max clique
	 */

    bool findMaxClique(std::vector<Vertex> &clq, bool uncolored, bool inNoOtherClique);
    
	/**
	 * \brief find independent cliques
	 */
	
	bool findIndepCliques(std::vector<std::vector<Vertex> > &indClq, bool uncolored, bool inNoOtherClique);
    
	/**
	 * \brief remove a specific vertex from the independent cliques
	 */
	
	bool removeVinIndClq(Vertex &v, std::vector<std::vector<Vertex> > &indClq);
    
	/**
	 * \brief use a specific, hard programmed independent cliques
	 */
	
	void useIndepClique(std::vector< std::vector< Vertex > > &indClq);

	/**
	 * \brief check if the actual coloring is allowed
	 */

    bool checkColoring() const;

	/**
	 * \brief color the clique, starting with the a specific color
	 * \param[in] clq clique
	 * \param[in] startColor start with this color
	 * \return TRUE, if no error occurs
	 */

    bool colorClique(std::vector<Vertex> &clq, int startColor);
    
	/**
	 * \brief check if the clique is independent
	 * \return TRUE, if its independent
	 */
	
	bool checkIndepClique(std::vector<std::vector<Vertex> > &indClq);

	/**
	 * \brief calculate the init. upper bound for the graph
	 * \return the upper bound
	 */

    int calcUB();
    
	/**
	 * \brief calculate the init. lower bound for the graph
	 * \return the lower bound
	 */
	
	int calcLB();

	/**
	 * \brief use a greedy algorithm to color the graph
	 */

    bool greedyColoring(Graph &g);
    
	/**
	 * \brief calculate the smallest possible color for vertex v
	 */
	
	int smallestPosColor(Vertex v) const;

	/**
	 * \brief color vertex v with a color
	 */

    bool colorVertex(Vertex v, int color);
    
	/**
	 * \brief uncolor vertex v
	 */
	
	bool uncolorVertex(Vertex v);
    
	/**
	 * \brief update the forbidden colors map for vertex v
	 */

    bool addFBC(Vertex v, int color);
    
	/**
	 * \brief update the forbidden colors map for vertex v
	 */
	
	bool removeFBC(Vertex v, int color);
    
	/**
	 * \brief update the color class
	 */

    bool incColorClass(int color);
    
	/**
	 * \brief update the color class
	 */

	bool decColorClass(int color);

	/**
	 * \brief function to calculate the upper gauss of a number
	 */

    int upperGauss(double x);
};

#endif
