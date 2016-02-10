/**
 * @file Coloring.hpp
 * @brief coloring a graph
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#ifndef COLORING_HPP
#define COLORING_HPP

#include "structs.hpp"
#include "Input.hpp"
#include "Heuristic.hpp"

/**
 * @brief color a graph
 * @class Coloring
 */
class Coloring {
  protected:
	  /**
	   * @brief graph
	   */
    Graph g;
	/**
	 * @brief bounds
	 */
    Bounds b;
	/**
	 * @brief property map
	 */
    PropertyMap pm;
	/**
	 * @brief time
	 */
    Time t;
	/**
	 * @brief counts
	 */
    Count c;
	/**
	 * @brief vector containing all vertices in the start clique
	 */
    std::vector<Vertex> startClique;
	/**
	 * @brief 2D vector which contains all cliques
	 */
    std::vector<std::vector<Vertex> > indClq;
	/**
	 * @brief current informations
	 */
    Current curr;
	/**
	 * @brief input parameters
	 */
    Parameters parm;
	/**
	 * @brief backtrack informations
	 */
    Backtracking bt;
	/**
	 * @brief color classes
	 */
    Colors cc;
	/**
	 * @brief clique informations
	 */
    Cliques cl;

	/**
	 * @brief initialize property maps, i.e. connect it with the specific graph
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initPropMap();
	/**
	 * @brief initialize neighbours for everey node
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initNeighbours();
	/**
	 * @brief resize vectors (upper bound)
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initResize();
	/**
	 * @brief initialize cliques
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initCliques();
	/**
	 * @brief initialize some useful counting informations
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initCounts();

	/**
	 * @brief main initialisation funct. that calls all other init. functions
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool initVar();

	/**
	 * @brief initialize time values
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool setTime();
	/**
	 * @brief set upper and lower bounds for the chromatic (equitable) number
	 * @param LB lower bound
	 * @param UB upper bound
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool setBounds(int LB, int UB);
	/**
	 * @brief set variables for the backup graphs
	 * @param nNodesSameCl how much nodes are colored
	 * @param nInCliqueBp how much nodes are in the cliques
	 * @return TRUE, if it works, otherwise FALSE
	 */
	bool setBackupVar(int nNodesSameCl, int nInCliqueBp);
	/**
	 * @brief set clique informations
	 * @param nodesInClique how much nodes are in the cliques
	 * @param nCliques how much cliques
	 * @param newClique TRUE, if it is a new clique, otherwise FALSE
	 * @return 
	 */
    bool setClique(long nodesInClique, long nCliques, bool newClique);
	/**
	 * @brief set current informations
	 * @param c color
	 * @param r rank
	 * @param node vertex
	 * @param uncoloredVertices uncolored vertices
	 * @param T how much max. color classes
	 * @param M index of the max. color class
	 * @param nColors actual used colors
	 * @param createNewGraphs TRUE, if new graphs should be created, otherwise FALSE
	 * @return 
	 */
    bool setCurr(int c, int r, Vertex node, int uncoloredVertices, int T, int M, long nColors, bool createNewGraphs);
	/**
	 * @brief set backtracking information
	 * @param status TRUE, if the algorithm backtracks, otherwise FALSE
	 * @param node to which node
	 * @param toRank to which rank
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool setBacktracking(bool status, Vertex node, int toRank);
	/**
	 * @brief set some starting parametery given by the input
	 * @param n number of vertices
	 * @param p density
	 * @param npr pruning rule
	 * @param tl time limit
	 * @param th treshold
	 * @param res ressource file
	 * @param nrg number of random graphs
	 * @param variant variant 
	 * @param pNewCl percentage to search for a better clique
	 * @return 
	 */
    bool setParm(long n = 40, double p = 0.5, long npr = 1, long tl = 3600, long th = 2, std::string res = "res/queen7_7.col", long nrg = 200, char variant = 'R', double pNewCl = 0.5);
	/**
	 * @brief set some informations for color classes
	 * @param T amount of max. color classes
	 * @param M index of max. color class
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool setTandM(int T, int M);

	/**
	 * @brief find max clique with vertex v in it
	 * @param clq vector containing the clique
	 * @param v vertex contained in the clique
	 * @param uncolored TRUE, if only uncolored vertices are allowed, otherwise FALSE
	 * @param inNoOtherClique TRUE, if only vertices in no other cliques are allowd, otherwise FALSE
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool findMaxClique(std::vector<Vertex> &clq, Vertex v, bool uncolored, bool inNoOtherClique);
	/**
	 * @brief calculate the max. saturation degree
	 * @return the max. saturation degree
	 */
    int findMaxSatDeg();
	/**
	 * @brief find vertices with the given saturation degree
	 * @param satDeg saturation degree
	 * @return vector with vertices
	 */
    std::vector<Vertex> findVertexSatDeg(int satDeg);

	/**
	 * @brief calculate the time
	 */
    void calcTime();
	/**
	 * @brief compare the degree for the given vertices
	 * @param v vertex 1
	 * @param w vertex 2
	 * @return TRUE, if the degree of vertex 1 is lower than the degree of vertex 2, otherwise FALSE
	 */
    bool compareDegree(Vertex v, Vertex w);

	/**
	 * @brief put node with specific parameters in the vector
	 * @param v vertex
	 * @param tmp vector with nodes
	 * @param uncolored TRUE, if the vertex should be uncolored, otherwise FALSE
	 * @param inOtherClique TRUE, if the vertex should not be in a clique, otherwise FALSE
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool putNodeWithParm(Vertex v, std::vector<Vertex> &tmp, bool uncolored, bool inOtherClique);
	/**
	 * @brief change the graph property (clique property) for the vertices in the clique
	 * @param clq clique
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool putInClique(std::vector<Vertex> &clq);

	/**
	 * @brief print neighbours
	 * @param v vertex
	 */
    void printNeighbours(const Vertex &v) const;
	/**
	 * @brief print vertex infos
	 * @param v vertex
	 */
    void printVertexInfo(const Vertex &v) const;
	/**
	 * @brief print forbidden colors for vertex
	 * @param v vertex
	 */
    void printFBC(const Vertex &v) const;
	/**
	 * @brief print counting informations
	 */
    void printCounts() const;

	/**
	 * @brief check if its a coloring for vertex and his neighbours 
	 * @param v vertex
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool checkColoring(const Vertex &v) const;
	/**
	 * @brief check if its a clique
	 * @param clq clique
	 * @return TRUE, if it works, otherwise FALSE
	 */
    bool checkClique(const std::vector<Vertex> &clq) const;
	/**
	 * @brief check if we need to backtrack from vertex v
	 * @param v vertex
	 * @return TRUE, if we should backtrack, otherwise FALSE
	 */
    bool checkForBacktracking(Vertex v);
	/**
	 * @brief check propertys for vertex
	 * @param v vertex
	 * @param uncolored 
	 * @param inNoOtherClique
	 * @return 
	 */
    bool checkNodeParm(Vertex v, bool uncolored, bool inNoOtherClique);
	/**
	 * @brief check for overtime
	 * @return TRUE, if its a timeout, otherwise FALSE
	 */
    bool checkOvertime();

	/**
	 * @brief update infos T and M for color classes
	 * @param lastColor last used color
	 * @param inc TRUE, if its incremented, otherwise FALSE
	 * @return 
	 */
    bool updateTandM(int lastColor, bool inc);
	/**
	 * @brief helper function
	 * @param M
	 * @return 
	 */
    int helpUpdateTandM(int M);

	/**
	 * @brief increment saturation degree
	 * @param v vertex 
	 * @param color color
	 * @return 
	 */
    bool incSatDeg(Vertex v, int color);
	/**
	 * @brief decrease saturation degree
	 * @param v vertex 
	 * @param color color
	 * @return 
	 */
    bool decSatDeg(Vertex v, int color);

	/**
	 * @brief calculate the backtrack rank
	 * @param v vertex
	 * @return 
	 */
    int backtrackToRank(Vertex v);
	/**
	 * @brief backtracking for leafs (search tree)
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
	 * @brief algorithm (heuristic) to find an upper bound
	 * @param g graph
	 * @return upper bound
	 */
    long naiveUB(const Graph &g);
	/**
	 * @brief algorithm (heuristic) to find a lower bound
	 * @param g graph
	 * @return lower bound
	 */
    long eqlLB(const Graph &g);

	/**
	 * @brief vertex selection strategy
	 * @return selected vertex
	 */
    Vertex passVSS();
    int helpPassVSS(Vertex v, int maxSatDeg);

	/**
	 * @brief main function for the normal dsatur algorithm
	 * @return 
	 */
    bool node();

	/**
	 * @brief print some vertex infos
	 */
    void printVertexInfo() const;
	/**
	 * @brief print adjacency list
	 */
    void printAdjMatrix() const;
	/**
	 * @brief print infos for the cliques
	 */
    void printCliqueInfo() const;
	/**
	 * @brief print forbidden color classes
	 */
    void printFBC() const;
	/**
	 * @brief print the given clique
	 * @param clq clique
	 */
    void printClique(const std::vector<Vertex> &clq) const;
	/**
	 * @brief print independent cliques
	 * @param indClq independent cliques (2D vector)
	 */
    void printIndepCliques(const std::vector<std::vector<Vertex> > &indClq) const;
	/**
	 * @brief print graph header
	 */
    void printGraphHeaders() const;
	/**
	 * @brief print ALL informations
	 */
    void printAll() const;
	/**
	 * @brief print upper and lower bound
	 */
    void printBounds() const;
	/**
	 * @brief print informations for the current step
	 */
    void printCurrent() const;
	/**
	 * @brief print specific color class
	 * @param i color = i+1
	 */
    void printColorClass(int i) const;
	/**
	 * @brief print all color classes
	 */
    void printColorClasses() const;
	/**
	 * @brief print formatted infos for the random graph studies
	 */
    void printStudyRandomGraphs() const;
    
	/**
	 * @brief find a max. clique w.r.t. to the properties given by the args
	 * @param clq clique
	 * @param uncolored TRUE, if they should be uncolored, otherwise FALSE
	 * @param inNoOtherClique TRUE, if they should be in no clique, otherwise FALSE
	 * @return 
	 */
    bool findMaxClique(std::vector<Vertex> &clq, bool uncolored, bool inNoOtherClique);
	/**
	 * @brief find independent cliques w.r.t. to the properties given by the args
	 * @param indClq 2D vector containing the cliques
	 * @param uncolored TRUE, if they should be uncolored, otherwise FALSE
	 * @param inNoOtherClique TRUE, if they should be in no clique, otherwise FALSE
	 * @return 
	 */
    bool findIndepCliques(std::vector<std::vector<Vertex> > &indClq, bool uncolored, bool inNoOtherClique);
	/**
	 * @brief remove a vertex from the 2D vector containg the cliques
	 * @param v vertex
	 * @param indClq 2D vector containing the cliques
	 * @return 
	 */
    bool removeVinIndClq(Vertex v, std::vector<std::vector<Vertex> > &indClq);
	/**
	 * @brief function that update the indep. cliques (searching for better indep. cliques) if desired
	 * @param indClq cliques
	 */
    void useIndepClique(std::vector< std::vector< Vertex > > &indClq);

	/**
	 * @brief check ifs a graph coloring
	 * @return 
	 */
    bool checkColoring() const;

	/**
	 * @brief color the clique with start color, start color + 1, start color + 2, ... and so on
	 * @param clq clique
	 * @param startColor start color
	 * @return 
	 */
    bool colorClique(std::vector<Vertex> &clq, int startColor);
	/**
	 * @brief check if the cliques are "independent"
	 * @param indClq 2D vector containing the cliques
	 * @return 
	 */
    bool checkIndepClique(std::vector<std::vector<Vertex> > &indClq);

	/**
	 * @brief calculate an upper bound
	 * @return upper bound
	 */
    int calcUB();
	/**
	 * @brief calculate a lower bound
	 * @return lower bound
	 */
    int calcLB();

	/**
	 * @brief color greedy a graph
	 * @param g graph
	 * @return 
	 */
    bool greedyColoring(Graph &g);
	/**
	 * @brief calculate the smallest possible color for vertex v
	 * @param v vertex
	 * @return 
	 */
    int smallestPosColor(Vertex v) const;

	/**
	 * @brief color a vertex
	 * @param v vertex
	 * @param color color
	 * @return 
	 */
    bool colorVertex(Vertex v, int color);
	/**
	 * @brief uncolor a vertex
	 * @param v vertex
	 * @return 
	 */
    bool uncolorVertex(Vertex v);
    
	/**
	 * @brief add a forbidden color for a vertex
	 * @param v vertex
	 * @param color color
	 * @return 
	 */
    bool addFBC(Vertex v, int color);
	/**
	 * @brief remove a forbidden color for a vertex
	 * @param v vertex 
	 * @param color color
	 * @return 
	 */
    bool removeFBC(Vertex v, int color);
    
	/**
	 * @brief increment color class
	 * @param color color
	 * @return 
	 */
    bool incColorClass(int color);
	/**
	 * @brief decrease color class
	 * @param color color
	 * @return 
	 */
    bool decColorClass(int color);

	/**
	 * @brief upper gauss for x
	 * @param x value
	 * @return 
	 */
    int upperGauss(double x);
};

#endif
