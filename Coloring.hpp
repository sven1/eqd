#ifndef COLORING_HPP
#define COLORING_HPP

#include "structs.hpp"
#include "Input.hpp"
#include "Heuristic.hpp"

class Coloring {
  protected:
    Graph g;
    Bounds b;
    PropertyMap pm;
    Time t;
    Count c;
    std::vector<Vertex> startClique;
    std::vector<std::vector<Vertex> > indClq;
    Current curr;
    Parameters parm;
    Backtracking bt;
    Colors cc;
    Cliques cl;

    bool initPropMap();
    bool initNeighbours();
    bool initResize();
    bool initCliques();
    bool initCounts();

    bool initVar();

    bool setTime();
    bool setBounds(int LB, int UB);
	bool setBackupVar(int nNodesSameCl, int nInCliqueBp);
    bool setClique(long nodesInClique, long nCliques, bool newClique);
    bool setCurr(int c, int r, Vertex node, int uncoloredVertices, int T, int M, long nColors, bool createNewGraphs);
    bool setBacktracking(bool status, Vertex node, int toRank);
    bool setParm(long n = 40, double p = 0.5, long npr = 1, long tl = 3600, long th = 2, std::string res = "res/queen7_7.col", long nrg = 200, char variant = 'R', double pNewCl = 0.5);
    bool setTandM(int T, int M);

    bool findMaxClique(std::vector<Vertex> &clq, Vertex v, bool uncolored, bool inNoOtherClique);
    int findMaxSatDeg();
    std::vector<Vertex> findVertexSatDeg(int satDeg);

    void calcTime();
    bool compareDegree(Vertex v, Vertex w);

    bool putNodeWithParm(Vertex v, std::vector<Vertex> &tmp, bool uncolored, bool inOtherClique);
    bool putInClique(std::vector<Vertex> &clq);

    void printNeighbours(const Vertex &v) const;
    void printVertexInfo(const Vertex &v) const;
    void printFBC(const Vertex &v) const;
    void printCounts() const;

    bool checkColoring(const Vertex &v) const;
    bool checkClique(const std::vector<Vertex> &clq) const;
    bool checkForBacktracking(Vertex v);
    bool checkNodeParm(Vertex v, bool uncolored, bool inNoOtherClique);
    bool checkOvertime();

    bool updateTandM(int lastColor, bool inc);
    int helpUpdateTandM(int M);

    bool incSatDeg(Vertex v, int color);
    bool decSatDeg(Vertex v, int color);

    int backtrackToRank(Vertex v);
    void newUBBacktracking();

    Colors cls;

  public:
    Coloring();
    Coloring(const Parameters &parm);
    Coloring(const Graph &g);
    Coloring(const Graph &g, const Parameters &parm);
  
    ~Coloring();

    long naiveUB(const Graph &g);
    long eqlLB(const Graph &g);

    Vertex passVSS();
    int helpPassVSS(Vertex v, int maxSatDeg);

    bool node();

    void printVertexInfo() const;
    void printAdjMatrix() const;
    void printCliqueInfo() const;
    void printFBC() const;
    void printClique(const std::vector<Vertex> &clq) const;
    void printIndepCliques(const std::vector<std::vector<Vertex> > &indClq) const;
    void printGraphHeaders() const;
    void printAll() const;
    void printBounds() const;
    void printCurrent() const;
    void printColorClass(int i) const;
    void printColorClasses() const;
    void printStudyRandomGraphs() const;
    
    bool findMaxClique(std::vector<Vertex> &clq, bool uncolored, bool inNoOtherClique);
    bool findIndepCliques(std::vector<std::vector<Vertex> > &indClq, bool uncolored, bool inNoOtherClique);
    bool removeVinIndClq(Vertex v, std::vector<std::vector<Vertex> > &indClq);
    void useIndepClique(std::vector< std::vector< Vertex > > &indClq);

    bool checkColoring() const;

    bool colorClique(std::vector<Vertex> &clq, int startColor);
    bool checkIndepClique(std::vector<std::vector<Vertex> > &indClq);

    int calcUB();
    int calcLB();

    bool greedyColoring(Graph &g);
    int smallestPosColor(Vertex v) const;

    bool colorVertex(Vertex v, int color);
    bool uncolorVertex(Vertex v);
    
    bool addFBC(Vertex v, int color);
    bool removeFBC(Vertex v, int color);
    
    bool incColorClass(int color);
    bool decColorClass(int color);

    int upperGauss(double x);
};

#endif
