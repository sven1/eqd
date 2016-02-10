/**
 * @file Heuristic.cpp
 * @brief Source File
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#include "Heuristic.hpp"

bool Heuristic::findIndepCliques(Graph &g, PropertyMap &pm){
  return true;
}

double Heuristic::rn(const double LB, const double UB){
  return ((double) rand() / RAND_MAX) * (UB - LB) + LB;
}

bool Heuristic::constructRandomGraph(Graph &g, int n, double p, int INIT_SRAND){
  double rNumber;
  bool error;
  vertexIter vIt1, vIt2, tmpIt;

  srand(INIT_SRAND);
  
  for(int i=1; i <= n; i++){
    add_vertex(g);
  }

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    tmpIt = vIt1 + 1;

    while(tmpIt != vIt2){
      rNumber = rn(0,1);

      if(rNumber <= p){
        error = add_edge(*vIt1, *tmpIt, g).second;

        if(!error){
          std::cout << "edge already exists" << std::endl;

          return false;
        }
      }

      tmpIt++;
    }
  }

  return true;
}
