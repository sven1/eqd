#include "boost.hpp"
#include "structs.hpp"

class Heuristic{
  private:
    static double rn(const double LB, const double UB);
  
  public:
    static bool findIndepCliques(Graph &g, PropertyMap &pm);

    static bool constructRandomGraph(Graph &g, int n, double p, int INIT_SRAND);
};
