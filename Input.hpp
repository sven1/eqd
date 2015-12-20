#include "boost.hpp"
#include "structs.hpp"
#include <string>

class Input{
  private:
    static bool readGraph(const std::string &filename, const std::string &sType , Graph &g);
    static bool readGraphCol(const std::string &filename, Graph &g);

  public:
    static bool readGraph(const std::string &filename, Graph &g);
    static bool readClique(Graph &g);
    static bool readInputArgs(int argc, char** argv, Parameters &p);
};
