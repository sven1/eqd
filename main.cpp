#include "Coloring.hpp"
#include "EqColoring.hpp"

void usage(){
  std::cout << "Usage: ./eqDsatur [N | R] [N | C] nPruningRule timeLimit threshold [FILE ClPerc | nRandomGraphs n p INIT_SRAND ClPerc]" << std::endl;
}

int main(int argc, char** argv){
  Parameters p;

  if(!Input::readInputArgs(argc, argv, p)){
    usage();

    exit(-1);
  }

  EqColoring e(p);

  return 0;
}
