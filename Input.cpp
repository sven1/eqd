#include "Input.hpp"
#include <cstdlib>
#include <string>
#include <fstream>

bool Input::readInputArgs(int argc, char** argv, Parameters &p){
  if(argc < 2){
    std::cout << "no mode" << std::endl;

    return false;
  }else if(argc == 2){
	std::cout << "no variant" << std::endl;

	return false;
  }else{
    if(std::string(argv[1]) == "N"){
      if(argc != 9){
        std::cout << "Missing some arguments in (N)ormal mode" << std::endl;

        return false;
      }else{
        p.eqDsatur = *argv[2];
        p.nPruningRule = atol(argv[3]);
        p.timeLimit = atol(argv[4]);
        p.threshold = atol(argv[5]);
        p.ressource = argv[6];
        p.pNewCl = atof(argv[7]);
		p.wStartCl = atoi(argv[8]); 
        p.variant = 'N';
        p.n = 0;
        p.p = 0;
	    p.INIT_SRAND = 0;
      }
    }else if(std::string(argv[1]) == "R"){
      if(argc != 12){
        std::cout << "Missing some arguments in (R)andom mode" << std::endl;

        return false;
      }else{
        p.eqDsatur = *argv[2];
        p.nPruningRule = atol(argv[3]);
        p.timeLimit = atol(argv[4]);
        p.threshold = atol(argv[5]);
        p.nRandomGraphs = atol(argv[6]);
        p.n = atol(argv[7]);
        p.p = atof(argv[8]);
	    p.INIT_SRAND = atoi(argv[9]);
        p.pNewCl = atof(argv[10]);
		p.wStartCl = atoi(argv[11]); 
        p.variant = 'R';
      }
    }else{
      std::cout << "Not a valid Variant" << std::endl;

      return false;
    }
  }

  return true;
}

bool Input::readGraph(const std::string &filename, Graph &g){
  std::string sType;
  std::size_t index = filename.find(".");

  if(index != std::string::npos){
    sType = filename.substr(index+1);
        
    if(readGraph(filename, sType, g)){
      return true;
    }else{
      return false;
    }
  }else{
    std::cout << "Not a expected *.* file format!" << std::endl;

    return false;
  }
}

bool Input::readClique(Graph &g){
  return true;
}

bool Input::readGraph(const std::string &filename, const std::string &sType, Graph &g){
  if(sType.compare("col") == 0){
    if(readGraphCol(filename, g)){
      return true;
    }else{
      std::cout << "Error: Can not read graph." << std::endl;

      return false;
    }
  }else{
    std::cout << "No suitable parser avalaible." << std::endl;

    return false;
  }
}

bool Input::readGraphCol(const std::string &filename, Graph &g){
  std::ifstream file(filename.c_str(), std::ios::in);

  char control = 0;
  std::string tmp;
  int nVertices, nEdges, fromV, toV;
  bool error;

  if(file.good()){
    while(file >> control){
      switch(control){
        case 'c':{
          getline(file, tmp);
        } break;

        case 'p':{
          file >> tmp;

          if(tmp.compare("edge") != 0){
            std::cout << "Not \"edge\" format." << std::endl;

            return false;
          }

          file >> nVertices;
          file >> nEdges;

        } break;

        case 'e':{
          file >> fromV;
          file >> toV;

          if(fromV < toV){
            error = add_edge(fromV, toV, g).second;
          }

          if(!error){
            std::cout << "edge already exists" << std::endl;

            return false;
          }
        } break;

        default:{
          std::cout << "Unknown control character \"" << control << "\" in DIMACS format." << std::endl;

          file.close();
          
          return false;
        }break;
      }
    }
  }else{
    std::cout << "Can't open file." << std::endl;

    file.clear();
    file.close();
    
    return false;
  }

  remove_vertex(0, g);

  file.close();
  
  return true;
}
