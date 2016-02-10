/**
 * @file EqColoring.cpp
 * @brief Source File
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());
  b.startUB = b.UB;
  b.startLB = b.LB;


  // same bounds -> quit
  if(b.LB != b.UB){
	// start normal dsatur
	if(parm.eqDsatur == 'N'){
		dsatur();
	// start dsatur with clique
	}else if(parm.eqDsatur == 'C'){
		//init property maps
		initPrevGraphsFF();

		dsaturClique();
	}
  }
  
  // calculate the exact total time
  calcTime();

  // print relevant data
  printStudyRandomGraphs();
}


/**
 * init property maps
 */

bool EqColoring::initPrevGraphsFF(){
  pmf.c = get(edge_capacity, gf.g);
  pmf.re = get(edge_reverse, gf.g);
  pmf.rc = get(edge_residual_capacity, gf.g);
  pmf.rf = get(ref_vertex_t(), gf.g);

  return true;
}

EqColoring::EqColoring(const Graph &g) : Coloring(g){

}

EqColoring::EqColoring(const Graph &g, const Parameters &parm) : Coloring(g, parm){

}

EqColoring::~EqColoring(){

}

/**
 * main function to traverse the searching tree (normal dsatur)
 */

bool EqColoring::node(){
	// check if timeout is already reacher
  if(checkOvertime()){
    return true;
  }

  // check if all nodes are already colored
  if(curr.uncoloredVertices == 0){
	// check if the current coloring uses (weniger) nodes than the upper bound is
    if(curr.nColors < b.UB){
		// check if the actual coloring is a equitable coloring
	    if(checkEqColoring()){
		    b.UB = curr.nColors;

			//std::cout << "NEW UB = " << b.UB << std::endl;
			newUBBacktracking();
	    }
    }

    return true;
  }
  
  //select a node
  Vertex v = passVSS();

  // pruning rule from the paper
  if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
	//loop through all colors
	for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
      //check if the color can be applied to this node (i.e. check if the color is available for this node)
	  if(pm.fbc[v][i - 1] == 0){
		// variable to count visited nodes in the searching tree
        c.visitedNodes++;

		// color vertex with color i
        colorVertex(v, i);

		//std::cout << c.visitedNodes << " Faerbe Knoten = " << v << " mit Farbe = " << i <<  std::endl;

		// do the same again (recursive fnct.)
        node(); 
       
		// check if we reached already the timeout
	      if(t.timeout){
	        return true;
	      }
   
		  //check if we do a backtrack at the moment
        if(bt.status){
			//check if the rank we suppose to backtrack is the current rank of the node
          if(bt.toRank == curr.rank){
			  // then stop backtracking
            bt.status = false;
          }else if(bt.toRank < curr.rank){
			  // continue to backtrack
            uncolorVertex(v);

            return true;
          }
        }

        uncolorVertex(v);
      }
    }
  }

  // check if a backtrack is possible, if yes then do backtrack
  checkForBacktracking(v);

  //std::cout << "BACKTRACK to = " << bt.toRank << std::endl;

  return false;
}

bool EqColoring::pruneFF(){
  t.startFF = std::clock();
  
  for(unsigned int i = std::max(b.LB, curr.nColors); i < b.UB; i++){
	t.startBuildNetwork = std::clock();

	// check pruning rule
	if(curr.M > ceil(parm.n*1./i)){
		continue;
	}
    
	// initialize backup graph
	if(!initBackupGraphs(i)){
		t.timeBuildNetwork += (std::clock() - t.startBuildNetwork) / (double) CLOCKS_PER_SEC;
		continue;
	}

    t.timeBuildNetwork += (std::clock() - t.startBuildNetwork) / (double) CLOCKS_PER_SEC;
	//std::cout << "test for color = " << i << std::endl;

    c.nFF++;

	// check for pruning rule with color i
    if(!pruneFF(i)){
      t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

	  //std::cout << "work for color = " << i << std::endl;
      return false;
    }
  }

  t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

	//std::cout << "no color suitable " << std::endl;

	// cant prune
  return true;
}

bool EqColoring::initNetwork(int colors){
  EdgeFord eF1, eF2;

  //std::cout << "source <-> clique colors" << std::endl;
  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < gf.vertCliques[i].size(); j++){
		eF1 = add_edge(gf.source1, gf.vertCliques[i][j], gf.g).first;
		eF2 = add_edge(gf.vertCliques[i][j], gf.source1, gf.g).first;

		pmf.re[eF1] = eF2;
		pmf.re[eF2] = eF1;

		pmf.c[eF1] = 1;
		pmf.c[eF2] = 0;

		//std::cout << "//clique nodes <-> clique colors" << std::endl;
		for(int k=0; k<colors; k++){
			if(pm.fbc[pmf.rf[gf.vertCliques[i][j]]][k] == 0){
				eF1 = add_edge(gf.vertCliques[i][j], gf.vertCliquesColors[i][k], gf.g).first;
				eF2 = add_edge(gf.vertCliquesColors[i][k], gf.vertCliques[i][j], gf.g).first;

				pmf.re[eF1] = eF2;
				pmf.re[eF2] = eF1;

				pmf.c[eF1] = 1;
				pmf.c[eF2] = 0;
			}
		}
	}
  }

  //std::cout << "clique colors <-> colors" << std::endl;
  // clique colors <-> colors
  for(unsigned int i=0; i<indClq.size(); i++){
	  if(gf.vertCliques[i].size() == 0){
		  std::cout << "A1: Clique leer" << std::endl;
		  continue;
	  }
	for(int k=0; k<colors; k++){
		eF1 = add_edge(gf.vertCliquesColors[i][k], gf.vertColors[k], gf.g).first;
		eF2 = add_edge(gf.vertColors[k], gf.vertCliquesColors[i][k], gf.g).first;

		pmf.re[eF1] = eF2;
		pmf.re[eF2] = eF1;

		pmf.c[eF1] = 1;
		pmf.c[eF2] = 0;
	}
  }

  //std::cout << "source <-> remaining vert" << std::endl;
	// source <-> remaining vert
	for(unsigned int i = 0; i<gf.vertRemaining.size(); i++){
		eF1 = add_edge(gf.source1, gf.vertRemaining[i], gf.g).first;
		eF2 = add_edge(gf.vertRemaining[i], gf.source1, gf.g).first;

		pmf.re[eF1] = eF2;
		pmf.re[eF2] = eF1;

		pmf.c[eF1] = 1;
		pmf.c[eF2] = 0;
	}

  //std::cout << "remaining vert <-> colors" << std::endl;
	// remaining vert <-> colors
	for(unsigned int i = 0; i<gf.vertRemaining.size(); i++){
		for(int k=0; k<colors; k++){
			if(pm.fbc[pmf.rf[gf.vertRemaining[i]]][k] == 0){
				eF1 = add_edge(gf.vertRemaining[i], gf.vertColors[k], gf.g).first;
				eF2 = add_edge(gf.vertColors[k], gf.vertRemaining[i], gf.g).first;

				pmf.re[eF1] = eF2;
				pmf.re[eF2] = eF1;

				pmf.c[eF1] = 1;
				pmf.c[eF2] = 0;
			}
		}
	}
	
  //std::cout << "colors <-> target" << std::endl;
  for(int k=0; k<colors; k++){
	int rU = ceil(parm.n*1. / colors);
    int rL = floor(parm.n*1. / colors);

    rU -= cc.n[k];
    rL -= cc.n[k];

	if(rL < 0){
		rL = 0;
	}

	gf.sumLB += rL;

		eF1 = add_edge(gf.vertColors[k], gf.target1, gf.g).first;
		eF2 = add_edge(gf.target1, gf.vertColors[k], gf.g).first;

		pmf.re[eF1] = eF2;
		pmf.re[eF2] = eF1;

		pmf.c[eF1] = rU - rL;
		pmf.c[eF2] = 0;

		eF1 = add_edge(gf.vertColors[k], gf.target2, gf.g).first;
		eF2 = add_edge(gf.target2, gf.vertColors[k], gf.g).first;

		pmf.re[eF1] = eF2;
		pmf.re[eF2] = eF1;

		pmf.c[eF1] = rL;
		pmf.c[eF2] = 0;
  }

  // tartget 1 <-> source 1
	eF1 = add_edge(gf.target1, gf.source1, gf.g).first;
	eF2 = add_edge(gf.source1, gf.target1, gf.g).first;

	pmf.re[eF1] = eF2;
	pmf.re[eF2] = eF1;

	pmf.c[eF1] = gf.uncoloredVertices;
	pmf.c[eF2] = 0;

	// source 2 <-> target 1
	eF1 = add_edge(gf.source2, gf.target1, gf.g).first;
	eF2 = add_edge(gf.target1, gf.source2, gf.g).first;

	pmf.re[eF1] = eF2;
	pmf.re[eF2] = eF1;

	pmf.c[eF1] = gf.sumLB;
	pmf.c[eF2] = 0;

	return true;
}

/**
 * remove dependencies for the second ford fulkerson
 */

bool EqColoring::initNetwork2(int colors){
  EdgeFord eF1, eF2;

  // source 2 <-> target 1
	eF1 = edge(gf.source2, gf.target1, gf.g).first;
	eF2 = edge(gf.target1, gf.source2, gf.g).first;

	pmf.c[eF1] = 0;
	pmf.c[eF2] = 0;

	pmf.rc[eF1] = 0;
	pmf.rc[eF2] = 0;

	// source 1 <-> target 1
	eF1 = edge(gf.source1, gf.target1, gf.g).first;
	eF2 = edge(gf.target1, gf.source1, gf.g).first;

	pmf.c[eF1] = 0;
	pmf.c[eF2] = 0;

	pmf.rc[eF1] = 0;
	pmf.rc[eF2] = 0;

	// colors <-> target 2, target 1
	for(int k=0; k<colors; k++){
		int rU = ceil(parm.n*1. / colors);
		int rL = floor(parm.n*1. / colors);

		rU -= cc.n[k];
		rL -= cc.n[k];

		if(rL < 0){
			rL = 0;
		}

		eF1 = edge(gf.vertColors[k], gf.target2, gf.g).first;
		eF2 = edge(gf.target2, gf.vertColors[k], gf.g).first;

		pmf.c[eF1] = 0;
		pmf.c[eF2] = 0;

		pmf.rc[eF1] = 0;
		pmf.rc[eF2] = 0;

		eF1 = edge(gf.vertColors[k], gf.target1, gf.g).first;
		eF2 = edge(gf.target1, gf.vertColors[k], gf.g).first;

		pmf.c[eF1] = rU;
		pmf.c[eF2] = 0;
	}

	return true;
}

/**
 * initialize the backup graph for a specific color which is given by the argument
 */

bool EqColoring::initBackupGraphs(int color){
	// clear graph
  gf.g.clear();
  gf.vertCliques.clear();
  gf.vertCliquesColors.clear();
  gf.vertColors.clear();
  gf.vertRemaining.clear();

  gf.uncoloredVertices = curr.uncoloredVertices;

  //if we use a start clique
  if(parm.wStartCl){
	gf.n = 1 + curr.uncoloredVertices + (cl.nCliques-1) * color + color + 1;
  }else{
	gf.n = 1 + curr.uncoloredVertices + cl.nCliques * color + color + 1;
  }
  gf.nNew = gf.n + 2;

  // add source node for phase 1
  gf.source1 = add_vertex(gf.g);
  int counterTmp = 0;

  //add nodes in cliques for phase one
  gf.vertCliques.resize(indClq.size());
  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < indClq[i].size(); j++){
		if(pm.c[indClq[i][j]] == 0){
			gf.vertCliques[i].push_back(add_vertex(gf.g));
			//set relation between nodes for the graph in FF with the nodes from the main graph given by the input
			pmf.rf[gf.vertCliques[i].back()] = indClq[i][j];
			counterTmp++;
		}else{
			std::cout << "ERROR WITH INDCLQ" << std::endl;
		}
	}
  }

  //add remaining nodes (not in cliques)
  vertexIter vIt1, vIt2;
  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
	if(pm.c[*vIt1] == 0 && pm.cl[*vIt1] == 0){
		gf.vertRemaining.push_back(add_vertex(gf.g));
		//set relation between nodes for the graph in FF with the nodes from the main graph given by the input
		pmf.rf[gf.vertRemaining.back()] = *vIt1;
		counterTmp++;
	}
  }
  
  if(counterTmp != gf.uncoloredVertices){
	  std::cout << "error uncolored vertices != counter" << std::endl;
  }

  //add nodes for colors for cliques for phase one
  gf.vertCliquesColors.resize(indClq.size());
  for(unsigned int i = 0; i < indClq.size(); i++){
	  if(gf.vertCliques[i].size() == 0){
		  std::cout << "Groesse war NULL" << std::endl;
		continue;
	  }
    for(unsigned int j = 0; j < (unsigned) color; j++){
		gf.vertCliquesColors[i].push_back(add_vertex(gf.g));
	}
  }

  // add color nodes
  for(unsigned int i = 0; i < (unsigned) color; i++){
	gf.vertColors.push_back(add_vertex(gf.g));
  }
  
  // add target node for phase 1
  gf.target1 = add_vertex(gf.g);

  // add source node for phase 2
  gf.source2 = add_vertex(gf.g);
  // add target node for phase 2
  gf.target2 = add_vertex(gf.g);

  gf.sumLB = 0;

  if(!initNetwork(color)){
	std::cout << "initNetwork error" << std::endl;
	return false;
  }

  return true;
}

void EqColoring::resetCap(){
  EdgesOutFordIter ei1, ei2;
  VertexFordIter it1, it2;

  for(tie(it1, it2) = vertices(gf.g); it1 != it2; it1++){
    for(tie(ei1, ei2) = out_edges(*it1, gf.g); ei1 != ei2; ei1++){
      pmf.rc[*ei1] = 0; 
    }
  }
}

bool EqColoring::pruneFF(int color){
  long flow;

  //start ford fulkerson from source2 as source to target2 as target
  flow = performEKMF(gf.source2, gf.target2);

  //std::cout << "flow = " << flow << " sumLB = " << gf.sumLB << std::endl;
  if(!(flow == gf.sumLB)){
	//there exist no a equitable coloring
	return true;
  }else{
	initNetwork2(color);
    
    flow = performEKMF(gf.source1, gf.target1);

	//std::cout << "flow = " << flow << " uncoloredVertices = " << gf.uncoloredVertices << std::endl;
    if(flow == gf.uncoloredVertices){
	  //there exist a equitable coloring
      return false;
    }else{
	  //there exist no a equitable coloring
      return true;
    }
  }
}

long EqColoring::performEKMF(VertexFord &vs, VertexFord &vt){
  std::vector<default_color_type> col(num_vertices(gf.g));
  std::vector<Traits::edge_descriptor> pred(num_vertices(gf.g));

  Traits::vertex_descriptor s = vs, t = vt;
  
  //boost function to perform a ford fulkerson
  long flow = edmonds_karp_max_flow(gf.g, s, t, pmf.c, pmf.rc, pmf.re, &col[0], &pred[0]);

  return flow;
}

bool EqColoring::useNewIndepCliques(bool sBetterClique){
	//save current data in backup variables
  vertexIter vIt1, vIt2;
  Graph tmpG;
  Cliques tmpCl = cl;
  std::vector< std::vector<Vertex> > tmpIndClq = indClq;

  if(sBetterClique == true){
    copy_graph(g, tmpG);
  }

  //delete all nodes in clique and search for better cliques
  if(parm.wStartCl){
	for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
		if(pm.cl[*vIt1] > 1){
			pm.cl[*vIt1] = 0;
		}
	}

	setClique(startClique.size(), 1, false);
  }else{
	for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
		if(pm.cl[*vIt1] >= 1){
			pm.cl[*vIt1] = 0;
		}
	}

	setClique(0, 0, false);
  }

  findIndepCliques(indClq, true, true);

  //if we do search for better cliques (its always the case)
  if(sBetterClique == true){
	  //check if in the new cliques are more nodes than in the old one
    if(cl.nodesInClique < tmpCl.nodesInClique){
      cl = tmpCl;
      g = tmpG;
      indClq = tmpIndClq;

      return false;
    }else if(cl.nodesInClique == tmpCl.nodesInClique){
	//if the same amount of nodes are in both cliques check which one has the lowest clique
		if(cl.vertInMinCl < tmpCl.vertInMinCl){
			cl = tmpCl;
			g = tmpG;
			indClq = tmpIndClq;

			return false;
		}else{
			return true;
		}
	}else{
		return true;
    }
  }

  std::cout << "ERROR!!!! error updating clique" << std::endl;
  return false;
}


bool EqColoring::updateIndepCliques(Vertex &v){
	int tmp1=1;

	if(!parm.wStartCl){
		tmp1=0;
	}

	//if the node v is in a clique
  if(pm.cl[v] > tmp1 ){
    cl.nNodesSameCl++;

	//check if a specific amount of nodes from the cliques are already visited (given by the probability pNewCl)
    if(cl.nNodesSameCl / (cl.nInCliqueBp * 1.0) > parm.pNewCl){
      if(useNewIndepCliques(true)){
		  //check for better cliques
        cl.nNodesSameCl = 0;
      }else{
		  //no better cliques are found
		if(!removeVinIndClq(v, indClq)){
			//remove vertex from the cliques
			std::cout << "error while removing specific vertex from indepent cliques" << std::endl;
		}

		return false;
      }
    }else{
      if(!removeVinIndClq(v, indClq)){
        std::cout << "error while removing specific vertex from indepent cliques" << std::endl;
      }

      return false;
    }
  }else if(pm.cl[v] == 1 && parm.wStartCl){
    std::cout << "UpdateIndepClique Knoten v aus Startclique!" << std::endl;
  }else{
    return false;
  }

  curr.rankNC = curr.rank;
  curr.createNewGraphs = true;

  return true;
}

bool EqColoring::nodeClique(){
	//check if the algorithm has reached the timeout limit
  if(checkOvertime()){
    return true;
  }

  //if all vertices are colored check if we found a better UB
  if(curr.uncoloredVertices == 0){
	  //check if actual used colors is lower than the current upper bound
    if(curr.nColors < b.UB){
		//check if the current coloring is an equitable one
	    if(checkEqColoring()){
		    b.UB = curr.nColors;

			//std::cout << "NEW UB = " << b.UB << std::endl;
        
			newUBBacktracking();
      }
    }
    
    return true;
  }

  //select a vertex
  Vertex v = passVSS();
  
  //std::cout << "curr.uncoloredVertices = " << curr.uncoloredVertices << std::endl;
  //printClique(startClique);
  //printIndepCliques(indClq);

  //std::cout << "n = " << parm.n << " M = " << curr.M << " T = " << curr.T << " max(nColors, b.LB) = " << std::max(curr.nColors,b.LB) << std::endl;
  //check for pruning rule (paper)
  if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
	  //check for pruning rule (with clique)
	  if(!pruneFF()){	
		  //loop through the colors
		for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
		  //std::cout << "i = " << i << " b.UB = " << b.UB << " curr.nColors = " << curr.nColors << " min = " << std::min(b.UB-1,curr.nColors+1) << std::endl;
		  //check if the color is available for this node
			if(pm.fbc[v][i - 1] == 0){

			c.visitedNodes++; 

			//color this vertex
			colorVertex(v, i);

			//std::cout << "Faerbe Knoten = " << v << " mit Farbe = " << i << " visited node = " << c.visitedNodes <<  std::endl;

			t.startUIC = std::clock();
	        
			//update the independent cliques
			updateIndepCliques(v);
	  
			t.timeUIC += (std::clock() - t.startUIC) / (double) CLOCKS_PER_SEC;

			//curr.node = v;

			//recursive function
			nodeClique(); 

			if(t.timeout){
				return true;
			}

			//check if the algorithm is actual backtracking
          if(bt.status){
			  //if it should backtrack to the actual node then check for better cliques
            if(bt.toRank == curr.rank){
              bt.status = false;

			  if(useNewIndepCliques(true)){
				  cl.nNodesSameCl = 0;
				  c.newCliques++;
			  }
            }else if(bt.toRank < curr.rank){
				//otherwise continue the backtracking
              uncolorVertex(v);

              return true;
            }
          }

          uncolorVertex(v);
        }
      }
    }else{
    
	}
  }else{

  }

  //check for backtracking, if its possible then backtrack
  checkForBacktracking(v);

  //std::cout << "BACKTRACK to = " << bt.toRank << std::endl;
  return false;
}

bool EqColoring::dsatur(){
  if(node()){
    return true;
  }else{
    return false;
  }
}

bool EqColoring::dsaturClique(){
  if(nodeClique()){
    return true;
  }else{
    return false;
  }
}

bool EqColoring::pruningRulePaper(){
  return true;
}

bool EqColoring::pruningRuleFF(){
  return true;
}

/**
 * check if the actual coloring is equitable
 */

bool EqColoring::checkEquitability() const{
  for(unsigned int i = 0; i < curr.nColors; i++){
    for(unsigned int j = i+1; j < curr.nColors; j++){
      if(std::abs(cc.n[i] - cc.n[j]) > 1){
        return false; 
      }  
    }
  }

  return true;
}

/**
 * check if it is a coloring and equitable
 */

bool EqColoring::checkEqColoring() const{
  if(!Coloring::checkColoring()){
    return false;
  }

  if(!checkEquitability()){
    return false;
  }

  return true;
}

/**
 * calculate an initial upper bound
 */

int EqColoring::naiveUB(){
  vertexIter vIt1, vIt2;
  Vertex v;
  bool haveSwapped = false;
  std::pair<int, int> cMinMax;

  //use an inital coloring and swap every time a node from a low color class (with a low amount of nodes in it) with one in a high color class (with a high amount of nodes in it)
  //is this not possible introduce a new color class (i.e. add a new one)
  while(!checkEquitability()){
    cMinMax = findMinMaxColorClass(INT_MAX, INT_MIN);

    for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
      v = *vIt1;

      if(cc.n[pm.c[*vIt1]-1] == cMinMax.second){
        for(int j = 0; j < curr.nColors; j++){
          if(cc.n[j] == cMinMax.first){
            haveSwapped = swapNodeToColor(v, j + 1);
            
            if(haveSwapped){
              break;
            }
          }
        }

        if(haveSwapped){
          break;
        }
      }
    }

    if(!haveSwapped){
      curr.nColors++;
      cc.n[pm.c[v]-1]--;
      pm.c[v] = curr.nColors;
      cc.n[curr.nColors-1]++;
    }
  }

  return curr.nColors;
}

bool EqColoring::swapNodeToColor(Vertex v, int color){
  int tmpColor;
  
  tmpColor = pm.c[v];
  pm.c[v] = color;

  if(!checkColoring()){
    pm.c[v] = tmpColor;

    return false;
  }else{
    cc.n[tmpColor - 1]--;
    cc.n[color - 1]++;

    return true;
  }
}

/**
 * find a the indecis for color class with a lowest and highest amount of colors in it
 */

std::pair<int, int> EqColoring::findMinMaxColorClass(int cMin, int cMax){
  std::pair<int, int> cMinMax = std::make_pair(cMin, cMax);

  for(int i = 0; i < curr.nColors; i++){
    if(cc.n[i] > cMinMax.second){
      cMinMax.second = cc.n[i];
    }
    
    if(cc.n[i] < cMinMax.first){
      cMinMax.first = cc.n[i];
    }
  }

  return cMinMax;
}

/**
 * calculate the upper bound, make use of the naiveUB fnct which is the main task for this
 */

int EqColoring::calcUB(){
  Graph tmpG;
  Current tmp;
  Colors cc_tmp;
  int UB;

  copy_graph(g, tmpG);
  tmp = curr;
  cc_tmp = cc;

  //use an initial coloring
  greedyColoring(g);
  //make it equitable by swapping nodes from the lowest and highest color classes or introduce a new color class
  UB = naiveUB();

  g = tmpG;
  curr = tmp;
  cc = cc_tmp;

  return UB;
}
