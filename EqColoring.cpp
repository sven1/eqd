#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());
  //b.LB = 3;
  //b.UB = 5;

  if(parm.eqDsatur == 'N'){
    dsatur();
  }else if(parm.eqDsatur == 'C'){
    initPrevGraphsFF();

    dsaturClique();
  }
  
  calcTime();

  printStudyRandomGraphs();
}

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

bool EqColoring::node(){
  if(checkOvertime()){
    return true;
  }

  if(curr.uncoloredVertices == 0){
    if(curr.nColors < b.UB){
	    if(checkEqColoring()){
		    b.UB = curr.nColors;

			//std::cout << "NEW UB = " << b.UB << std::endl;
	    }
    }

    return true;
  }
  
  Vertex v = passVSS();

  for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
    if(pm.fbc[v][i - 1] == 0){
      if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
        c.visitedNodes++;

        colorVertex(v, i);

		//std::cout << c.visitedNodes << " Faerbe Knoten = " << v << " mit Farbe = " << i <<  std::endl;

        node(); 
       
	      if(t.timeout){
	        return true;
	      }
   
        if(bt.status){
          if(bt.toRank == curr.rank){
            bt.status = false;
          }else if(bt.toRank < curr.rank){
            uncolorVertex(v);

            return true;
          }
        }

        uncolorVertex(v);
      }
    }
  }

  checkForBacktracking(v);

  //std::cout << "BACKTRACK to = " << bt.toRank << std::endl;

  return false;
}

bool EqColoring::pruneFF(){
  t.startFF = std::clock();
  
  for(unsigned int i = std::max(b.LB, curr.nColors); i < b.UB; i++){
	t.startBuildNetwork = std::clock();

	if(curr.M > ceil(parm.n*1./i)){
		continue;
	}
    
	if(!initBackupGraphs(i)){
		t.timeBuildNetwork += (std::clock() - t.startBuildNetwork) / (double) CLOCKS_PER_SEC;
		continue;
	}

    t.timeBuildNetwork += (std::clock() - t.startBuildNetwork) / (double) CLOCKS_PER_SEC;
	//std::cout << "test for color = " << i << std::endl;

    c.nFF++;

    if(!pruneFF(i)){
      t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

	  //std::cout << "work for color = " << i << std::endl;
      return false;
    }
  }

  t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

	//std::cout << "no color suitable " << std::endl;

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

	eF1 = add_edge(gf.target1, gf.source1, gf.g).first;
	eF2 = add_edge(gf.source1, gf.target1, gf.g).first;

	pmf.re[eF1] = eF2;
	pmf.re[eF2] = eF1;

	pmf.c[eF1] = gf.uncoloredVertices;
	pmf.c[eF2] = 0;

	eF1 = add_edge(gf.source2, gf.target1, gf.g).first;
	eF2 = add_edge(gf.target1, gf.source2, gf.g).first;

	pmf.re[eF1] = eF2;
	pmf.re[eF2] = eF1;

	pmf.c[eF1] = gf.sumLB;
	pmf.c[eF2] = 0;

	return true;
}

bool EqColoring::initNetwork2(int colors){
  EdgeFord eF1, eF2;

	eF1 = edge(gf.source2, gf.target1, gf.g).first;
	eF2 = edge(gf.target1, gf.source2, gf.g).first;

	pmf.c[eF1] = 0;
	pmf.c[eF2] = 0;

	pmf.rc[eF1] = 0;
	pmf.rc[eF2] = 0;

	eF1 = edge(gf.source1, gf.target1, gf.g).first;
	eF2 = edge(gf.target1, gf.source1, gf.g).first;

	pmf.c[eF1] = 0;
	pmf.c[eF2] = 0;

	pmf.rc[eF1] = 0;
	pmf.rc[eF2] = 0;


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

bool EqColoring::initBackupGraphs(int color){
  gf.g.clear();
  gf.vertCliques.clear();
  gf.vertCliquesColors.clear();
  gf.vertColors.clear();
  gf.vertRemaining.clear();

  gf.uncoloredVertices = curr.uncoloredVertices;
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

  flow = performEKMF(gf.source2, gf.target2);

  //std::cout << "flow = " << flow << " sumLB = " << gf.sumLB << std::endl;
  if(!(flow == gf.sumLB)){
	return true;
  }else{
	initNetwork2(color);
    
    flow = performEKMF(gf.source1, gf.target1);

	//std::cout << "flow = " << flow << " uncoloredVertices = " << gf.uncoloredVertices << std::endl;
    if(flow == gf.uncoloredVertices){
      return false;
    }else{
      return true;
    }
  }
}

long EqColoring::performEKMF(VertexFord &vs, VertexFord &vt){
  std::vector<default_color_type> col(num_vertices(gf.g));
  std::vector<Traits::edge_descriptor> pred(num_vertices(gf.g));

  Traits::vertex_descriptor s = vs, t = vt;
  
  long flow = edmonds_karp_max_flow(gf.g, s, t, pmf.c, pmf.rc, pmf.re, &col[0], &pred[0]);

  return flow;
}

bool EqColoring::useNewIndepCliques(bool sBetterClique){
  vertexIter vIt1, vIt2;
  Graph tmpG;
  Cliques tmpCl = cl;
  std::vector< std::vector<Vertex> > tmpIndClq = indClq;

  if(sBetterClique == true){
    copy_graph(g, tmpG);
  }

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

  if(sBetterClique == true){
    if(cl.nodesInClique < tmpCl.nodesInClique){
      cl = tmpCl;
      g = tmpG;
      indClq = tmpIndClq;

      return false;
    }else if(cl.nodesInClique == tmpCl.nodesInClique){
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
	 //if(cl.nodesInClique - startClique.size() == 0){
      //if(useNewIndepCliques(true)){
        //cl.nNodesSameCl = 0;
        //c.newCliques++;

        //return true;
      //}else{
		  //std::cout << "uncolored V = " << curr.uncoloredVertices << std::endl;
		//std::cout << "KEINE BESSERE CLIQUE" << std::endl;
	  //}
  //}

	int tmp1=1;

	if(!parm.wStartCl){
		tmp1=0;
	}

  if(pm.cl[v] > tmp1 ){
    cl.nNodesSameCl++;

    if(cl.nNodesSameCl / (cl.nInCliqueBp * 1.0) > parm.pNewCl){
      if(useNewIndepCliques(true)){
        cl.nNodesSameCl = 0;
      }else{
		if(!removeVinIndClq(v, indClq)){
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
  if(checkOvertime()){
    return true;
  }

  if(curr.uncoloredVertices == 0){
    if(curr.nColors < b.UB){
	    if(checkEqColoring()){
		    b.UB = curr.nColors;

			//std::cout << "NEW UB = " << b.UB << std::endl;
        
			newUBBacktracking();
      }
    }
	//else{
		//checkForBacktracking(curr.node);
	//}
    
    return true;
  }

  Vertex v = passVSS();
  
  std::cout << "curr.uncoloredVertices = " << curr.uncoloredVertices << std::endl;
  printClique(startClique);
  printIndepCliques(indClq);

  //std::cout << "n = " << parm.n << " M = " << curr.M << " T = " << curr.T << " max(nColors, b.LB) = " << std::max(curr.nColors,b.LB) << std::endl;
  if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
	  if(!pruneFF()){	
		for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
			if(pm.fbc[v][i - 1] == 0){

			c.visitedNodes++; 

			colorVertex(v, i);

			//std::cout << "Faerbe Knoten = " << v << " mit Farbe = " << i << " visited node = " << c.visitedNodes <<  std::endl;

			t.startUIC = std::clock();
	        
			updateIndepCliques(v);
	  
			t.timeUIC += (std::clock() - t.startUIC) / (double) CLOCKS_PER_SEC;

			//curr.node = v;

			nodeClique(); 

			if(t.timeout){
				return true;
			}

          if(bt.status){
            if(bt.toRank == curr.rank){
              bt.status = false;

			  if(useNewIndepCliques(true)){
				  cl.nNodesSameCl = 0;
				  c.newCliques++;
			  }
            }else if(bt.toRank < curr.rank){
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

bool EqColoring::checkEqColoring() const{
  if(!Coloring::checkColoring()){
    return false;
  }

  if(!checkEquitability()){
    return false;
  }

  return true;
}

int EqColoring::naiveUB(){
  vertexIter vIt1, vIt2;
  Vertex v;
  bool haveSwapped = false;
  std::pair<int, int> cMinMax;

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

int EqColoring::calcUB(){
  Graph tmpG;
  Current tmp;
  Colors cc_tmp;
  int UB;

  copy_graph(g, tmpG);
  tmp = curr;
  cc_tmp = cc;

  greedyColoring(g);
  UB = naiveUB();

  g = tmpG;
  curr = tmp;
  cc = cc_tmp;

  return UB;
}
