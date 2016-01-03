#include "EqColoring.hpp"

EqColoring::EqColoring() : Coloring(){
}

EqColoring::EqColoring(const Parameters &parm) : Coloring(parm){
  setBounds(b.LB, calcUB());

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

  return false;
}

bool EqColoring::pruneFF(){
  t.startFF = std::clock();
  
  for(unsigned int i = std::max(b.LB,curr.nColors); i < b.UB; i++){
    if(!initBackupGraphs(i)){
		continue;
	}

    c.nFF++;

    if(!pruneFF(i)){
      t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

      return false;
    }
  }

  t.timeFF += (std::clock() - t.startFF) / (double) CLOCKS_PER_SEC;

  return true;
}

bool EqColoring::initA1(){
  EdgeFord eF1, eF2;

  for(int i = 1; i <= gf.uncoloredVertices; i++){
    eF1 = add_edge(gf.vert[0], gf.vert[i], gf.g).first;
    eF2 = add_edge(gf.vert[i], gf.vert[0], gf.g).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = 1;
    pmf.c[eF2] = 0;
  }

  return true;
}

bool EqColoring::initA2andA3(int l){
  EdgeFord eF1, eF2;
  VertexFordIter vIt1, vIt2;

  int aUVPos = 1;
  int tmpIndex, colorPos;
  int color = l;

  for(unsigned int i = 0; i < indClq.size(); i++){
    for(unsigned int j = 0; j < indClq[i].size(); j++){
      pmf.rf[gf.vert[aUVPos]] = indClq[i][j];

      for(int k = 0; k < color; k++){
        if(pm.fbc[indClq[i][j]][k] == 0){

	        //Vorknoten von Farbe k der i-ten Clique
          tmpIndex = gf.uncoloredVertices + (k + 1) + i * color; 

          eF1 = add_edge(gf.vert[aUVPos], gf.vert[tmpIndex], gf.g).first;
          eF2 = add_edge(gf.vert[tmpIndex], gf.vert[aUVPos], gf.g).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          //zur Farbe
          colorPos = gf.uncoloredVertices + (k + 1) + cl.nCliques * color; 
			if(edge(gf.vert[tmpIndex], gf.vert[colorPos], gf.g).second == false){
          eF1 = add_edge(gf.vert[tmpIndex], gf.vert[colorPos], gf.g).first;
          eF2 = add_edge(gf.vert[colorPos], gf.vert[tmpIndex], gf.g).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;
			}
        }
      }

      aUVPos++;
    }
  }

  for(tie(vIt1,vIt2) = vertices(g); vIt1 != vIt2; vIt1++){
    if(pm.c[*vIt1] == 0 && pm.cl[*vIt1] == 0){
      pmf.rf[gf.vert[aUVPos]] = *vIt1;

      for(int k = 0; k < color; k++){
        if(pm.fbc[*vIt1][k] == 0){
          tmpIndex = gf.uncoloredVertices + (k + 1) + (cl.nCliques - 1) * color; 

          eF1 = add_edge(gf.vert[aUVPos], gf.vert[tmpIndex], gf.g).first;
          eF2 = add_edge(gf.vert[tmpIndex], gf.vert[aUVPos], gf.g).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = 1;
          pmf.c[eF2] = 0;

          colorPos = gf.uncoloredVertices + (k + 1) + cl.nCliques * color; 

          eF1 = add_edge(gf.vert[tmpIndex], gf.vert[colorPos], gf.g).first;
          eF2 = add_edge(gf.vert[colorPos], gf.vert[tmpIndex], gf.g).first;

          pmf.re[eF1] = eF2;
          pmf.re[eF2] = eF1;

          pmf.c[eF1] = gf.uncoloredVertices - (cl.nodesInClique - startClique.size());
          pmf.c[eF2] = 0;
        }
      }

      aUVPos++;
    }
  }

  if(aUVPos != gf.uncoloredVertices + 1){
    std::cout << "aUVPos = " << aUVPos << " and uncoloredVertices = " << curr.uncoloredVertices << std::endl;
    std::cout << "error while constructing network" << std::endl;
    std::cout << "not using all uncolored Vertices" << std::endl;

    return false;
  }

  return true;
}

int EqColoring::initA4(int l){
  int sumLB = 0, rU, rL, colorPos, n = gf.n, nNew = gf.nNew, color = l;
  EdgeFord eF1, eF2;

  for(int i = 1; i <= color; i++){
    rU = ceil(parm.n / color);
    rL = floor(parm.n / color);

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;

	  return 0;
    }

    sumLB += rL;

    colorPos = gf.uncoloredVertices + i + cl.nCliques * color; 

    //von Farben zur alten Senke
    eF1 = add_edge(gf.vert[colorPos], gf.vert[n - 1], gf.g).first;
    eF2 = add_edge(gf.vert[n - 1], gf.vert[colorPos], gf.g).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = rU - rL;
    pmf.c[eF2] = 0;

    //von Farben zur neuen Senke
    eF1 = add_edge(gf.vert[colorPos], gf.vert[nNew - 1], gf.g).first;
    eF2 = add_edge(gf.vert[nNew - 1], gf.vert[colorPos], gf.g).first;

    pmf.re[eF1] = eF2;
    pmf.re[eF2] = eF1;

    pmf.c[eF1] = rL;
    pmf.c[eF2] = 0;
  }

  return sumLB;
}

bool EqColoring::initRespectLB(int sumLB){
  EdgeFord eF1, eF2;
  int n = gf.n, nNew = gf.nNew;

  //von neuer Quelle zur alten Senke
  eF1 = add_edge(gf.vert[nNew - 2], gf.vert[n - 1], gf.g).first;
  eF2 = add_edge(gf.vert[n - 1], gf.vert[nNew - 2], gf.g).first;

  pmf.re[eF1] = eF2;
  pmf.re[eF2] = eF1;

  pmf.c[eF1] = sumLB;
  pmf.c[eF2] = 0;

  //von alter Senke zur alten Quelle
  eF1 = add_edge(gf.vert[n - 1], gf.vert[0], gf.g).first;
  eF2 = add_edge(gf.vert[0], gf.vert[n - 1], gf.g).first;

  pmf.re[eF1] = eF2;
  pmf.re[eF2] = eF1;

  pmf.c[eF1] = INT_MAX;
  pmf.c[eF2] = 0;

  return true;
}

bool EqColoring::removeRespectLB(int l, int sumLB){
  int rL, rU, colorPos, n = gf.n, nNew = gf.nNew, color = l;

  EdgeFord eF1, eF2;

  for(int i = 1; i <= color; i++){
    rU = ceil(parm.n / color);
    rL = floor(parm.n / color);

    rU -= cc.n[i-1];
    rL -= cc.n[i-1];

    if(rL < 0){
      rL = 0;
    }

    colorPos = gf.uncoloredVertices + i + cl.nCliques * color; 

    //von Farbe zur alten Senke
    eF1 = edge(gf.vert[colorPos], gf.vert[n - 1], gf.g).first;
    eF2 = edge(gf.vert[n - 1], gf.vert[colorPos], gf.g).first;

    pmf.c[eF1] = rU;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = rU - rL;
    pmf.rc[eF2] = 0;

    //von Farbe zur neuen Senke
    eF1 = edge(gf.vert[colorPos], gf.vert[nNew - 1], gf.g).first;
    eF2 = edge(gf.vert[nNew - 1], gf.vert[colorPos], gf.g).first;

    pmf.c[eF1] = 0;
    pmf.c[eF2] = 0;

    pmf.rc[eF1] = 0;
    pmf.rc[eF2] = 0;
  }

  //von alter Senke zur alten Quelle
  eF1 = edge(gf.vert[n - 1], gf.vert[0], gf.g).first;
  eF2 = edge(gf.vert[0], gf.vert[n - 1], gf.g).first;

  pmf.c[eF1] = 0;
  pmf.c[eF2] = 0;

  pmf.rc[eF1] = 0;
  pmf.rc[eF2] = 0;

  //von neuer Quelle zur alten Senke
  eF1 = add_edge(gf.vert[nNew - 2], gf.vert[n - 1], gf.g).first;
  eF2 = add_edge(gf.vert[n - 1], gf.vert[nNew - 2], gf.g).first;

  pmf.c[eF1] = 0;
  pmf.c[eF2] = 0;

  pmf.rc[eF1] = 0;
  pmf.rc[eF2] = 0;

  return true;
}

bool EqColoring::initBackupGraphs(int color){
  gf.g.clear();
  gf.vert.clear();

  gf.uncoloredVertices = curr.uncoloredVertices;
  gf.n = 1 + curr.uncoloredVertices + cl.nCliques * color + color + 1;
  gf.nNew = gf.n + 2;

  for(int j = 0; j < gf.nNew; j++){
    gf.vert.push_back(add_vertex(gf.g));
  }
    
  initA1();
  
  initA2andA3(color);
    
  int sumLB = initA4(color);

  if(sumLB == 0){
	return false;
  }

  initRespectLB(sumLB);

  gf.sumLB = sumLB;
 
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

  flow = performEKMF(gf.vert[gf.nNew - 2], gf.vert[gf.nNew - 1]);

  if(!(flow == gf.sumLB)){
	return true;
  }else{
    removeRespectLB(color, gf.sumLB);
    
    flow = performEKMF(gf.vert[0], gf.vert[gf.n - 1]);

    if(flow == curr.uncoloredVertices){
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
    if(cl.nodesInClique <= tmpCl.nodesInClique){
      cl = tmpCl;
      g = tmpG;
      indClq = tmpIndClq;

      return false;
    }else{
		return true;
    }
  }

  std::cout << "ERROR!!!! error updating clique" << std::endl;
  return false;
}


bool EqColoring::updateIndepCliques(Vertex &v){
     if(cl.nodesInClique - startClique.size() == 0){
      if(useNewIndepCliques(true)){
        cl.nNodesSameCl = 0;
        c.newCliques++;

        return true;
      }else{
		  //std::cout << "uncolored V = " << curr.uncoloredVertices << std::endl;
		//std::cout << "KEINE BESSERE CLIQUE" << std::endl;
	  }
  }

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
        
			newUBBacktracking();
      }
    }

    
    return true;
  }

  Vertex v = passVSS();

  if(parm.n >= (curr.M - 1) * std::max(curr.nColors, b.LB) + curr.T){
	  if(!pruneFF()){	
		for(int i = 1; i <= std::min(b.UB - 1, curr.nColors + 1); i++){
			if(pm.fbc[v][i - 1] == 0){

			c.visitedNodes++; 

			colorVertex(v, i);

			t.startUIC = std::clock();
	        
			updateIndepCliques(v);
	  
			t.timeUIC += (std::clock() - t.startUIC) / (double) CLOCKS_PER_SEC;

			nodeClique(); 

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
    }else{
    }
  }else{
  }

  checkForBacktracking(v);

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
