/**
 * @file boost.hpp
 * @brief boost typedefs for an easy handeling with boost lib
 * @author Sven Förster
 * @version 1.1.0
 * @date 2016-02-10
 */

#ifndef BOOST_HPP
#define BOOST_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/bind.hpp>

#include <vector>
#include <iostream>

using namespace boost;

/**
 * \brief property for a vertex that represents to which clique the vertex belongs
 */

struct clique_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex that stores all neighbours for the vertex
 */

struct nachbarn_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex that stores all forbidden colors for this vertex. I.e. colors that a vertex could not use
 */

struct fb_colors_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex that stores all available colors for this vertex. I.e. colors that this vertex can use
 */

struct av_color_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex that stores the rank, i.e. the time (step) when he was colored
 */

struct rank_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex that stores the saturation degree for this vertex
 */

struct grad_t {
    typedef vertex_property_tag kind;
};

/**
 * \brief property for a vertex in a backup graph that stores the corresponding vertex in the original graph (if existant)
 */

struct ref_vertex_t {
    typedef vertex_property_tag kind;
};

typedef property<vertex_index_t, int, property<vertex_color_t, int, property<clique_t, int, property<nachbarn_t, std::vector<int>, property<fb_colors_t, std::vector<int>, property<av_color_t, int, property<rank_t, int, property<grad_t , int > > >  > > > > > VertexProperty;
typedef property<edge_weight_t, std::pair<int,int> > EdgeProperty;
typedef adjacency_list<vecS, vecS, undirectedS, VertexProperty, EdgeProperty> Graph;
typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
typedef adjacency_list<listS, vecS, directedS, property<vertex_name_t, std::string, property<ref_vertex_t, int> >, property<edge_capacity_t, long, property<edge_residual_capacity_t, long, property<edge_reverse_t, Traits::edge_descriptor> > > > GraphFord;
typedef graph_traits<GraphFord>::edge_descriptor EdgeFord;
typedef graph_traits<GraphFord>::vertex_descriptor VertexFord;
typedef graph_traits<GraphFord>::vertex_iterator VertexFordIter;
typedef graph_traits<GraphFord>::out_edge_iterator EdgesOutFordIter;
typedef graph_traits<GraphFord>::edge_iterator EdgeFordIter;
typedef property_map<GraphFord, ref_vertex_t>::type RefVertexMap;
typedef property_map<GraphFord, edge_capacity_t>::type CapacityMap;
typedef property_map<GraphFord, edge_reverse_t>::type ReverseEdgeMap;
typedef property_map<GraphFord, edge_residual_capacity_t>::type ResidualCapacityMap;
typedef property_map<Graph, vertex_index_t>::type IndexMap;
typedef property_map<Graph, vertex_color_t>::type ColorMap;
typedef property_map<Graph, edge_weight_t>::type WeightMap;
typedef property_map<Graph, clique_t>::type CliqueMap;
typedef property_map<Graph, nachbarn_t>::type NachbarMap;
typedef property_map<Graph, fb_colors_t>::type FBCMap;
typedef property_map<Graph, av_color_t>::type AvColorMap;
typedef property_map<Graph, grad_t>::type GradMap;
typedef property_map<Graph, rank_t>::type RankMap;
typedef graph_traits<Graph>::vertex_iterator vertexIter;
typedef graph_traits<Graph>::edge_iterator edgeIter;
typedef graph_traits<Graph>::out_edge_iterator outEdgeIter;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor Edge;
typedef graph_traits<Graph>::adjacency_iterator adjaIter;
typedef graph_traits<Graph>::degree_size_type degree_type;

#endif	
