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
 * \brief property for a vertex in a FF graph that stores the corresponding vertex in the original graph (if existant)
 */

struct ref_vertex_t {
    typedef vertex_property_tag kind;
};

/**
 * @brief properties for a vertex
 */
typedef property<vertex_index_t, int, property<vertex_color_t, int, property<clique_t, int, property<nachbarn_t, std::vector<int>, property<fb_colors_t, std::vector<int>, property<av_color_t, int, property<rank_t, int, property<grad_t , int > > >  > > > > > VertexProperty;
/**
 * @brief properties for an edge
 */
typedef property<edge_weight_t, std::pair<int,int> > EdgeProperty;
/**
 * @brief graph type
 */
typedef adjacency_list<vecS, vecS, undirectedS, VertexProperty, EdgeProperty> Graph;
/**
 * @brief type to define reverse edges for FF graphes
 */
typedef adjacency_list_traits<vecS, vecS, directedS> Traits;
/**
 * @brief FF graph type
 */
typedef adjacency_list<listS, vecS, directedS, property<vertex_name_t, std::string, property<ref_vertex_t, int> >, property<edge_capacity_t, long, property<edge_residual_capacity_t, long, property<edge_reverse_t, Traits::edge_descriptor> > > > GraphFord;
/**
 * @brief edge type in FF graphs
 */
typedef graph_traits<GraphFord>::edge_descriptor EdgeFord;
/**
 * @brief vertex type in FF graphs
 */
typedef graph_traits<GraphFord>::vertex_descriptor VertexFord;
/**
 * @brief vertex iterator in FF graphs
 */
typedef graph_traits<GraphFord>::vertex_iterator VertexFordIter;
/**
 * @brief iterator for outgoing edges in FF graphs
 */
typedef graph_traits<GraphFord>::out_edge_iterator EdgesOutFordIter;
/**
 * @brief iterator for edges in FF graphs
 */
typedef graph_traits<GraphFord>::edge_iterator EdgeFordIter;
/**
 * @brief property for a vertex in FF graphs, describes the connection to a vertex in the "original graph"
 */
typedef property_map<GraphFord, ref_vertex_t>::type RefVertexMap;
/**
 * @brief property for an edge in FF graph, describes the capacity
 */
typedef property_map<GraphFord, edge_capacity_t>::type CapacityMap;
/**
 * @brief property for an edge in FF graph, describes the connection to the reverse edge (preperation for input args for FF from boost lib)
 */
typedef property_map<GraphFord, edge_reverse_t>::type ReverseEdgeMap;
/**
 * @brief property for an edge in FF graph, describes the residual capacity
 */
typedef property_map<GraphFord, edge_residual_capacity_t>::type ResidualCapacityMap;
/**
 * @brief property for a vertex in a graph, describes the index
 */
typedef property_map<Graph, vertex_index_t>::type IndexMap;
/**
 * @brief property for a vertex in a graph, describes the color
 */
typedef property_map<Graph, vertex_color_t>::type ColorMap;
/**
 * @brief property for an edge in a graph, describes the weight
 */
typedef property_map<Graph, edge_weight_t>::type WeightMap;
/**
 * @brief property for a vertex in a graph, describes the clique
 */
typedef property_map<Graph, clique_t>::type CliqueMap;
/**
 * @brief property for a vertex in a graph, that stores all neighbours
 */
typedef property_map<Graph, nachbarn_t>::type NachbarMap;
/**
 * @brief property for a vertex in a graph, that stores all forbidden colors
 */
typedef property_map<Graph, fb_colors_t>::type FBCMap;
/**
 * @brief property for a vertex in a graph, that stores all available colors
 */
typedef property_map<Graph, av_color_t>::type AvColorMap;
/**
 * @brief property for a vertex in a graph, that stores the saturation grad
 */
typedef property_map<Graph, grad_t>::type GradMap;
/**
 * @brief property for a vertex in a graph, that stores the rank (step when its colored by the algorithm)
 */
typedef property_map<Graph, rank_t>::type RankMap;
/**
 * @brief iterator for a vertex
 */
typedef graph_traits<Graph>::vertex_iterator vertexIter;
/**
 * @brief iterator for an edge
 */
typedef graph_traits<Graph>::edge_iterator edgeIter;
/**
 * @brief iterator for the outgoing edges
 */
typedef graph_traits<Graph>::out_edge_iterator outEdgeIter;
/**
 * @brief vertex type
 */
typedef graph_traits<Graph>::vertex_descriptor Vertex;
/**
 * @brief edge type
 */
typedef graph_traits<Graph>::edge_descriptor Edge;
/**
 * @brief iter for adjacent vertices for a given vertex
 */
typedef graph_traits<Graph>::adjacency_iterator adjaIter;
/**
 * @brief type for the degree of a vertex
 */
typedef graph_traits<Graph>::degree_size_type degree_type;

#endif	
