#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <string>


namespace GraphSegment {
// -- Original Boost-based definitions --
typedef boost::property<boost::vertex_name_t, int64_t> vertex_p;
typedef boost::property<boost::edge_weight_t, double> edge_p;

typedef boost::adjacency_list<
    boost::vecS, boost::vecS,
    boost::bidirectionalS,
    vertex_p,
    edge_p,
    boost::no_property
> Graph;

typedef boost::adjacency_list<
    boost::vecS, boost::vecS,
    boost::undirectedS,
    vertex_p,
    boost::no_property,
    boost::no_property
> UndirectedGraph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
using vertex_t = int32_t;

struct NeighborInfo {
    int neighbor;    // neighbor vertex index in newG
    double weight;   // edge weight
};

std::vector<std::vector<int>> get_simple_path(const UndirectedGraph& G);

std::vector<int> find_next_node(
    int current_vertex,
    double th_min,
    double th_add,
    const std::vector<std::vector<NeighborInfo>>& adjacency_list,
    const std::vector<bool> &used_vertex
);

std::vector<std::vector<int>> build_roads(
    int starting_vertex,
    std::function<std::vector<int>(int)> next_node_fn,
    std::vector<bool> &used_vertex
);

Graph cleanup_graph(const Graph& G, double cc_cut);

std::vector<std::vector<int>> get_tracks(
    const Graph &G, double cc_cut, double th_min, double th_add);

void write_tracks(
    const std::vector<std::vector<int>>& tracks, const std::string& filename);

} // namespace GraphSegment