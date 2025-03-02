#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

#include "graph_segment/cc_and_walk.h"

// Helper to convert NetworkX graph to Boost Graph
GraphSegment::Graph convert_nx_to_bgl(pybind11::object nx_graph) {
    GraphSegment::Graph G;

    // Add nodes
    pybind11::list nodes = nx_graph.attr("nodes")();
    for (auto node : nodes) {
        boost::add_vertex(G);  // Assume 0-indexed node IDs
    }

    // Add edges
    pybind11::list edges = nx_graph.attr("edges")();
    for (auto edge : edges) {
        auto u = edge.cast<std::pair<int, int>>().first;
        auto v = edge.cast<std::pair<int, int>>().second;
        boost::add_edge(u, v, G);
    }

    return G;
}

// Pybind11 module
PYBIND11_MODULE(_core, m) {
    m.doc() = "C++ module for Connected Components and Walk Through";

    m.def("get_tracks", [](pybind11::object nx_graph, double cc_cut, double th_min, double th_add) {
        auto G = convert_nx_to_bgl(nx_graph);
        return GraphSegment::get_tracks(G, cc_cut, th_min, th_add);
    }, "Process tracks in a graph",
       pybind11::arg("nx_graph"), pybind11::arg("cc_cut"), pybind11::arg("th_min"), pybind11::arg("th_add"));
}
