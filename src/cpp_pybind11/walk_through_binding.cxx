#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <boost/graph/adjacency_list.hpp>
#include <vector>

#include "graph_segment/cc_and_walk.h"
#include "graph_segment/bestpred_walk.h"

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

    using IndexArray = pybind11::array_t<int64_t, pybind11::array::c_style | pybind11::array::forcecast>;
    using ScoreArray = pybind11::array_t<float, pybind11::array::c_style | pybind11::array::forcecast>;
    m.def("cc_and_bestpred_walk",
          [](int num_nodes, IndexArray edge_index, ScoreArray scores, float cc_cut, float th_min,
             float th_add, float min_track_score, int min_track_hits) {
              if (edge_index.ndim() != 2 || edge_index.shape(0) != 2)
                  throw std::invalid_argument("edge_index must have shape (2, num_edges)");
              const auto m = static_cast<size_t>(edge_index.shape(1));
              const int64_t* e = edge_index.data();
              std::vector<int64_t> src(e, e + m), dst(e + m, e + 2 * m);
              std::vector<float> s(scores.data(), scores.data() + scores.size());
              GraphSegment::Tracks t;
              {
                  pybind11::gil_scoped_release release;
                  t = GraphSegment::cc_and_bestpred_walk(num_nodes, src, dst, s, cc_cut, th_min, th_add,
                                                         min_track_score, min_track_hits);
              }
              return pybind11::make_tuple(t.cc, t.walk);
          },
          "CC + best-predecessor walk. edge_index (2, E) must point outward (acyclic).\n"
          "Returns (cc_tracks, walk_tracks), each a list of node-index lists.",
          pybind11::arg("num_nodes"), pybind11::arg("edge_index"), pybind11::arg("scores"),
          pybind11::arg("cc_cut") = 0.01f, pybind11::arg("th_min") = 0.3f, pybind11::arg("th_add") = 0.8f,
          pybind11::arg("min_track_score") = 0.f, pybind11::arg("min_track_hits") = 0);
}
