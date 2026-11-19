#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <stdexcept>
#include <vector>

#include "graph_segment/bestpred_walk.h"

namespace py = pybind11;
using IndexArray = py::array_t<int64_t, py::array::c_style | py::array::forcecast>;
using ScoreArray = py::array_t<float, py::array::c_style | py::array::forcecast>;

// (2, E) edge_index and (E,) scores -> src, dst, scores; then run fn without the GIL
template <typename Fn>
py::tuple run(IndexArray edge_index, ScoreArray scores, Fn fn) {
    if (edge_index.ndim() != 2 || edge_index.shape(0) != 2)
        throw std::invalid_argument("edge_index must have shape (2, num_edges)");
    const auto m = static_cast<size_t>(edge_index.shape(1));
    const int64_t* e = edge_index.data();
    std::vector<int64_t> src(e, e + m), dst(e + m, e + 2 * m);
    std::vector<float> s(scores.data(), scores.data() + scores.size());
    GraphSegment::Tracks t;
    {
        py::gil_scoped_release release;
        t = fn(src, dst, s);
    }
    return py::make_tuple(t.cc, t.walk);
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "Graph segmentation for GNN tracking; use graph_segment.get_tracks";

    m.def("cc_and_walk",
          [](int n, IndexArray ei, ScoreArray sc, float cc_cut, float th_min, float th_add,
             bool walk, bool best_predecessor, float min_track_score, int min_track_hits) {
              return run(ei, sc, [&](auto& s, auto& d, auto& w) {
                  return GraphSegment::cc_and_walk(n, s, d, w, cc_cut, th_min, th_add, walk,
                                                   best_predecessor, min_track_score, min_track_hits);
              });
          },
          "CC (+ walk). edge_index (2, E) must point outward (acyclic).\n"
          "Returns (cc_tracks, walk_tracks), each a list of node-index lists.",
          py::arg("num_nodes"), py::arg("edge_index"), py::arg("scores"), py::arg("cc_cut"),
          py::arg("th_min"), py::arg("th_add"), py::arg("walk"), py::arg("best_predecessor"),
          py::arg("min_track_score") = 0.f, py::arg("min_track_hits") = 0);

    m.def("cc_and_bestpred_walk",
          [](int n, IndexArray ei, ScoreArray sc, float cc_cut, float th_min, float th_add,
             float min_track_score, int min_track_hits) {
              return run(ei, sc, [&](auto& s, auto& d, auto& w) {
                  return GraphSegment::cc_and_bestpred_walk(n, s, d, w, cc_cut, th_min, th_add,
                                                            min_track_score, min_track_hits);
              });
          },
          "CC + best-predecessor walk. edge_index (2, E) must point outward (acyclic).\n"
          "Returns (cc_tracks, walk_tracks), each a list of node-index lists.",
          py::arg("num_nodes"), py::arg("edge_index"), py::arg("scores"),
          py::arg("cc_cut") = 0.01f, py::arg("th_min") = 0.3f, py::arg("th_add") = 0.8f,
          py::arg("min_track_score") = 0.f, py::arg("min_track_hits") = 0);

    m.def("cc_and_junction_removal",
          [](int n, IndexArray ei, ScoreArray sc, float cc_cut, float junction_cut,
             int min_chain_length, float min_track_score, int min_track_hits) {
              return run(ei, sc, [&](auto& s, auto& d, auto& w) {
                  return GraphSegment::cc_and_junction_removal(n, s, d, w, cc_cut, junction_cut,
                                                               min_chain_length, min_track_score,
                                                               min_track_hits);
              });
          },
          "CC + junction removal. Returns (chain_tracks, junction_tracks).",
          py::arg("num_nodes"), py::arg("edge_index"), py::arg("scores"), py::arg("cc_cut"),
          py::arg("junction_cut"), py::arg("min_chain_length"),
          py::arg("min_track_score") = 0.f, py::arg("min_track_hits") = 0);
}
