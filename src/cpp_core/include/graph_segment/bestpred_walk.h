#pragma once

#include <cstdint>
#include <vector>

namespace GraphSegment {

struct Tracks {
    std::vector<std::vector<int>> cc;    // components that are simple paths
    std::vector<std::vector<int>> walk;  // paths walked through the other components
};

/// CC + best-predecessor walk on a directed graph given as edge arrays.
///
/// The edges must already point outward (from the hit nearer the interaction point),
/// so the graph is acyclic. Edges scoring <= cc_cut are dropped first. A weakly
/// connected component in which every node has at most one incoming and one
/// outgoing edge, and more than two nodes, is a track. The other components are walked:
/// from each unused node in topological order, follow every neighbour scoring
/// above th_add, or only the best one if none does and it scores above th_min, and
/// keep the longest path. The best-predecessor rule removes a step u -> v when another
/// unused hit, not on the current path, has a higher-scoring edge into v.
///
/// Post-filters: tracks with fewer than min_track_hits hits, or whose mean edge
/// score (over the edges between their own hits) is below min_track_score, are
/// dropped. Their hits stay used.
///
/// The output is the same as acorn's cc_and_walk_utils.cc_and_bestpred_walk (the
/// same tie-breaks), with tracks as node indices.
Tracks cc_and_bestpred_walk(
    int num_nodes, const std::vector<int64_t>& src, const std::vector<int64_t>& dst,
    const std::vector<float>& scores, float cc_cut, float th_min, float th_add,
    float min_track_score = 0.f, int min_track_hits = 0);

}  // namespace GraphSegment
