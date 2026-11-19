#include "graph_segment/bestpred_walk.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace GraphSegment {

namespace {

struct Edge {
    int src;
    int dst;
    float score;
};

void check_edges(int num_nodes, const std::vector<int64_t>& src, const std::vector<int64_t>& dst,
                 const std::vector<float>& scores) {
    if (src.size() != dst.size() || src.size() != scores.size())
        throw std::invalid_argument("src, dst and scores must have the same length");
    for (size_t i = 0; i < src.size(); ++i)
        if (src[i] < 0 || dst[i] < 0 || src[i] >= num_nodes || dst[i] >= num_nodes)
            throw std::invalid_argument("edge index out of range");
}

int find_root(std::vector<int>& parent, int x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

// Drops tracks shorter than min_hits or whose mean internal edge score is below min_score
std::vector<std::vector<int>> post_filter(std::vector<std::vector<int>> tracks,
                                          const std::vector<Edge>& edges, int num_nodes,
                                          float min_score, int min_hits) {
    if (min_hits > 0) {
        tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                                    [&](const auto& t) { return static_cast<int>(t.size()) < min_hits; }),
                     tracks.end());
    }
    if (min_score <= 0.f || tracks.empty()) return tracks;

    std::vector<int> label(num_nodes, -1);
    for (size_t i = 0; i < tracks.size(); ++i)
        for (int h : tracks[i]) label[h] = static_cast<int>(i);
    std::vector<double> sum(tracks.size(), 0.);
    std::vector<int> count(tracks.size(), 0);
    for (const auto& e : edges) {
        if (label[e.src] >= 0 && label[e.src] == label[e.dst]) {
            sum[label[e.src]] += e.score;
            ++count[label[e.src]];
        }
    }
    std::vector<std::vector<int>> kept;
    for (size_t i = 0; i < tracks.size(); ++i) {
        const double mean = count[i] > 0 ? sum[i] / count[i] : 0.;
        if (mean >= min_score) kept.push_back(std::move(tracks[i]));
    }
    return kept;
}

}  // namespace

Tracks cc_and_walk(int num_nodes, const std::vector<int64_t>& src,
                   const std::vector<int64_t>& dst, const std::vector<float>& scores,
                   float cc_cut, float th_min, float th_add, bool walk, bool best_predecessor,
                   float min_track_score, int min_track_hits) {
    check_edges(num_nodes, src, dst, scores);

    // A repeated (u, v) keeps its first position and its last score (as networkx's
    // to_networkx does in acorn); then the edges at or below cc_cut go.
    std::vector<Edge> edges;
    {
        std::unordered_map<int64_t, size_t> seen;
        seen.reserve(src.size());
        for (size_t i = 0; i < src.size(); ++i) {
            const int64_t key = src[i] * num_nodes + dst[i];
            auto [it, inserted] = seen.try_emplace(key, edges.size());
            if (inserted)
                edges.push_back({static_cast<int>(src[i]), static_cast<int>(dst[i]), scores[i]});
            else
                edges[it->second].score = scores[i];
        }
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                   [&](const Edge& e) { return !(e.score > cc_cut); }),
                    edges.end());
    }

    // Weakly connected components; a component is "bad" if any node branches
    std::vector<int> indeg(num_nodes, 0), outdeg(num_nodes, 0), parent(num_nodes);
    std::iota(parent.begin(), parent.end(), 0);
    for (const auto& e : edges) {
        ++outdeg[e.src];
        ++indeg[e.dst];
        parent[find_root(parent, e.src)] = find_root(parent, e.dst);
    }
    std::vector<int> comp(num_nodes);
    for (int v = 0; v < num_nodes; ++v) comp[v] = find_root(parent, v);
    std::vector<char> bad(num_nodes, 0);
    for (int v = 0; v < num_nodes; ++v)
        if (indeg[v] > 1 || outdeg[v] > 1) bad[comp[v]] = 1;
    auto has_edge = [&](int v) { return indeg[v] + outdeg[v] > 0; };

    Tracks out;

    // CC: simple-path components with more than two hits, hits in index order
    {
        std::unordered_map<int, size_t> index;
        std::vector<std::vector<int>> groups;
        for (int v = 0; v < num_nodes; ++v) {
            if (!has_edge(v) || bad[comp[v]]) continue;
            auto [it, inserted] = index.try_emplace(comp[v], groups.size());
            if (inserted) groups.emplace_back();
            groups[it->second].push_back(v);
        }
        for (auto& g : groups)
            if (g.size() > 2) out.cc.push_back(std::move(g));
    }

    // Walk graph: the edges of the branching components, in edge order
    std::vector<std::vector<std::pair<int, float>>> succ(num_nodes), pred(num_nodes);
    for (const auto& e : edges) {
        if (!bad[comp[e.src]]) continue;
        succ[e.src].push_back({e.dst, e.score});
        pred[e.dst].push_back({e.src, e.score});
    }

    // Topological order as networkx gives it: Kahn generations, the first in index order
    std::vector<int> order, remaining(num_nodes, 0), gen;
    for (int v = 0; v < num_nodes; ++v) {
        if (!has_edge(v) || !bad[comp[v]]) continue;
        remaining[v] = static_cast<int>(pred[v].size());
        if (remaining[v] == 0) gen.push_back(v);
    }
    while (!gen.empty()) {
        order.insert(order.end(), gen.begin(), gen.end());
        std::vector<int> next;
        for (int u : gen)
            for (const auto& [v, s] : succ[u])
                if (--remaining[v] == 0) next.push_back(v);
        gen = std::move(next);
    }

    std::vector<char> used(num_nodes, 0);
    std::vector<int> path, best;
    auto free_hit = [&](int v) {
        return !used[v] && std::find(path.begin(), path.end(), v) == path.end();
    };

    // Next hits from path.back(): the best neighbour, plus every other one above th_add;
    // none if the best is at or below th_min. Then the best-predecessor rule.
    auto next_hits = [&](int cur) {
        std::vector<std::pair<int, float>> nb;
        for (const auto& p : succ[cur])
            if (free_hit(p.first)) nb.push_back(p);
        std::vector<int> hits;
        if (nb.empty()) return hits;
        std::stable_sort(nb.begin(), nb.end(),
                         [](const auto& a, const auto& b) { return a.second > b.second; });
        if (nb[0].second <= th_min) return hits;
        for (size_t i = 0; i < nb.size(); ++i) {
            if (i > 0 && nb[i].second <= th_add) break;
            const auto [v, s] = nb[i];
            const bool beaten = best_predecessor && std::any_of(pred[v].begin(), pred[v].end(), [&](const auto& q) {
                return q.first != cur && free_hit(q.first) && q.second > s;
            });
            if (!beaten) hits.push_back(v);
        }
        return hits;
    };

    // Depth-first over the branches in score order; among the longest complete paths
    // the last one wins, as acorn's choose_longest_road picks it.
    auto walk_from = [&](auto& self) -> void {
        const auto hits = next_hits(path.back());
        if (hits.empty()) {
            if (path.size() >= best.size()) best = path;
            return;
        }
        for (int h : hits) {
            path.push_back(h);
            self(self);
            path.pop_back();
        }
    };

    for (int node : walk ? order : std::vector<int>{}) {
        if (used[node]) continue;
        path = {node};
        if (next_hits(node).empty()) {
            used[node] = 1;
            continue;
        }
        best.clear();
        walk_from(walk_from);
        for (int h : best) used[h] = 1;
        out.walk.push_back(best);
    }
    path.clear();

    out.cc = post_filter(std::move(out.cc), edges, num_nodes, min_track_score, min_track_hits);
    out.walk = post_filter(std::move(out.walk), edges, num_nodes, min_track_score, min_track_hits);
    return out;
}

Tracks cc_and_junction_removal(int num_nodes, const std::vector<int64_t>& src,
                               const std::vector<int64_t>& dst, const std::vector<float>& scores,
                               float cc_cut, float junction_cut, int min_chain_length,
                               float min_track_score, int min_track_hits) {
    check_edges(num_nodes, src, dst, scores);
    const size_t m = src.size();

    // components of the edges above cc_cut; repeated edges count in the degrees
    std::vector<int> indeg(num_nodes, 0), outdeg(num_nodes, 0), parent(num_nodes);
    std::iota(parent.begin(), parent.end(), 0);
    for (size_t i = 0; i < m; ++i) {
        if (!(scores[i] > cc_cut)) continue;
        ++outdeg[src[i]];
        ++indeg[dst[i]];
        parent[find_root(parent, src[i])] = find_root(parent, dst[i]);
    }
    std::vector<int> comp(num_nodes), size(num_nodes, 0);
    std::vector<char> bad(num_nodes, 0);
    for (int v = 0; v < num_nodes; ++v) ++size[comp[v] = find_root(parent, v)];
    for (int v = 0; v < num_nodes; ++v)
        if (indeg[v] > 1 || outdeg[v] > 1 || size[comp[v]] < min_chain_length) bad[comp[v]] = 1;

    auto group = [&](const std::vector<int>& label, size_t min_size) {
        std::unordered_map<int, size_t> index;
        std::vector<std::vector<int>> groups;
        for (int v = 0; v < num_nodes; ++v) {
            if (label[v] < 0) continue;
            auto [it, inserted] = index.try_emplace(label[v], groups.size());
            if (inserted) groups.emplace_back();
            groups[it->second].push_back(v);
        }
        std::vector<std::vector<int>> out;
        for (auto& g : groups)
            if (g.size() >= min_size) out.push_back(std::move(g));
        return out;
    };

    std::vector<int> chain(num_nodes);
    for (int v = 0; v < num_nodes; ++v) chain[v] = bad[comp[v]] ? -1 : comp[v];

    // on the other hits: edges above junction_cut that are the best in and best out
    std::vector<int64_t> best_in(num_nodes, -1), best_out(num_nodes, -1);
    for (size_t i = 0; i < m; ++i) {
        if (!(scores[i] > junction_cut) || chain[src[i]] >= 0 || chain[dst[i]] >= 0) continue;
        if (best_in[dst[i]] < 0 || scores[i] > scores[best_in[dst[i]]]) best_in[dst[i]] = i;
        if (best_out[src[i]] < 0 || scores[i] > scores[best_out[src[i]]]) best_out[src[i]] = i;
    }
    std::iota(parent.begin(), parent.end(), 0);
    for (int v = 0; v < num_nodes; ++v) {
        const int64_t i = best_in[v];
        if (i >= 0 && best_out[src[i]] == i)
            parent[find_root(parent, src[i])] = find_root(parent, dst[i]);
    }
    std::vector<int> junction(num_nodes);
    for (int v = 0; v < num_nodes; ++v) junction[v] = chain[v] >= 0 ? -1 : find_root(parent, v);

    // edges for the post-filter: above the lower cut, first position, last score
    std::vector<Edge> edges;
    if (min_track_score > 0.f) {
        const float cut = std::min(cc_cut, junction_cut);
        std::unordered_map<int64_t, size_t> seen;
        for (size_t i = 0; i < m; ++i) {
            auto [it, inserted] = seen.try_emplace(src[i] * num_nodes + dst[i], edges.size());
            if (inserted)
                edges.push_back({static_cast<int>(src[i]), static_cast<int>(dst[i]), scores[i]});
            else
                edges[it->second].score = scores[i];
        }
        edges.erase(std::remove_if(edges.begin(), edges.end(),
                                   [&](const Edge& e) { return !(e.score > cut); }),
                    edges.end());
    }
    Tracks out;
    out.cc = post_filter(group(chain, 2), edges, num_nodes, min_track_score, min_track_hits);
    out.walk = post_filter(group(junction, 2), edges, num_nodes, min_track_score, min_track_hits);
    return out;
}

}  // namespace GraphSegment
