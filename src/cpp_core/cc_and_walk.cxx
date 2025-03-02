#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <fstream>
#include <string>
#include <chrono>

#include "graph_segment/cc_and_walk.h"

namespace GraphSegment {    
//------------------------------------------------------------------------------
// 1) Functions handling undirected connected components to find "simple paths"
//------------------------------------------------------------------------------
std::vector<std::vector<int>> get_simple_path(const UndirectedGraph& G)
{
    // Weakly connected components (well, 'connected' because G is undirected)
    std::vector<vertex_t> component(num_vertices(G));
    size_t num_components = boost::connected_components(G, &component[0]);

    // Group vertices by component ID
    std::vector<std::vector<Vertex>> component_groups(num_components);
    for(size_t i = 0; i < component.size(); ++i) {
        component_groups[component[i]].push_back(static_cast<Vertex>(i));
    }

    std::vector<std::vector<int>> final_tracks;
    final_tracks.reserve(num_components);

    // For each component, check if it's a simple chain (all nodes deg <= 2)
    for(const auto& sub_graph : component_groups) {
        if (sub_graph.size() < 3) continue;

        bool is_signal_path = true;
        for (auto node : sub_graph) {
            if (degree(node, G) > 2) {
                is_signal_path = false;
                break;
            }
        }
        if (!is_signal_path) continue;

        // Collect hit IDs
        std::vector<int> track;
        track.reserve(sub_graph.size());
        for (auto node : sub_graph) {
            int hit_id = boost::get(boost::vertex_name, G, node);
            track.push_back(hit_id);
        }
        final_tracks.push_back(std::move(track));
    }
    return final_tracks;
}

//------------------------------------------------------------------------------
// 3) Updated find_next_node using the adjacency list
//------------------------------------------------------------------------------
std::vector<int> find_next_node(
    int current_vertex,
    double th_min,
    double th_add,
    const std::vector<std::vector<NeighborInfo>>& adjacency_list,
    const std::vector<bool> &used_vertex
)
{
    // Gather candidate neighbors & scores
    const auto& neighbors = adjacency_list[current_vertex];
    std::vector<std::pair<int, double>> neighbors_scores;
    neighbors_scores.reserve(neighbors.size());

    for (auto &nb : neighbors) {
        if (nb.weight <= th_min) continue;         // skip low-score edges
        if (nb.neighbor == current_vertex) continue; // skip self-loop
        // We do NOT remove used hits here yet; do that after scoring
        neighbors_scores.push_back({nb.neighbor, nb.weight});
    }
    if (neighbors_scores.empty()) return {};

    // Find the single "best neighbor" by max score
    auto best_neighbor_it = std::max_element(neighbors_scores.begin(), neighbors_scores.end(),
                                                [](auto &a, auto &b){
                                                    return a.second < b.second;
                                                });

    // Collect neighbors above th_add
    std::vector<int> next_hits;
    next_hits.reserve(neighbors_scores.size());
    for (auto &p : neighbors_scores) {
        if (p.second > th_add) {
            next_hits.push_back(p.first);
        }
    }
    // If none above th_add, pick just the best
    if (next_hits.empty()) {
        next_hits.push_back(best_neighbor_it->first);
    }

    // Finally remove any that are already used
    // (You could also do this before picking best neighbor, depending on logic.)
    next_hits.erase(std::remove_if(next_hits.begin(), next_hits.end(),
                    [&](int v){ return used_vertex[v]; }),
                    next_hits.end());
    return next_hits;
}

//------------------------------------------------------------------------------
// 4) Updated build_roads to use adjacency_list & used_vertex
//------------------------------------------------------------------------------
std::vector<std::vector<int>> build_roads(
    int starting_vertex,
    std::function<std::vector<int>(int)> next_node_fn,
    std::vector<bool> &used_vertex
)
{
    // 'path' is a list of partial paths under expansion
    std::vector<std::vector<int>> path = {{starting_vertex}};

    while (true) {
        std::vector<std::vector<int>> new_path;
        new_path.reserve(path.size()); // optional but helpful
        bool all_done = true;

        for (const auto &p : path) {
            int last_v = p.back();
            if (last_v == -1) {
                // already ended path
                new_path.push_back(p);
                continue;
            }

            std::vector<int> candidates = next_node_fn(last_v);
            if (candidates.empty()) {
                new_path.push_back(p);
            }
            else {
                // We branch out for each candidate
                all_done = false;
                for (auto c : candidates) {
                    std::vector<int> extended = p;
                    extended.push_back(c);
                    new_path.push_back(std::move(extended));
                }
            }
        }

        path = std::move(new_path);
        if (all_done) break;
    }
    return path;
}

//------------------------------------------------------------------------------
// 5) Cleanup graph: remove isolated vertices & edges below a threshold
//------------------------------------------------------------------------------
Graph cleanup_graph(const Graph& G, double cc_cut)
{
    Graph newG;
    // Map old vertex -> new vertex
    std::vector<int> old_to_new(num_vertices(G), -1);

    int old_id = 0;
    int new_id = 0;
    for (auto v : boost::make_iterator_range(vertices(G))) {
        // skip isolated
        if (in_degree(v, G) == 0 && out_degree(v, G) == 0) {
            old_id++;
            continue;
        }
        // create new vertex
        auto name = boost::get(boost::vertex_name, G, v);
        add_vertex(name, newG);
        old_to_new[old_id++] = new_id++;
    }

    // add edges with weight > cc_cut
    auto edges_range = boost::edges(G);
    for (auto it = edges_range.first; it != edges_range.second; ++it) {
        int s = static_cast<int>(boost::source(*it, G));
        int t = static_cast<int>(boost::target(*it, G));
        double w = boost::get(boost::edge_weight, G, *it);
        if (w <= cc_cut) continue;

        int ns = old_to_new[s];
        int nt = old_to_new[t];
        if (ns < 0 || nt < 0) continue; // was isolated
        add_edge(ns, nt, w, newG);
    }
    return newG;
}

//------------------------------------------------------------------------------
// 6) get_tracks: orchestrates the pipeline and uses local adjacency structures
//------------------------------------------------------------------------------
std::vector<std::vector<int>> get_tracks(const Graph &G, double cc_cut, double th_min, double th_add)
{
    // 1) Prune the original graph
    Graph newG = cleanup_graph(G, cc_cut);
    size_t nV = boost::num_vertices(newG);

    // 2) Build a local adjacency list: adjacency_list[v] => list of (neighbor, weight)
    std::vector<std::vector<NeighborInfo>> adjacency_list(nV);
    {
        auto edges_range = boost::edges(newG);
        adjacency_list.reserve(nV); // optional
        for (auto it = edges_range.first; it != edges_range.second; ++it) {
            auto s = boost::source(*it, newG);
            auto t = boost::target(*it, newG);
            double w = boost::get(boost::edge_weight, newG, *it);

            // s, t are Vertex descriptors [0..nV-1], because newG has vecS
            adjacency_list[s].push_back({static_cast<int>(t), w});
        }
    }

    // 3) Store "hit ID" for each vertex, so we don't call get(vertex_name, ...)
    std::vector<int> vertex_hit_id(nV);
    for (size_t v = 0; v < nV; ++v) {
        vertex_hit_id[v] = boost::get(boost::vertex_name, newG, static_cast<Vertex>(v));
    }

    // 4) Build an undirected copy for the "simple path" step
    UndirectedGraph ugraph;
    ugraph.m_vertices.reserve(nV); // vector-based optimization

    for (size_t v = 0; v < nV; ++v) {
        add_vertex(vertex_hit_id[v], ugraph);
    }
    {
        auto edges_range = boost::edges(newG);
        for (auto it = edges_range.first; it != edges_range.second; ++it) {
            auto s = boost::source(*it, newG);
            auto t = boost::target(*it, newG);
            add_edge(s, t, ugraph);
        }
    }

    // 5) Find "simple paths" in the undirected graph
    std::vector<std::vector<int>> sub_graphs = get_simple_path(ugraph);

    // 6) Mark used vertices in newG that appear in those simple paths
    //    We'll track usage by "vertex index," so we need to invert: hit_id -> vertex
    //    but we only have vertex -> hit_id. Let's build a simple map for that:
    //    (Alternatively, you could do this right inside get_simple_path if you stored
    //     the newly assigned vertex index instead of the original name.)
    std::unordered_map<int, int> hit_id_to_vertex; // for quick usage check
    hit_id_to_vertex.reserve(nV);
    for (size_t v = 0; v < nV; ++v) {
        hit_id_to_vertex[ vertex_hit_id[v] ] = static_cast<int>(v);
    }

    std::vector<bool> used_vertex(nV, false);
    for (const auto& track : sub_graphs) {
        for (int hid : track) {
            if (hit_id_to_vertex.count(hid)) {
                used_vertex[ hit_id_to_vertex[hid] ] = true;
            }
        }
    }
    int num_simple_paths = (int)sub_graphs.size();

    // 7) Topological sort the pruned graph
    std::vector<Vertex> topo_order;
    topo_order.reserve(nV);
    boost::topological_sort(newG, std::back_inserter(topo_order));

    // We'll create a small lambda that, given a vertex index, returns next hits
    auto next_node_fn = [&](int current_vertex) {
        return find_next_node(current_vertex, th_min, th_add, adjacency_list, used_vertex);
    };

    // 8) Walk in reverse topo order and expand roads
    for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it) {
        int v = static_cast<int>(*it);
        if (used_vertex[v]) continue; // skip if already used

        // Build roads from this vertex
        auto roads = build_roads(v, next_node_fn, used_vertex);
        used_vertex[v] = true; // mark it used

        if (roads.empty()) continue;

        // Find the "longest" road among them
        auto longest_it = std::max_element(roads.begin(), roads.end(),
            [](auto &a, auto &b){ return a.size() < b.size(); });

        if (longest_it->size() >= 3) {
            // Copy to a track of "hit IDs"
            std::vector<int> track;
            track.reserve(longest_it->size());
            for (int vtx : *longest_it) {
                used_vertex[vtx] = true;
                track.push_back(vertex_hit_id[vtx]);
            }
            sub_graphs.push_back(std::move(track));
        }
    }

    // Optional: print or return how many new tracks we found
    std::cout << "From CC&&Walk: Number of tracks found by Walkthrough: "
                << sub_graphs.size() - num_simple_paths << std::endl;

    return sub_graphs;
}

//------------------------------------------------------------------------------
// 7) Write tracks to file
//------------------------------------------------------------------------------
void write_tracks(const std::vector<std::vector<int>>& tracks, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open file: " << filename << std::endl;
        return;
    }
    // Each track: space-separated hit IDs, then -1 at end
    for (const auto &track : tracks) {
        for (int hid : track) {
            file << hid << " ";
        }
        file << "-1 ";
    }
    file << std::endl;
}    
} // namespace GraphSegment