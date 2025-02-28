#include <cugraph/graph.hpp>
#include <cugraph/graph_functions.hpp>               // for weakly_connected_components, etc.
#include <cugraph/algorithms.hpp>                    // for BFS
#include <rmm/device_uvector.hpp>
#include <raft/handle.hpp>
#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <chrono>

// Boost for reading DOT
#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

// -----------------------------------------------------------------------------
// Original "Boost Graph" definitions for CPU side only (parsing DOT):
// -----------------------------------------------------------------------------
typedef boost::property<boost::vertex_name_t, int64_t> vertex_p;
typedef boost::property<boost::edge_weight_t, double>  edge_p;

typedef boost::adjacency_list<
    boost::vecS, boost::vecS,
    boost::bidirectionalS,
    vertex_p,
    edge_p
> CPU_Graph;

// -----------------------------------------------------------------------------
// We'll define a struct for edges we keep after filtering (like "cleanup_graph")
// -----------------------------------------------------------------------------
struct EdgeData {
    int32_t src;
    int32_t dst;
    float   w;
};

// -----------------------------------------------------------------------------
// Example function to replicate logic in "get_tracks", but via cuGraph
// -----------------------------------------------------------------------------
std::vector<std::vector<int>> run_cugraph_pipeline(const CPU_Graph& G,
                                                   double cc_cut,
                                                   double th_min,
                                                   double th_add)
{
    // ---------------------------
    // 1) Build host-side adjacency from G, skipping edges <= cc_cut
    // ---------------------------
    int32_t nV = static_cast<int32_t>(boost::num_vertices(G));
    if (nV == 0) {
        return {};
    }
    int32_t nE_original = static_cast<int32_t>(boost::num_edges(G));

    // We'll store "hit ID" in a host array
    std::vector<int64_t> host_hit_id(nV);
    for (auto v_it = boost::vertices(G).first; v_it != boost::vertices(G).second; ++v_it) {
        int32_t v = static_cast<int32_t>(*v_it);
        int64_t name = boost::get(boost::vertex_name, G, *v_it);
        host_hit_id[v] = name;
    }

    // Collect edges, skipping low-weight
    std::vector<EdgeData> edges;
    edges.reserve(nE_original);

    for (auto e_it = boost::edges(G).first; e_it != boost::edges(G).second; ++e_it) {
        int32_t s = static_cast<int32_t>(boost::source(*e_it, G));
        int32_t t = static_cast<int32_t>(boost::target(*e_it, G));
        double w  = boost::get(boost::edge_weight, G, *e_it);
        if (w > cc_cut) {    // keep edge
            edges.push_back({ s, t, static_cast<float>(w) });
        }
    }
    // Edges are now a "cleaned" version.

    // ---------------------------
    // 2) Create cuGraph device structures
    //    We'll build an SG (single-GPU) graph for simplicity
    // ---------------------------
    raft::handle_t handle; // The RAFT handle for all rapids libs
    // We store source/dest in device vectors
    rmm::device_uvector<int32_t> d_src(edges.size(), handle.get_stream());
    rmm::device_uvector<int32_t> d_dst(edges.size(), handle.get_stream());
    rmm::device_uvector<float>   d_weights(edges.size(), handle.get_stream());

    // Copy edges to device
    for (size_t i = 0; i < edges.size(); ++i) {
        d_src.element(i) = edges[i].src;
        d_dst.element(i) = edges[i].dst;
        d_weights.element(i) = edges[i].w;
    }

    // We assume no multi-GPU, no renumbering needed if vertex IDs are 0..nV-1
    bool has_data = (edges.size() > 0);
    bool sorted_by_source = false;   // We didn't necessarily sort edges by source
    bool sorted_by_dest   = false;
    bool do_renumber      = false;   // If G was not guaranteed to have contiguous IDs, set true

    // Create a GraphCSR or GraphCOO if you prefer. 
    // We'll build an "edgelist_t" then call cugraph::create_graph.
    cugraph::edgelist_t<int32_t, int32_t, float> edgelist{
        d_src.data(),
        d_dst.data(),
        d_weights.data(),
        static_cast<int32_t>(edges.size())
    };

    // Graph properties
    cugraph::graph_properties_t graph_props{ /* is_symmetric = ???, etc.*/ };
    // For a bidirectionalS (like original) we do a directed = false approach if you want
    // to treat it as undirected. But let's say we want it directed:
    bool store_transposed = false; // Typically BFS on a forward graph

    auto graph_tuple = cugraph::create_graph<int32_t, int32_t, float, false, false>(
        handle,
        edgelist,
        cugraph::graph_meta_t<int32_t, int32_t>{ do_renumber, sorted_by_source, sorted_by_dest },
        nV,
        graph_props,
        store_transposed
    );

    // The returned tuple is: (unique_ptr<graph_t>, unique_ptr<resource_handle>, optional renumber_map)
    auto& graph_ptr = std::get<0>(graph_tuple); // graph_t
    auto& edge_weights_ptr = std::get<1>(graph_tuple); // optional edge_weights_t 
    // no renumber_map in get<2> if do_renumber=false

    auto graph_view = graph_ptr->view();
    auto weight_view = edge_weights_ptr ? edge_weights_ptr->view() : std::optional<rmm::device_uvector<float>>{};

    // ---------------------------
    // 3) Weakly Connected Components
    //    If you want an "undirected" sense of components
    // ---------------------------
    // For a directed graph, you might do strongly_connected_components. 
    // Let's do a "weakly" approach:
    rmm::device_uvector<int32_t> d_components(nV, handle.get_stream());

    bool connect_undirected = true; // treat edges as undirected
    cugraph::weakly_connected_components<int32_t, int32_t>(
        handle,
        graph_view,
        d_components.data(),
        connect_undirected
    );

    // bring them back to CPU
    std::vector<int32_t> h_components(nV);
    raft::update_host(h_components.data(), d_components.data(), nV, handle.get_stream());
    handle.sync_stream();

    // Group vertices by component
    std::unordered_map<int32_t, std::vector<int32_t>> comp_map;
    for (int32_t v = 0; v < nV; ++v) {
        comp_map[h_components[v]].push_back(v);
    }

    // We'll create a "sub_graphs" for the "simple paths" if deg <= 2, size >=3
    // We'll do degree check by calling cugraph::degree
    auto out_degs = cugraph::degree<int32_t>(handle, graph_view, cugraph::degree_direction::OUT); 
    // If you want in_degs + out_degs, you'd compute them or do an undirected approach:
    // out_degs has length nV. Let's also do an in_degs if needed:
    auto in_degs  = cugraph::degree<int32_t>(handle, graph_view, cugraph::degree_direction::IN);
    // We'll combine them for total deg in CPU
    std::vector<int32_t> h_out_degs(nV), h_in_degs(nV);
    raft::update_host(h_out_degs.data(), out_degs.data(), nV, handle.get_stream());
    raft::update_host(h_in_degs.data(), in_degs.data(), nV, handle.get_stream());
    handle.sync_stream();

    // Build "simple path" sub-graphs
    std::vector<std::vector<int>> sub_graphs;
    for (auto &kv : comp_map) {
        auto &compVec = kv.second;
        if ((int)compVec.size() < 3) continue;
        bool is_simple = true;
        for (auto vv : compVec) {
            // total deg = in + out
            int total_deg = h_out_degs[vv] + h_in_degs[vv];
            if (total_deg > 2) {
                is_simple = false;
                break;
            }
        }
        if (!is_simple) continue;
        // build final "track" of hit IDs
        std::vector<int> track;
        track.reserve(compVec.size());
        for (auto vv : compVec) {
            track.push_back(static_cast<int>(host_hit_id[vv]));
        }
        sub_graphs.push_back(std::move(track));
    }
    int num_simple_paths = (int) sub_graphs.size();

    // ---------------------------
    // 4) BFS example: a stand-in for "walkthrough"
    //    cugraph has BFS from a single source. For multiple, you'd loop or do multi-source BFS
    // ---------------------------
    // We'll do BFS from an arbitrary starting vertex. 
    // In your code, you might do multiple starts in reverse topological order, 
    // but cugraph doesn't have a built-in topological sort. 
    // We'll do a single BFS from vertex 0 for demonstration.
    rmm::device_uvector<int32_t> d_distances(nV, handle.get_stream());
    rmm::device_uvector<int32_t> d_predecessors(nV, handle.get_stream());

    // init BFS
    cugraph::bfs<int32_t, int32_t, float>(
        handle,
        graph_view,
        weight_view,        // optional edge weights
        d_distances.data(),
        d_predecessors.data(),
        0,                  // start from vertex 0
        false,              // direction_optimizing
        std::numeric_limits<int32_t>::max()
    );
    // BFS results are in d_distances, d_predecessors. We can build a path from 0 to each reachable vertex
    // We'll do a trivial example: any path of length >= 3 => add it to sub_graphs
    std::vector<int32_t> h_pred(nV);
    raft::update_host(h_pred.data(), d_predecessors.data(), nV, handle.get_stream());
    handle.sync_stream();

    // For demonstration: pick a vertex X, rebuild path back to 0
    int X = nV - 1; 
    std::vector<int> path;
    while (X != -1 && X < nV) {
        path.push_back((int)host_hit_id[X]);
        X = h_pred[X];
        if (X < 0) break;
    }
    if ((int)path.size() >= 3) {
        std::reverse(path.begin(), path.end()); 
        sub_graphs.push_back(std::move(path));
    }

    std::cout << "From CC&&Walk: BFS tracks found: " 
              << (sub_graphs.size() - num_simple_paths) << std::endl;

    return sub_graphs;
}

// -----------------------------------------------------------------------------
// Simple I/O for final track writing
// -----------------------------------------------------------------------------
void write_tracks(const std::vector<std::vector<int>>& tracks, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Unable to open " << filename << std::endl;
        return;
    }
    for (auto &t : tracks) {
        for (auto hid : t) {
            file << hid << " ";
        }
        file << "-1 ";
    }
    file << "\n";
}

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dot_file>\n";
        return 1;
    }
    std::string dot_file_name(argv[1]);

    // -- CPU side: read the DOT file using Boost
    std::ifstream dot_file(dot_file_name);
    if (!dot_file) {
        std::cerr << "Error: Unable to open " << dot_file_name << std::endl;
        return 1;
    }

    CPU_Graph G;
    boost::dynamic_properties dp;
    dp.property("hit_id", boost::get(boost::vertex_name, G));
    dp.property("edge_scores", boost::get(boost::edge_weight, G));

    if (!boost::read_graphviz(dot_file, G, dp, "hit_id")) {
        std::cerr << "Error: Unable to parse graph from DOT file.\n";
        return 1;
    }
    double cc_cut = 0.01, th_min = 0.1, th_add = 0.6;

    // Time the pipeline
    auto start_t = std::chrono::high_resolution_clock::now();
    auto final_tracks = run_cugraph_pipeline(G, cc_cut, th_min, th_add);
    auto end_t = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double,std::milli>(end_t - start_t).count();

    std::cout << "Found " << final_tracks.size() << " total tracks.\n";
    std::cout << "Time (cuGraph) = " << elapsed_ms << " ms\n";

    // Write out results
    write_tracks(final_tracks, "tracks.txt");
    return 0;
}
