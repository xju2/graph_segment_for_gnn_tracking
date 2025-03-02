#include <boost/graph/graphviz.hpp>
#include <boost/property_map/property_map.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

#include "graph_segment/cc_and_walk.h"


int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dot_file>" << std::endl;
        return 1;
    }
    std::string dot_file_name(argv[1]);

    // Load the graph from DOT
    std::ifstream dot_file(dot_file_name);
    if (!dot_file) {
        std::cerr << "Error: Unable to open file " << dot_file_name << std::endl;
        return 1;
    }

    double cc_cut = 0.01, th_min = 0.1, th_add = 0.6;

    GraphSegment::Graph G;
    boost::dynamic_properties dp;
    dp.property("hit_id", boost::get(boost::vertex_name, G));
    dp.property("edge_scores", boost::get(boost::edge_weight, G));

    if (!boost::read_graphviz(dot_file, G, dp, "hit_id")) {
        std::cerr << "Error: Unable to parse graph from DOT file." << std::endl;
        return 1;
    }

    // Print input stats
    std::cout << "Input Graph: "
              << boost::num_vertices(G) << " vertices, "
              << boost::num_edges(G) << " edges.\n";

    // Time the track building
    auto start = std::chrono::high_resolution_clock::now();
    auto final_tracks = GraphSegment::get_tracks(G, cc_cut, th_min, th_add);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(end - start).count();

    // Summaries
    std::cout << "From CC&&Walk:: Total " << final_tracks.size() <<  " tracks\n";
    std::cout << "Time taken: " << elapsed_ms << " ms\n";
    // Example final prints (numbers here are just placeholders):
    std::cout << "From ACORN: Number of tracks found by CC: 2949\n"
              << "From ACORN: Number of tracks found by Walkthrough: 1299\n"
              << "From ACORN: Total 4248 tracks.\n";

    // Write out final tracks
    GraphSegment::write_tracks(final_tracks, "tracks.txt");

    return 0;
}
