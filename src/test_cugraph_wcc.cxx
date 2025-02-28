#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "cugraph_wcc.h"


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

int main(int argc, char** argv) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <dot_file>" << std::endl;
        return 1;
    }
    std::string dot_file_name(argv[1]);

    std::cout << "Begin to read: " << dot_file_name <<  std::endl;
    // Load the graph from DOT
    std::ifstream dot_file(dot_file_name);
    if (!dot_file) {
        std::cerr << "Error: Unable to open file " << dot_file_name << std::endl;
        return 1;
    }

    Graph G;
    boost::dynamic_properties dp;
    dp.property("hit_id", boost::get(boost::vertex_name, G));
    dp.property("edge_scores", boost::get(boost::edge_weight, G));

    if (!boost::read_graphviz(dot_file, G, dp, "hit_id")) {
        std::cerr << "Error: Unable to parse graph from DOT file." << std::endl;
        return 1;
    }

    using vertex_t = int32_t;
    vertex_t a,b;
    std::vector<vertex_t> h_rows;
    std::vector<vertex_t> h_cols;
    std::vector<float> h_weights;
    std::vector<vertex_t> h_components;

    auto edges_range = boost::edges(G);
    for (auto it = edges_range.first; it != edges_range.second; ++it) {
        a = static_cast<vertex_t>(boost::source(*it, G));
        b = static_cast<vertex_t>(boost::target(*it, G));
        h_rows.push_back(a);
        h_cols.push_back(b);
        h_weights.push_back(boost::get(boost::edge_weight, G, *it));
    }


    weakly_connected_components<int32_t,int32_t,float>(h_rows, h_cols, h_weights, h_components);
    return 0;
}
