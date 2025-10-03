add_library(gs_core STATIC cc_and_walk.cxx)
set_target_properties(
  gs_core
  PROPERTIES LINKER_LANGUAGE CXX
             CXX_STANDARD 20
             PUBLIC_HEADER include/graph_segment/cc_and_walk.h)
target_include_directories(gs_core PUBLIC include)
target_link_libraries(gs_core PUBLIC Boost::graph)
