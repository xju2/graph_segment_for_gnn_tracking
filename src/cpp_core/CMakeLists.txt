add_library(gs_core STATIC bestpred_walk.cxx)
set_target_properties(
  gs_core
  PROPERTIES LINKER_LANGUAGE CXX
             CXX_STANDARD 20
             PUBLIC_HEADER "include/graph_segment/bestpred_walk.h")
target_include_directories(gs_core PUBLIC include)
