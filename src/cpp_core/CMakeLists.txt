find_package(Boost REQUIRED COMPONENTS graph)
message(STATUS "Boost Libaries: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")

add_library(gs_core SHARED
  cc_and_walk.cxx
)
set_target_properties(
  gs_core PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 20
  PUBLIC_HEADER include/graph_segment/cc_and_walk.h
)
target_include_directories(gs_core PUBLIC include)
target_link_libraries(gs_core PUBLIC Boost::graph)


install(TARGETS gs_core
  LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME}/lib64
  PUBLIC_HEADER DESTINATION ${SKBUILD_HEADERS_DIR}/graph_segment
)
