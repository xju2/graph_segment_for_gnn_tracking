find_package(Boost REQUIRED COMPONENTS graph)
message(STATUS "Boost Libaries: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")

add_library(gs_core SHARED
  cc_and_walk.cxx
)
set_target_properties(
  gs_core PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 20
)
target_include_directories(gs_core PUBLIC include)
target_link_libraries(gs_core PUBLIC Boost::graph)
target_sources(gs_core
  PUBLIC
    "include/graph_segment/cc_and_walk.h"
)

install(TARGETS gs_core
  LIBRARY DESTINATION lib
)
install(DIRECTORY include/ # Replace "include" with your header directory
  DESTINATION include
  FILES_MATCHING PATTERN "*.h*" # Install all header files
)
