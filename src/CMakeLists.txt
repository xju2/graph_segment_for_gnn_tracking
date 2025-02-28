if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_compile_options(-g)
else()
  add_compile_options(-O4)
endif()

# cpp compile options
add_compile_options(-Wfatal-errors)
add_compile_options(-Wno-deprecated)
add_compile_options(-std=c++17)
# add_compile_options(-Wc++11-extensions)
add_compile_options(-Wall)
add_compile_options(-fPIC)
add_compile_options(-pedantic)
# add_compile_options(-Wgnu-zero-variadic-macro-arguments)

find_package(Boost REQUIRED COMPONENTS graph)
message(STATUS "Boost Libaries: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")

add_executable(walk_through walk_through.cxx)
set_target_properties(
  walk_through PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 20
)
target_link_libraries(walk_through PUBLIC Boost::graph)

add_subdirectory(pybinding)

