set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

# add_library(graph_segment MODULE walk_through_binding.cxx)
# target_link_libraries(graph_segment PRIVATE pybind11::module)
# target_include_directories(graph_segment PRIVATE ${Python3_INCLUDE_DIRS})
# set_target_properties(graph_segment PROPERTIES PREFIX "" SUFFIX ".so")
# target_compile_definitions(graph_segment PRIVATE VERSION_INFO=${PROJECT_VERSION})

pybind11_add_module(graph_segment walk_through_binding.cxx)

install(TARGETS graph_segment LIBRARY DESTINATION ${SKBUILD_PLATLIB_DIR}/graph_segment/lib64)
