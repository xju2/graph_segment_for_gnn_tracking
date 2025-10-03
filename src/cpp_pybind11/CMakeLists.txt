pybind11_add_module(_core walk_through_binding.cxx)
target_link_libraries(_core PRIVATE gs_core)

install(TARGETS _core LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
