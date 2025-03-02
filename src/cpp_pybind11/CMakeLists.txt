set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core walk_through_binding.cxx)
target_link_libraries(_core PUBLIC gs_core)
# set_target_properties(
#   _core PROPERTIES
#   BUILD_WITH_INSTALL_RPATH TRUE
#   INSTALL_RPATH_DIR lib64
# )

install(TARGETS _core LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
