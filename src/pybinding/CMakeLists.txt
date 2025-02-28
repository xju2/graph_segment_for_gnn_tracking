set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core walk_through_binding.cxx)

install(TARGETS _core LIBRARY DESTINATION graph_segment)
