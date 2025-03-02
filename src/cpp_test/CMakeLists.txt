add_executable(eval_walk_through walk_through.cxx)
set_target_properties(
  eval_walk_through PROPERTIES
  LINKER_LANGUAGE CXX
  CXX_STANDARD 20
)
target_link_libraries(eval_walk_through PUBLIC gs_core)


install(TARGETS eval_walk_through
  RUNTIME DESTINATION ${SKBUILD_SCRIPTS_DIR}
)
