
##################
# BUILD C EXAMPLES
##################

if(METIS_FOUND)
   add_executable(colored_grains_partitioning colored_grains_partitioning.c ${PROJECT_SOURCE_DIR}/utilities/lplib4_helpers.c)
   target_link_libraries(colored_grains_partitioning LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
   install (TARGETS colored_grains_partitioning DESTINATION share/LPlib/examples COMPONENT examples)
endif()

add_executable(basic_test basic_test.c)
target_link_libraries(basic_test LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS basic_test DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(tetrahedra_neighbours tetrahedra_neighbours.c)
target_link_libraries(tetrahedra_neighbours LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS tetrahedra_neighbours DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(indirect_writes indirect_writes.c)
target_link_libraries(indirect_writes LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS indirect_writes DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(adaptive_blocking adaptive_blocking.c)
target_link_libraries(adaptive_blocking LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS adaptive_blocking DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(build_edges build_edges.c)
target_link_libraries(build_edges LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS build_edges DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(atomic_locks atomic_locks.c)
target_link_libraries(atomic_locks LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS atomic_locks DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(mem_clear_copy mem_clear_copy.c)
target_link_libraries(mem_clear_copy LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS mem_clear_copy DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(compare_sorts compare_sorts.c)
target_link_libraries(compare_sorts LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${METIS_LIBRARIES})
install (TARGETS compare_sorts DESTINATION share/LPlib/examples COMPONENT examples)
