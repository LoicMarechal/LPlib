
##################
# BUILD C EXAMPLES
##################

add_executable(basic_test basic_test.c)
target_link_libraries(basic_test LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS basic_test DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(tetrahedra_neighbours tetrahedra_neighbours.c)
target_link_libraries(tetrahedra_neighbours LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS tetrahedra_neighbours DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(indirect_writes indirect_writes.c)
target_link_libraries(indirect_writes LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS indirect_writes DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(adaptive_blocking adaptive_blocking.c)
target_link_libraries(adaptive_blocking LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS adaptive_blocking DESTINATION share/LPlib/examples COMPONENT examples)
