
##################
# BUILD C EXAMPLES
##################

add_executable(basic_test basic_test.c)
target_link_libraries(basic_test LP.3 ${libMeshb_LIBRARIES} ${LINK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS basic_test DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(tetrahedra_neighbours tetrahedra_neighbours.c)
target_link_libraries(tetrahedra_neighbours LP.3 ${libMeshb_LIBRARIES} ${LINK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS tetrahedra_neighbours DESTINATION share/LPlib/examples COMPONENT examples)

add_executable(indirect_writes indirect_writes.c)
target_link_libraries(indirect_writes LP.3 ${libMeshb_LIBRARIES} ${LINK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS indirect_writes DESTINATION share/LPlib/examples COMPONENT examples)


########################
# BUILD FORTRAN EXAMPLES
########################

if (CMAKE_Fortran_COMPILER)
   add_executable(basic_test_f77 basic_test.f)
   target_link_directories(basic_test_f77 PRIVATE ${libMeshb_LINK_DIRS})
   target_link_libraries(basic_test_f77 LP.3 gfortran ${libMeshb_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
   install (TARGETS basic_test_f77 DESTINATION share/LPlib/examples COMPONENT examples)

   add_executable(indirect_writes_f77 indirect_writes.f)
   target_link_directories(indirect_writes_f77 PRIVATE ${libMeshb_LINK_DIRS})
   target_link_libraries(indirect_writes_f77 LP.3 gfortran ${libMeshb_Fortran_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
   install (TARGETS indirect_writes_f77 DESTINATION share/LPlib/examples COMPONENT examples)
endif ()
