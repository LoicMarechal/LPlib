
#################
# BUILD UTILITIES
#################

add_executable(hilbert hilbert.c ${libMeshb_HELPERS_DIRS}/libmeshb7_helpers.c)
target_link_libraries(hilbert LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (FILES lplib3_helpers.c DESTINATION helpers COMPONENT sources)
install (FILES lplib3_helpers.h DESTINATION include COMPONENT headers)
install (TARGETS hilbert DESTINATION bin COMPONENT applications)
