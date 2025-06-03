
#################
# BUILD UTILITIES
#################

add_executable(hilbert hilbert.c ${libMeshb_HELPERS_DIRS}/libmeshb8_helpers.c)
target_link_libraries(hilbert LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (FILES lplib4_helpers.c DESTINATION helpers COMPONENT sources)
install (FILES lplib4_helpers.h DESTINATION include COMPONENT headers)
install (TARGETS hilbert DESTINATION bin COMPONENT applications)

add_executable(cpu_bandwidth cpu_bandwidth.c)
target_link_libraries(cpu_bandwidth LP.4 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS cpu_bandwidth DESTINATION bin COMPONENT applications)
