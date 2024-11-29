
#################
# BUILD UTILITIES
#################

add_executable(hilbert hilbert.c)
target_link_libraries(hilbert LP.3 ${libMeshb_LIBRARIES} ${math_LIBRARIES} ${LINK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (FILES lplib3_helpers.h DESTINATION include COMPONENT headers)
install (TARGETS hilbert DESTINATION bin COMPONENT applications)
