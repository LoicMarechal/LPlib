
#################
# BUILD UTILITIES
#################

add_executable(hilbert hilbert.c)
target_link_libraries(hilbert LP.3 ${libMeshb_C_LIBRARY} ${LINK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
install (TARGETS hilbert DESTINATION bin COMPONENT applications)
