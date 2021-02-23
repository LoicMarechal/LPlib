
##########################
# BUILD THE STATIC LIBRARY
##########################

add_library(LP.3 lplib3.c)
install (FILES lplib3.h lplib3.ins DESTINATION include COMPONENT headers)
install (TARGETS LP.3 EXPORT LPlib-target DESTINATION lib COMPONENT libraries)
install (EXPORT LPlib-target DESTINATION lib/cmake/${PROJECT_NAME})
export  (PACKAGE LPlib)
