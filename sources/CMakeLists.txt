
##########################
# BUILD THE STATIC LIBRARY
##########################

add_library(LP.4 lplib4.c)
install (FILES lplib4.h DESTINATION include COMPONENT headers)
install (TARGETS LP.4 EXPORT LPlib-target DESTINATION lib COMPONENT libraries)
install (EXPORT LPlib-target DESTINATION lib/cmake/${PROJECT_NAME})
export  (PACKAGE LPlib)
