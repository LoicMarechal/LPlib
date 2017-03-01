
################################
# GLOBAL CONFIGURATION VARIABLES
################################

cmake_minimum_required (VERSION 3.0.1)
set (LPlib_VERSION_MAJOR 3)
set (LPlib_VERSION_MINOR 5)

if (WIN32)
   set (CMAKE_C_FLAGS "-DWIN32")
   set (CMAKE_CONFIGURATION_TYPES "Release")
   if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
     set (CMAKE_INSTALL_PREFIX "$ENV{LP3}" CACHE PATH "Default LPlib install path" FORCE)
   endif()
else ()
   set (CMAKE_C_FLAGS "-O3")
   set (CMAKE_Fortran_FLAGS "-O3")
   if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
     set (CMAKE_INSTALL_PREFIX "/opt/LPlib" CACHE PATH "Default LPlib install path" FORCE)
   endif()
endif()

project (LPlib VERSION ${LPlib_VERSION_MAJOR}.${LPlib_VERSION_MINOR} LANGUAGES C)

include (CheckLanguage)
check_language (Fortran)
if(CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()


#######################################
# SET FILES AND DIRECTORIES TO BE BUILT
#######################################

include_directories ("${PROJECT_SOURCE_DIR}/sources")
add_subdirectory (sources)

if (WIN32)
   find_library (MESHB_LIB Meshb.7 PATHS $ENV{LM7}/lib)
   find_library (MESHBF_LIB Meshbf.7 PATHS $ENV{LM7}/lib)
   find_path (MESHB_INC libmeshb7.h PATHS $ENV{LM7}/include)
else ()
   find_library (MESHB_LIB Meshb.7 PATHS /opt/libMeshb/lib)
   find_library (MESHBF_LIB Meshbf.7 PATHS /opt/libMeshb/lib)
   find_path (MESHB_INC libmeshb7.h PATHS /opt/libMeshb/include)
endif ()

if(MESHB_LIB AND MESHBF_LIB AND MESHB_INC)
   include_directories ("${MESHB_INC}")
   add_subdirectory (examples)
   add_subdirectory (utilities)
else()
   message (WARNING "
   libMeshb.7 could not be found.
   Consequently, examples and utilities will not be built.
   You can download it from GitHUB at https://github.com/LoicMarechal/libMeshb")
endif()

install (FILES LICENSE copyright.txt DESTINATION share/LPlib)
install (DIRECTORY sample_meshes DESTINATION share/LPlib)


######################################
# SET PACKAGE AND DEPLOYMENT VARIABLES
######################################

include (InstallRequiredSystemLibraries)
set (CPACK_GENERATOR "TGZ")
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/copyright.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${LPlib_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${LPlib_VERSION_MINOR}")
set(CPACK_COMPONENTS_ALL applications examples libraries headers)
include (CPack)
