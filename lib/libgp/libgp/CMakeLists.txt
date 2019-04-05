# Created by Manuel Blum on 2011-05-25.
# Copyright 2013 University of Freiburg.

cmake_minimum_required( VERSION 2.8 )

#project( libgp LANGUAGES CXX C)

# if no option is given, standard is release
if( NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE Release )
endif( )

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_FLAGS_RELEASE "-DCLSQUARE -Wall -O2 -fPIC" )
set( CMAKE_CXX_FLAGS_DEBUG   "-DCLSQUARE -Wall -g -O0 -fPIC" )


set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

find_package( Eigen3 3.0.1 REQUIRED )

# ----- Add include directories -----
include_directories( SYSTEM ${EIGEN3_INCLUDE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include/core )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include/wrapper )




# ----- Add sources -----
include( "Sources.cmake" )

# ----- Add gplib library target -----
add_library( gp  ${LIBGP_SRC}  ${LIBGP_INTERFACES} )

# ----- Add gplib wrapper library target -----
add_library( gpwrap  ${LIBGP_WRAPPER_SRC}  ${LIBGP_WRAPPER_INTERFACES} )
target_link_libraries( gpwrap gp )

# ----- Generate PKgconfig file ----
get_property( LIB64  GLOBAL  PROPERTY  FIND_LIBRARY_USE_LIB64_PATHS )
if( "${LIB64}" STREQUAL "TRUE" )
    set(LIBSUFFIX 64)
else()
    set(LIBSUFFIX "")
endif()

set( INSTALL_LIB_DIR  lib${LIBSUFFIX}  CACHE  PATH  "Installation directory for libraries" )
mark_as_advanced( INSTALL_LIB_DIR )
configure_file( ${CMAKE_SOURCE_DIR}/contrib/libgp.pc.in      ${CMAKE_BINARY_DIR}/libgp.pc      @ONLY )
configure_file( ${CMAKE_SOURCE_DIR}/contrib/libgpwrap.pc.in  ${CMAKE_BINARY_DIR}/libgpwrap.pc  @ONLY )





# ----- Add example targets -----
option( BUILD_EXAMPLES "Build tests" ON )
if( BUILD_EXAMPLES )
  add_subdirectory( examples )
endif()




# ----- Testing -----
option( BUILD_TESTS "Build examples" ON )
if( BUILD_TESTS )
  add_subdirectory( tests )
endif()





# ----- Install -----
install( FILES ${LIBGP_INTERFACES} DESTINATION include/gp )

install(  TARGETS gp  gpwrap
          RUNTIME DESTINATION bin
          LIBRARY DESTINATION ${INSTALL_LIB_DIR}
          ARCHIVE DESTINATION ${INSTALL_LIB_DIR} )

install( FILES ${CMAKE_BINARY_DIR}/libgp.pc      DESTINATION ${INSTALL_LIB_DIR}/pkgconfig )
install( FILES ${CMAKE_BINARY_DIR}/libgpwrap.pc  DESTINATION ${INSTALL_LIB_DIR}/pkgconfig )
