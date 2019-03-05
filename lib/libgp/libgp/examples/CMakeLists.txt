# libgp - Gaussian process library for Machine Learning
# Copyright (c) 2013, Manuel Blum <mblum@informatik.uni-freiburg.de>
# All rights reserved.

add_executable( gpdense gp_example_dense.cc )
target_link_libraries( gpdense gp )

add_executable( gpwrap_example  gp_wrapper_example.cpp )
target_link_libraries( gpwrap_example  gp  gpwrap )


file( COPY ${CMAKE_SOURCE_DIR}/examples/data DESTINATION .)
