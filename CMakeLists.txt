cmake_minimum_required( VERSION 2.8.7 )
project( GraphSpectralEmbedding )

find_package( VTK )
include( ${VTK_USE_FILE} )

find_package( Armadillo )
include_directories( ${PROJECT_BINARY_DIR} ${ARMADILLO_INCLUDE_DIR} )

add_executable( GraphSpectralEmbedding src/GraphSpectralEmbedding.cpp )
target_link_libraries( GraphSpectralEmbedding ${ARMADILLO_LIBRARIES} ${VTK_LIBRARIES} )


