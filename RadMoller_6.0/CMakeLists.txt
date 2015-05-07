#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(RadMollerGen)

#----------------------------------------------------------------------------
# Setup CERN root 
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
FIND_PACKAGE(ROOT REQUIRED)

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
			${PROJECT_SOURCE_DIR}
		    ${ROOT_INCLUDE_DIR})
file(GLOB sourcesSIM ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headersSIM ${PROJECT_SOURCE_DIR}/include/*.h  )

file(GLOB sourcesCOM ${PROJECT_SOURCE_DIR}/src/COM/*.cc)
file(GLOB headersCOM ${PROJECT_SOURCE_DIR}/include/COM/*.h  )


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(RadMoller RadMoller.cc ${sourcesSIM} ${headersSIM} )

target_link_libraries(RadMoller "-L${ROOT_LIBRARY_DIR}" Core Cint Tree RIO Net FTGL MathCore MathMore Thread Hist Matrix Physics Gpad Minuit Foam)

target_link_libraries(RadMoller ${cadmesh_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build errorpropagation. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#


foreach(_script ${errorpropagation_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(errorpropagation DEPENDS errprop)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS RadMoller DESTINATION bin)

 
