######################################################################
## CMakeLists.txt --- gsElasticity
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

## gsElasticity module
project(gsElasticity)

option(GISMO_ELASTICITY       "With Elasticity module"            false  )
if  (${GISMO_ELASTICITY})
message ("  GISMO_ELASTICITY        ${GISMO_ELASTICITY}")
endif()

# Apply G+Smo config
include(gsConfig)

## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules" )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" )

# add filedata folder
add_definitions(-DELAST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

# add example files
if(GISMO_BUILD_EXAMPLES)
    add_custom_target(${PROJECT_NAME})
    add_custom_target(${PROJECT_NAME}-examples)
    add_subdirectory(examples)
else()
    add_subdirectory(examples EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

# Print the contents of {PROJECT_NAME}-examples
get_property(dirs TARGET ${PROJECT_NAME}-examples PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "dirs: ${dirs}")


# # add unittests
# aux_gs_cpp_directory(${PROJECT_SOURCE_DIR}/unittests unittests_SRCS)
# set(gismo_UNITTESTS ${gismo_UNITTESTS} ${unittests_SRCS}
#   CACHE INTERNAL "gismo list of unittests")
