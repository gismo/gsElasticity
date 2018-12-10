### CMakeLists.txt ---
##
## Author: Angelos Mantzaflaris
## Copyright (C) 2018
######################################################################

## gsElasticity module
project(gsElasticity)

## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# Apply G+Smo config
include( gsConfig)

# Add module
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP}
  )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules"
  #EXCLUDE_FROM_ALL 1
  #EXCLUDE_FROM_DEFAULT_BUILD 1
  )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" )
