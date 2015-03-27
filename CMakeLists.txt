######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

get_filename_component(GM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#project(${GM_NAME}_dev)

## #################################################################
## Collect files
## #################################################################

aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_H  )
aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_HPP)

## #################################################################
## Add module
## #################################################################

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_INS)
  if(${GM_NAME}_INS)
    LIST( REMOVE_ITEM ${GM_NAME}_CPP ${${GM_NAME}_INS})
  endif()
endif()

include_directories(${GISMO_DEV_INCLUDE_DIRS})

add_library(${GM_NAME}_dev OBJECT
  ${${GM_NAME}_H}
  ${${GM_NAME}_HPP}
  ${${GM_NAME}_CPP}
  )

set_target_properties(${GM_NAME}_dev PROPERTIES
  COMPILE_DEFINITIONS gismo_dev_EXPORTS 
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules"
  )

set(gismo_dev_MODULES ${gismo_dev_MODULES} $<TARGET_OBJECTS:${GM_NAME}_dev>
    CACHE INTERNAL "G+Smo-dev modules" )
