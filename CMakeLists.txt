######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library. 
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012 - 2015 RICAM-Linz.
######################################################################

get_filename_component(GM_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#project(${GM_NAME})

## #################################################################
## Collect files
## #################################################################

aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_H  )
aux_source_directory     (${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR} ${GM_NAME}_HPP)

## #################################################################
## Add module
## #################################################################

if(GISMO_BUILD_LIB)

  #include_directories(${GISMO_INCLUDE_DIRS})

  add_library(${GM_NAME}_dev OBJECT
  ${${GM_NAME}_H}
  ${${GM_NAME}_HPP}
  #${${GM_NAME}_CPP}
  )

  set_target_properties(${GM_NAME}_dev PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  EXCLUDE_FROM_ALL 1
  EXCLUDE_FROM_DEFAULT_BUILD 1)

set(gismo_dev_MODULES ${gismo_dev_MODULES} $<TARGET_OBJECTS:${GM_NAME}_dev>
    CACHE INTERNAL "G+Smo modules" )
set (gismo_dev_SOURCES ${gismo_dev_SOURCES} ${${GM_NAME}_CPP}
     CACHE INTERNAL "gismo_dev source files")

endif()
