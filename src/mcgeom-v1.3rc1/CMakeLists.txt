project (MCGEOM)

cmake_minimum_required(VERSION 2.8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})


file (GLOB FILES src/*.c)
foreach (FILE ${FILES})
    set (mcgeomfiles ${mcgeomfiles} ${FILE})

#  to exclude a particular file, do this

#  get_filename_component(RELFILE ${FILE} NAME)
#  if (NOT (${RELFILE} MATCHES main) )
#    set (mcgeomfiles ${mcgeomfiles} ${FILE})
#  endif ()

endforeach ()


add_library(mcgeom ${mcgeomfiles})
include_directories(${PROJECT_SOURCE_DIR}/include)

if (CMAKE_INSTALL_PREFIX)
  if (APPEND_ARCHOS_TO_INSTALLPATH)
    install(TARGETS mcgeom DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${ARCHOS})
  else ()
    install(TARGETS mcgeom DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  endif ()
else ()
  if (APPEND_ARCHOS_TO_INSTALLPATH)
    install(TARGETS mcgeom DESTINATION ${PROJECT_SOURCE_DIR}/lib/${ARCHOS})
  else ()
    install(TARGETS mcgeom DESTINATION ${PROJECT_SOURCE_DIR}/lib)
  endif ()
endif ()


