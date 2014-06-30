project (MKDUAL)

cmake_minimum_required(VERSION 2.8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})


#
# MSTK Libraries and includes
#


if (NOT MSTK_INCLUDE_DIR)
  find_path(MSTK_INCLUDE_DIR MSTK.h
    HINTS ${MSTK_DIR}
          ${MSTK_DIR}/include
        )
    
   if (NOT MSTK_INCLUDE_DIR)
     message(FATAL_ERROR "Error: could not find MSTK include files\n"
     "Try defining MSTK_DIR: (-D MSTK_DIR:FILEPATH=<MSTK install prefix>)\n"
     "Or MSTK_INCLUDE_DIR")
   endif()
endif()

if (NOT USE_DEBUG_MSTK_LIBS)
  find_library(MSTK_LIBRARY mstk
    HINTS ${MSTK_LIB_DIR}
          ${MSTK_DIR}
          ${MSTK_DIR}/lib
          ${MSTK_DIR}/lib/${ARCHOS}
          )
else ()
  find_library(MSTK_LIBRARY mstk-d
    HINTS ${MSTK_LIB_DIR}
          ${MSTK_DIR}
          ${MSTK_DIR}/lib
          ${MSTK_DIR}/lib/${ARCHOS}
          )
endif ()

if (NOT MSTK_LIBRARY)
  message(FATAL_ERROR "Error: could not find MSTK library\n"
     "Try defining MSTK_DIR: (-D MSTK_DIR:FILEPATH=<MSTK install prefix>)\n"
     "Or MSTK_LIB_DIR")
endif()



#
# MCGEOM Libraries and includes
#


if (NOT MCGEOM_INCLUDE_DIR)
  find_path(MCGEOM_INCLUDE_DIR MCGeom.h
    HINTS ${MCGEOM_DIR}
          ${MCGEOM_DIR}/include
        )
    
   if (NOT MCGEOM_INCLUDE_DIR)
     message(FATAL_ERROR "Error: could not find MCGeom include files\n"
     "Try defining MCGeom_DIR: (-D MCGeom_DIR:FILEPATH=<MCGeom install prefix>)\n"
     "Or MCGeom_INCLUDE_DIR")
   endif()
endif()

if (NOT USE_DEBUG_MCGeom_LIBS)
  find_library(MCGEOM_LIBRARY mcgeom
    HINTS ${MCGEOM_LIB_DIR}
          ${MCGEOM_DIR}
          ${MCGEOM_DIR}/lib
          ${MCGEOM_DIR}/lib/${ARCHOS}
          )
else ()
  find_library(MCGEOM_LIBRARY mcgeom-d
    HINTS ${MCGEOM_LIB_DIR}
          ${MCGEOM_DIR}
          ${MCGEOM_DIR}/lib
          ${MCGEOM_DIR}/lib/${ARCHOS}
          )
endif ()

if (NOT MCGEOM_LIBRARY)
  message(FATAL_ERROR "Error: could not find MCGEOM library\n"
     "Try defining MCGEOM_DIR: (-D MCGEOM_DIR:FILEPATH=<MCGEOM install prefix>)\n"
     "Or MCGEOM_LIB_DIR")
endif()


if (ENABLE_ExodusII)

# Is there a way to look at the MSTK installation being pointed to 
# and tell if it is compiled with Exodus support

#
# Exodus II libraries
#

find_library(ExodusII_LIBRARY exoIIv2c
  HINTS ${ExodusII_LIB_DIR}
        ${ExodusII_DIR}
        ${ExodusII_DIR}/lib
        ${ExodusII_DIR}/lib/${ARCHOS}
        )

if (NOT ExodusII_LIBRARY)
  message(FATAL_ERROR "Error: could not find ExodusII library\n"
     "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII install prefix>)\n"
     "Or ExodusII_LIB_DIR")
endif()


#
# NetCDF libraries
#

find_library(NetCDF_LIBRARY netcdf
  HINTS ${NetCDF_LIB_DIR}
        ${NetCDF_DIR}
        ${NetCDF_DIR}/lib
        ${NetCDF_DIR}/lib/${ARCHOS}
        )

if (NOT NetCDF_LIBRARY)
  message(FATAL_ERROR "Error: could not find NetCDF library\n"
     "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF install prefix>)\n"
     "Or NetCDF_LIB_DIR")
endif()

endif()



file (GLOB FILES src/*.c)
foreach (FILE ${FILES})
  get_filename_component(RELFILE ${FILE} NAME)
  if (NOT (${RELFILE} MATCHES main) )
    set (mkduallibfiles ${mkduallibfiles} ${FILE})
  endif ()
endforeach ()


add_library(mkdual ${mkduallibfiles})

add_executable(run_mkdual src/main.c)
include_directories(${PROJECT_SOURCE_DIR}/include ${MSTK_INCLUDE_DIR} ${MCGEOM_INCLUDE_DIR})
target_link_libraries(run_mkdual mkdual ${MSTK_LIBRARY} ${MCGEOM_LIBRARY} ${ExodusII_LIBRARY} ${NetCDF_LIBRARY} m)

if (CMAKE_INSTALL_PREFIX)
  if (APPEND_ARCHOS_TO_INSTALLPATH)
    install(TARGETS run_mkdual DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/${ARCHOS})
    install(TARGETS mkdual DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${ARCHOS})
  else ()
    install(TARGETS run_mkdual DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
    install(TARGETS mkdual DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  endif ()
else ()
  if (APPEND_ARCHOS_TO_INSTALLPATH)
    install(TARGETS run_mkdual DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS})
    install(TARGETS mkdual DESTINATION ${PROJECT_SOURCE_DIR}/lib/${ARCHOS})
  else ()
    install(TARGETS run_mkdual DESTINATION ${PROJECT_SOURCE_DIR}/bin)
    install(TARGETS mkdual DESTINATION ${PROJECT_SOURCE_DIR}/lib)
  endif ()
endif ()


