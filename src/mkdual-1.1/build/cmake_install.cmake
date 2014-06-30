# Install script for directory: /home/quanb/arctic_polyGen/mkdual-1.1

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/quanb/arctic_polyGen/mkdual-1.1")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux" TYPE EXECUTABLE FILES "/home/quanb/arctic_polyGen/mkdual-1.1/build/run_mkdual")
  IF(EXISTS "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/n/local_linux/bin/strip" "$ENV{DESTDIR}/home/quanb/arctic_polyGen/mkdual-1.1/bin/x86_64_Linux/run_mkdual")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "/home/quanb/arctic_polyGen/mkdual-1.1/lib/x86_64_Linux" TYPE STATIC_LIBRARY FILES "/home/quanb/arctic_polyGen/mkdual-1.1/build/libmkdual.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/home/quanb/arctic_polyGen/mkdual-1.1/build/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/home/quanb/arctic_polyGen/mkdual-1.1/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
