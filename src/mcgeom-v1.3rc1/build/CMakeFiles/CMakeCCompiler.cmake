SET(CMAKE_C_COMPILER "/scratch/ymp/packages/gcc/gcc-4.5.1/bin/gcc")
SET(CMAKE_C_COMPILER_ARG1 "")
SET(CMAKE_C_COMPILER_ID "GNU")
SET(CMAKE_C_PLATFORM_ID "Linux")

SET(CMAKE_AR "/n/local_linux/bin/ar")
SET(CMAKE_RANLIB "/n/local_linux/bin/ranlib")
SET(CMAKE_LINKER "/n/local_linux/bin/ld")
SET(CMAKE_COMPILER_IS_GNUCC 1)
SET(CMAKE_C_COMPILER_LOADED 1)
SET(CMAKE_COMPILER_IS_MINGW )
SET(CMAKE_COMPILER_IS_CYGWIN )
IF(CMAKE_COMPILER_IS_CYGWIN)
  SET(CYGWIN 1)
  SET(UNIX 1)
ENDIF(CMAKE_COMPILER_IS_CYGWIN)

SET(CMAKE_C_COMPILER_ENV_VAR "CC")

IF(CMAKE_COMPILER_IS_MINGW)
  SET(MINGW 1)
ENDIF(CMAKE_COMPILER_IS_MINGW)
SET(CMAKE_C_COMPILER_ID_RUN 1)
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS c)
SET(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
SET(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
SET(CMAKE_C_SIZEOF_DATA_PTR "8")
SET(CMAKE_C_COMPILER_ABI "ELF")

IF(CMAKE_C_SIZEOF_DATA_PTR)
  SET(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
ENDIF(CMAKE_C_SIZEOF_DATA_PTR)

IF(CMAKE_C_COMPILER_ABI)
  SET(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
ENDIF(CMAKE_C_COMPILER_ABI)

SET(CMAKE_C_HAS_ISYSROOT "")


SET(CMAKE_C_IMPLICIT_LINK_LIBRARIES "c")
SET(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/home/quanb/local/lib64;/usr/lib64;/scratch/ymp/packages/gcc/gcc-4.5.1/lib/gcc/x86_64-unknown-linux-gnu/4.5.1;/scratch/ymp/packages/gcc/gcc-4.5.1/lib64;/lib64;/home/quanb/local/lib;/usr/lanl/lib;/scratch/ymp/packages/gcc/gcc-4.5.1/lib")
