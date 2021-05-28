cmake_minimum_required (VERSION 2.8)
project (dsmc)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
   set (CMAKE_BUILD_TYPE RELEASE)
endif ()

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT) 
	SET(CMAKE_INSTALL_PREFIX "/usr/local" ) 
ENDIF() 

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -fopenmp ")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -Wall -Wextra -Warray-temporaries \
   -fimplicit-none -fbacktrace -ffree-line-length-0 -Wtabs \
  -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan")

elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2  -qopenmp ")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -check all -fpe0 -warn -traceback \
  -debug extended")
endif ()


file(GLOB_RECURSE sources  
    src/allocateKernel.f90
    src/dataTypes.f90
    src/collisionKernel.f90
    src/coreKernel.f90
    src/particleListModule.f90
    src/sampleKernel.f90
    src/inputKernel.f90
    src/main.f90)
add_executable(dsmc ${sources})
