# Cmake version, project name, language
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(fortran-metis-interface Fortran)

# Set output paths for modules, archives, and executables
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release
# If not specified default to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

# Compiler flags for gfortran.
if (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(common "-std=f2018 -Wall")
    set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3 -march=native")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -O0 -g -C -fbacktrace -fdefault-real-8")
endif()

# Compiler flags for ifort.
if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(common "-std18 -warn all")
    set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${common} -O0 -g -C -traceback")
endif ()

# Compiler flags for Cray ftn
if(CMAKE_Fortran_COMPILER_ID MATCHES Cray)
  set(common " -h noomp")
  set(CMAKE_Fortran_FLAGS_RELEASE "${common} -O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "${common} -O0 -g")
endif()

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
endif ()
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")


# Set the METIS include and lib path (assuming standard install path)
# Edit the following two lines for a local METIS folder 
set(METIS_INCLUDE /usr/local/include)
set(METIS_LIB /usr/local/lib)

include_directories(${METIS_INCLUDE})
link_directories (${METIS_LIB})

# Library to archive (libfmetis.a)
add_library(fmetis src/lib/metis_interface.f90 src/lib/metis_enum.f90 src/lib/metis_io.f90 src/lib/metis_oo_interface.f90)
target_link_libraries(fmetis metis)

# Tests
enable_testing()
foreach(execid PartMeshNodal1 PartMeshNodal2 PartGraphRecursive1 PartGraphRecursive2 PartGraphKway NodeND )
    add_executable(test_${execid} src/tests/test_${execid}.f90)
    target_link_libraries(test_${execid} fmetis)
    add_test(test_${execid} bin/test_${execid})
endforeach()

# Examples
add_subdirectory(examples)

install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")
install(DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" DESTINATION "${CMAKE_INSTALL_PREFIX}")