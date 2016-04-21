# molec - Molecular Dynamics Framework
#
#  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
#                      Michel Breyer  (mbreyer@student.ethz.ch)
#                      Florian Frei   (flofrei@student.ethz.ch)
#                      Fabian Thuring (thfabian@student.ethz.ch)
#
#  This file is distributed under the MIT Open Source License.
#  See LICENSE.txt for details.

cmake_minimum_required(VERSION 2.8)
project(molec C)

if(NOT(CMAKE_MAJOR_VERSION LESS 3) AND NOT(CMAKE_MINOR_VERSION LESS 2))
    cmake_policy(SET CMP0054 OLD)
endif()

### Compiler options

# Set optimization flags (Release mode automatically adds -O3)
if(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC")
    message(STATUS "Using Visual C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "Using GNU C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
    message(STATUS "Using Clang C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "Using Intel C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
    message(STATUS "Using PGI C compiler -- Broken")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")

endif()

### CMake options
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release" FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
endif(NOT CMAKE_BUILD_TYPE)

option(MOLEC_UNITTEST "Build the UnitTests" ON)

option(MOLEC_TIMER "Enable timing framework" OFF)
option(MOLEC_FORCE_TIMER "Enable force timer" ${MOLEC_TIMER})
option(MOLEC_INTEGRATOR_TIMER "Enable integrator timer" ${MOLEC_TIMER})
option(MOLEC_PERIODIC_TIMER "Enable periodic timer" ${MOLEC_TIMER})

if(MOLEC_TIMER)
    if(MOLEC_FORCE_TIMER)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_FORCE)
    endif()

    if(MOLEC_INTEGRATOR_TIMER)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_INTEGRATOR)
    endif()

    if(MOLEC_PERIODIC_TIMER)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_PERIODIC)
    endif()
endif(MOLEC_TIMER)

option(MOLEC_SINGLE_PRECISION "Compile for singe precision floating point (4 bytes)" ON)
if(MOLEC_SINGLE_PRECISION)
    add_definitions(-DMOLEC_SINGLE_PRECISION)
endif(MOLEC_SINGLE_PRECISION)

option(MOLEC_CELLLIST_COUNT "Count missrate of cell-list interactions" OFF)
if(MOLEC_CELLLIST_COUNT)
    add_definitions(-DMOLEC_CELLLIST_COUNT_INTERACTION=1)
endif(MOLEC_CELLLIST_COUNT)

option(MOLEC_DUMP_COORDINATES "Writes the particle coordinates to a file" OFF)
if(MOLEC_DUMP_COORDINATES)
    add_definitions(-DMOLEC_DUMP_COORDINATES=1)
endif(MOLEC_DUMP_COORDINATES)

### Doxygen Documentation
find_package(Doxygen)
option(MOLEC_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)"
       ${DOXYGEN_FOUND})
set(DOC_DIR ${PROJECT_SOURCE_DIR}/doc)

### Compilation
file(GLOB_RECURSE MOLEC_HEADERS
     ${CMAKE_CURRENT_SOURCE_DIR}/include
     ${CMAKE_CURRENT_SOURCE_DIR}/include/molec/*.h)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

if(MOLEC_UNITTEST)
    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif(MOLEC_UNITTEST)

if(DOXYGEN_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
endif(DOXYGEN_FOUND)
