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

include(ExternalProject)
include(CMakeParseArguments)

enable_testing()

if(NOT(CMAKE_MAJOR_VERSION LESS 3) AND NOT(CMAKE_MINOR_VERSION LESS 2))
    cmake_policy(SET CMP0054 OLD)
endif()

### Compiler options

# Set optimization flags (Release mode automatically adds -O3)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

if(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC")
    message(STATUS "Using Visual C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    add_definitions(-D_CRT_SECURE_NO_WARNINGS=1)

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "Using GNU C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
    
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
    message(STATUS "Using Clang C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")
    
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-implicit-function-declaration")
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "Using Intel C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")

elseif(${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
    message(STATUS "Using PGI C compiler")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")

endif()

### CMake options
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(MOLEC_PROFILING "Enable flags for proflinig (-fprofile-arcs)" OFF)
if(MOLEC_PROFILING)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
endif(MOLEC_PROFILING)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release" FORCE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
endif(NOT CMAKE_BUILD_TYPE)

option(MOLEC_UNITTEST "Build the UnitTests" ON)

option(MOLEC_TIMER_ENABLE "Enable timing framework" OFF)
option(MOLEC_TIMER_FORCE "Enable force timer" ${MOLEC_TIMER})
option(MOLEC_TIMER_INTEGRATOR "Enable integrator timer" ${MOLEC_TIMER})
option(MOLEC_TIMER_PERIODIC "Enable periodic timer" ${MOLEC_TIMER})
option(MOLEC_TIMER_SIMULATION "Enable simulation timer" ${MOLEC_TIMER})

if(MOLEC_TIMER_ENABLE)
    if(MOLEC_TIMER_FORCE)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_FORCE)
    endif()

    if(MOLEC_TIMER_INTEGRATOR)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_INTEGRATOR)
    endif()

    if(MOLEC_TIMER_PERIODIC)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_PERIODIC)
    endif()

    if(MOLEC_TIMER_SIMULATION)
        add_definitions(-DMOLEC_TIME)
        add_definitions(-DMOLEC_TIME_SIMULATION)
    endif()
else(MOLEC_TIMER_ENABLE)
    if(MOLEC_TIMER_FORCE OR MOLEC_TIMER_INTEGRATOR OR MOLEC_TIMER_PERIODIC OR MOLEC_TIMER_SIMULAITON)
        message(WARNING "MOLEC_TIMER_ENABLE is disabled")
    endif()
endif(MOLEC_TIMER_ENABLE)

option(MOLEC_CELLLIST_COUNT "Count missrate of cell-list interactions" OFF)
if(MOLEC_CELLLIST_COUNT)
    add_definitions(-DMOLEC_CELLLIST_COUNT_INTERACTION=1)
endif(MOLEC_CELLLIST_COUNT)

option(MOLEC_DUMP_COORDINATES "Writes the particle coordinates to a file" OFF)
if(MOLEC_DUMP_COORDINATES)
    add_definitions(-DMOLEC_DUMP_COORDINATES=1)
endif(MOLEC_DUMP_COORDINATES)


### Build argtable
set(ARGTABLE_INSTALL_DIR "${CMAKE_BINARY_DIR}/external/argtable2/install")
set(CMAKE_EXTERNAL_ARGUMENTS  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                              -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
                              -Wno-dev
                              -DCMAKE_INSTALL_PREFIX=${ARGTABLE_INSTALL_DIR}
                              -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
)

ExternalProject_Add(libargtable2
     URL             "${CMAKE_CURRENT_SOURCE_DIR}/external/argtable2"
     PREFIX          "external/argtable2"
     CMAKE_ARGS      ${CMAKE_EXTERNAL_ARGUMENTS}
 )
 message(STATUS "Building argtable2 from: ${CMAKE_CURRENT_SOURCE_DIR}/external/argtable2")

# Set the link and include directories
link_directories(${ARGTABLE_INSTALL_DIR}/lib)
include_directories(${ARGTABLE_INSTALL_DIR}/include)

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
include_directories(${PROJECT_SOURCE_DIR}/external)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)

if(MOLEC_UNITTEST)
    add_subdirectory(${PROJECT_SOURCE_DIR}/test)
endif(MOLEC_UNITTEST)

if(DOXYGEN_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
endif(DOXYGEN_FOUND)
