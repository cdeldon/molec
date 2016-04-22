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

set(UNITTEST_SOURCE
    Unittest.c
    UnittestMain.c
    UnittestPeriodic.c
    UnittestTimer.c
    UnittestParser.c
    UnittestCellVector.c
    UnittestSort.c
    UnittestForce.c
    UnittestIntegrator.c
)

set(UNITTEST_HEADER
    Unittest.h
)

add_executable(molecTest ${UNITTEST_SOURCE} ${UNITTEST_HEADER} ${MOLEC_HEADERS})
target_link_libraries(molecTest molecCore)

# Link against math library (-lm)
if(NOT(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC"))
    target_link_libraries(molecTest m)
endif(NOT(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC"))

# Add CTest
if(WIN32)
    add_test(molecRunTests "${CMAKE_BINARY_DIR}/molecTest.exe")
else()
    add_test(molecRunTests "${CMAKE_BINARY_DIR}/molecTest")
endif()