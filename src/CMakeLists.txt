# molec - Molecular Dynamics Framework
#
#  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
#                      Michel Breyer  (mbreyer@student.ethz.ch)
#                      Florian Frei   (frofrei@student.ethz.ch)
#                      Fabian Thuring (thfabian@student.ethz.ch)
# 
#  This file is distributed under the MIT Open Source License. 
#  See LICENSE.txt for details.

cmake_minimum_required(VERSION 2.8)

set(MOLEC_SOURCE
    Common.c
    Force.c
    Integrator.c
    InitialCondition.c
    main.c
    Parameter.c
    Periodic.c
    Simulation.c
    Timer.c
)

add_executable(molec ${MOLEC_SOURCE} ${MOLEC_HEADERS})

# Link against math library (-lm)
if(NOT(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC"))
    target_link_libraries(molec m)
endif(NOT(${CMAKE_C_COMPILER_ID} STREQUAL "MSVC"))

