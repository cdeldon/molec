/*                   _
 *   _ __ ___   ___ | | ___  ___
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (flofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 *
 *  This file is distributed under the MIT Open Source License.
 *  See LICENSE.txt for details.
 */

#ifndef MOLEC_MAIN_H
#define MOLEC_MAIN_H

#include <molec/Simulation.h>
#include <molec/Force.h>
#include <molec/Timer.h>
#include <molec/Integrator.h>
#include <molec/LoadConfig.h>
#include <stdlib.h>

#ifdef MOLEC_MEASURE_TIME

    // maximal number of timers (has to be >= than number of defined (0 or 1) timers)
    #define NUM_TIMERS               3

    #define MOLEC_MEASUREMENT_INIT   molec_measurement_init(NUM_TIMERS);
    #define MOLEC_MEASUREMENT_FINISH molec_measurement_finish();


    // timer name, if defined as 1, timing will happen, else nothing will be performed
    #define MOLEC_MEASURE_FORCE      1
    #define MOLEC_MEASURE_INTEGRATOR 1
    #define MOLEC_MEASURE_PERIODIC   1

    #if MOLEC_MEASURE_FORCE
        #define MOLEC_MEASUREMENT_FORCE_START molec_measurement_start(0);
        #define MOLEC_MEASUREMENT_FORCE_STOP  molec_measurement_stop(0);
    #else
        #define MOLEC_MEASUREMENT_FORCE_START
        #define MOLEC_MEASUREMENT_FORCE_STOP
    #endif

    #if MOLEC_MEASURE_INTEGRATOR
        #define MOLEC_MEASUREMENT_INTEGRATOR_START molec_measurement_start(1);
        #define MOLEC_MEASUREMENT_INTEGRATOR_STOP  molec_measurement_stop(1);
    #else
        #define MOLEC_MEASUREMENT_INTEGRATOR_START
        #define MOLEC_MEASUREMENT_INTEGRATOR_STOP
    #endif

    #if MOLEC_MEASURE_PERIODIC
        #define MOLEC_MEASUREMENT_PERIODIC_START molec_measurement_start(2);
        #define MOLEC_MEASUREMENT_PERIODIC_STOP  molec_measurement_stop(2);
    #else
        #define MOLEC_MEASUREMENT_PERIODIC_START
        #define MOLEC_MEASUREMENT_PERIODIC_STOP
    #endif

#else
    #define MOLEC_MEASUREMENT_INIT
    #define MOLEC_MEASUREMENT_FINISH

    #define MOLEC_MEASURE_FORCE                     0
    #define MOLEC_MEASURE_INTEGRATOR                0
    #define MOLEC_MEASURE_PERIODIC                  0

    #define MOLEC_MEASUREMENT_FORCE_START
    #define MOLEC_MEASUREMENT_FORCE_STOP

    #define MOLEC_MEASUREMENT_INTEGRATOR_START
    #define MOLEC_MEASUREMENT_INTEGRATOR_STOP

    #define MOLEC_MEASUREMENT_PERIODIC_START
    #define MOLEC_MEASUREMENT_PERIODIC_STOP

#endif // MOLEC_MEASURE_TIME


#endif // MOLEC_MAIN_H
