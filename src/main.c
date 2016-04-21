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

#include <molec/Force.h>
#include <molec/Integrator.h>
#include <molec/LoadConfig.h>
#include <molec/Simulation.h>
#include <molec/Timer.h>

int main(int argc, const char* argv[])
{
    srand(42);

    molec_load_parameters(argc, argv, 1);

    //molec_force_calculation force_calculation = &molec_force_N2_refrence;
    molec_force_calculation force_calculation = &molec_force_cellList_double_pointer;
    molec_force_integration force_integration = &molec_integrator_leapfrog_unroll_2;

    MOLEC_MEASUREMENT_INIT;

    MOLEC_MEASUREMENT_SIMULATION_START();
    molec_run_simulation(force_calculation, force_integration);
    MOLEC_MEASUREMENT_SIMULATION_STOP();


    molec_uint64_t force_cycles = 1;
    molec_uint64_t integrator_cycles = 1;

    force_cycles = MOLEC_MEASUREMENT_FORCE_GET_MEDIAN();
    integrator_cycles = MOLEC_MEASUREMENT_INTEGRATOR_GET_MEDIAN();

    printf("Simulation cycles: %llu\n", MOLEC_MEASUREMENT_SIMULATION_GET_MEDIAN());

    MOLEC_MEASUREMENT_FINISH;
}
;
