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

#include <molec/main.h>

int main(int argc, const char* argv[])
{
    srand(42);

    molec_load_parameters(argc, argv, 1);

    //molec_force_calculation force_calculation = &molec_force_N2_refrence;
    molec_force_calculation force_calculation = &molec_force_cellList_for_swap;
    molec_force_integration force_integration = &molec_integrator_leapfrog_refrence;

    MOLEC_MEASUREMENT_INIT

    molec_run_simulation(force_calculation, force_integration);

    molec_uint64_t force_cycles = 1;
    molec_uint64_t integrator_cycles = 1;
    
#if MOLEC_MEASURE_FORCE
    force_cycles = molec_measurement_get_median(0);
    printf("Force cycles: %llu\n", force_cycles);
#endif
#if MOLEC_MEASURE_INTEGRATOR
    integrator_cycles = molec_measurement_get_median(1);
    printf("Integrator cycles: %llu\n", integrator_cycles);
#endif
    

    printf("Ratio: %f3.2\n", 1.0 * force_cycles / integrator_cycles);

    MOLEC_MEASUREMENT_FINISH
}
;
