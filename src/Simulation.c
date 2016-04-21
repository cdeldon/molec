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
#include <molec/InitialCondition.h>
#include <molec/Parameter.h>
#include <molec/Periodic.h>
#include <molec/Simulation.h>
#include <molec/Dump.h>
#include <molec/Timer.h>
#include <molec/main.h>

#include <stdlib.h>
#include <stdio.h>

molec_Simulation_SOA_t* molec_init_simulation_SOA()
{
    if(molec_parameter == NULL)
        molec_error("molec_parameter is nullptr\n");

    const int N = molec_parameter->N;

    molec_Simulation_SOA_t* simulation = malloc(sizeof(molec_Simulation_SOA_t));

    MOLEC_MALLOC(simulation->x, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->y, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->z, sizeof(Real) * N);

    MOLEC_MALLOC(simulation->v_x, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->v_y, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->v_z, sizeof(Real) * N);

    MOLEC_MALLOC(simulation->f_x, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->f_y, sizeof(Real) * N);
    MOLEC_MALLOC(simulation->f_z, sizeof(Real) * N);

    return simulation;
}

void molec_free_simulation_SOA(molec_Simulation_SOA_t* simulation)
{
    MOLEC_FREE(simulation->x);
    MOLEC_FREE(simulation->y);
    MOLEC_FREE(simulation->z);

    MOLEC_FREE(simulation->v_x);
    MOLEC_FREE(simulation->v_y);
    MOLEC_FREE(simulation->v_z);

    MOLEC_FREE(simulation->f_x);
    MOLEC_FREE(simulation->f_y);
    MOLEC_FREE(simulation->f_z);

    MOLEC_FREE(simulation);
}

void molec_run_simulation(void (*molec_compute_force)( molec_Simulation_SOA_t*, Real*, int),
                          void (*molec_force_integration)(Real*, Real*, const Real*, Real*, const int))
{
    if(MOLEC_DUMP_COORDINATES)
    {
        dump_file_name = "molec_simulation_dump.xyz";
        molec_dump_file = fopen(dump_file_name, "w");
        if(molec_dump_file == NULL)
            molec_error("Unable to open file %s to dump particle coordinates\n", dump_file_name);
    }
    // Set parameters
    if(molec_parameter == NULL)
        molec_error("molec_parameter is nullptr\n");

    // Local alias
    const int N = molec_parameter->N;
    const int Nstep = molec_parameter->Nstep;

    molec_Simulation_SOA_t* sim = molec_init_simulation_SOA();;

    // Set initial conditions
    molec_set_initial_condition(sim);

    if(MOLEC_DUMP_COORDINATES)
    {
        // print the number of atoms
        fprintf(molec_dump_file, "%d\n", N);
    }

    // Run sim
    Real Ekin_x = 0.0, Ekin_y = 0.0, Ekin_z = 0.0;
    Real Epot = 0.0;

    printf("\n      ================ MOLEC - Simulation steps ================\n\n");
    printf("%10s\t%15s\t%15s\t%15s\n", "Step", "Ekin", "Epot", "Etot");
    for(int n = 1; n <= Nstep; ++n)
    {
        if(MOLEC_DUMP_COORDINATES)
        {
            molec_dump_coordinates(sim, N);
        }
        Ekin_x = Ekin_y = Ekin_z = 0.0;
        Epot = 0.0;

        // 1. Compute force

        MOLEC_MEASUREMENT_FORCE_START

        molec_compute_force(sim, &Epot, N);
        
        MOLEC_MEASUREMENT_FORCE_STOP

        // 2. Integrate ODE
        MOLEC_MEASUREMENT_INTEGRATOR_START
        molec_force_integration(sim->x, sim->v_x, sim->f_x, &Ekin_x, N);
        molec_force_integration(sim->y, sim->v_y, sim->f_y, &Ekin_y, N);
        molec_force_integration(sim->z, sim->v_z, sim->f_z, &Ekin_z, N);
        MOLEC_MEASUREMENT_INTEGRATOR_STOP

        // 3. Apply periodic boundary conditions
        molec_periodic_refrence(sim->x, N);
        molec_periodic_refrence(sim->y, N);
        molec_periodic_refrence(sim->z, N);

        // 4. Report result
        Real Ekin = Ekin_x + Ekin_y + Ekin_z;
        Real Etot = Ekin + Epot;
        printf("%10i\t%15.6f\t%15.6f\t%15.6f\n", n, Ekin, Epot, Etot);

        /* molec_print_simulation_SOA(sim); */
    }

    // Free memory
    molec_free_simulation_SOA(sim);
    MOLEC_FREE(molec_parameter);

    if(MOLEC_DUMP_COORDINATES)
    {
        fclose(molec_dump_file);
    }
}

void molec_print_simulation_SOA(const molec_Simulation_SOA_t* sim)
{
    const int N = molec_parameter->N;
    for(int i = 0; i < N; ++i)
        printf(" (%f, %f, %f)\t(%f, %f, %f)\n", sim->x[i], sim->y[i],
               sim->z[i], sim->v_x[i], sim->v_y[i], sim->v_z[i]);
}
