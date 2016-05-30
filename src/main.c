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

#include <argtable2.h>
#include <molec/Force.h>
#include <molec/Integrator.h>
#include <molec/LoadConfig.h>
#include <molec/Periodic.h>
#include <molec/Simulation.h>
#include <molec/Timer.h>
#include <string.h>

molec_force_calculation arg_get_force_routine(const char* key)
{
    if (strcmp(key, "N2") == 0)
        return &molec_force_N2_refrence;
    else if(strcmp(key, "cell_ref") == 0)
        return &molec_force_cellList_reference;
    else if(strcmp(key, "cell_v1") == 0)
        return &molec_force_cellList_v1;
    else if(strcmp(key, "cell_v2") == 0)
        return &molec_force_cellList_v2;
    else if(strcmp(key, "knuth") == 0)
        return &molec_force_cellList_knuth;
    else if(strcmp(key, "q") == 0)
        return &molec_force_quadrant;
    else if(strcmp(key, "q_g") == 0)
        return &molec_force_quadrant_ghost;
    else if(strcmp(key, "q_g_u") == 0)
        return &molec_force_quadrant_ghost_unroll;
    else if(strcmp(key, "q_g_avx") == 0)
        return &molec_force_quadrant_ghost_avx;
    else if(strcmp(key, "q_g_fma") == 0)
        return &molec_force_quadrant_ghost_fma;
    else        
        molec_error("invalid parameter '%s' for option \"--force\"\n", key);
    return NULL;
}

molec_force_integration arg_get_integration_routine(const char* key)
{
    if(strcmp(key, "lf") == 0)
        return &molec_integrator_leapfrog_refrence;
    else if(strcmp(key, "lf2") == 0)
        return &molec_integrator_leapfrog_unroll_2;
    else if(strcmp(key, "lf4") == 0)
        return &molec_integrator_leapfrog_unroll_4;
    else if(strcmp(key, "lf8") == 0)
        return &molec_integrator_leapfrog_unroll_8;
    else if(strcmp(key, "lf_avx") == 0)
        return &molec_integrator_leapfrog_avx;
    else
        molec_error("invalid parameter '%s' for option \"--integrator\"\n", key);
    return NULL;
}


molec_periodic arg_get_periodic_routine(const char* key)
{
    if(strcmp(key, "ref") == 0)
        return &molec_periodic_refrence;
    else if(strcmp(key, "c4") == 0)
        return &molec_periodic_close4;
    else if (strcmp(key, "c") == 0)
        return &molec_periodic_close;
    else
        molec_error("invalid parameter '%s' for option \"--periodic\"\n", key);
    return NULL;
}

int main(int argc, char** argv)
{
    // force routine can appear at most once --> arg_str0
    struct arg_str* arg_force_routine
        = arg_str0("f", "force", "<string>",
                   "Specify the force subroutine.\n"
                   "                             - N2         N2 refrence\n"
                   "                             - cell_ref   Cell-list refrence\n"
                   "                             - cell_v1    Cell-list improvement 1\n"
                   "                             - cell_v2    Cell-list improvement 2\n"
                   "                             - knuth      Cell-list (Knuth)\n"
                   "                             - q          Quadrant\n"
                   "                             - q_g        Quadrant (ghost)\n"
                   "                             - q_g_u      Quadrant (ghost unrolled)\n"
                   "                             - q_g_avx    Quadrant (ghost-avx)\n"
                   "                             - q_g_fma    Quadrant (ghost-fma)");

    // integrator routine can appear at most once --> arg_str0
    struct arg_str* arg_integrator_routine
        = arg_str0("i", "integrator", "<string>",
                   "Specify the integrator subroutine.\n"
                   "                             - lf         Leap-frog (refrence)\n"
                   "                             - lf2        Leap-frog (unroll x2)\n"
                   "                             - lf4        Leap-frog (unroll x4)\n"
                   "                             - lf8        Leap-frog (unroll x8)\n"
                   "                             - lf_avx     Leap-frog (avx)");
    // periodic routine
    struct arg_str* arg_periodic_routine
        = arg_str0("p", "periodic", "<string>",
                   "Specify the periodic subroutine.\n"
                   "                             - ref        Refrence implementation\n"
                   "                             - c          With assumption\n"
                   "                             - c4         With assumption (unroll x4)");
    // parameter can appear at most once --> arg_file0
    struct arg_file* arg_parameters
        = arg_file0("c", "config", "<file>", "Path to to the configuration file.");

    // contains the number of desired particles
    struct arg_int* arg_desired_particles
        = arg_int0("n", "N", "<int>", "Set the number of particles.");

    // contains the number of desired timesteps
    struct arg_int* arg_desired_steps
        = arg_int0("s", "step", "<int>", "Set the number of timesteps.");

    // contains the number of desired timesteps
    struct arg_dbl* arg_desired_density
        = arg_dbl0("r", "rho", "<flaot>", "Set the particle density.");

    // help
    struct arg_lit* arg_help = arg_lit0("h", "help", "Print this help statement and exit.");
    // verbosity
    struct arg_int* arg_verb
        = arg_int0("v", "verbose", "<int>",
                   "Set verbosity level.\n"
                   "                             - 0:  Only print timers in the end\n"
                   "                             - 1:  Print settings (default)\n"
                   "                             - 2:  Full output (print energies every step)");

    // maximal number of errors = 20
    struct arg_end* end_struct = arg_end(20);

    void* argtable[] = {arg_force_routine,
                        arg_integrator_routine,
                        arg_periodic_routine,
                        arg_parameters,
                        arg_desired_particles,
                        arg_desired_steps,
                        arg_desired_density,
                        arg_help,
                        arg_verb,
                        end_struct};

    char* progname = "molec";
    // verify the argtable[] entries were allocated sucessfully
    if(arg_nullcheck(argtable) != 0)
    {
        // NULL entries were detected, some allocations must have failed
        printf("%s: insufficient memory\n", progname);
        exit(1);
    }

    // set any command line default values prior to parsing
    arg_force_routine->sval[0] = "cell_ref";
    arg_integrator_routine->sval[0] = "lf";
    arg_periodic_routine->sval[0] = "ref";
    arg_parameters->filename[0] = "";
    arg_desired_particles->ival[0] = 1000;
    arg_desired_steps->ival[0] = 100;
    arg_desired_density->dval[0] = 1.25f;
    arg_verb->ival[0] = 1;

    // parse argtable
    int nerrors = arg_parse(argc, argv, argtable);
    // special case: '--help | -h' takes precedence over error reporting
    if(arg_help->count > 0)
    {
        printf("Usage: %s", progname);
        arg_print_syntax(stdout, argtable, "\n");
        arg_print_glossary(stdout, argtable, "  %-25s %s\n");
        exit(0);
    }
    if(nerrors > 0)
    {
        // Display the error details contained in the arg_end struct
        arg_print_errors(stdout, end_struct, progname);
        printf("Try '%s --help' for more information.\n", progname);
        exit(1);
    }

    // normal case: take the command line options at face value
    molec_force_calculation force_calculation = arg_get_force_routine(arg_force_routine->sval[0]);
    molec_force_integration force_integration
        = arg_get_integration_routine(arg_integrator_routine->sval[0]);
    molec_periodic periodic = arg_get_periodic_routine(arg_periodic_routine->sval[0]);
    const char* config_file_name = arg_parameters->filename[0];
    const int desired_N = arg_desired_particles->ival[0];
    const float desired_rho = arg_desired_density->dval[0];
    molec_verbose = arg_verb->ival[0];

    srand(42);
    molec_load_parameters(config_file_name, 1, desired_N, desired_rho);

    // set the number of steps
    molec_parameter->Nstep = arg_desired_steps->ival[0];

    // print the used routines passed as argument
    if(molec_verbose)
    {
        printf("\n      ================ MOLEC - Settings ================\n\n");
        printf("      %-20s %10s\n", "Force routine:", arg_force_routine->sval[0]);
        printf("      %-20s %10s\n", "Integrator routine:", arg_integrator_routine->sval[0]);
        printf("      %-20s %10s\n", "Periodic routine:", arg_periodic_routine->sval[0]);
    }

    MOLEC_MEASUREMENT_INIT;

    MOLEC_MEASUREMENT_SIMULATION_START();
    molec_run_simulation(force_calculation, force_integration, periodic);
    MOLEC_MEASUREMENT_SIMULATION_STOP();

    printf("\n");
    MOLEC_MEASUREMENT_PRINT;

    // free memory 
    MOLEC_FREE(molec_parameter);
    MOLEC_MEASUREMENT_FINISH;
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
}
