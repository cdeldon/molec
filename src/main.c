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
#include <molec/Periodic.h>
#include <molec/LoadConfig.h>
#include <molec/Simulation.h>
#include <molec/Timer.h>
#include <argtable2.h>
#include <string.h>

molec_force_calculation arg_get_force_routine(const char* key)
{
    if(strcmp(key, "N2") == 0)
        return &molec_force_N2_refrence;
    else if(strcmp(key, "cell") == 0)
        return &molec_force_cellList;
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
    else
        molec_error("invalid parameter '%s' for option \"--integrator\"\n", key);
    return NULL;
}


molec_periodic arg_get_periodic_routine(const char* key)
{
    if(strcmp(key, "ref") == 0)
        return &molec_periodic_refrence;
    else
        molec_error("invalid parameter '%s' for option \"--periodic\"\n", key);
    return NULL;
}

int main(int argc, char** argv)
{
    // force routine can appear at most once --> arg_str0
    struct arg_str* arg_force_routine = arg_str0(
        "f", "force", "string", "force routine specification: "
                                "[ N2 | cell | cellfor | cellswap | celldp1 | celldp2 ]");
    // integrator routine can appear at most once --> arg_str0
    struct arg_str* arg_integrator_routine
        = arg_str0("i", "integrator", "string", "integrator routine specification: "
                                                "[ lf | lf2 ]");
    // periodic routine
    struct arg_str* arg_periodic_routine
        = arg_str0("p", "periodic", "string", "periodic routine specification: "
                                              "[ ref ]");
    // parameter can appear at most once --> arg_file0
    struct arg_file* arg_parameters
        = arg_file0("c", "config", "string", "path to configuration file");

    // contains the number of desired particles
    struct arg_int* arg_desired_particles = arg_int0("n", "N", "int", "number of particles");

    // help
    struct arg_lit* arg_help = arg_lit0("h", "help", "prints help");
    // verbosity
    struct arg_int* arg_verb
        = arg_int0("v", "verbose", "int", "verbose output: "
                                          "[ 0: silent, 1: settings, 2: full ]");

    // maximal number of errors = 20
    struct arg_end* end_struct = arg_end(20);

    void* argtable[] = {arg_force_routine,
                        arg_integrator_routine,
                        arg_periodic_routine,
                        arg_parameters,
                        arg_desired_particles,
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
    arg_force_routine->sval[0] = "cell";
    arg_integrator_routine->sval[0] = "lf";
    arg_periodic_routine->sval[0] = "ref";
    arg_parameters->filename[0] = "";
    arg_desired_particles->ival[0] = -1;
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
    molec_verbose = arg_verb->ival[0];

    srand(42);
    molec_load_parameters(config_file_name, 1, desired_N);

    // print the used routines passed as argument
    if(molec_verbose)
    {
        printf("\n      ================ MOLEC - Settings ================\n\n");
        printf("      %-20s %10s\n", "Force routine:",arg_force_routine->sval[0]);
        printf("      %-20s %10s\n", "Integrator routine:", arg_integrator_routine->sval[0]);
        printf("      %-20s %10s\n", "Periodic routine:", arg_periodic_routine->sval[0]);
    }

    MOLEC_MEASUREMENT_INIT;

    MOLEC_MEASUREMENT_SIMULATION_START();
    molec_run_simulation(force_calculation, force_integration, periodic);
    MOLEC_MEASUREMENT_SIMULATION_STOP();


    MOLEC_MEASUREMENT_PRINT;

    MOLEC_FREE(molec_parameter);

    MOLEC_MEASUREMENT_FINISH;

    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
}
