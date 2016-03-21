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

#include <molec/LoadConfig.h>
#include <stdio.h>
#include <string.h>

molec_Loader_t* molec_loader = NULL;

void molec_load_parameters(const int argc, const char* argv[])
{
    molec_parameter_init(1000);
    molec_loader = malloc(sizeof(molec_Loader_t));

    // Read the simulation parameters from the file passed as first argument
    // to the executable
    if(argc < 2)
        printf("Running simulation with default parameters\n");

    printf("Running simulation with parameters specified in <%s>\n", argv[1]);

    MOLEC_MALLOC(molec_loader->filename, MOLEC_FILENAME_MAX_LENGTH * sizeof(char));
    memcpy(molec_loader->filename, argv[1], MOLEC_FILENAME_MAX_LENGTH * sizeof(char));

    FILE* parametersFile = fopen(molec_loader->filename, "r");
    if(parametersFile == NULL)
        molec_error("File %s does not exist\n", molec_loader->filename);

    char tag[256];
    char value[256];

    int tokens;
    while(!feof(parametersFile))
    {
        tokens = fscanf(parametersFile, "%255s = %255s", tag, value);

        if(tokens == 2 && tag[0] != '#')
        {
            printf("tag: <%s>, value: <%s>\n", tag, value);
            // Store value in 'value' char array into molec_parameter
            if(tag[0] == 'N' && strlen(tag) == 1)
            {
                molec_parameter->N = atoi(value);
            }
            else if(strcmp(tag, "dt") == 0)
            {
                molec_parameter->dt = atof(value);
            }
            else if(strcmp(tag, "Nstep") == 0)
            {
                molec_parameter->Nstep = atoi(value);
            }
            else if(strcmp(tag, "L") == 0)
            {
                molec_parameter->L = atof(value);
            }
            else if(strcmp(tag, "mass") == 0)
            {
                molec_parameter->mass = atof(value);
            }
            else if(strcmp(tag, "Rcut") == 0)
            {
                molec_parameter->Rcut = atof(value);
                molec_parameter->Rcut2 = molec_parameter->Rcut * molec_parameter->Rcut;
            }
            else if(strcmp(tag, "epsLJ") == 0)
            {
                molec_parameter->epsLJ = atof(value);
            }
            else if(strcmp(tag, "sigLJ") == 0)
            {
                molec_parameter->sigLJ = atof(value);
            }
            else if(strcmp(tag, "scaling") == 0)
            {
                molec_parameter->scaling = atof(value);
            }
            else
            {
                printf("Unrecongized parameter : \"%s\"\n", tag);
            }
        }
    }
    molec_cell_init();
}
