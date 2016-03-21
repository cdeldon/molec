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
#include <string.h>
#include <stdio.h>

molec_Loader_t* molec_loader = NULL;

void molec_load_parameters(const int argc, const char *argv[])
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

    FILE *parametersFile = fopen(molec_loader->filename, "r");
    if(parametersFile == NULL)
        molec_error("File %s does not exist\n", molec_loader->filename);

    char tag[256];
    char value[256];

    int tokens;
    while (!feof(parametersFile))
    {
        tokens = fscanf(parametersFile, "%s = %s",tag, value);
        printf("tag: <%s>, value: <%s>\n", tag, value);

        if(tokens == 2)
        {
            // Store value in 'value' char array into molec_parameter
            if (strcmp(tag,'N')==0)
            {
                molec_parameter->N = atoi(value);
            }
            if (strcmp(tag,"Nstep")==0)
            {
                molec_parameter->Nstep = atoi(value);
            }
            else
            {
                printf("Unrecongized parameter : \"%s\"\n", tag);
            }
        }
    }
    molec_cell_init();
    else
    {
        molec_error("Failed to read parameters file\n");
    }
}
