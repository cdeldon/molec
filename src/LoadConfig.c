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

    char *buffer = malloc(MOLEC_LINE_MAX_LENGTH * sizeof(char));
    char *tag = malloc(MOLEC_LINE_MAX_LENGTH * sizeof(char));
    char *value = malloc(MOLEC_LINE_MAX_LENGTH * sizeof(char));

    char *delimPtr = malloc(MOLEC_LINE_MAX_LENGTH * sizeof(char));
    char *endPtr = malloc(MOLEC_LINE_MAX_LENGTH * sizeof(char));

    while (fscanf(parametersFile, "%s", buffer) == 1)
    {
        // Process line buffer
        memcpy(tag, buffer, MOLEC_LINE_MAX_LENGTH);
        if (buffer[0]!='#' && fgets(buffer, MOLEC_LINE_MAX_LENGTH, parametersFile) != NULL)
        {
            delimPtr = strstr(buffer,"=");
            endPtr = strstr(delimPtr,";");

            memcpy(value,delimPtr+2 ,endPtr-delimPtr-2);

            // Store value in 'value' char array into molec_parameter
            if (strcmp(tag,"N")==0)
            {
                molec_parameter->N = atoi(value);
            }
            else
            {
                printf("Unrecongized parameter : \"%s\"\n", tag);
            }
        }
    }
    if (feof(parametersFile))
    {
        // Regenerate cell list parameters
        molec_cell_init();
    }
    else
    {
        molec_error("Failed to read parameters file\n");
    }
}
