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

#include "Unittest.h"
#include <molec/Timer.h>
#include <molec/LoadConfig.h>
#include <string.h>

/**
 * Test the parser
 */
TEST_CASE(molec_UnittestParser)
{
    const char unittestFile[] = "molec_UnittestParser.txt";

    // Create dummy file
    FILE * fp;
    fp = fopen (unittestFile, "w+");
    CHECK(fp != NULL)

    // Fill file
    fputs("N = 1234\n", fp);
    fputs("dt = 0.123\n", fp);
    fputs("Nstep = 123\n", fp);
    fputs("rho = 1.05\n", fp);
    fputs("mass = 3.0\n", fp);
    fputs("Rcut = 13.1\n", fp);
    fputs("epsLJ = 14.5\n", fp);
    fputs("sigLJ = 12.5\n", fp);
    fputs("scaling = 0.25\n", fp);

    fclose(fp);

    // Setup file
    int argc = 2;

    const char** argv;
    MOLEC_MALLOC(argv, argc * sizeof(char*));
    MOLEC_MALLOC(argv[1], sizeof(unittestFile));
    memcpy((void*) argv[1], unittestFile, sizeof(unittestFile));

    // Parse file
    molec_load_parameters(unittestFile , 0, 1234);

    CHECK_EQ_INTEGER(molec_parameter->N, 1234);
    CHECK_EQ_DOUBLE(molec_parameter->dt, (Real) 0.123);
    CHECK_EQ_INTEGER(molec_parameter->Nstep, 123);
    CHECK_EQ_DOUBLE(molec_parameter->rho, (Real) 1.05);
    CHECK_EQ_DOUBLE(molec_parameter->mass, (Real) 3.0);
    CHECK_EQ_DOUBLE(molec_parameter->Rcut, (Real) 13.1);
    CHECK_EQ_DOUBLE(molec_parameter->epsLJ, (Real) 14.5);
    CHECK_EQ_DOUBLE(molec_parameter->sigLJ, (Real) 12.5);
    CHECK_EQ_DOUBLE(molec_parameter->scaling, (Real) 0.25);

    remove("molec_UnittestParser.txt");
}
