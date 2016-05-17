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

#define TINYTEST_MAIN
#include "Unittest.h"
#include <stdlib.h>

/**
 * Setup the Unittest environment
 *
 * For documentaion see https://github.com/thfabian/tinytest
 */
int main(int argc, char* argv[])
{
    srand(42);
    tinytest_init(argc, argv);

    // Register new testcases here
    REGISTER_TEST_CASE(molec_UnittestAlign);
    REGISTER_TEST_CASE(molec_UnittestTimer);
    REGISTER_TEST_CASE(molec_UnittestParser);
    REGISTER_TEST_CASE(molec_UnittestPeriodic);
    REGISTER_TEST_CASE(molec_UnittestCellVectorDirections);
    REGISTER_TEST_CASE(molec_UnittestSort);
    REGISTER_TEST_CASE(molec_UnittestCompare);
    REGISTER_TEST_CASE(molec_UnittestForce);
    REGISTER_TEST_CASE(molec_UnittestIntegrator);
    REGISTER_TEST_CASE(molec_UnittestGhost);

    int ret = tinytest_run();
    tinytest_free();

    return ret;
}
