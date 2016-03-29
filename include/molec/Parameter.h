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

#ifndef MOLEC_PARAMETER_H
#define MOLEC_PARAMETER_H

#include <molec/Common.h>
#include <molec/CellListParam.h>

/**
 * @brief Parameter of the simulation
 *
 * To globally access the paramters use the pointer @c molec_parameter.
 * @code{.c}
 *     double dt = molec_parameter->dt; // local alias
 * @endcode
 */
typedef struct molec_Parameter
{
    /** Number of atoms */
    int N;

    /** Number of simulation steps */
    int Nstep;

    /** Density of particles per unit volume */
    Real rho;

    /** Extend of the bounding box */
    Real L;

    /** Time step */
    Real dt;

    /** Mass of the atoms */
    Real mass;

    /** Cut-off radius */
    Real Rcut;

    /** Cut-off radius squared */
    Real Rcut2;

    /** Perturbation [0, 1) of the regular grid in the initial positions */
    Real scaling;
    
    /** Lennard-Jones parameter: epsilon */
    Real epsLJ;

    /** Lennard-Jones parameter: sigma */
    Real sigLJ;

    /** Parameters of the cell list */
    molec_CellList_Parameter_t cellList;

} molec_Parameter_t;

/**
 * Global access to the parameters
 */
extern molec_Parameter_t* molec_parameter;

/**
 * @brief Initialize the parameter struct
 *
 * Allocate the paramter pointer @c molec_parameter and set default values
 */
void molec_parameter_init(int N);

/**
 * @brief Prints the simulation paramteres on the terminal
 */
void molec_print_parameters();

#endif

