/*                   _           
 *   _ __ ___   ___ | | ___  ___ 
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__ 
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *          
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (frofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 * 
 *  This file is distributed under the MIT Open Source License. 
 *  See LICENSE.txt for details.
 */

#ifndef MOLEC_INITIAL_CONDITION_H
#define MOLEC_INITIAL_CONDITION_H
 
#include <molec/Simulation.h>

/**
 * Initialize the simulation (position, velocity and force) by applying initial condition
 *
 * The atoms will be placed on a regular grid and randomly offset. The parameter 
 * @c molec_parameter->scaling will make the particles do not overlap.
 * 
 * @param simulation    Simulation struct holding the position, velocity and force arrays
 */
void molec_set_initial_condition(molec_Simulation_SOA_t* simulation);

#endif
