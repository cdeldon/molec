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
#include <molec/Parameter.h>
#include <molec/CellVector.h>
#include <molec/Inline.h>

void molec_force_gonnet(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    //local alias for L
    //const Real L = molec_parameter->L;
    const Real Rcut = molec_parameter->Rcut;

    molec_CellList_Parameter_t cellList_parameter = molec_parameter->cellList;
    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    // for each particle compute cell index
    int* c_idx;
    MOLEC_MALLOC(c_idx, sizeof(int) * N);

    /* Note that the cell list will be traversed in the following form:
     *    for(z=0;z<cellList_parameter.N_z;++z)
     *      for(y=0;y<cellList_parameter.N_y;++y)
     *          for(x=0;x<cellList_parameter.N_x;++x)
     *
     * the fastest running index is x, while the slowest is z
     */
    for(int i = 0; i < N; ++i)
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i] / cellList_parameter.c_x;
        int idx_y = y[i] / cellList_parameter.c_y;
        int idx_z = z[i] / cellList_parameter.c_z;

        // linear index of cell
        int idx = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
        c_idx[i] = idx;
    }


    //list of indices of the particles in cell k
    int* cellList[cellList_parameter.N];
    //indices for numbering;
    int* index_array;
    int k;
    MOLEC_MALLOC(index_array,sizeof(int)*cellList_parameter.N);


    //allocate space for all pointers
    for(k=0;k<cellList_parameter.N;++k)
    {
        int* vec;
        MOLEC_MALLOC(vec,sizeof(int)*2*(N/cellList_parameter.N));
        cellList[k]=vec;
        index_array[k]=0;
    }

    int p;
    //add particles to indivitual cellList
    for(p=0;p<N;++p)
    {
        cellList[c_idx[p]][index_array[c_idx[p]]]=p;
        ++index_array[c_idx[p]];
    }

    // Loop over the cells
    for(int idx_z = 0; idx_z < cellList_parameter.N_z; ++idx_z)
    for(int idx_y = 0; idx_y < cellList_parameter.N_y; ++idx_y)
    for(int idx_x = 0; idx_x < cellList_parameter.N_x; ++idx_x)
    {
        // compute scalar cell index
        const int idx = idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);

        //index for the static lookup
        int index_for_norm_vector=0;

        // loop over neighbour cells
        for(int d_z = -1; d_z <= 1; ++d_z)
        for(int d_y = -1; d_y <= 1; ++d_y)
        for(int d_x = -1; d_x <= 1; ++d_x)
        {
            // compute cell index considering periodic BC
            int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
            int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
            int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);

            // linear index
            int n_idx = n_idx_x
                        + cellList_parameter.N_x
                              * (n_idx_y + cellList_parameter.N_y * n_idx_z);

            //start algorithmus of paper
            //get index array of the two cells we look at
            int* cAlpha=cellList[idx];
            int* cBeta=cellList[n_idx];

            //number of atoms to look at
            int number_of_relevant_atoms=index_array[idx]+index_array[n_idx]+2;

            //allocate space for the method
            Real* diff;
            MOLEC_MALLOC(diff,sizeof(Real)*number_of_relevant_atoms);
            Real* parts;
            MOLEC_MALLOC(parts,sizeof(Real)*number_of_relevant_atoms);

            //get vector
            const Real* norm_vec=molec_CellLookupTable[index_for_norm_vector];

            //index for diff and parts
            int nparts=0;
            int q;
            //pAlpha ind cAlpha
            for(q=0;q<index_array[idx]+1;++q)
            {
            diff[q]=norm_vec[0]*x[cAlpha[q]]+norm_vec[1]*y[cAlpha[q]]+norm_vec[2]*z[cAlpha[q]];
            parts[q]=cAlpha[q];
            }
            //pBeta ind cBeta
            for(q=0;q<index_array[n_idx]+1;++q)
            {
            diff[q+index_array[idx]+1]=norm_vec[0]*x[cBeta[q]]+norm_vec[1]*y[cBeta[q]]+norm_vec[2]*z[cBeta[q]]-Rcut;
            parts[q]=cBeta[q];
            }

            //sort diff and parts according to d


            //mainloop over sortet atoms
            for(q=0;q<number_of_relevant_atoms;++q)
            {

            }




            //next vector
            ++index_for_norm_vector;

        }//end loop over neighbors

    }//end loop over all cells


    // free memory
    MOLEC_FREE(c_idx);
    MOLEC_FREE(index_array);
    for(k=0;k<cellList_parameter.N;++k)
    {
        MOLEC_FREE(cellList[k]);
    }

    *Epot = Epot_;
}
