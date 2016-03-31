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

static int neighbors[27];
extern molec_CellList_Parameter_t cellList_parameter;
/**
 * Calculate the positive modulo between two integers, used for periodic BC
 */
MOLEC_INLINE int mod(int b, int m)
{
    return (b % m + m) % m;
}

/**
  * Transform cellNumber to respective 1D indices
  */
MOLEC_INLINE void molec_1to3trans(int cellNumber,int idx_x,int idx_y,int idx_z)
{
    idx_z=cellNumber / (cellList_parameter.N_x*cellList_parameter.N_y);
    idx_y=cellNumber / cellList_parameter.N_x - cellList_parameter.N_y*idx_z;
    idx_x=cellNumber - cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
}

/**
  * Transform 1D indices to respective cellNumber
  */
MOLEC_INLINE void molec_3to1trans(int idx_x,int idx_y,int idx_z,int cellNumber)
{
    cellNumber=idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
}

/**
  * Calculate neighbors
  */
MOLEC_INLINE void molec_get_neighbors(int cellNumber)
{
    int idx_x,idx_y,idx_z;
    int index=0;
    molec_1to3trans(cellNumber,idx_x,idx_y,idx_z);
    for(int d_z = -1; d_z <= 1; ++d_z)
    for(int d_y = -1; d_y <= 1; ++d_y)
    for(int d_x = -1; d_x <= 1; ++d_x)
    {
        int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
        int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
        int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);
        molec_3to1trans(n_idx_x,n_idx_y,n_idx_z,neighbors[index]);
        ++index;
    }

}


void molec_force_reference_dp(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
     cellList_parameter = molec_parameter->cellList;
    int molec_cellList[cellList_parameter.N][2*(molec_parameter->N/cellList_parameter.N)];

    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    int* index_array; //number of atoms in the cell
    MOLEC_MALLOC(index_array,sizeof(int)*cellList_parameter.N);
    memset(index_array, 0, sizeof(int)*cellList_parameter.N); //set number to zero in each cell

    for(int i = 0; i < N; ++i) //structure the celllist
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i] / cellList_parameter.c_x;
        int idx_y = y[i] / cellList_parameter.c_y;
        int idx_z = z[i] / cellList_parameter.c_z;

        // linear index of cell
        int cellNmbr;
        molec_3to1trans(idx_x,idx_y,idx_z,cellNmbr);

        molec_cellList[cellNmbr][index_array[cellNmbr]]=i;
        ++index_array[cellNmbr];
        //possibly copy the real values?
        //molec_cellList[idx][index_array[idx]][0]=x[i];
        //molec_cellList[idx][index_array[idx]][1]=y[i];
        //molec_cellList[idx][index_array[idx]][2]=z[i];
    }

    int cellNumber;
    //loop over all cells
    for(cellNumber=0;cellNumber<cellList_parameter.N;++cellNumber)
    {

        //get neighbors
        molec_get_neighbors(cellNumber);

        //get my index list
        int* myIndexList=molec_cellList[cellNumber];

        int nNumber;
        //loop over all neighbors
        for(nNumber=0;nNumber<27;++nNumber)
        {
        int mneighbor=neighbors[nNumber];

        //get corresponding list
        int* neighborIndexList=molec_cellList[mneighbor];

        //case one for different cells
        if(cellNumber<mneighbor)
        {
            //calculate interaction with loop over the two lists
            int k,i;
            for(k=0;k<index_array[cellNumber];++k)//first list
            {
                for(i=0;i<index_array[mneighbor];++i)//second list
                {


                }
             }

        }
        else if(cellNumber==mneighbor) //case two for the same cell
        {
            //calculate interaction with loop over the two lists
            int k,i;
            for(k=0;k<index_array[cellNumber];++k)//first list
            {
                for(i=k+1;i<index_array[mneighbor];++i)//second list
                {


                }
             }

        }

        }//end loop over neighbors


    }//end loop over cells

    *Epot = Epot_;
}
