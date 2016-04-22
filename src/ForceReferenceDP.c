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
#include <string.h>

//#define NofCells (molec_parameter->cellList.N)
//#define NofAtoms (molec_parameter->rho*2)

static int neighbors[27];
molec_CellList_Parameter_t cellList_parameter;
static int rho;
static int MaxNumber;
static int** molec_cellList;//[NofCells][(int)NofAtoms];

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
MOLEC_INLINE void molec_1to3trans(int cellNumber,int* idx_x,int* idx_y,int* idx_z)
{
    *idx_z=cellNumber / (cellList_parameter.N_x*cellList_parameter.N_y);
    *idx_y=cellNumber / cellList_parameter.N_x - cellList_parameter.N_y*(*idx_z);
    *idx_x=cellNumber - cellList_parameter.N_x * (*idx_y + cellList_parameter.N_y * *idx_z);
}

/**
  * Transform 1D indices to respective cellNumber
  */
MOLEC_INLINE void molec_3to1trans(int idx_x,int idx_y,int idx_z,int* cellNumber)
{
    *cellNumber=idx_x + cellList_parameter.N_x * (idx_y + cellList_parameter.N_y * idx_z);
}

/**
  * Calculate neighbors
  */
MOLEC_INLINE void molec_get_neighbors(int cellNumber)
{
    int idx_x,idx_y,idx_z;
    int index=0;
    molec_1to3trans(cellNumber,&idx_x,&idx_y,&idx_z);
    for(int d_z = -1; d_z <= 1; ++d_z)
    for(int d_y = -1; d_y <= 1; ++d_y)
    for(int d_x = -1; d_x <= 1; ++d_x)
    {
        int n_idx_z = mod(idx_z + d_z, cellList_parameter.N_z);
        int n_idx_y = mod(idx_y + d_y, cellList_parameter.N_y);
        int n_idx_x = mod(idx_x + d_x, cellList_parameter.N_x);
        molec_3to1trans(n_idx_x,n_idx_y,n_idx_z,&neighbors[index]);
        ++index;
    }

}

/**
 * Calculate distance between x and y taking periodic boundaries into account
 */
MOLEC_INLINE Real dist(Real x, Real y, Real L)
{
    Real r = x - y;
    if(r < -L / 2)
        r += L;
    else if(r > L / 2)
        r -= L;
    return r;
}

void molec_force_celllist_dp(molec_Simulation_SOA_t* sim, Real* Epot, const int N)
{
    //get molec parameter
    cellList_parameter = molec_parameter->cellList;
    //get fraction
    rho=molec_parameter->N/cellList_parameter.N;
    //assert(rho>0);
    //printf("%i \n",rho);
    MaxNumber=10*rho;
    int index;
    assert(molec_parameter);

    if(molec_cellList==NULL)
    {
        MOLEC_MALLOC(molec_cellList,cellList_parameter.N*sizeof(int*));
        for(index = 0;index<cellList_parameter.N;++index)
        {
        MOLEC_MALLOC(molec_cellList[index],MaxNumber*sizeof(int));
        }
    }
    //int** molec_cellList;
    //MOLEC_MALLOC(molec_cellList,cellList_parameter.N*2*rho*sizeof(int*));

    const Real sigLJ = molec_parameter->sigLJ;
    const Real epsLJ = molec_parameter->epsLJ;
    const Real L = molec_parameter->L;
    const Real Rcut2 = molec_parameter->Rcut2;

    // Local aliases
    const Real* x = sim->x;
    const Real* y = sim->y;
    const Real* z = sim->z;
    Real* f_x = sim->f_x;
    Real* f_y = sim->f_y;
    Real* f_z = sim->f_z;

    Real Epot_ = 0;

    //number of atoms in the cell[cellnumber]
    int* NumberOfAtomsInCell;
    MOLEC_MALLOC(NumberOfAtomsInCell,sizeof(int)*cellList_parameter.N);
    //set number to zero in each cell
    memset(NumberOfAtomsInCell, 0, sizeof(int)*cellList_parameter.N);

    int* CellIndexOfAtom;//storing corresponding cellnumber
    MOLEC_MALLOC(CellIndexOfAtom,sizeof(int)*molec_parameter->N);


    for(int i = 0; i < N; ++i) //structure the celllist
    {
        // linear one dimensional index of cell associated to i-th particle
        int idx_x = x[i] / cellList_parameter.c_x;
        int idx_y = y[i] / cellList_parameter.c_y;
        int idx_z = z[i] / cellList_parameter.c_z;

        //set forces to zero
        f_x[i] = f_y[i] = f_z[i] = 0.0;

        // linear index of cell
        int cellNmbr;
        molec_3to1trans(idx_x,idx_y,idx_z,&cellNmbr);
        CellIndexOfAtom[i]=cellNmbr;

        molec_cellList[cellNmbr][NumberOfAtomsInCell[cellNmbr]]=i; //segfault?
        NumberOfAtomsInCell[cellNmbr]+=1;
        //possibly copy the real values?
        //molec_cellList[idx][index_array[idx]][0]=x[i];
        //molec_cellList[idx][index_array[idx]][1]=y[i];
        //molec_cellList[idx][index_array[idx]][2]=z[i];
        assert(NumberOfAtomsInCell[cellNmbr]<=MaxNumber);
    }

    int CellNumberOne;
    //loop over all cells
    for(CellNumberOne=0;CellNumberOne<cellList_parameter.N;++CellNumberOne)
    {

        //get neighbors
        molec_get_neighbors(CellNumberOne);

        //get my index list
        int* IndexArrayOfAtomsInCellOne=molec_cellList[CellNumberOne];

        int NeighborNumber;
        //loop over all neighbors
        for(NeighborNumber=0;NeighborNumber<27;++NeighborNumber)
        {
        int CellNumberTwo=neighbors[NeighborNumber];

        //get corresponding list
        int* IndexArrayOfAtomsInCellTwo=molec_cellList[CellNumberTwo];

        //case one for different cells
        if(CellNumberOne<CellNumberTwo)
        {
            //calculate interaction with loop over the two lists
            int k,i;
            for(k=0;k<NumberOfAtomsInCell[CellNumberOne];++k)//first list
            {

                //index one
                int IndexOfAtomOne=IndexArrayOfAtomsInCellOne[k];

                const Real xi = x[IndexOfAtomOne];
                const Real yi = y[IndexOfAtomOne];
                const Real zi = z[IndexOfAtomOne];

                Real f_xi = f_x[IndexOfAtomOne];
                Real f_yi = f_y[IndexOfAtomOne];
                Real f_zi = f_z[IndexOfAtomOne];

                for(i=0;i<NumberOfAtomsInCell[CellNumberTwo];++i)//second list
                {

                    //index two
                    int IndexOfAtomTwo=IndexArrayOfAtomsInCellTwo[i];

                    const Real xij = dist(xi, x[IndexOfAtomTwo], L);
                    const Real yij = dist(yi, y[IndexOfAtomTwo], L);
                    const Real zij = dist(zi, z[IndexOfAtomTwo], L);

                    const Real r2 = xij * xij + yij * yij + zij * zij;

                    if(r2 < Rcut2)
                    {
                        // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                        const Real s2 = (sigLJ * sigLJ) / r2;
                        const Real s6 = s2 * s2 * s2;

                        Epot_ += 4 * epsLJ * (s6 * s6 - s6);

                        const Real fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                        f_xi += fr * xij;
                        f_yi += fr * yij;
                        f_zi += fr * zij;

                        f_x[IndexOfAtomTwo] -= fr * xij;
                        f_y[IndexOfAtomTwo] -= fr * yij;
                        f_z[IndexOfAtomTwo] -= fr * zij;
                    }
                }//end second list

                f_x[IndexOfAtomOne] = f_xi;
                f_y[IndexOfAtomOne] = f_yi;
                f_z[IndexOfAtomOne] = f_zi;

             }//end first list

        }
        else if(CellNumberOne==CellNumberTwo) //case two for the same cell
        {
            //calculate interaction with loop over the two lists
            int k,i;
            for(k=0;k<NumberOfAtomsInCell[CellNumberOne];++k)//first list
            {
                //index one
                int IndexOfAtomOne=IndexArrayOfAtomsInCellOne[k];

                const Real xi = x[IndexOfAtomOne];
                const Real yi = y[IndexOfAtomOne];
                const Real zi = z[IndexOfAtomOne];

                Real f_xi = f_x[IndexOfAtomOne];
                Real f_yi = f_y[IndexOfAtomOne];
                Real f_zi = f_z[IndexOfAtomOne];

                for(i=k+1;i<NumberOfAtomsInCell[CellNumberTwo];++i)//second list
                {
                    //index two
                    int IndexOfAtomTwo=IndexArrayOfAtomsInCellTwo[i];

                    const Real xij = dist(xi, x[IndexOfAtomTwo], L);
                    const Real yij = dist(yi, y[IndexOfAtomTwo], L);
                    const Real zij = dist(zi, z[IndexOfAtomTwo], L);

                    const Real r2 = xij * xij + yij * yij + zij * zij;

                    if(r2 < Rcut2)
                    {
                        // V(s) = 4 * eps * (s^12 - s^6) with  s = sig/r
                        const Real s2 = (sigLJ * sigLJ) / r2;
                        const Real s6 = s2 * s2 * s2;

                        Epot_ += 4 * epsLJ * (s6 * s6 - s6);

                        const Real fr = 24 * epsLJ / r2 * (2 * s6 * s6 - s6);

                        f_xi += fr * xij;
                        f_yi += fr * yij;
                        f_zi += fr * zij;

                        f_x[IndexOfAtomTwo] -= fr * xij;
                        f_y[IndexOfAtomTwo] -= fr * yij;
                        f_z[IndexOfAtomTwo] -= fr * zij;

                    }//end if

                }//end second list

                f_x[IndexOfAtomOne] = f_xi;
                f_y[IndexOfAtomOne] = f_yi;
                f_z[IndexOfAtomOne] = f_zi;

                }//end first list

            }//end else

        }//end neighbors

        }//end loop cells

    *Epot = Epot_;
}
