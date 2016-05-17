<p align="center">
    <img src="https://github.com/thfabian/molec/blob/master/doc/logo/logo.png">
</p>
# Force implementations
The following list contains a description for each force routine implementation
## N^2 algorithm - *N2*
Base implementation, computational complexity of ![equation](https://latex.codecogs.com/png.latex?%5Cmathcal%7BO%7D%28N%5E2%29)

## CellList reference - *cell_ref*
Base implementation of celllist algorithm.

### Memory
* No aligned allocation
* No reuse of datastructure

### Algorithm
* iterate over each cell, compute cell neighbors
* for each particle *i* in cell, perform the following
    * for each cell neighbor, iterate over particles *j*
    * only if *i*<*j* compute the following
        * compute distances and forces
        * update force for particle *i* and *j*

### Optimizations
No explicit optimizations

# CellList v1 - *cell_v1*
Celllist algorithm with scalar replacement, arrays are memory aligned, discriminate to avoid double interaction on cell level (not on particle level like in *cell_ref*)
### Memory
* aligned allocation
* No reuse of datastructure

### Algorithm
* iterate over each cell *idx*, compute cell neighbors
* iterate over neighbor cell *n_idx* only if *idx* > *n_idx*
* for each particle *i* in cell *idx* and particle *j* in cell *n_idx*, perform the following
    * compute distances and forces
    * update force for particle *i* and *j*

### Optimizations
* Full scalar replacement

# CellList v2 - *cell_v2*
Celllist algorithm with scalar replacement, arrays are memory aligned, discriminate to avoid double interaction on cell level (not on particle level like in *cell_ref*), cell-neighbors are computed only once, distance computation improved
### Memory
* aligned allocation
* Some reuse of datastructure

### Algorithm
* iterate over each cell *idx*, compute cell neighbors
* iterate over neighbor cell *n_idx* only if *idx* > *n_idx*
* for each particle *i* in cell *idx* and particle *j* in cell *n_idx*, perform the following
    * compute distances and forces
    * update force for particle *i* and *j*

### Optimizations
* Full scalar replacement
* New *dist* function which does not need to compute *L/2*


# Quadrants - *q*
The datastructure over which the force calculation works is completely redesigned, in order to have contiguous memory access to particle *position*, *velocities* and *forces* inside each cell.

### Memory
The particles are stored *per-quadrant*, i.e. all particles lying inside the same cell are stored in the same quadrant struct. The computation of forces is performed between entire quadrants.

### Algorithm
* iterate over each cell *idx*
* iterate over neighbor cell *n_idx* only if *idx* > *n_idx*
* for each particle *i* in cell *idx* and particle *j* in cell *n_idx* perform force computation

### Optimizations
* Full scalar replacement
* Aligned and contiguous memory access
* Particles are likely to be in the correct order for next timestep as after force computation, particles are written back to *SOA* datastructure following the cell order

# Quadrants Ghost - *q_g*
One of the points wich had a big negative performance impact in the previous implementation is the branching inside the *dist* computation which needs to be performed in order to deal correctly with periodic boundary conditions.
The *Quadrant Ghost* implementation trades memory for efficientcy by copying particle coordinates lying in boundary cells to *ghost* cells with shifted coordinates.
Using this method allows to compute the distance between particles without branches.

### Memory
Using *Ghost Quadrants* has an impact in the memory used by the program, as some particle coordinates need to be copied multiple times. Memory aliasing is also introduced to allow superposition of force arrays between different (ghost and mirror) quadrants.

### Algorithm
* Generate data structure based on ghost quadrants (with internal cross references for coordinates that do not need to be shifted)
* Loop over all internal cells *idx*
* For each cell *idx* loop over negihbor cells *n_idx*
* Compute interaction between two quadrants without needing to check boundary conditions

### Optimizations
* Full scalar replacement
* Memory management allows saving branching

