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
