<p align="center">
    <img src="https://github.com/thfabian/molec/blob/master/doc/logo/logo.png">
</p>
## Molecular Dynamics Framework [![Build Status](https://travis-ci.org/thfabian/molec.svg?branch=master)](https://travis-ci.org/thfabian/molec)
> * Group Name: **Yellow Lobster**
> * Group participants names: Michel Breyer, Carlo Del Don, Florian Frei, Fabian Thüring
> * Project Title: **MD Simulation**

# Group Project
Molecular Dynamics project for the course *How to Write Fast Numerical Code* (263-2300).

## TODO list
- [ ] sorting routine (at the end of force calculation, exploiting the already processed elements in the same order
- [ ] dummy implementation of cell list (no scalar replacement, arrays of arrays)

## Molecular Dynamics (MD)
Molecular dynamics simulation consists of the numerical, step-by-step, solution of the
classical equations of motion, which for a simple atomic system may be written:

![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cbegin%7Balign*%7D%20%26m_i%20%5Cddot%7B%5Cvec%7Bx%7D%7D_i%20%3D%20%5Cvec%7Bf%7D_i%5C%5C%26%5Cvec%7Bf%7D_i%20%3D%20-%5Cnabla_%7B%5Cvec%7Bx%7D_i%7D%5Cmathcal%7BU%7D%28%5Cvec%7Bx_1%7D%2C%5Cldots%2C%5Cvec%7Bx_N%7D%29%20%5Cend%7Balign*%7D)

For this purpose we need to be able to compute the forces ![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cvec%7Bf%7D_i) for each particle of the system, which is derived by the negative gradient of a potential field  ![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cmathcal%7BU%7D%28%5Cvec%7Bx%7D_1%2C%5Cldots%2C%5Cvec%7Bx%7D_N%29).

## Running the code
A configuration file can be passed to the executable in order to set some default simulation parameters. The configuration file has to have the following structure:
```{sh}
# Number of particles
N = 10000
# time step
dt = 0.005 
# number of simulation steps
Nstep = 100
# particle density (particles per unit volume)
rho = 1.2
# mass of each particle
mass = 1
# cutoff radius
Rcut = 2.5
# Lennard-Jones parameters
epsLJ = 1
sigLJ = 1
# particles initial disturbance (in [0,1) )
scaling = 0.05
```

## References 
 * [1] http://udel.edu/~arthij/MD.pdf
 * [2] http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf
 
