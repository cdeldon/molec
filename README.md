<p align="center">
    <img src="/doc/logo/logo.png">
</p>
## Molecular Dynamics Framework [![Build Status](https://travis-ci.org/thfabian/molec.svg?branch=master)](https://travis-ci.org/thfabian/molec) [![Build Status](https://ci.appveyor.com/api/projects/status/qw7cek485ppr003w?svg=true)](https://ci.appveyor.com/project/thfabian/molec/branch/master)

> * Group Name: **Yellow Lobster**
> * Group participants names: Michel Breyer, Carlo Del Don, Florian Frei, Fabian Thüring
> * Project Title: **MD Simulation**

# Group Project
Molecular Dynamics project for the course *How to Write Fast Numerical Code* (263-2300).

## Molecular Dynamics (MD)
The molecular dynamics simulation consists of the numerical, step-by-step, solution of the
classical equations of motion, which for a simple atomic system may be written as:

![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cbegin%7Balign*%7D%20%26m_i%20%5Cddot%7B%5Cvec%7Bx%7D%7D_i%20%3D%20%5Cvec%7Bf%7D_i%5C%5C%26%5Cvec%7Bf%7D_i%20%3D%20-%5Cnabla_%7B%5Cvec%7Bx%7D_i%7D%5Cmathcal%7BU%7D%28%5Cvec%7Bx_1%7D%2C%5Cldots%2C%5Cvec%7Bx_N%7D%29%20%5Cend%7Balign*%7D)

For this purpose we need to be able to compute the forces ![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cvec%7Bf%7D_i) for each particle in the system. The force can be derived by the negative gradient of a potential field  ![equation](https://latex.codecogs.com/png.latex?%5Clarge%20%5Cmathcal%7BU%7D%28%5Cvec%7Bx%7D_1%2C%5Cldots%2C%5Cvec%7Bx%7D_N%29).

In general, this kind of problems are computationally expensive for large number of particles. The naive algorithm would require the computation of ![equation](https://latex.codecogs.com/png.latex?N%28N-1%29/2) interactions, where ![equation](https://latex.codecogs.com/png.latex?N) is the number of particles in the system.

Note, that the asymptotic complexity of such a problem is quadratic in the number of particles ![equation](https://latex.codecogs.com/png.latex?N).

Introducing a cut-off radius ![equation](https://latex.codecogs.com/png.latex?r_c) will drastically reduce the computational complexity of the problem as particles which lie *far* away from each other are treated as if there was no interactions between them. A particle only interacts with the particles whose distance is at most ![equation](https://latex.codecogs.com/png.latex?r_c).

### Celllist
In order to exploit the short-range nature of the particle-particle interaction, a spatial subdivision structure is introduced. This allows performing fast neighbor queries to find out which particles interact with each other.

Below is a small animation which shows the difference between the naive algorithm and the celllist based algorithm:

<p align="center">
  <img src="https://github.com/thfabian/molec/blob/master/doc/video/no-cell-gif.gif" width="350"/>
  <img src="https://github.com/thfabian/molec/blob/master/doc/video/cell-gif.gif" width="350"/>
</p>

### Molec in action
A **very** small demo of one simulation performed with *molec*. In the video only ![equation](https://latex.codecogs.com/png.latex?7%5E3%20%3D%20343) particles are shown, note that with *molec* it is possible to compute the interaction for thousands of particles.

<p align="center">
<a href="https://www.youtube.com/watch?v=RcpJUXjaxks" target="_blank"><img src="http://img.youtube.com/vi/RcpJUXjaxks/hqdefault.jpg" 
alt="Molec in action" width="400" height="320" border="10" /></a>
</p>

## Running the code
A configuration file can be passed to the executable in order to set some default simulation parameters. The configuration file has to have the following structure:
```{sh}
# time step
dt = 0.005
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

The flow of molec is shown in the following image:
<p align="center">
    <img src="https://github.com/thfabian/molec/blob/master/doc/readme/program-flow.png">
</p>

## References 
 * [[1]](http://udel.edu/~arthij/MD.pdf "Gonnet paper") Gonnet paper
 * [[2]](http://cacs.usc.edu/education/cs596/01-1LinkedListCell.pdf "Cell list implementation") Cell list implementation
 
