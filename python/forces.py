#!usr/bin/env python3
#                   _
#   _ __ ___   ___ | | ___  ___
#  | '_ ` _ \ / _ \| |/ _ \/ __|
#  | | | | | | (_) | |  __/ (__
#  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
#
#  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
#                      Michel Breyer  (mbreyer@student.ethz.ch)
#                      Florian Frei   (flofrei@student.ethz.ch)
#                      Fabian Thuring (thfabian@student.ethz.ch)
#
#  This file is distributed under the MIT Open Source License.
#  See LICENSE.txt for details.

from pymolec import *

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os.path

# seaborn formatting
sns.set_context("notebook", font_scale=1.1)
sns.set_style("darkgrid")
sns.set_palette('deep')
deep = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]

def measure_performance():

#    forces = ['cell_ref', 'q', 'q_g']
#    N = np.logspace(3, 7.5, 15, base=10).astype(np.int32)
#    steps = np.array([1000, 800, 800, 600, 600, 600, 600, 400, 400, 400, 300, 200, 100, 50, 30]);

    forces = ['cell_ref', 'q', 'q_g']
    N = np.logspace(3, 7.5, 15, base=10).astype(np.int32)
    steps = np.array([1000, 800, 800, 600, 600, 600, 600, 400, 400, 400, 300, 200, 100, 50, 30]);


    rho = 1.25
    rc  = 2.5

    flops =  N * rc**3 * rho * (18 * np.pi + 283.5)
    
    if os.path.isfile("performances-forces.npy"):
        print("Loading data from <performances-forces.npy>")
        performances = np.load("performances-forces.npy")
        return performances, N, forces
    else:
        performances = np.zeros((len(forces), len(N)))
        
        for force_idx,force in enumerate(forces):
            p = pymolec(N=N, force=force, steps=steps, rho=rho)
            times = p.run()
            
            # store the performance in the array
            perf = flops / times[0,:]
            performances[force_idx, :] = perf
            
        print("Saving performance data to <performances-forces.npy>")    
        np.save("performances-forces", performances)
        
        return performances, N ,forces

def plot_performance(performances, N, forces):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1);

    #iterate over the forces
    for force_idx, force in enumerate(forces):
        perf = performances[force_idx, :]
        ax.semilogx(N, perf, 'o-')

    ax.set_xlim([np.min(N)*0.9, np.max(N)*1.1])
    ax.set_ylim([0, 2.5])

    ax.set_xlabel('Number of particles')
    ax.set_ylabel('Performance [Flops/Cycle]',
                  rotation=0,
                  horizontalalignment = 'left')
    ax.yaxis.set_label_coords(-0.055, 1.05)

    plt.legend(forces)

    filename = 'forces.pdf'
    print("saving '%s'" % filename )
    plt.savefig(filename)


if __name__ == '__main__':
    perf, N, forces = measure_performance()
    plot_performance(perf, N, forces)
