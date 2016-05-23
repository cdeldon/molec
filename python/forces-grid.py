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

    forces = ['q'];
    
    N     = np.logspace(4,7,8).astype(np.int32)
    steps = np.array([100, 100, 90, 80, 65, 50, 35, 20])
    rhos  = np.array([0.5, 1., 2., 4., 6.,8.,10.])


    rc  = 2.5

    if os.path.isfile("performances-grid-forces-density.npy"):
        print("Loading data from <performances-grid-forces-density.npy")
        performances = np.load("performances-grid-forces-density.npy")
        return performances, N, rhos
    else:

        performances = np.zeros((len(rhos), len(N)))

        for rho_idx, rho in enumerate(rhos):
            flops =  N * rc**3 * rho * (18 * np.pi + 283.5)

            p = pymolec(N=N, rho=rho, force=forces, steps=steps, integrator='lf8', periodic='c4')
            output = p.run()

            perf = flops / output['force']
            performances[len(rhos)-1-rho_idx, :] = perf

        print("Saving performance data to <performances-grid-forces-density.npy>")
        np.save("performances-grid-forces-density", performances)

    return performances, N, rhos

def plot_performance(performances, N, rhos):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1);

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(10, 133, n = 256, as_cmap=True)

    ax = sns.heatmap(performances, linewidths=1,
                      yticklabels=rhos[::-1], xticklabels=N,
                      vmax=0.2*np.round(np.max(np.max(performances))*5),
                      vmin=0.2*np.round(np.min(np.min(performances))*5),
                      cmap=cmap, annot=False
                      )


    cax = plt.gcf().axes[-1]
    pos_old = cax.get_position()
    pos_new = [pos_old.x0 - 0.01, pos_old.y0 + 0,  pos_old.width, pos_old.height*((len(rhos)-1)*1./len(rhos))]
    cax.set_position(pos_new)
    cax.tick_params(labelleft=False, labelright=True)

    ax.text(len(N)+0.35, len(rhos), 'Performance\n[flops/cycle]', ha='left', va='top')


    rho_labels_short = ['%.2f' % a for a in rhos]
    ax.set_yticklabels(rho_labels_short)
    
    N_labels_short = ['10$^{%1.2f}$' % a for a in np.array(np.log10(N))]
    ax.set_xticklabels(N_labels_short)

    ax.set_xlabel('Number of particles $N$')
    ax.set_ylabel('Particle density',
                    rotation=0, horizontalalignment = 'left')
    ax.yaxis.set_label_coords(0., 1.01)
    plt.yticks(rotation=0)

    filename = 'forces-grid.pdf'
    print("saving '%s'" % filename )
    plt.savefig(filename)


if __name__ == '__main__':
    perf, N, rhos = measure_performance()
    plot_performance(perf, N, rhos)
