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

# seaborn formatting
sns.set_context("notebook", font_scale=1.1)
sns.set_style("darkgrid")
sns.set_palette('deep')
deep = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]

def main():

    integrators = ['lf', 'lf2', 'lf4', 'lf8']
    N = np.array([100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 10000])

    flops = 9 * N

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1);

    for integrator in integrators:
        p = pymolec(N=N, force='cell_ref', integrator=integrator)
        times = p.run()

        perf = flops / times[1,:]
        ax.plot(N, perf, 'o-')


    ax.set_ylim([0, 2])

    ax.set_xlabel('Number of particles')
    ax.set_ylabel('Performance [Flops/Cycle]',
                  rotation=0,
                  horizontalalignment = 'left')
    ax.yaxis.set_label_coords(-0.055, 1.05)

    plt.legend(integrators)

    filename = 'integrators.pdf'
    print("saving '%s'" % filename )
    plt.savefig('integrators.pdf')


if __name__ == '__main__':
    main()
