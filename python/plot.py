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

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import json

# seaborn formatting
sns.set_context("notebook", font_scale=1.1)
sns.set_style("darkgrid")
sns.set_palette('deep')
deep = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]

try:
    filename = sys.argv[1]
except IndexError as ie:
    print('usage: plot results.txt')
    sys.exit(1)

# load results from json object
with open(filename, 'r') as infile:
    results = json.load(infile)

N   = np.array(results['N'])
rho = np.array(results['rho'])

del results['N']
del results['rho']

#----- plot runtime ------

fig = plt.figure()
ax = fig.add_subplot(1,1,1);

for k in sorted(results):
    if 'cell_ref' in results:
        ax.semilogx(N, np.array(results['cell_ref']) / np.array(results[k]), 'o-', label=k)
    elif 'lf' in results:
        ax.semilogx(N, np.array(results['lf']) / np.array(results[k]), 'o-', label=k)


ax.set_xlabel('Number of particles $N$')
ax.set_ylabel('Runtime Speedup',
              rotation=0,
              horizontalalignment = 'left')
ax.yaxis.set_label_coords(-0.055, 1.05)

ax.set_xlim([np.min(N)*0.9, np.max(N)*1.1])
ax.set_ylim([0.0, 1.2 * ax.get_ylim()[1]])

ax.legend(loc='upper right')

plt.savefig(filename[:filename.rfind('.')]+'-runtime.pdf')

#----- plot performance -----

flops = dict()
flops['cell_ref'] = lambda N, rho : 301 * N * rho * 2.5**3
flops['q']        = lambda N, rho : 301 * N * rho * 2.5**3
flops['q_g']      = lambda N, rho : 180 * N * rho * 2.5**3
flops['q_g_avx']  = lambda N, rho : N * (205 * rho * 2.5**3 + 24)
flops['lf']     = lambda N, rho : 9 * N
flops['lf2']    = lambda N, rho : 9 * N
flops['lf4']    = lambda N, rho : 9 * N
flops['lf8']    = lambda N, rho : 9 * N
flops['lf_avx'] = lambda N, rho : 9 * N

fig = plt.figure()
ax = fig.add_subplot(1,1,1);

for k in sorted(results):
    ax.semilogx(N, flops[k](N,rho) / np.array(results[k]), 'o-', label=k)

ax.set_xlabel('Number of particles $N$')
ax.set_ylabel('Performance [Flops/Cycles]',
              rotation=0,
              horizontalalignment = 'left')
ax.yaxis.set_label_coords(-0.055, 1.05)

ax.set_xlim([np.min(N)*0.9, np.max(N)*1.1])
ax.set_ylim([-0.1, 1.4 * ax.get_ylim()[1]])

ax.legend(loc='upper right')

plt.savefig(filename[:filename.rfind('.')]+'-performance.pdf')
