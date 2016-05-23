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

for force in results:
    ax.semilogx(N, results[force], label=force)

ax.set_xlim([np.min(N)*0.9, np.max(N)*1.1])
ax.set_xlabel('Number of particles $N$')
ax.set_ylabel('Runtime [Cycles]',
              rotation=0,
              horizontalalignment = 'left')
ax.yaxis.set_label_coords(-0.055, 1.05)
ax.legend()
plt.savefig(filename[:filename.rfind('.')]+'-runtime.pdf')

#----- plot performance -----

flops = dict()
flops['cell_ref'] = (lambda n, r : 301 * n * r * 2.5**3)
flops['q'] = lambda n, r : 301 * n * r * 2.5**3
flops['q_g'] = lambda n, r : 180 * n * r * 2.5**3
flops['q_g_avx'] = lambda n, r :  n * (205 * r * 2.5**3 + 24)

fig = plt.figure()
ax = fig.add_subplot(1,1,1);

for force in results:
    ax.semilogx(N, flops[force](N,rho) / np.array(results[force]), label=force)

ax.set_xlim([np.min(N)*0.9, np.max(N)*1.1])
ax.set_xlabel('Number of particles $N$')
ax.set_ylabel('Performance [Flops/Cycles]',
              rotation=0,
              horizontalalignment = 'left')
ax.yaxis.set_label_coords(-0.055, 1.05)
ax.legend()
plt.savefig(filename[:filename.rfind('.')]+'-performance.pdf')
