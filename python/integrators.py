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
import json
import sys

#------------------------------------------------------------------------------

integrators = ['lf', 'lf2', 'lf4', 'lf8', 'lf_avx']

N = np.logspace(2, 5, 12, base=10).astype(np.int32)
steps = np.array([25])

rho = 1.0
rc  = 2.5

#------------------------------------------------------------------------------

filename = sys.argv[1]

results = {}

for integrator in integrators:
    p = pymolec(N=N, rho=rho, steps=steps, force='q_g_avx', integrator=integrator)
    output = p.run()

    results['N'] = output['N'].tolist()
    results['rho'] = output['rho'].tolist()
    results[integrator] = output['integrator'].tolist()

print('Saving performance data to ' + filename)

with open(filename, 'w') as outfile:
    json.dump(results, outfile, indent=4)
