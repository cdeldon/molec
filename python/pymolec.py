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
import os, subprocess

class pymolec:

    def __init__(self, N=np.array([1000]), steps=100, force="N2", integrator="lf",
             periodic="ref"):

        self.N = N
        self.steps = steps

        self.force = force
        self.integrator = integrator
        self.periodic = periodic
        
        print("DONE")

    def run(self, path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'build', 'molec')):

        times = np.zeros((4, len(self.N)))

        for i in range(len(self.N)):
            cmd = [path]
            cmd += ["--N=" + str(self.N[i])]
            cmd += ["--force=" + self.force]
            cmd += ["--integrator=" + self.integrator]
            cmd += ["--periodic=" + self.periodic]
            cmd += ["--step=" + str(self.steps)]
            cmd += ["--verbose=0"]

            print(cmd[1:])

            out = subprocess.check_output(cmd).decode(encoding='utf-8').split('\t')

            times[0,i] = int(out[2])
            times[1,i] = int(out[4])
            times[2,i] = int(out[6])
            times[3,i] = int(out[8])

        return times
    
