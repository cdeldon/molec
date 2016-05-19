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
import time, sys, os, subprocess

class pymolec:

    def __init__(self, N=np.array([1000]), rho=1.25, steps=np.array([100]),
                 force="cell_ref", integrator="lf", periodic="ref"):

        self.N = N
        self.rho = rho

        
        if hasattr(steps, "__len__"):
            if len(N) != len(steps):
                self.steps = np.full(len(N), steps[0], dtype=np.int)
        else:
            self.steps = np.full(len(N), steps, dtype=np.int)
        

        self.force = force
        self.integrator = integrator
        self.periodic = periodic

    def run(self, path = None):

        # Use default path
        if not path:
            script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)))
            if os.name == 'nt':
                path = os.path.join(script_path, '..', 'build', 'molec.exe')
            else:
                path = os.path.join(script_path, '..', 'build', 'molec')

        # Check if molec exists
        if not os.path.exists(path):
            raise IOError("no such file or directory: %s" % path)

        times = np.zeros((4, len(self.N)))

        print ("Running molec: %s" % path)
        print ("rho = {0}, force = {1}, integrator = {2}, periodic = {3}".format(
            self.rho, self.force, self.integrator, self.periodic))

        for i in range(len(self.N)):
            cmd = [path]
            cmd += ["--N=" + str(self.N[i])]
            cmd += ["--rho=" + str(self.rho)]
            cmd += ["--step=" + str(self.steps[i])]
            cmd += ["--force=" + self.force]
            cmd += ["--integrator=" + self.integrator]
            cmd += ["--periodic=" + self.periodic]
            cmd += ["--verbose=0"]

            # Print status
            start = time.time()
            print(" - N = %9i ..." % self.N[i], end='')
            sys.stdout.flush()

            try:
                out = subprocess.check_output(cmd).decode(encoding='utf-8').split('\t')

                print(" %20f s" % (time.time() - start))

                times[0,i] = int(out[2]) # force
                times[1,i] = int(out[4]) # integrator
                times[2,i] = int(out[6]) # periodic
                times[3,i] = int(out[8]) # simulation

            except subprocess.CalledProcessError as e:
                print(e.output)

        return times

def main():
    p = pymolec()
    print(p.run())

if __name__ == '__main__':
    main()
