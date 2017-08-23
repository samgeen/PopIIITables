'''
Test the fortran code for reading the tables
Sam Geen, August 2017
'''

import matplotlib.pyplot as plt
import glob, os, sys
import numpy as np
import scipy.io

def run(var = "radius"):
    filenames = glob.glob("PopIII/track*.dat")
    filenames.sort()
    plt.clf()
    for fn in filenames:
        # Read
        print "Running for file", fn
        f = scipy.io.FortranFile(fn,"r")
        a = f.read_reals(dtype=np.float64)
        ms = f.read_reals(dtype=np.float64)
        rs = f.read_reals(dtype=np.float64)
        ts = f.read_reals(dtype=np.float64)
        f.close()
        # Plot
        if var == "radius":
            plt.plot(ms,rs,label=a)
        else:
            plt.ylim([0,5])
            plt.plot(ms,ts,label=a)
    plt.legend(frameon=False,fontsize="x-small")
    filename = "dump/testfortran_"+var+".pdf"
    plt.savefig(filename)

if __name__=="__main__":
    run("radius")
    run("temperature")
