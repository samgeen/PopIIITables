'''
Read Population III star tables
Sam Geen, June 2017
'''

import numpy as np
import scipy.interpolate
import scipy.io

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from collections import OrderedDict

interporder = "linear"

Msun = 1.98855e33
year = 31556926.0
Rsun = 6.95700e10

class PopIIIReader(object):

    def __init__(self):
        self._header = None
        self._tables = None
        self._filestubs = ["1em5","3em5","1em4","3em4","1em3"]
        self._accs = [float(x.replace("m","-")) for x in self._filestubs]

    def Header(self):
        '''
        Header for each file, list of variables
        '''
        if self._header is None:
            raw = open("PopIII/header.txt").read()
            self._header = raw.replace("\n","").split(" ")
        return self._header

    def Tables(self):
        if self._tables is None:
            filenames = ["PopIII/tab"+str(x)+".dat" \
                for x in self._filestubs]
            self._tables = OrderedDict()
            for acc, filename in zip(self._accs,filenames):
                self._tables[acc] = np.loadtxt(filename)
        return self._tables

    def Make2DGrid(self,varname):
        '''
        Make a 2D grid from the input tables
        '''
        naccs = 50
        nmasses = 100
        accs = []
        masses = []
        values = []
        ivar = self.Header().index(varname)
        imass = self.Header().index("mass")
        for acc, table in self.Tables().iteritems():
            accs = np.concatenate((accs,table[:,imass]*0.0+acc))
            masses = np.concatenate((masses,table[:,imass]))
            values = np.concatenate((values,table[:,ivar]))
        # Make grid
        points = np.zeros((len(values),2))
        points[:,0] = np.log10(accs)
        points[:,1] = masses
        #points = np.array([np.log10(accs),masses])
        amin = np.min(self._accs)
        amax = np.max(self._accs)
        mmin = masses.min()
        mmax = masses.max()
        newaccs = np.linspace(np.log10(amin),np.log10(amax),naccs)
        newmasses = np.linspace(mmin,mmax,nmasses)
        grida, gridm = np.meshgrid(newaccs,newmasses)
        newgrid = scipy.interpolate.griddata(points,values,(grida,gridm),
                                             method=interporder)
        return 10.0**newaccs, newmasses, newgrid

    def VarFromTables(self,varname):
        '''
        Return 2D table of a single variable
        '''
        ivar = self.Header().index(varname)
        imass = self.Header().index("mass")
        nt = len(self.Tables())
        # Prepare to interpolate
        nv = 1000
        masses = np.array([])
        for ml, table in self.Tables().iteritems():
            masses = np.concatenate((masses,table[:,imass]))
        mlow = np.min(masses)
        mhigh = np.max(masses)
        mgrid = np.arange(mlow,mhigh*1.00001,(mhigh-mlow)/float(nv-1))
        grid = np.zeros((len(self._accs),nv))
        for acc, table in self.Tables().iteritems():
            iacc = self._accs.index(acc)
            vals = table[:,ivar]
            f = scipy.interpolate.interp1d(table[:,imass], table[:,ivar],
                                           bounds_error=False,fill_value=0.0)
            grid[iacc,:] = f(mgrid)
        return self._accs, mgrid, grid
            
    def PlotVals(self,varname):
        accs,ms,v = self.VarFromTables(varname)
        plt.clf()
        for i in range(0,len(accs)):
            plt.plot(ms,v[i,:],label=str(accs[i])+" Msun/yr")
        plt.legend()
        plt.xlabel("Mass / Msun")
        plt.ylabel(varname)
        #plt.yscale("log")
        plt.savefig("dump/popiiitest_"+varname+".pdf")

    def Make2DTables(self,varname):
        acc,mass,vals = self.Make2DGrid(varname)
        # Convert units to linear cgs
        acc *= Msun/year
        mass *= Msun
        if varname == "radius":
            vals = 10.0**vals * Rsun # log(R/Rsun) -> cm
        if varname == "T":
            vals = 10.0**vals # K
        # Write Fortran file
        # Write as 4D table used by lookup_table_module.f90
        filename = "PopIII/popiii_grid_"+varname+".dat"
        ndim = np.array([2])
        sizes = np.array([len(acc),len(mass),1,1,])
        dum = np.array([1.0])
        f = scipy.io.FortranFile(filename,"w")
        # Dimensions
        f.write_record(ndim.astype(np.int32))
        # Sizes of arrays
        f.write_record(sizes.astype(np.int32))
        # Array axes
        f.write_record(acc.astype(np.float64))
        f.write_record(mass.astype(np.float64))
        f.write_record(dum.astype(np.float64))
        f.write_record(dum.astype(np.float64))
        # Values in ND table
        f.write_record(vals.astype(np.float64))
        f.close()
        
    def Plot2DGrid(self,varname):
        acc,mass,vals = self.Make2DGrid(varname)
        vals[np.isnan(vals)] = vals[np.isreal(vals)].min()
        #print np.min(acc), np.min(mass), np.min(vals)
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        nacc = len(acc)
        for iacc in range(0,nacc):
            label = None
            if (iacc % (nacc//5)) == 0 or iacc == nacc-1:
                label = r"$log(\dot{M})=$"+str(np.log10(acc[iacc]))
            r = float(iacc)/float(nacc)
            col = (r,0.0,0.0)
            ax.plot(mass,vals[:,iacc],color=col,label=label)
        ax.set_xlabel('Mass / Msun')
        ax.set_ylabel(varname)
        #ax.set_xscale("log")
        ax.set_yscale("log")
        ax.legend(frameon=False,fontsize="x-small")
        #plt.show()
        fig.savefig("dump/sideview_1_"+interporder+"_"+varname+".pdf")

if __name__=="__main__":
    reader = PopIIIReader()
    for v in ["radius","T"]:
        #reader.PlotVals(v)
        #reader.Plot2DGrid(v)
        reader.Make2DTables(v)
