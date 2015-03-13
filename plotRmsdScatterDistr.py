#!/usr/bin/python

# usage: rmsdEnergy.py rmsd_file energy_file xlim1 xlim2 ylim1 ylim2

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import savefig
from sklearn.neighbors.kde import KernelDensity
import os
import sys
#import csv
#from scipy import stats

rmsd =[]
energy = []
my_folders = []
rmsdFile = sys.argv[1]
energyFile = sys.argv[2]
my_xlim1 = float(sys.argv[3])
my_xlim2 = float(sys.argv[4])
my_ylim1 = float(sys.argv[5])
my_ylim2 = float(sys.argv[6])

my_rmsd = np.genfromtxt(rmsdFile,usecols=(0),delimiter='  ',dtype=None) 
my_energy = np.genfromtxt(energyFile,usecols=(0),delimiter='  ',dtype=None)

#### Density Overlay ########
my_rmsd_reshaped = my_rmsd.reshape(-1,1)
#
kde1 = KernelDensity(bandwidth=0.05).fit(my_rmsd_reshaped)

datamin = my_rmsd_reshaped.min()  ############# changed back to dist1
datamax = my_rmsd_reshaped.max()
datarange = datamax - datamin
#
my_x = np.linspace(datamin - (datarange*0.1), datamax + (datamax*0.01), 300).reshape(-1, 1)
#
my_density1 = np.exp(kde1.score_samples(my_x))
#
fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
ax.set_xlim(my_xlim1,my_xlim2)
ax.set_ylim(my_ylim1,my_ylim2)
ax.scatter(my_rmsd,my_energy,alpha=0.2)
#
ax2.plot(my_x, my_density1, color = "red")
ax2.set_ylim(0,5)
#plt.show()
plt.savefig("rmsdEnergy.png")

