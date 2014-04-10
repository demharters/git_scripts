#!/usr/bin/python

from urllib import urlopen
import brewer2mpl
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import rcParams
import matplotlib.gridspec as gridspec
from pylab import savefig
import os
import sys

my_file = sys.argv[1]

file = open(my_file)

my_data = pd.read_csv(file, sep="  ", header=None)

for i in range(0,len(my_data.ix[0,:])):

        Y = my_data.ix[:,i].dropna().values

        x1 = np.linspace(Y.min(), Y.max(), 100)

        plt.hist(Y,bins=x1)

        plt.show()

