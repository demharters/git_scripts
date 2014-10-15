#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pylab import savefig
from sklearn.neighbors.kde import KernelDensity
import os
import sys
import csv
from scipy import stats

my_folders = []

def import_text(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator,
                           skipinitialspace=True):
        if line:
            yield line


my_if = sys.argv[1]
my_path = os.getcwd()
my_file = my_path + "/" + my_if
my_base = my_path.split('/')[-1]

my_data = np.loadtxt(my_if)
my_rmsd = my_data[0]
my_energy = my_data[1]

slope, intercept, r_value, p_value, std_err = stats.linregress(my_rmsd,my_energy)

print("slope = %s"%slope)
