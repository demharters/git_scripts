#! /usr/bin/python


import MDAnalysis
import sys
from pylab import *

my_traj = sys.argv[1]
myTitle = sys.argv[2]

u = MDAnalysis.Universe(my_traj,my_traj)
OH = u.selectAtoms("segid B and resid 18 and name HO")

xArr = []
yArr = []
zArr = []
data = []

for ts in u.trajectory:
    xArr.append(OH.coordinates()[0,0])
    yArr.append(OH.coordinates()[0,1])
    zArr.append(OH.coordinates()[0,2])


def normalise(myArray):

    newArray = []
    mean = sum(myArray)/float(len(myArray))
    
    for i in range(0,len(myArray)):
        newArray.append(myArray[i]-mean)

    return newArray


for myArray in (xArr,yArr,zArr):
    newArray = normalise(myArray)
    data.append(newArray)

boxplot(data)
xticks([1,2,3],["x","y","z"])
ylim(-0.4,0.4)
title("%s_HO"%myTitle)
savefig("%s_HO.png"%myTitle)

