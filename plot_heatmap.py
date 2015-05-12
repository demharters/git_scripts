#! /usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import operator
from collections import Counter

myPath = os.path.abspath(".") + "/"
myFiles = os.listdir(myPath + "/chm_notors_cat/")
myFolders = ["nochm_notors_cat/","chm_notors_cat/"]
#my_params=('shift',"slide","tilt","roll","rise","twist","propeller","buckle","opening","shear","stagger","stretch",'major_gw_refined','major_gw_pp','minor_gw_refined','minor_gw_pp')
my_params=('shift',"slide","tilt","roll","rise","twist","propeller","buckle","opening","shear","stagger","stretch")
myDict = {}
myTotalDict = {}


for folder in myFolders:
    myTotalDict[folder] = {}
    for myFile in myFiles:
        
        print myFile

        myData = []

        myFilePath = myPath + folder + myFile
        myData = np.genfromtxt(myFilePath,delimiter="\t",dtype=None)
        
        
        # initialise dictionary with 0s
        for i in range(1,len(myData[0])):
                myDict[i]=0

        print len(myData[0])

        # fill dictionary with sums of columns
        for line in myData:
            for j in range(1,len(myData[0])):
                try:    
                    myDict[j]+=line[j]
                except:
                    pass
            break

        print len(myData[0])
        print len(myDict.values())

        # calculate average for each value
        for key,value in zip(myDict,myDict.values()):
            if value != 0:
                myDict[key]=value/len(myData)
            else:
                pass

        myTotalDict[folder][myFile] = myDict.values()
        #print len(myDict.values())


nochmDict = myTotalDict[myTotalDict.keys()[0]]
chmDict = myTotalDict[myTotalDict.keys()[1]]




#diffDict = {key: map(operator.sub,nochmDict[key],chmDict.get(key,0)) for key in nochmDict.keys()}
#diffDict = {key: Counter(nochmDict[key])-Counter(chmDict[key]) for key in nochmDict.keys()}
#    pass
#for key in nochmDict.keys():
#    print key
#    #print chmDict(key,0) 
#    print chmDict[key] 
#
diffDict = {}
for key in nochmDict:
    try:
        diffDict[key] = map(operator.sub,nochmDict[key],chmDict[key])
    except:
        diffDict[key] = map(operator.sub,nochmDict[key][:-1],chmDict[key])


print diffDict
