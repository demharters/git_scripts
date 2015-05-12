import numpy as np
from MDAnalysis import *
import sys, os
import csv
from collections import defaultdict

# my_path = os.getcwd() + "/"

mySetups = ["init_nochm_full","init_chm_full","init_chm_fullfixedCHM","init_chm_fullfixedGC"]
my_path = "/home/scratch/mosaics/arm/dna/schofield/fixed_tors/"
os.chdir(my_path)
myBigFile = sys.argv[1] + "_cat.dat"
mySmallFile= sys.argv[1] + ".dat"
myColumn = int(sys.argv[2])

def extract_pdbs(myFrames):
    u = Universe("sampled.pos.pdb","sampled.pos.pdb")
    protein = u.selectAtoms("all")
    with Writer("sampled_extracted.pos.pdb", protein.numberOfAtoms()) as W:
        for ts in u.trajectory:
            if ts.frame in myFrames:
                W.write(protein)

for setup in mySetups:

    myFile = my_path + setup + "_cat/" + myBigFile

    myVals = []

    with open(myFile) as f:
        reader = csv.reader(f,delimiter="\t")
        myList = list(reader)

    for i in range(0,len(myList)):
        myVals.append(myList[i][myColumn])

    myMedian = np.median(map(float,myVals))

    for dir,subdirs,files in os.walk("."):
        for myDir in subdirs:
            if setup in myDir and "rep" in myDir:

                os.chdir(myDir)
                myIndexes = []
                myVals2 = []

                with open(mySmallFile) as f2:
                    reader = csv.reader(f2,delimiter="\t")
                    myList = list(reader)

                for i in range(0,len(myList)):
                    myVals2.append(myList[i][myColumn])

                for i, j in enumerate(map(float,myVals2)):
                    if j == myMedian:
                        myIndexes.append(i+1)

                extract_pdbs(myIndexes)
                os.chdir(my_path)

            else:
                pass





    # myFrames = defaultdict(list)
    #
    # for i in myIndexes:
    #     if i <= 1000:
    #         myFrames["1"] = i
    #     if i > 1000 and len(str(i)) == 4:
    #         myRep = str(i)[0]
    #         myFrames[myRep] = int(str(i)[1:])
    #     if i > 1000 and len(str(i)) == 5:
    #         myRep = str(i)[0:2]
    #         myFrames[myRep] = int(str(i)[2:])
    #
    # print myFrames
    # print len(myFrames)
    # print len(myIndexes)