#!/usr/bin/python

import sys
import os
import subprocess
import re
import time

def replaceEntry(fileIn, fileOut, oldWord, newWord):
    with open(fileIn) as oldFile, open(fileOut, 'w') as newFile:
        for line in oldFile:
            if oldWord not in line:
                newFile.write(line)
            else:
                newLine = line.replace(oldWord, newWord)
                newFile.write(newLine)

print os.getcwd()
myFile = sys.argv[1]
idx = myFile.find(".pdb")
outName = myFile[0:idx] + "_edited.pdb"

remove_H = "remove_H.pl %s" % myFile
process1 = subprocess.Popen(remove_H.split(), stdout=subprocess.PIPE)
output1 = process1.communicate()[0]
with open("1.pdb","w") as f1:
    f1.write(output1)

h_add = "pymol_hadd.py 1.pdb"
process2 = subprocess.Popen(h_add.split(), stdout=subprocess.PIPE)
# pymol takes a moment to open up. Wait 5s.
time.sleep(5)

badWords = [ "CONECT", "TER" ]
with open('1_hadd.pdb') as oldFile, open('2.pdb', 'w') as newFile:
    for line in oldFile:
        if not any(badWord in line for badWord in badWords) and "HETATM" not in line:
            newFile.write(line)
        if "HETATM" in line:
            newLine = line.replace("HETATM", "ATOM  ")
            newFile.write(newLine)

replaceEntry("2.pdb","3.pdb","DG","G ")
replaceEntry("3.pdb","2.pdb","DT","T ")
replaceEntry("2.pdb","3.pdb","DC","C ")
replaceEntry("3.pdb","4.pdb","DU","U ")
replaceEntry("4.pdb","2.pdb","DA","A ")

replaceEntry("2.pdb","3.pdb","RSQ","CFO")
replaceEntry("3.pdb","2.pdb","HOH","TIP")
replaceEntry("2.pdb","3.pdb","OP1","O1P")
replaceEntry("3.pdb","2.pdb","OP2","O2P")
replaceEntry("2.pdb","3.pdb","C10","C  ")
replaceEntry("3.pdb","7.pdb","O30","O  ")

rename_H = "rename_RNA_H.pl 7.pdb"
process8 = subprocess.Popen(rename_H.split(), stdout=subprocess.PIPE)
output8 = process8.communicate()[0]
with open("8.pdb","w") as f8:
    f8.write(output8)

#regex = re.compile(r"H0?", re.IGNORECASE)
with open("8.pdb") as oldFile, open(outName, "w") as newFile:
    for line in oldFile:
        if "TIP" in line:
            #newLine = regex.sub(r"H", line)
            newLine = line.replace("H01","H1 ")
            newLine = newLine.replace("H02","H2 ")
            newFile.write(newLine)
        else:
            newFile.write(line)

delete_tmpFiles = "rm 1.pdb 2.pdb 3.pdb 4.pdb 5.pdb 6.pdb 7.pdb 8.pdb 1_hadd.pdb"
processEnd = subprocess.Popen(delete_tmpFiles.split(), stdout=subprocess.PIPE)
