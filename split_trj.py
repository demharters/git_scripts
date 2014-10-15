#!/usr/bin/python

import sys
import os

my_file = sys.argv[1]
#my_file = "sampled.pos.pdb"

my_frames = []
trj_1 = []
trj_2 = []

with open(my_file) as f:

    for line in f:
        if "MODEL" in line:
            my_frames.append(int(line[6:9].strip()))

        else:
            pass

trj_1_status = True
trj_1.append(my_frames[0])

for frame in range(1,len(my_frames)):

    if abs(my_frames[frame] - my_frames[frame-1]) < 10:
        mswitch = 0

    elif abs(my_frames[frame] - my_frames[frame-1]) > 10:
        mswitch = 1

    else:
        pass


    if mswitch == True and trj_1_status == True:
        trj_2.append(my_frames[frame])
        trj_1_status = False

    elif mswitch == True and trj_1_status == False:
        trj_1.append(my_frames[frame])
        trj_1_status = True

    elif mswitch == False and trj_1_status == True:
        trj_1.append(my_frames[frame])
        trj_1_status = True

    elif mswitch == False and trj_1_status == False:
        trj_2.append(my_frames[frame])
        trj_1_status = False

    else:
        pass


print my_frames
print (trj_1,trj_2)
print (my_frames[-1],my_frames[-2])

my_trj_file = "trj_frames.ndx"

with open(my_trj_file,'a') as f:

    f.write("[ Trajectory 1 ]\n")
    for iter in trj_1:
        f.write("%s " % iter)

    f.write("\n[ Trajectory 2 ]\n")
    for iter in trj_2:
        f.write("%s " % iter)
