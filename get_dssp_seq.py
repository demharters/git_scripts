#! /usr/bin/python3

import DSSPData as dd
import sys

dd_ob = dd.DSSPData()

dd_ob.parseDSSP(sys.argv[1])

my_ss = dd_ob.getSecStruc()
my_res = dd_ob.getResnums()

seg_counter = 0

segment = []
segment_idx = []

start_idx = 0
end_idx = 0

prev_index = 0

for i in range(1,len(my_ss)):

    if my_ss[i] == my_ss[prev_index]:
        segment.append(my_ss[i])
        segment_idx.append(my_res[i])

    else:
        segment_no = "segment_%s" % seg_counter
        seg_counter += 1
        
        if len(segment_idx) >= 1:

            start_idx = segment_idx[0]
            end_idx = segment_idx[-1]

            print(segment_no,"\n",start_idx,"-",end_idx,"\n",segment)

        else:
            pass
        
        
        # reinitialise 
        segment = []
        segment_idx = []
   
    prev_index += 1

