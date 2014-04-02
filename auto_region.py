#! /usr/bin/python3

import math
import os
from sys import argv
import params

my_dir = os.getcwd()
my_file = "%s/%s" % (my_dir,params.my_pdb)

my_chains = params.chains
flex_res = params.flex

#Get number of chains
no_chains = len(my_chains)

#Instantiate lists
my_inputs = []
in_first_res = []
in_last_res = []
segments = []

#Instantiate lists for first_res and last_res of all segments
all_first_res = []
all_last_res = []
all_centre = []
stride_list =[]

header = "CBLC ~"
for i in params.chains:
    header += i
header += "\n"

continue_var = 'y'
cter = 0


for j in my_chains:

    f = open(my_file)
    
    i = 0
    prev_res = ""

    stride = "STRIDE >"

    for line in f:
        if "ATOM" in line and line[21] in j and line[23:26] != prev_res:
            
            current_res = line[23:26].strip()

            if i == 0:
                in_first_res.append(int(line[23:26]))

            if current_res in flex_res[cter]:
                stride += "C"

            else:
                stride += "R"

            i += 1
            prev_res = line[23:26]

        else:
            pass
    
    in_last_res.append(int(in_first_res[cter]) + i-1)
    
    stride += "\n"
    stride_list.append(stride)
    f.close()
    cter += 1

my_tmp = "temp.pdb"
f = open(my_file)
f_tmp = open(my_tmp,"w")

for line in f:
    if "STRIDE" in line or "CBLC" in line or line in ['\n', '\r\n',"END"]:
        pass
    else:
        f_tmp.write(line)

f.close()
f_tmp.close()
os.rename(my_tmp,my_file)

for i in reversed(range(len(my_chains))):
    
    with open(my_file, 'r') as original: data = original.read()
    with open(my_file, 'w') as modified: modified.write(stride_list[i] + data + "\n")

with open(my_file, 'r') as original: data = original.read()
with open(my_file, 'w') as modified: modified.write(header + data)

# A function to calculate centres between flexible residues i.e. segments
def calc_centres(start,end,flex):

    start = int(start)-1
    end = int(end)+1
    
    centre = []
    seg_start = []
    seg_end= []
    
    len_seg = int()
    temp_centre = int()

    # Add start and end residues to flexible residues list for calculation of segment lengths
    print(flex)
    flex.insert(0,start)
    flex.append(end)

    # loop through flexible residues list to calculate segment lengths and centres
    for i in range (0,len(flex)-1):

        # make lists of last and first res
        seg_start.append(int(flex[i]) + 1)
        seg_end.append(int(flex[i + 1]) - 1)

        # Calculate length of segments
        len_seg = (seg_end[i] - seg_start[i]) + 1       

        # If length is even then divide by 2 and add to start residue
        
        if len_seg % 2 == 0:

            temp_centre = int((len_seg/2))
            centre.append(temp_centre + int(seg_start[i])-1)
        
        # If length is odd then divide by 2, round up (ceiling) and add to start residue
        else:

            temp_centre = int(math.ceil(len_seg/2))
            centre.append(temp_centre + int(seg_start[i])-1)

    return(seg_start,seg_end,centre)

# Store ouput from calc_centres in separate lists
for i in range(int(no_chains)):

    segments.append(calc_centres(in_first_res[i],in_last_res[i],flex_res[i][:]))

    all_first_res.append(segments[i][0])
    all_last_res.append(segments[i][1])
    all_centre.append(segments[i][2])
    

#Print region template
def print_region(start,end,centre,no_segments):

    in_nseg = no_segments
    in_ncenter = no_segments
    in_first = start
    in_last = end
    in_centre = centre
    rot_setting = params.rot
    trans_setting = params.trans
    rot_free_setting = '0'
    prop_free_setting = '0'

    element = "{segment}"
    dependency = "{independent}"
    nseg = "{%s}" % in_nseg
    ncenter = "{%s}" % in_ncenter
    firstres = "{%s}" % in_first
    lastres = "{%s}" % in_last
    baseres = "{%s}" % centre
    centers = "{%s}" % centre
    prop_trans_sig = "{%s}" % trans_setting
    prop_rot_sig = "{%s}" % rot_setting
    prop_trans_sig_freeres = "{%s}" % rot_free_setting
    prop_rot_sig_freeres = "{%s}" % prop_free_setting

    template = """~region[\element_top_type{0}
    
    \dependency_type{1}
    
    \\nseg{2}
    \\ncenter{3}
    \segments_firstres{4}
    \segments_lastres{5}
    
    \segments_baseres{6}
    
    \centers{7}
    
    \prop_trans_sig{8}
    \prop_rot_sig{9}
    \prop_trans_sig_freeres{10}
    \prop_rot_sig_freeres{11}
]
    """

    print(template.format(element, dependency, nseg, ncenter, firstres, lastres, baseres, centers, prop_trans_sig,
    prop_rot_sig, prop_trans_sig_freeres, prop_rot_sig_freeres),file=f)

def print_options(no_chains,all_first_res,all_last_res,all_centre):
    cter = int()
    for i in range(int(no_chains)):
        print("\nchain %d" % i)
        cter = 0
        for j in all_centre[i]:
            print("segment %d" % cter, all_first_res[i][cter],all_last_res[i][cter],all_centre[i][cter])
            cter +=1
            
    return


def make_regions(all_first_res,all_last_res,all_centre,seg_comb):

    seg_comb_list = seg_comb.split()
    seg_com_list_split = [i.split(':') for i in seg_comb_list]

    custom_start = []
    custom_end = []
    custom_centre = []
    
    no_segments = int()

    for i in range(len(seg_com_list_split)):

        chain = int(seg_com_list_split[i][0])
        segment = int(seg_com_list_split[i][1])

        custom_start.append(my_chains[chain]+":"+str(all_first_res[chain][segment]))
        custom_end.append(my_chains[chain]+":"+str(all_last_res[chain][segment]))
        custom_centre.append(my_chains[chain]+":"+str(all_centre[chain][segment]))
        
    no_segments = len(custom_centre)

    return(custom_start,custom_end,custom_centre,no_segments)


seg_comb = params.segment_combination 


f = open('region.data', 'w')

for i in range(len(seg_comb)):

    #print_options(no_chains,all_first_res,all_last_res,all_centre)
    first,last,centre,no_segments = make_regions(all_first_res,all_last_res,all_centre,seg_comb[i])
    first = ",".join(map(str, first))
    last = ",".join(map(str, last))
    centre = ",".join(map(str, centre))
    print_region(first,last,centre,no_segments)

f.close()
