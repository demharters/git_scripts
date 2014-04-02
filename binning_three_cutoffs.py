#! /usr/bin/python3

# Usage: binning.py filter_data cutoff1 cutoff2 cutoff3 


import sys
import os

def main(*argv):

    fn = sys.argv[1]
    glob_cutoff = float(sys.argv[2])
    bg_cutoff = float(sys.argv[3])
    glob_cutoff2 = float(sys.argv[4])
    
    fullbind_vals = []
    halfbind_vals =[]

    f = open(fn)

    i = 0

    for line in f:

        i += 1
        
        if float(line.split()[1]) < bg_cutoff and float(line.split()[0]) < glob_cutoff and float(line.split()[2]) < glob_cutoff2:
            
            fullbind_vals.append(i)

        else:
        
            halfbind_vals.append(i)

    f.close()


    fout_fullbind = open("fullfullbind.pdb",'w')
    fout_halfbind = open("halfbind.pdb","w")

    with open("sampled.pos.pdb") as f:

        parsing = False
        i = 1
        j = 1

        for line in f:

            if "MODEL" in line and int(line.split()[1]) in fullbind_vals:

                parsing = True 

                print(str(line.split()[0]),str(i),file = fout_fullbind)

                i += 1

            elif "MODEL" in line and int(line.split()[1]) not in fullbind_vals:
             
                parsing = False

                print(str(line.split()[0]),str(j),file = fout_halfbind)
            
                j += 1

            else:
                pass

            if "MODEL" not in line and parsing:

                print(line.rstrip("\n"),file = fout_fullbind)

            elif "MODEL" not in line and not parsing:

                print(line.rstrip("\n"),file = fout_halfbind)

            else:
                pass

    fout_fullbind.close()
    fout_halfbind.close()



if __name__ == "__main__":
	main(sys.argv)
