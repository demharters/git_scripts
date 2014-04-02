#! /usr/bin/python3

# Usage: binning.py filter_data cutoff1 cutoff2 


import sys
import os

def main(*argv):

    fn = sys.argv[1]
    glob_cutoff = float(sys.argv[2])
    bg_cutoff = float(sys.argv[3])
    
    fullbind_vals = []
    topbind_vals = []
    lowbind_vals = []
    nobind_vals = []

    f = open(fn)

    i = 0

    for line in f:

        i += 1
        
        if float(line.split()[1]) > bg_cutoff and float(line.split()[0]) > glob_cutoff:
            
            nobind_vals.append(i)

        elif float(line.split()[1]) < bg_cutoff and float(line.split()[0]) > glob_cutoff:
            
            topbind_vals.append(i)

        elif float(line.split()[1]) > bg_cutoff and float(line.split()[0]) < glob_cutoff:
            
            lowbind_vals.append(i)

        else:
        
            fullbind_vals.append(i)

    f.close()


    fout_fullbind = open("fullbind.pdb",'w')
    fout_topbind = open("topbind.pdb",'w')
    fout_lowbind = open("lowbind.pdb","w")
    fout_nobind = open("nobind.pdb","w")

    with open("sampled.pos.pdb") as f:

        parsing = False
        i = 1
        j = 1
        k = 1
        l = 1

        for line in f:

            if "MODEL" in line and int(line.split()[1]) in fullbind_vals:

                parsing = "fullbind"

                print(str(line.split()[0]),str(i),file = fout_fullbind)

                i += 1

            elif "MODEL" in line and int(line.split()[1]) in topbind_vals:

                parsing = "topbind"

                print(str(line.split()[0]),str(j),file = fout_topbind)

                j += 1

            elif "MODEL" in line and int(line.split()[1]) in lowbind_vals:

                parsing = "lowbind" 

                print(str(line.split()[0]),str(k),file= fout_lowbind)

                k += 1

            elif "MODEL" in line and int(line.split()[1]) in nobind_vals:

                parsing = "nobind" 

                print(str(line.split()[0]),str(l),file= fout_nobind)
                
                l += 1

            else:
                pass


            if "MODEL" not in line and parsing == "fullbind":

                print(line.rstrip("\n"),file = fout_fullbind)

            elif "MODEL" not in line and parsing == "topbind":

                print(line.rstrip("\n"),file = fout_topbind)

            elif "MODEL" not in line and parsing == "lowbind":

                print(line.rstrip("\n"),file = fout_lowbind)

            elif "MODEL" not in line and parsing == "nobind":

                print(line.rstrip("\n"),file = fout_nobind)
            
            else:
                pass

    fout_fullbind.close()
    fout_topbind.close()
    fout_lowbind.close()
    fout_nobind.close()



if __name__ == "__main__":
	main(sys.argv)
