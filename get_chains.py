#! /usr/bin/python3

# Usage: get_chains.py [file_id] [chain 1] ... [chain n] 


import sys
import os

def main(*argv):

    my_dir = os.getcwd()
    aa3 = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR"
    files_id = sys.argv[1]


    for fn in os.listdir(my_dir):

        if files_id in fn:

            tmp1 = fn[0:4]
            tmp2 = "%s_clean" % tmp1

            fout = open(tmp2,"w")
            f = open(fn)

            for line in f:
                
                if "ATOM" in line and line[17:20] in aa3 and line[21] in sys.argv[2:]:
					
                    fout.write(line)

                else:
                    pass

            fout.close()
            f.close()

        else:
            pass

if __name__ == "__main__":
	main(sys.argv)
