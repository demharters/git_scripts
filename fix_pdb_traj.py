#! /usr/bin/python3

# Usage: binning.py traj_file 


import sys
import os

def main(*argv):

    fn = sys.argv[1]

    fout_name = 'fixed_' + str(fn)

    fout = open(fout_name,'w')

    my_vals = ["ATOM","TER","ENDMDL"]

    with open(fn) as f:

        i = 1

        for line in f:

            if "MODEL" in line:

                print(str(line.split()[0]),str(i),file= fout)
                i += 1

            elif line.split()[0] in my_vals:

                print(line.rstrip("\n"),file = fout)

            else:
                pass

    fout.close()


if __name__ == "__main__":
    main(sys.argv)

