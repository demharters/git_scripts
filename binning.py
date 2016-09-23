#! /usr/local/bin/python3

# Usage: binning.py filter_data trajectory column cutoff 


import sys
import os

def main(*argv):

    fn = sys.argv[1]
    my_traj = sys.argv[2]
   
    # generate filename for output 
    fout_tmp = my_traj
    if fout_tmp.endswith('.pdb'):
        fout_tmp = fout_tmp[:-4]

    my_column = int(sys.argv[3])-1
    my_cutoff = float(sys.argv[4])

    higher_vals = []
    lower_vals = []

    f = open(fn)

    i = 0

    for line in f:

        i += 1
        
        if float(line.split()[my_column]) > my_cutoff:
			
            higher_vals.append(i)

        elif float(line.split()[my_column]) < my_cutoff:

            lower_vals.append(i)

        else:
            pass

    f.close()

    fout_higher_name = fout_tmp + '_higher_' + str(my_cutoff) + '.pdb'
    fout_lower_name = fout_tmp + '_lower_' + str(my_cutoff) + '.pdb'


    fout_higher = open(fout_higher_name,'w')
    fout_lower = open(fout_lower_name,'w')

    with open(my_traj) as f:

        parsing = False
        i = 1
        j = 1

        for line in f:

            if "MODEL" in line and int(line.split()[1]) in higher_vals:

                parsing = True

                print(str(line.split()[0]),str(i),file= fout_higher)

                i += 1

            elif "MODEL" in line and int(line.split()[1]) not in higher_vals:

                parsing = False

                print(str(line.split()[0]),str(j),file= fout_lower)
            
                j += 1
            
            else:
                pass


            if "MODEL" not in line and parsing:

                print(line.rstrip("\n"),file = fout_higher)

            elif "MODEL" not in line and not parsing:

                print(line.rstrip("\n"),file = fout_lower)

            else:
                pass

    fout_higher.close()
    fout_lower.close()



if __name__ == "__main__":
	main(sys.argv)
