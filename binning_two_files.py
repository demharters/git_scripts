#! /usr/local/bin/python3

# Usage: binning.py filter_data trajectory column cutoff 


import sys
import os

def main(*argv):

    fn = sys.argv[1]
    my_data = sys.argv[2]
    my_column = int(sys.argv[3])-1
    my_cutoff = float(sys.argv[4])
   
    # generate filename for output 
    fout_tmp = my_data
    if fout_tmp.endswith('.dat'):
        fout_tmp = fout_tmp[:-4]
    
    fout_higher_name = fout_tmp + '_higher_' + str(my_cutoff) + "_col" + str(my_column+1) + '.dat'
    fout_lower_name = fout_tmp + '_lower_' + str(my_cutoff) + "_col" + str(my_column+1) + '.dat'

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

    fout_higher = open(fout_higher_name,'w')
    fout_lower = open(fout_lower_name,'w')

    with open(my_data) as f:

        parsing = False
        i = 1
        j = 1
        lineCter = 0

        for line in f:

            lineCter += 1

            if lineCter in higher_vals:
                print(line.strip("\n"),file= fout_higher)

                i += 1

            elif lineCter in lower_vals:

                print(line.strip("\n"),file= fout_lower)
            
                j += 1
            
            else:
                pass


    fout_higher.close()
    fout_lower.close()



if __name__ == "__main__":
	main(sys.argv)
