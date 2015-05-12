#! /usr/bin/python3

# Usage: binning.py filter_data data_to_extract_from column1 column2 cutoff 


import sys
import os

def main(*argv):

    fn = sys.argv[1]
    my_data = sys.argv[2]
    my_column1 = int(sys.argv[3])-1
    my_column2 = int(sys.argv[4])-1
    my_cutoff = float(sys.argv[5])
   
    # generate filename for output 
    fout_tmp = my_data
    if fout_tmp.endswith('.dat'):
        fout_tmp = fout_tmp[:-4]
    
    fout_higher_name = fout_tmp + '_higher2col_' + str(my_cutoff) + '.dat'
    fout_lower_name = fout_tmp + '_lower2col_' + str(my_cutoff) + '.dat'

    higher_vals = []
    lower_vals = []

    f = open(fn)

    i = 0

    for line in f:

        i += 1
       

        if (float(line.split()[my_column1]) > my_cutoff) and (float(line.split()[my_column2]) > my_cutoff):
            #print("Higher than %s: %s %s" % (my_cutoff,line.split()[my_column1],line.split()[my_column2]))
            higher_vals.append(i)

        elif (float(line.split()[my_column1]) < my_cutoff) and (float(line.split()[my_column2]) < my_cutoff):
            #print("Lower than %s: %s %s" % (my_cutoff,line.split()[my_column1],line.split()[my_column2]))
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
