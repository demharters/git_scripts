#! /usr/bin/python

import sys
import os

def main(*argv):

    my_dir = os.getcwd() 
    files_id = sys.argv[1]
    chains = sys.argv[2:]


    for fn in os.listdir(my_dir):

        if files_id in fn:

            for j in chains:

                f = open(fn)
                i = 0
                prev_res = ""

                for line in f:

                    if "ATOM" in line and line[21] in j and line[23:26] != prev_res:
		       
                        if i == 0:
                            first_res = int(line[23:26])

                        i += 1
                        prev_res = line[23:26]

                    else:
                        pass

                f.close()

                last_res = first_res + i-1

                print "%s chain %s first_res = %d last_res = %d length = %d" % (fn,j,first_res,last_res,i)	

if __name__ == "__main__":
	main(sys.argv)
