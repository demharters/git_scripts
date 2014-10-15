#! /usr/bin/python3



import sys
import os

def main(*argv):

    my_file = sys.argv[1]
    my_chothia = sys.argv[2]
    end = my_file.find(".pdb")
    fout_name = my_file[0:end] + '_renum.pdb'
    fout = open(fout_name,'w')
    my_vals = ["TER","ENDMDL"]
    chothia_h = []
    chothia_l = []
    resiter_h = 0
    resiter_l = 0

    my_open_chothia = open(my_chothia)
    
    for line in my_open_chothia:

        if line.split()[0] == "ATOM" and line.split()[4] == "H":
        
            chothia_h.append(line[23:27])

        elif line.split()[0] == "ATOM" and line.split()[4] == "L":

            chothia_l.append(line[23:27])

        else:
            pass

    my_open_chothia.close()


    my_open_file = open(my_file)

    for line in my_open_file:

        i = 1

        if line.split()[0] == "MODEL":
        
            print(str(line.split()[0]),str(i),file= fout)
            i += 1
            resiter_h = 0
            resiter_l = 0

        elif line.split()[0] == "ATOM" and line.split()[4] == "H" and resiter_h < len(chothia_h):
        
            print("%s%s%s" % (line[0:23],chothia_h[resiter_h],line[27:].rstrip("\n")), file = fout)
            resiter_h += 1
        
        elif line.split()[0] == "ATOM" and line.split()[4] == "L" and resiter_l < len(chothia_l):
        
            print("%s%s%s" % (line[0:23],chothia_l[resiter_l],line[27:].rstrip("\n")), file = fout)
            resiter_l += 1
        
        else:
           
           print(line.rstrip("\n"),file = fout)


    fout.close()



if __name__ == "__main__":
	main(sys.argv)
