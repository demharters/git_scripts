#!/usr/local/bin/python3.5

import sys
import os

def main(*argv):

    fn = sys.argv[1]
    f = open(fn)

    my_residue_no = []
    my_chains = []
    my_lines = []
    i = 0

    parsing = False
    for line in f:

        if "ATOM" == line[0:4] and line[21:22] not in my_chains:
            parsing = True

        elif parsing and "TER" not in line and "ENDMDL" not in line:

            # my_residue_no.append(int(line[23:26]))
            my_chains.append(line[21:22])

        if "ENDMDL" in line:
            parsing = False

        # if not parsing
        else:
            pass

    for chain in my_chains:

        for oldLine in f:

            i += 1

            if line[21:22] == chain & i % 3 == False:

                my_lines.append(oldLine)


    # add final entry so parser does not crash when looking for i+1 in last iteration
    my_chains.append("END")
    
    f.close()


    for i in range(len(my_atom_name)):

        # reset temporary list
        my_tmp = []
        my_tmp2 = []

        if my_atom_name[i].strip() == "CA":
            
            my_tmp.append(my_atom_no[i])
            my_tmp.append(my_atom_no[i+1])
            
            my_tmp2.append(my_atom_no[i])
            my_tmp2.append(my_atom_no[i+2])
            
            my_conect.append(my_tmp)
            my_conect.append(my_tmp2)

        # connect oxygens with CAs of consecutive residues. Stop at last residue.
        elif my_atom_name[i].strip() == "O" and my_residue_no[i] != my_residue_no[-1] and my_chains[i] == my_chains[i+2]:
            
            my_tmp.append(my_atom_no[i])
            my_tmp.append(my_atom_no[i+2])
            
            my_conect.append(my_tmp)
            

        else:
            pass


    with open(fn, "a") as f:

        for conect_pair in my_conect:
       
            new_line = "CONECT" + str(conect_pair[0]).rjust(5) + str(conect_pair[1]).rjust(5) + "\n"
            f.write(new_line)


if __name__ == "__main__":
        main(sys.argv)
