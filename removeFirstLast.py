#! /usr/bin/python

import sys

my_file = sys.argv[1]
myIdx = my_file.find(".pdb")
fout_name = my_file[:myIdx] + "_trunc.pdb"

my_dict = {}
myMax = []
myMin = []

# Create dictionary with chains as keys and lists of residues as values
with open(my_file,"r") as f:
    for line in f:
        if "ATOM" in line:
            chain = line[21]

            if chain in my_dict:
                my_dict[chain].append(line[23:26].strip())
            else:
                my_dict[chain] = [line[23:26].strip()]

        else:
            pass

# Create lists for Min and Max residues for each chain
for i in my_dict:
    myMax.append(max(my_dict[i], key=lambda x:int(x)))
    myMin.append(min(my_dict[i], key=lambda x:int(x)))

# Copy input file without Max and Min residues
with open(my_file,"r") as f:
    with open(fout_name,"w") as fout:

        for line in f:

            if "ATOM" in line:
                k = 0
                for i in my_dict:



                    #print "if %s in %s and (%s == %s or %s == %s):"%(i,line[21],line[23:26].strip(),myMax[k],line[23:26].strip(),myMin[k])

                    if i in line[21] and (line[23:26].strip() == myMax[k] or line[23:26].strip() == myMin[k]):
                        break
                    
                    elif i in line[21] and (line[23:26].strip() != myMax[k] or line[23:26].strip() != myMin[k]):
                        fout.write(line)
                    
                    else:
                        pass
        
                    k += 1
            
            else:
                fout.write(line)
