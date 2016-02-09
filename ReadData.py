#! /usr/bin/python

import sys
from os.path import join
import pickle

if __name__=='__main__':
    input_file = sys.argv[1]
    data = pickle.load(open(join('pickles',input_file),'rb'))
   
    #myOutFile = "../dists.dat"

    dists1 = []
    frames = []
    cter = 0
    
    for elem in data:
        for elem_2 in data[elem]:
            try:
                dists1.append(elem_2[2][1][1])
                frames.append(elem_2[0])
                cter+=1
            except IndexError:
                pass

        print cter
        break

    #with open(myOutFile,"w") as f:
    #    for i in range(0,len(dists1)):
    #        f.write("%s\t%s\n"%(frames[i],dists1[i]))
    
    #Output format of this example:
	#o_waters (model->100, atom pair-->(('A', 11, 'HO'), ('A', 11, 'O2P')), interacting water with distances -->[((81, ''), 4.023108623937465, 3.8371495931224784)])
