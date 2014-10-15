#! /usr/bin/python3

import sys

my_file = sys.argv[1]

f = open(my_file)
a = []
a.append(0)
b = []
c = []
i=0

for line in f:
    if "ATOM" in line[0:4]:
        a.append(line[23:26])
        b.append(line[13:16])
        c.append(line)
    else:
        pass

b.append("filler")
b.append("filler")
b.append("filler")
b.append("filler")

for i in range(1,len(a)):
    
    if int(a[i])-int(a[i-1]) > 1:
        print(a[i])

    #elif b[i] != b[i+3]:
        #print(c[i])

print(a[-1])

f.close()
