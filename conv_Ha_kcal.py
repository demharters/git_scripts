#! /usr/bin/python

import sys

my_if = sys.argv[1]
my_open_if = open(my_if)

my_of = my_if + "_conv"
my_open_of = open(my_of,"w")

conv_factor = 627.509607998

for line in my_open_if:

    my_energy = float(line.split()[0])
    my_new_energy = my_energy * conv_factor
    my_open_of.write(str(my_new_energy) + '\n')


