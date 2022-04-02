#!/usr/bin/python3

#add some  contamination for meta-genomic simulation

import sys
import os
import random
from os import listdir
from os.path import isfile, join
binning = sys.argv[1]
color = {}
color['0'] = "grey"
color['a'] = "violet"
color['m'] = "red"
color['p'] = "blue"
print ("Name,Colour")
for line in open(binning):
    arr = line.split()
    print (f'{arr[0]},{color[arr[1]]}')
