#!/usr/bin/python3


import sys
import os
import random
from os import listdir
from os.path import isfile, join
binning = sys.argv[1]
kmer = 31
for line in open(binning):
    arr = line.split()
    hap = arr[1]
    pat = int(arr[2])
    mat = int(arr[3])
    if pat > 2*kmer and mat > 2* kmer:
        if pat < 10*mat and mat < 10 * pat:
            print (line.strip())
