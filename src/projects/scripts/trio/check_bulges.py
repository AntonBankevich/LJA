#!/usr/bin/python3

#add some  contamination for meta-genomic simulation

import sys
import os
import random
from os import listdir
from os.path import isfile, join
binning = sys.argv[1]
bulges = sys.argv[2]
hapl = {}
defined = set({'m', 'p'})
for line in open(binning, 'r'):
    arr = line.split()
    hapl[arr[0]] = arr[1]
good = 0
fixable = 0
unfixable = 0
contradictive = 0
total = 0
for line in open (bulges, 'r'):
    arr= line.split()
    if len(arr) <=2:
        continue
    if (hapl[arr[0]] == 'm' and hapl[arr[1]] == 'p') or (hapl[arr[0]] == 'p' and hapl[arr[1]] == 'm'):
        good +=1
    elif (hapl[arr[0]] in defined or hapl[arr[1]] in defined):
        if hapl[arr[0]]!= hapl[arr[1]]:
            fixable +=1
        else:
            contradictive += 1
    else:
        unfixable +=1
        print (line.strip())

    total +=1
print (f'of total {total} bulges: ')
print (f'good: {good} fixable: {fixable} unfixable: {unfixable} contradictive: {contradictive}')
