#!/usr/bin/python3
import sys
aligned = {}
same = 0
notsame = 0
chr = 0
ctg = 0
amounts = {}
for line in open (sys.argv[1]):
    if not int(line) in amounts:
        amounts[int(line)] = 0
    amounts[int(line)] +=1
for i in sorted(amounts.keys()):
    print (f'{i} : {amounts[i]}')
