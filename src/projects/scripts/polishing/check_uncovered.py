#!/usr/bin/python3
import sys
aligned = {}
snps = {}
starts = {}
finishes = {}
#<coord.filtered> <snps>
#50711672 181660018 | 130948308 10 | 130948347 130948299 | 100.0 | chr5 ptg000001l |
for line in open (sys.argv[1]):
    arr = line.split()
    chr = arr[11]
    if chr not in starts.keys():
        starts[chr] = []
        finishes[chr] = []
    starts[chr].append(int(arr[0]))
    finishes[chr].append(int(arr[1]))
for chr in starts.keys():
    starts[chr].sort()
    finishes[chr].sort()
#print (starts)
#print (finishes)
uncovered = 0
total = 0
for line in open (sys.argv[2]):
    arr = line.split()
    coord = int(arr[2])    
    chr = arr[0]
    ind = 0
    covered = False
    while ind < len (starts[chr]):
        if finishes[chr][ind] < coord:
            ind+=1
            continue
        elif  starts[chr][ind] <= coord:
            covered = True
            break
        else:
            break
    if not covered:
        uncovered +=1
    total +=1
print (f'Total {total} snps, uncovered: {uncovered}')
