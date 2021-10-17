#!/usr/bin/python3
import sys
aligned = {}
same = 0
notsame = 0
chr = 0
ctg = 0
for line in open (sys.argv[1]):
    arr = line.split()
    if len(arr) < 5:
        continue
    if arr[3] == ".":
        chr += 1
    elif arr[4] == ".":
        ctg += 1
    else:
        continue
print (f'Total {ctg + chr} snps, deleted in chr- {chr}, in contig - {ctg}')
