#!/usr/bin/python3
import sys
aligned = {}
snps1 = set()
#<snps1> <snps2> <diameter>
diameter = int(sys.argv[3])
same = 0
notsame = 0
for line in open (sys.argv[1]):
    arr = line.split()
    coord = int(arr[2])   
    for i in range(coord - diameter, coord + diameter + 1):
        snps1.add(i) 
for line in open (sys.argv[2]):
    arr = line.split()
    coord = int(arr[2])   
    suff = ""
    if coord in snps1:
        same +=1
        suff = " shared"
    else:
        notsame +=1
    print (line.strip() + suff)
print (f'Total {same + notsame} snps, shared(with eps {diameter}): {same}')
