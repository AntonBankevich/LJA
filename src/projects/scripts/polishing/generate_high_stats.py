#!/usr/bin/python3
import sys
import os
aligned = {}
same = 0
notsame = 0
chr = 0
ctg = 0
snps = {}
snps_file = sys.argv[1]
debug_prefix = sys.argv[2]
length_cutoff = int(sys.argv[3])

for line in open (snps_file):
    arr = line.split()
    if len (arr) < 10:
        continue
    if arr[9] != "high_mult":
        continue
    contig_id = arr[1]
    
    pos = int(arr[5])
    addition = 0
    
    if not contig_id in snps.keys():
       snps[contig_id] = {}
    if arr[3] == ".":
        addition = -len (arr[4])
    elif arr[4] == ".":
        addition = len (arr[3])
    snps[contig_id][pos] = addition
#    if contig_id == "50":
#        print (f'{pos} {addition}')
realmults = {}
#27 6 3 1 C    1 1 4

shift = 0
for id in snps.keys():
    
    if not (os.path.exists(debug_prefix + id)):
        continue
    print (f'processing contig {id}')
    f = open (debug_prefix +id)
    skip = 0
    for line in f:
        skip -=1
        if (skip > 0):
            continue
        if len(line.split()) <3:
            continue
        arr = line.split()
        pos = int(arr[0])
        mult = int(arr[3])
        addition = 0
        
        for pp in range(pos, pos +mult):
            if pp in snps[id].keys():
                addition =  snps[id][pp]
                shift += 1
                break
        real_mult = mult + addition
        skip = mult
        if real_mult >= length_cutoff:
            if addition != 0:
                print (f'processing pos {pos} with add {addition}')
            
            if real_mult not in realmults.keys():
                realmults[real_mult] = [0,0,0]
            eq = 0
            ls = 0
            mr =0
            for  j in range (5, len(arr)):
                if int(arr[j]) > real_mult:
                    mr+= 1
                elif int(arr[j]) == real_mult:
                    eq +=1
                else:
                    ls +=1
            realmults[real_mult][0] += ls
            realmults[real_mult][1] += eq
            realmults[real_mult][2] += mr
    print (f'shifted {shift}')
    print ("read_fraction:\tless\tequal\tmore\ttotal_reads")
    for real_mult in range(length_cutoff, 70):
        if real_mult in realmults.keys():
            total_s = realmults[real_mult][0] + realmults[real_mult][1] + realmults[real_mult][2]
            print (f'{real_mult}\t{realmults[real_mult][0]/total_s:.5f}\t{realmults[real_mult][1]/total_s:.5f}\t{realmults[real_mult][2]/total_s:.5f}\t{total_s}')
