#!/usr/bin/python3
import sys
aligned = {}
#<alignments> <reads.fasta>
for line in open (sys.argv[1]):
    id = line.split()[0]
    if id not in aligned.keys():
        aligned[id] = []
    aligned[id].append(line.strip())
for line in open (sys.argv[2]):
    if line[0] == ">" or line[0] == "@":
        id = line.strip().split()[0][1:]
        if id in aligned.keys():
            for outp in aligned[id]:
                print(outp)
