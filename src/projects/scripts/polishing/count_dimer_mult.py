#!/usr/bin/python3
import sys
#<contigs> <cutoff>
counts = {}
min_count= int(sys.argv[2])
for line in open (sys.argv[1]):
    lenl = len(line)
    if lenl < 1000:
        continue
    start = 0
    while (start < lenl):
        multiplicity = 1
        while ((start + 2 * multiplicity + 1 < lenl) and line[start] == line[start + 2 * multiplicity ]  and line [start + 1] == line[start + 2 * multiplicity + 1]):
            multiplicity+= 1
        if multiplicity >= 3:
            start += 1
            if not multiplicity in counts:
                counts[multiplicity] = 0
            counts[multiplicity] += 1
        start = start + 2 * multiplicity - 1;
res = 0
for key, value in sorted(counts.items()):
    print("{} : {}".format(key, value))
    if key >= min_count:
        res += value
print (f'Total {res} dinucleotide repeats longer than {min_count}')
