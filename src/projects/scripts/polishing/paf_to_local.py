#!/usr/bin/python3
import sys
#SRR10963010.15  7018    856     1640    +       contig_4631     309059  31294   31929   76      784     0       tp:A:P  cm:i:4  s1:i:47 s2:i:47 dv:f:0.0211     rl:i:0
#read_id contig_id alignment_start alignment_end
for line in open (sys.argv[1], 'r'):
    arr = line.split()
    if len(arr) < 9:
        break
    read_id = arr[0]
    rc = False
    contig_id = arr[5]
    if arr[4] == '-':
        rc = True
    aln_len = int(arr[6])
    aln_start = int(arr[7])
    aln_end = int(arr[8])
    if rc:
        tmp = aln_start
        aln_start = aln_len - aln_end
        aln_end = aln_len - tmp
    minus = ""
    if rc:
        minus = "-"        
    print (f'{read_id} {minus}{contig_id} {aln_start} {aln_end}')
