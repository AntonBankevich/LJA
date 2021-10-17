#!/usr/bin/python3
import sys
if len(sys.argv) != 3:
    print ('Usage: {sys.argv[0]} <alignments> <contigs>')
cont_len = {}
cur_len = 0
cur_id = ""
for line in open(sys.argv[2]):
    if line[0] == ">":
        if cur_id != "":
            cont_len[cur_id] = cur_len
        cur_id = line.strip()[1:]
        cur_len = 0
    else:
        cur_len += len(line.strip())
cont_len[cur_id] = cur_len

starts = {}
finishs = {}
for line in open(sys.argv[1]):
    arr = line.split()
    id = arr[-3]
    start = int(arr[-2])
    finish = int(arr[-1])
    if id[0] == '-':
        id = id[1:]
        contig_len = cont_len[id];
        tmp = contig_len - finish;
        finish = contig_len - start;
        start = tmp;
    if not id in starts.keys():
        starts[id] = []
        finishs[id] = []
   
    starts[id].append([start, 0])
    starts[id].append([finish, 1])
#    print (f'{id} {start} {finish}')
#    print (line)
empty = {}
for id in starts.keys():
    starts[id].sort()
    current = 0
    cur_start = 0
    empty[id] = []
    for pair in starts[id]:
#        print (f'{id}  {pair[0]} {pair[1]}')              

        if pair[1] == 0:
            if current == 0 and pair[0] != cur_start:
                empty[id].append([cur_start, pair[0]])
            current +=1
            
        else:
            current -=1
            if current == 0:
                cur_start = pair[0]
    for pair in empty[id]:
        if pair[1] - pair[0] > 1:
            middle = " ends "
            if pair[0] !=0 and pair[1]!= cont_len[id]:
                middle = " MIDDLE "
            print (f'{id}  {pair[0]} {pair[1]}  {pair[1] - pair[0]}   {middle} {cont_len[id]}')              
