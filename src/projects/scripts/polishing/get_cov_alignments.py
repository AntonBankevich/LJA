#!/usr/bin/python3
import sys
if len(sys.argv) != 3:
    print ('Usage: {sys.argv[0]} <alignments> <position>')
used_names = set()
#Get all alignments record (RC ignored) covering exact position
pos = int(sys.argv[2])
for line in open (sys.argv[1]):
    arr = line.split()
    id = arr[0]
    name = arr[1]
    if name[0] == '-':
        continue
    start = int(arr[2])
    end = int(arr[3])
    if start < pos and end > pos:
        print (line.strip())

'''
cur_pos = 0
to_use = False
id = ""
#grep -P 'ccs 21|ccs -21'
for line in open (sys.argv[2]):
    if line[0] == ">":
        id = line.strip()[1:]
        if id in used_names:
            print (">" + id)        
    else:        
        if id in used_names:
            print (line.strip())
'''
