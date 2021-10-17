#!/usr/bin/python3
import sys
s=""
for line in open(sys.argv[1], 'r'):
    line = line.strip()
    if line[0] == ">" or line[len(line) -1] == "s":
        line = "\n" + line + "\n"
    print (line, end='')
