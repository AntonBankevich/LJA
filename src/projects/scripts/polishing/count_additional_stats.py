#!/usr/bin/python3
import sys
file = open(sys.argv[1])
low_cover_high = 0
mults = {}
line = file.readline()
while True:
    line = file.readline()

    if line == '':
        break
    if line.find("low_cov") > 0:
        line = file.readline()
        mult = int(line.split()[3])
        if mult >= 10:
            low_cover_high += 1
    elif line.find("high_")> 0:
        line = file.readline()
        mult = int(line.split()[3])
        if not mult in mults.keys():
            mults[mult] = 0
        mults[mult] +=1
    else:
        line = file.readline()
print (low_cover_high)
for i in range (2, 16):
    sum = 0
    for j in range(0, 5):
        if 5*i + j in mults.keys():
            sum += mults[5*i+j]
    print (f'range {5*i} - {5*i+ 5} : {sum}')
