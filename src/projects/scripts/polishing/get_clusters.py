#!/usr/bin/python3
import sys
aligned = {}
snps1 = set()
#<snps> <cluster_dist> <min_cluster_size>
file = sys.argv[1]
dist = int(sys.argv[2])
min_size = int(sys.argv[3])
coords = []
all_snps = {}
for line in open (sys.argv[1]):
    arr = line.split()
    if (arr[0]) != "chrX":
        continue
    coord = int(arr[2])   
    coords.append(coord)
    all_snps[coord] = line.strip()
coords.sort()
total_size = 0
cur_clust = [coords[0]]
for i in range (1, len(coords)):
    if coords[i] > coords[i -1] + dist:
        if len(cur_clust) >= min_size:
            print(f"cluster of {len(cur_clust)} errors")
            for f in cur_clust:
                print (all_snps[f])
            total_size += len(cur_clust)
        cur_clust = []
    cur_clust.append(coords[i])
if len(cur_clust) >= min_size:
    print(f"cluster of {len(cur_clust)} errors")
    for f in cur_clust:
        print (all_snps[f])
    total_size += len(cur_clust)
print (f'Total in clusters {total_size} snps of {len(coords)}')
