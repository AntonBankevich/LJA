#!/usr/bin/python3
import sys
def get_compressed_snps(snp_file, contig_file):
    comp_coord = 0
    prev = "Z"
    compressed = []
    for line in open (contig_file):
        if line[0] != ">":
            for c in line.strip():
                if c != prev:
                    comp_coord+= 1
                prev = c           
                compressed.append(comp_coord)
    res ={}
    for line in open (snp_file):
        arr = line.split()
        coord = int(arr[2])
        res[coord] = compressed[coord]
        arr[2] = str(compressed[coord])
        print ("\t".join(arr))
    return res

#<snps1> <uncompressed_contigs>
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f'Usage {sys.argv[0]} <snps1> <uncompressed_ref>')
        exit()
    get_compressed_snps(sys.argv[1], sys.argv[2])
    
