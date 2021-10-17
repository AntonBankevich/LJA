#!/usr/bin/python3
import sys
aligned = {}
snps = {}
if len(sys.argv) != 3:
    print (f'Usage: {sys.argv[0]} <snps_file> <debug files prefix>')
    exit()
def get_chr_pos(line):
    arr = line.split()
    if len(arr) < 6:
        return -1
    else:
        return int(arr[2])

def get_contig_pos(line):
    arr = line.split()
    if len(arr) < 6:
        return -1
    else:
        return int(arr[5])

for line in open (sys.argv[1]):
    arr = line.split()
    contig_id = arr[1]
    pos = int(arr[5])
    if not contig_id in snps.keys():
       snps[contig_id] = {}
    snps[contig_id][pos] = line.strip()
print ("ref\tcontig\tpos_ref\tvalue_ref\tvalue_cont\tpos_cont\ttotal_mult\tcov\tres_mult\tres\tcheck_base")
ordered_output = {}
for id in snps.keys():
    f = open (sys.argv[2] +id)
    for line in f:
        if len(line.split()) <3:
            continue

        pos = int(line.split()[0])
        if (pos) in snps[id].keys():
            snps_line = snps[id][pos]
            arr = snps_line.split()
            seq = arr[3]
            res = "unknown"
            if seq == ".":
                seq = arr[4]
            count = 0
            nucls = set()
            for i in range (0, len(seq)):
                nucls.add(seq[i])
#            print (f'{seq} {count}')
            debug_arr = line.split()
#            if arr[3] != "." and arr[4] != ".":
#                res = "mismatch"
#            el
            if len(nucls) >= 3:
                res = "complex_indel"
            else:
                total_mult = int(debug_arr[1])
                cov = int(debug_arr[2])
                res_mult = int(debug_arr[3])
                updated_res = ""
                if (total_mult == 0):
                    if (cov < 5):
                        updated_res="low_covered MSA"
                    else:
                        updated_res="normal MSA"
                elif (cov < 5):
                    updated_res= "low_covered regular"
                elif (res_mult >= 10):
                    updated_res = " high_mult"                
                if updated_res != "" and res == "unknown":
                    res = updated_res
                if (cov > 0):
                    remainder = total_mult/cov - total_mult //cov
                    if remainder > 0.45 and remainder < 0.55:
                        res += " close_call"
            if get_chr_pos(snps_line) not in ordered_output.keys():
                ordered_output[get_chr_pos(snps_line)] = []
            ordered_output[get_chr_pos(snps_line)].append(f'{snps[id][pos]}\t{total_mult}\t{cov}\t{res_mult}\t{res}\t{debug_arr[4]}')
            ordered_output[get_chr_pos(snps_line)].append(line.strip())
for id in sorted(ordered_output.keys()):
    for l in ordered_output[id]:
        print (l)
                

'''
'''
