#!/usr/bin/python3
import sys
file = sys.argv[1]
num_lines = sum(1 for line in open(file))
total = (num_lines -1 )/2
inf = open(file, "r")
data = inf.read()
unknown = data.count("unknown")
unknown_close = data.count("unknown close")
all_close = data.count("close")
low_MSA = data.count("low_covered MSA")
low_reg = data.count("low_covered regular")
high_mult = data.count("high_mult")
reg_MSA = data.count("normal MSA")
complex = data.count("complex_indel")
print (f'unknown: {unknown* 100/total:.2f} ')
print (f'   of them: close call: {unknown_close* 100/total:.2f} ')
print (f'low_MSA: {low_MSA* 100/total:.2f} ')
print (f'low_reg: {low_reg* 100/total:.2f} ')
print (f'high_mult: {high_mult* 100/total:.2f} ')
print (f'reg_MSA: {reg_MSA* 100/total:.2f} ')
print (f'complex: {complex* 100/total:.2f} ')
print (f'  all close call: {all_close* 100/total:.2f} ')


print (f'total: {total } ')
print (f'unknown: {unknown } ')
print (f'   of them: close call: {unknown_close} ')
print (f'low_MSA: {low_MSA} ')
print (f'low_reg: {low_reg} ')
print (f'high_mult: {high_mult} ')
print (f'reg_MSA: {reg_MSA} ')
print (f'complex: {complex} ')
print (f'  all close call: {all_close} ')
