#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
from Bio import SeqIO
from Bio.Seq import Seq

start = time.clock()

alignment = list(SeqIO.parse('nad3.pir', 'pir'))

sequence_matrix = [list(s.seq) for s in alignment]

# Transposed sequence_matrix
positions = [list(i) for i in zip(*sequence_matrix)]

"""
positions[0]:
['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'M', '-', '-', '-']

i.e. the 0th character of each sequence

"""

nr_seqs = len(positions[0])
# User-adjustable parameters
_is     = int((nr_seqs * 0.5) + 1)
_fs     = int((nr_seqs * 0.85))
_cp     = 8
_bl_1   = 15
_bl_2   = 10

# -1 nonconserved, gap
# 0 nonconserved
# 1 conserved
# 2 highly conserved

def classifier(freq, most_common_freq):  
    if(most_common_freq < _is):
        return 0
    elif((most_common_freq >= _is) and (most_common_freq < _fs)):
        return 1
    elif(most_common_freq >= _fs):
        return 2

def set_conservation(pos):
    max_residue = max(pos,key=pos.count)
    max_freq = pos.count(max_residue)
    
    if(any(v == '-' for v in pos)):
        return(-1)
    else:
        return classifier(max_residue, max_freq)
    
counter = 0
cons_level = {}
for index, p in enumerate(positions):
    c = set_conservation(p)
    if(c <= 0):
        counter += 1
    else:
        counter = 0
    if(counter < 9):
        cons_level[index] = c
    elif(counter == 9):
        for k in range(index-8, index):
            del cons_level[k]



# cons_level    - a dict with all remaining positions and their conservation level
# blocks_1      - one dict for each block (i.e. split by consecutive nonconserved regions)
# hc_positions  - all highly conserved positions in blocks
# flanks        - flanking (first and last) highly conserved position within a block
# blocks_2      - blocks anchored by flanking positions

#def set_blocks(conservation):
#blocks = [{s:conservation[s] for s in set(v)} for k, v in groupby(conservation, key=lambda i,j=count(): i-next(j))]


def group_by_consecutives(list_of_dicts):
    return [{s:list_of_dicts[s] for s in set(v)} for k, v in groupby(list_of_dicts, key=lambda i,j=count(): i-next(j))]

from itertools import groupby, count
blocks_1 = group_by_consecutives(cons_level)

# Look for HC values (=2) and keep the indices of these in a list 
hc_positions = [[(k) for k,v in (d.items()) if v == 2] for d in blocks_1]
# Extract the first and last HC in each block
flanks = [(f[0], f[-1]) for f in hc_positions if len(f) > 1]

blocks_2 = [{key: value for key, value in cons_level.items() if (key >= start) and (key <= end)} for (start, end) in flanks]

    
all_gaps = []
for block in blocks_2:
    gaps = []
    is_gap = False
    before = 0
    after = 0
    for k,v in block.items():
        if(v == -1):
            is_gap = True
            gaps.extend([g for g in range(k-before,k+1)])
        if(v == 0 and is_gap == False):
            before += 1
        elif(v > 0 and is_gap == False):
            before = 0
        if(v == 0 and is_gap == True):
            after += 1
        elif(v > 0 and is_gap == True):
            gaps.extend([g for g in range(k-1,k-1+after)])
            before = 0
            is_gap = False
    if(len(gaps) != 0):
        all_gaps.extend(gaps)




blocks_3 = {}
for b in blocks_2:
    for k,v in b.items():
        if(k not in all_gaps):
            blocks_3[k] = v
    
blocks_4 = group_by_consecutives(blocks_3)

new_positions = []

for b in blocks_4:
    for k,v in b.items():
        #print(k,v)
        new_positions.append(positions[k])

new_alignment = [list(i) for i in zip(*new_positions)]



for i,s in enumerate(alignment):
    new_seq = Seq(''.join(new_alignment[i]))
    s.seq = new_seq

SeqIO.write(alignment, "output.fna", "fasta")

    
end = time.clock()
print(end-start)
