#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import time
from Bio import SeqIO
from Bio.Seq import Seq

start = time.clock()

alignment = list(SeqIO.parse('nad3.pir', 'pir'))

original_alignment = list(SeqIO.parse('nad3.pir', 'pir'))

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
 #   print(pos)
    if(any(v == '-' for v in pos)):
        return(-1)
    else:
        return classifier(max_residue, max_freq)
    
counter = 0
cons_level = {}
for index, p in enumerate(positions):
#    if(index > 339 and index < 351):
        c = set_conservation(p)
  #      print(index+1, c)
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
hc_positions = [[k for k,v in sorted(d.items()) if v == 2] for d in blocks_1]
# Extract the first and last HC in each block
flanks = [(f[0], f[-1]) for f in hc_positions if len(f) > 1]

blocks_2 = [{key: value for key, value in cons_level.items() if (key >= start) and (key <= end)} for (start, end) in flanks]

    
all_gaps = []
for block in blocks_2:
    gaps = []
    is_gap = False
    before = 0
    after = 0
    for k,v in sorted(block.items()):
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
            gaps.extend([g for g in range((k-after), k)])
            after = 0
            is_gap = False
    if(len(gaps) != 0):
        all_gaps.extend(gaps)

blocks_3 = {}
for b in blocks_2:
    for k,v in b.items():
        if(k not in all_gaps):
            blocks_3[k] = v

# Remove small remaining blocks <= BL2
blocks_to_remove = []
blocks_4 = group_by_consecutives(blocks_3)
for i, block in enumerate(blocks_4):
    if len(block) < _bl_2:
        blocks_to_remove.append(i)
blocks_5 = [b for i,b in enumerate(blocks_4) if i not in blocks_to_remove]

new_positions = []

for b in blocks_5:
    for k,v in sorted(b.items()):
        new_positions.append(positions[k])

new_alignment = [list(i) for i in zip(*new_positions)]

for i,s in enumerate(alignment):   
    new_seq = Seq(''.join(new_alignment[i]))
    s.seq = new_seq

SeqIO.write(alignment, "output5.fna", "fasta")


end = time.clock()
print(end-start)


######## HTML OUTPUT #########


import pandas as pd
#blosum = pd.read_csv('blosum62.csv')

df = pd.DataFrame(sequence_matrix)

from collections import defaultdict

def color_most_freq(series):
    count_dict = defaultdict(int)
    for s in series.values:
        count_dict[s] += 1
    max_val = max(count_dict, key=count_dict.get)
    
    if count_dict[max_val] < _is or max_val == '-':
        max_val = 0  
    
    # Color scheme based on http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
    
    colors = {0: 'black',
              'D': '#E60A0A', # ASP  - Bright Red
              'E': '#E60A0A', # GLU - Bright Red
              'C': '#E6E600', # CYS Yellow
              'M': '#E6E600', # MET Yellow
              'K': '#145AFF', # LYS Blue
              'R': '#145AFF', # ARG Blue
              'S': '#FA9600', # SER Orange
              'T': '#FA9600', # THR Orange
              'F': '#3232AA', # PHE Mid blue
              'Y': '#3232AA', # TYR Mid blue
              'N': '#00DCDC',
              'Q': '#00DCDC',
              'G': '#561010',#'#EBEBEB',
              'I': '#0F820F',
              'L': '#0F820F',
              'V': '#0F820F',
              'A': '#561010',#'#C8C8C8',
              'W': '#B45AB4',
              'H': '#8282D2',
              'P': '#DC9682',
              }
    
    color = colors[max_val]
    
    return ['color: %s' % color if v == max_val else '' for v in series.values]

def bold_most_freq(series):
    count_dict = defaultdict(int)
    for s in series.values:
        count_dict[s] += 1
    max_val = max(count_dict, key=count_dict.get)
    
    if count_dict[max_val] < _is or max_val == '-':
        max_val = 0  
    return ['font-weight: bold' if v == max_val else '' for v in series.values]

def hover(hover_color="#f1eeee"):
    return dict(selector="tr:hover",
                props=[("background-color", "%s" % hover_color)])
      
def highlight_cols(x):
    d = x.copy()
    hb = [item+1 for sublist in blocks_5 for item in sublist if item+1 in d.columns]
    d[hb] = 'background-color: #ded6d5'
    return d  

styles = [
    hover(),
    dict(selector="th",
                 props=[("font-size", "7pt"),
                        ("text-align", "left")]),
    dict(selector="thead",
                 props=[("font-size", "7pt"),
                        ("text-align", "center")]),
    dict(selector="tbody",
                 props=[("font-size", "9pt")])]

html = []

df.columns += 1
#df.index += 1

df.index = [alignment[i].id for i in range(len(df.index))]

nr_columns = len(df.columns)

col = [i for i in range(0,nr_columns) if (i % 40 == 0)]

for i,c in enumerate(col):
    if(c == col[-1]):
        html.append(df[df.columns[c:]])
    else:
        html.append(df[df.columns[c:col[i+1]]])




for part in range(len(html)):
    html[part] = html[part].style.set_properties(**{'text-align': 'center'})\
    .set_table_styles(styles)\
    .apply(highlight_cols, axis=None)\
    .apply(color_most_freq, axis=0)\
    .apply(bold_most_freq, axis=0)\
    .render()
    print(part)


wrapper = \
"""
<html>
    <head>
    </head>
    <body><style>table {width:640px; border-collapse: collapse; font-family: Verdana;}</style>
    <p><b>Parameters used</b><br>
    Minimum Number Of Sequences For A Conserved Position: %d<br>
    Minimum Number Of Sequences For A Flanking Position: %d<br>
    Maximum Number Of Contiguous Nonconserved Positions: %d<br>
    Minimum Length Of A Block: %d<br>
    Allowed Gap Positions: None<br>
    </p>
    <p><b>Flank positions of the %d selected block%s</b><br>
    Flanks: %s</p>
    <p><b>New number of positions in output: %d</b> (%s of the original %d positions)<br>
    <p>%s</p>

    </body>
    </html>
"""

fw = open('test.html','w')

flanks = [(list(k.keys())[0]+1, list(k.keys())[-1]+1) for k in blocks_5]

mult_flanks = ''
if len(flanks) > 1:
    mult_flanks = 's'

pos_perc = str(round((len(new_positions)/len(positions)) * 100,0)) + '%'


whole = wrapper % (_is, _fs, _cp, _bl_2, len(flanks), mult_flanks, flanks,
                   len(new_positions), pos_perc, len(positions), '<p> </p>'.join(html))

fw.write(whole)
fw.close()

end = time.clock()
print(end-start)



