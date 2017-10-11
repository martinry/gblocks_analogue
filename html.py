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
            #print(range(k-before,k+1))
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
    #new_seq = Seq(''.join([r if (ind == 0 or ind % 10 != 0) else (' ' + r) for ind,r in enumerate(new_alignment[i])]))
    new_seq = Seq(''.join(new_alignment[i]))
    s.seq = new_seq

SeqIO.write(alignment, "output5.fna", "fasta")


end = time.clock()
print(end-start)
    

import pandas as pd

#blosum = pd.read_csv('blosum62.csv')

df = pd.DataFrame(sequence_matrix)


def hover(hover_color="#ffff99"):
    return dict(selector="tr:hover",
                props=[("background-color", "%s" % hover_color)])


#def blocks(col='#E5D283'):
    
           

def highlight_cols(x):
    #copy df to new - original data are not changed
    d = x.copy()
    #select all values to default value - red color
    #print(d.columns)
    #d[d.columns] = 'background-color: red'
    #d.columns = 'background-color: red'
    #overwrite values grey color
    hb = [item for sublist in blocks_5 for item in sublist if item in d.columns]
    #print(hb)
    #test_cols = [1,2,3,4,5,6,7,22]
    #d[test_cols] = 'background-color: grey'
    d[hb] = 'background-color: #8ee5ee'
  #  d[[11,22]] = 'border-color: none'
    #return color df
    return d  

styles = [
    hover(),
    dict(selector="th",
                 props=[("font-size", "0pt")]),
    dict(selector="tbody",
                 props=[("font-size", "8pt")])]

#    dict(selector="th", props=[("font-size", "100%"),
#                               ("text-align", "center")])
#    dict(selector="caption", props=[("caption-side", "bottom")])


html = []


df.columns += 1
df.index += 1

df.index = [alignment[i].id for i in range(len(df.index))]


#app = df.style.set_table_styles(styles)




nr_columns = len(df.columns)

#df.columns = [i if (i%10==0) else '' for i in df.columns]

col = [i for i in range(0,nr_columns) if (i % 40 == 0)]

#html += df.style.apply(to_color, axis=0).render()

for i,c in enumerate(col):
    if(c == col[-1]):
        html.append(df[df.columns[c:]])
    else:
       # print(i)
        html.append(df[df.columns[c:col[i+1]]])


#html[0].style.bar(subset=[18], color='lightblue')
for part in range(len(html)):
    html[part] = html[part].style.set_properties(**{'width': '50px'}).set_table_styles(styles).apply(highlight_cols, axis=None).render()
    
#html[1] = html[1].style.apply(highlight_cols, axis=None).set_table_styles(styles).render()




#html = html.to_html()


#pd.set_option('display.width', '800')


#.highlight_null().render().split('\n')[:10]


"""
html = (
    df.style
    .set_properties(**{'font-size': '9pt', 'font-family': 'Calibri'})
    #.bar(subset=[18], color='lightblue')
    .render()
)
"""
#<style>table {width:840px;}</style>

wrapper = \
"""
<html>
    <head>
    </head>
    <body><style>table {max-width:840px;}</style>
    %s%s%s%s
    </body>
    </html>
"""

fw = open('test.html','w')

whole = wrapper % (html[0], html[1], html[2], html[3])

fw.write(whole)
fw.close()







