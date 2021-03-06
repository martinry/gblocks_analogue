#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

### Argparse ###
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.8'
)

parser.add_argument(
    '-i',
    dest='infile',
    type=argparse.FileType('r'),
    required=True,
    help='select an input file in any of the following formats: FASTA (.fasta, .fna, .fna), PIR (.pir)'
)


args = parser.parse_args()

######## PARSING #########

alignment = list(SeqIO.parse(args.infile, 'pir'))

######## FORMAT DATA STRUCTURE #########

sequence_matrix = [list(s.seq) for s in alignment]

# Transposed sequence_matrix
positions = [list(i) for i in zip(*sequence_matrix)]

"""
$ positions[0]
['-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', 'M', '-', '-', '-']

i.e. the 0th character of each sequence

"""

######## SIMILARITY MATRIX #########

# Exhaustive Matching of the Entire Protein Sequence Database
# http://www.jstor.org/stable/pdf/2876949.pdf?refreqid=excelsior%3Aa5ef7e218e7c3a6a572bc8e399158883
sim = pd.read_csv('gonnet120.csv', index_col='aa')
#sim = pd.read_csv('blosum62.csv', index_col='aa')

sim = sim.to_dict()

######## SET CONSERVATION #########

nr_seqs = len(positions[0])
# User-adjustable parameters

# Minimum Number Of Sequences For A Conserved Position
# (50% of the number of sequences + 1)
# Any integer bigger than half the number of sequences and smaller or equal than the total number of sequences
_b1     = int((nr_seqs * 0.5) + 1)

# Minimum Number Of Sequences For A Flank Position
# (85% of the number of sequences)
# Any integer equal or bigger than Minimum Number Of Sequences For A Conserved Position
_b2     = int((nr_seqs * 0.85))

# Maximum Number Of Contiguous Nonconserved Positions
# Any integer
_b3     = 8

# Minimum Length Of A Block
# Any integer equal or bigger than 2
_b4     = 10

# Allowed Gap Positions
# (None, With Half, All)
_b5     = None

# Use Similarity Matrices
# (Yes (1), No (0))
_b6     = 1


def classifier(most_common_freq, sim_score):  
    if most_common_freq < _b1:
        return 0
    elif most_common_freq >= _b1 and \
         most_common_freq < _b2:
             ### I haven't found a specific threshold for the similarity score, 
             ### 8 is what I have deducted from experimental runs
             if _b6 == 1 and sim_score > 8: 
                 return 2
             else:
                return 1
    elif most_common_freq >= _b2:
        return 2

# With a position-matrix as input,
# calculate the frequency of the most common element
# Returns:
# -1 if nonconserved, gap
#  0  if nonconserved, no gap
#  1  if conserved
#  2  if highly conserved

def set_conservation(pos):
    max_residue = max(pos,key=pos.count) # Find the most common residue
    max_freq = pos.count(max_residue) # 
    
    new_pos = [p for p in pos if p != max_residue]
    
    second_max_residue = ''
    second_max_freq = 0
    similarity_score = 0

    if len(new_pos) > 0:
        second_max_residue = max(new_pos,key=new_pos.count)
    
    if max_residue != '-' and second_max_residue != '-' and second_max_residue != '':
        second_max_freq = new_pos.count(second_max_residue)
        similarity_score = sim[max_residue][second_max_residue] * second_max_freq
    
    if(any(v == '-' for v in pos)):
        return(-1)
    else:
        return classifier(max_freq, similarity_score)


# Remove Contiguous Nonconserved Positions < _b3

def remove_contiguous_nonconserved(pos):
    counter = 0
    for index, p in enumerate(pos):
        conservation = set_conservation(p)
        if(conservation <= 0):  # If positions is nonconserved
            counter += 1        # Increment by 1
        else:
            counter = 0         # Reset if position is conserved
        if(counter <= _b3):     # if not too many contiguous nonconserved positions
            conservation_level[index] = conservation # Set {Index: Conservation level} 
        elif(counter == _b3+1): # If too many contiguous nonconserved positions
            for k in range(index-_b3, index):
                del conservation_level[k] # Remove the last _b3 indices

conservation_level = {}
remove_contiguous_nonconserved(positions)


######## DEFINE BLOCKS #########
from itertools import groupby, count

# conservation_level    - a dict with all remaining positions' indices and their conservation level
# blocks_1              - one dict for each block (i.e. split by consecutive nonconserved regions)
# hc_positions          - all highly conserved positions in blocks
# flanks                - flanking (first and last) highly conserved position within a block
# blocks_2              - blocks anchored by flanking positions

# Looks at the index of a dictionary and splits it into multiple dictionaries
# wherever a non-consecutive index occurs
# {1234789} -> {1234, 789}
def group_by_consecutives(d):
    return [{index:d[index] for index in set(v)} for k, v in groupby(d, key=lambda i,j=count(): i-next(j))]

def set_flanks():
    # The conservation_level dictionary split at occurrences of contiguous nonconserved positions
    _blocks = group_by_consecutives(conservation_level)
    
    # Look for highly conserved positions (where value is 2) and keep the indices of these in a list 
    hc_positions = [[k for k,v in sorted(d.items()) if v == 2] for d in _blocks]
    
    # Extract the first and last HC in each block
    flanks = [(f[0], f[-1]) for f in hc_positions if len(f) > 1]
    
    # Redefine blocks anchored by flanking positions
    _blocks = [{key: value for key, value in conservation_level.items() if 
              (key >= start) and (key <= end)} for (start, end) in flanks]
    return _blocks

flanked_blocks = set_flanks()

all_gaps = [] # All gap positions

# Find gaps and nonconserved positions adjacent to gaps from a block
def find_gaps(list_of_blocks):
    for block in list_of_blocks:
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

find_gaps(flanked_blocks)

blocks = {}
    
for block in flanked_blocks:
    for k,v in block.items():
        if(k not in all_gaps):
            blocks[k] = v

# Remove small remaining blocks <= _b4
blocks_to_remove = []
blocks = group_by_consecutives(blocks)
for i, block in enumerate(blocks):
    if len(block) < _b4:
        blocks_to_remove.append(i)
blocks = [b for i,b in enumerate(blocks) if i not in blocks_to_remove]

new_positions = []

for block in blocks:
    for k,v in sorted(block.items()):
        new_positions.append(positions[k])

# This is the final alignment as positions in a 2D-list
# We need to convert these into sequence records
new_alignment = [list(i) for i in zip(*new_positions)]

for i,s in enumerate(alignment): # Iterate through the original alignment
    new_seq = Seq(''.join(new_alignment[i])) # take the corresponding position of the new alignment
    s.seq = new_seq

cwd = os.getcwd()

filepath = "output/"

#os.makedirs(os.path.dirname(filepath), exist_ok=True)
os.makedirs(filepath, exist_ok = True)

os.chdir('output/')

print("Current working directory changed from %s to %s" % (cwd,os.getcwd()))
print()

output_fna = args.infile.name + '.fna'    
SeqIO.write(alignment, output_fna, "fasta")

######## HTML OUTPUT #########

def to_html(matrix):
    
    #blosum = pd.read_csv('blosum62.csv')
    
    df = pd.DataFrame(matrix)
    
    from collections import defaultdict
    
    def color_most_freq(series):
        count_dict = defaultdict(int)
        for s in series.values:
            count_dict[s] += 1
        max_val = max(count_dict, key=count_dict.get)
        
        if count_dict[max_val] < _b1 or '-' in series.values:
            max_val = 0  
        
        # Color scheme based on http://acces.ens-lyon.fr/biotic/rastop/help/colour.htm
        
        colors = {0: 'black',
                  'D': '#E60A0A', # ASP Bright Red
                  'E': '#E60A0A', # GLU Bright Red
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
        
        if count_dict[max_val] < _b1 or '-' in series.values:
            max_val = 0  
        return ['font-weight: bold' if v == max_val else '' for v in series.values]
    
    def hover(hover_color="#f1eeee"):
        return dict(selector="tr:hover",
                    props=[("background-color", "%s" % hover_color)])
          
    def highlight_cols(x):
        d = x.copy()
        hb = [item+1 for sublist in blocks for item in sublist if item+1 in d.columns]
        d[hb] = 'background-color: #ded6d5'
        return d  
    
    html = []
    
    df.columns += 1
    #df.index += 1
    
    df.index = [alignment[i].id for i in range(len(df.index))]
    
    nr_columns = len(df.columns)
    
    col = [i for i in range(0,nr_columns) if (i % 60 == 0)]
    
    for i,c in enumerate(col):
        if(c == col[-1]):
            html.append(df[df.columns[c:]])
        else:
            html.append(df[df.columns[c:col[i+1]]])
    
    #.style.set_properties(**{'text-align': 'center'})\
    
    def apply_styles(table):
        
        styles = [
        hover(),
        dict(selector="th",
                     props=[("font-size", "5pt"),
                            ("text-align", "left")]),
        dict(selector="tbody",
                     props=[("font-size", "8pt"),
                            ("font-weight", "normal")])]
        
        
        for part in range(len(table)):
            table[part] = table[part].style\
            .set_table_styles(styles)\
            .apply(highlight_cols, axis=None)\
            .apply(color_most_freq, axis=0)\
            .apply(bold_most_freq, axis=0)\
            .render()
            print("Generating HTML output part %s/%s" % (part+1, len(table)+1))
    
        return table

    apply_styles(html)
    
    wrapper = \
    """
    <html>
        <head>
        </head>
        <body><style>table {width: 883px; border-collapse: collapse; font-family: Verdana;}
        </style>
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
    output_html = args.infile.name + '.html'
    
    fw = open(output_html,'w')
    
    flanks = [(list(sorted(k.keys()))[0]+1, list(sorted(k.keys()))[-1]+1) for k in blocks]
    
    mult_flanks = ''
    if len(flanks) > 1:
        mult_flanks = 's'
    
    pos_perc = str(round((len(new_positions)/len(positions)) * 100,0)) + '%'
    
    
    whole = wrapper % (_b1, _b2, _b3, _b4, len(flanks), mult_flanks, flanks,
                       len(new_positions), pos_perc, len(positions), '<p> </p>'.join(html))
    
    fw.write(whole)
    fw.close()
    

to_html(sequence_matrix)
