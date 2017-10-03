#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 09:35:26 2017

@author: martin
"""
from collections import Counter
import pandas as pd

class Sequence(object):
    
    def __init__(self, name='', description='', sequence=''):
        self.name = name
        self.description = description
        self.sequence = sequence
        
# Input parser (pir)
"""
Assume our input data looks like the following:
    
SEQUENCES MUST BE OF EQUAL LENGTH

>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR 
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH 
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ 
  SDVPERSIPI TREEKPAIAG AQRKK*

>P1;CRAB_BOVIN
ALPHA CRYSTALLIN C CHAIN (ALPHA(C)-CRYSTALLIN).
  MDIAIHHPWI RRPFFPFHSP SRLFDQFFGE HLLESDLFPA STSLSPFYLR 
  PPSFLRAPSW IDTGLSEMRL EKDRFSVNLD VKHFSPEELK VKVLGDVIEV 
  HGKHEERQDE HGFISREFHR KYRIPADVDP LAITSSLSSD GVLTVNGPRK 
  QASGPERTIP ITREEKPAVT AAPKK*
...

"""

def parser(file_name, file_type):
    #accession_number = ''
    accession_id = 0
    tmp_seq = ''
    
    with open(file_name) as data:
        for i, line in enumerate(data):
            line = line.strip().replace(' ', '')
     
            if(line.startswith('>')):
                accession_id = i
                seq = Sequence(name=line)
                
            elif(i == (accession_id + 1)):
                seq.description = line
              
            else:
                tmp_seq += line # Combine fragmented sequences
                if('*' in line): # until *
                    tmp_seq = tmp_seq.replace('*', '')
                    seq.sequence = tmp_seq
                    tmp_seq = '' # Reset sequence var
                    
                    seq_objects.append(seq)

# --------------------------------------------------------------
                    

seq_objects = []

#parser('test.pir', '.pir')
parser('nad3.pir', '.pir')

seq_matrix = [seq_objects[n].sequence for n in range( len(seq_objects) )]

k = [[f for f in s] for s in seq_matrix]


####################################
############ DEFINITIONS ###########
####################################

nr_seqs = len(seq_objects)

# IS  - 50% of nr of sequences + 1
_is = int((nr_seqs * 0.5) + 1)
# FS  - 85% of nr of sequences
_fs = int((nr_seqs * 0.85))
_cp = 8
_bl1 = 15
_bl2 = 10

# Example:
# Nr seqs:  17
# IS:       9.5
# FS:       14.45


####################################


conservation_level = [0, 1, 2]
#nc = [0] * nr_seqs # The number of consecutive 'nonconserved' for a residue in the last row

conservation_index = []

#new_df = pd.DataFrame()

def classifier(freq, most_common_freq):
    """
    if(freq < _is):
        return conservation_level[0]
    
    elif((freq >= _is) and (freq < _fs) and (most_common_freq < _is)):
        return conservation_level[1]

    elif((freq >= _fs) and (most_common_freq >= _is)):
        return conservation_level[2]
    
    #else:
    #    return conservation_level[1]
    """
    
    if(most_common_freq < _is):
        return conservation_level[0]
    elif((most_common_freq >= _is) and (most_common_freq < _fs)):
        return conservation_level[1]
    elif(most_common_freq >= _fs):
        return conservation_level[2]
    
    

def deep_thought(index, row):
    most_common = row.value_counts().idxmax()
    most_common_freq = (row.value_counts().max())
#    print(row.name, most_common, freq)
#    print(row.name, row.values, most_common, freq, classifier(most_common, freq))
    #print(r)
    #c = row.value_counts()
    f = row.value_counts()
    freq = f[f > 1].sum()
    #print(most_common, most_common_freq)
    #print(row)
    #print(c)
    #print(f)
    #print(freq)
    
    #print('%s\t%s\t%s\t%s\t%s' % (index+1, most_common, most_common_freq, freq, row.values))
    
    
    if(any(v == '-' for v in row.values)):
        return(0)
    else:
        return classifier(freq, most_common_freq)
    


"""
def deep_thought(index, row):
    global nc

    for i, item in enumerate(row):
        freq = row.value_counts()[item] # Calculate frequency
        freq = (freq / len(row)) * 100  # Calculate frequency percentage
        
        cons = classifier(item, freq)   # Determine level of conservation
        #print(item, freq, cons)
        
        if(cons == 'n'):                # If current residue is nonconserved
            nc[i] += 1        # Increment counter for that position
        else:
            nc[i] = 0
        

    if(all(v <= _cp for v in nc)):
        print(nc)
        return temp.append(row.tolist())

#    if(any(v > _cp for v in nc)): # If ANY position has 9 nonconserved residues
#        nc = [0] * nr_seqs        # reset counter   
"""

# Untransposed dataframe
data_frame = pd.DataFrame(k)
# Transposed dataframe
df = data_frame.T #.transpose()

# Lookup dataframe, used to store information about each df entry
# in dictionaries
df_rows = df.shape[0]
df_cols = df.shape[1]

#print('%s\t%s\t%s\t%s\t%s' % ('Index', 'MC', 'MCFreq', 'Freq', 'Row values'))
for index,row in df.iterrows():
    conservation_index.append(deep_thought(index, row))
        

#new_df = new_df.append(temp)

"""
Rec for HC:
Most common amino acid must be > 13    
    
Req for C:
the most commmon amino acid must be >8
the combined amount of identical amino acids must be >13




"""


new = []
new2 = []

counter = 0

for i, x in enumerate(conservation_index):
    if(x == 0):
        counter += 1
    else:
        counter = 0
        
    if(counter < 9):
        new.append([i,x])
        new2.extend([i])
    elif(counter == 9):
        new = new[:-8]
        new2 = new2[:-8]

# Transform nested list into dictionary
n = {h[0] : h[1] for h in new}

# Transform dictionary into list of dictionaries based on blocks,
# new block for every non-consecutive index
from itertools import groupby, count
#blocks = [list(b) for a, b in groupby(n, key=lambda i,j=count(): i-next(j))]
blocks = [{s:n[s] for s in set(v)} for k, v in groupby(n, key=lambda i,j=count(): i-next(j))]

# Clean up this mess
# Look for HC values (=2) and keep the indices of these in a list 
r = [[(k) for k,v in (d.items()) if v == 2] for d in blocks]
# Extract the first and last HC in each block
values = [(val[0], val[-1]) for val in r if len(val) > 1]

u = [{key: value for key, value in n.items() if (key >= start) and (key <= end)} for (start, end) in values]



"""
tempnew = {}
for block in values:
    start = block[0]
    end = block[1]
    
    u = [{key: value} for key, value in n.items() if (key >= start) and (key <= end)]
    print(u)
"""

"""
#b = pd.DataFrame(blocks)
for d in blocks:
    tmp = []
    r = [(k) for k,v in (d.items()) if v == 2]
    print(r)
    if(len(r) > 1):
        print(r[0], r[-1])
    
"""



        

