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

# -------------------------------------------------------ME-------

seq_objects = []

#parser('test.pir', '.pir')
parser('nad3.pir', '.pir')

seq_matrix = [seq_objects[n].sequence for n in range( len(seq_objects) )]

k = [[f for f in s] for s in seq_matrix]

data_frame = pd.DataFrame(k)
df = data_frame.transpose()






t = []

seq_length = len(seq_objects[0].sequence)

for i in range(0, seq_length):
    residues = [n.sequence[i] for n in seq_objects]
    t.append(residues)
    
    


#for obj in seq_objects:
#    print(obj.name)------ME



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

conservation_level = ['n', 'c', 'h']

"""
def classifier(col, eq, most_abundant, max_freq, perc_freq, gaps):
    if(perc_freq < _is or gaps):
        return conservation_level[0]
    
    elif((perc_freq > _is) and (perc_freq < _fs)):
        return conservation_level[1]
    
    elif(perc_freq > _fs):
        return conservation_level[2]
"""

def classifier(is_gap, frequency_perc):
    if(frequency_perc < _is or is_gap):
        return conservation_level[0]
    
    elif((frequency_perc > _is) and (frequency_perc < _fs)):
        return conservation_level[1]
    
    elif(frequency_perc > _fs):
        return conservation_level[2]
    
    

####################################


#print('### STATS ###')
#print('Number of sequences: %d' % nr_seqs)
#print('\n')
#print('%s\t%s\t%s\t%s\t%s\t%s'%
#      ('Col.', 'All_eq', 'Most_ab', 'Max_frq', 'Frq %', 'Gaps'))



temp = []



for i, residue in enumerate(t):
    #for i in (len(residue)):
    cons = []
#    freq = []
    
    mc_count = Counter(residue)
    mc = mc_count.most_common(1)
    mc = mc[0][0]
    #mc = max(x for x in residue)
    eq = all(x==residue[0] for x in residue)
    rc = residue.count(mc)
    fp = (rc / len(residue)) * 100
    gap = ('-' in residue)
    
    for x in residue:
        frequency = residue.count(x)
        frequency_perc = (frequency / len(residue)) * 100
        is_gap = (x == '-')

        conservation = classifier(is_gap, frequency_perc)
        
           
        cons.extend(conservation)
    
    
        
#        freq.append(round(frequency_perc))
#    print(i+1, cons)
#    print(i+1, freq)
#    print(i+1, residue)
    
    
    
"""    

    #at = all(x==residue[0] for x in residue)
    print('%s\t%s\t%s\t%s\t%s\t%s\t%s'%
          (i,
          eq, # are all values equal?
          mc, # most common amino acid
          rc, # occurrence of most common amino acid
          fp, # occurrence %
          gap,# do gaps occur?
          conservation))
    print(residue)
"""    

    

####################################
########## CLASSIFICATION ##########
####################################
'''

IS, FS, CP, BL1, BL2
# nonconserved      - <IS identical residues or gap
# conserved         - ≥IS and <FS identical residues
# highly conservd   - ≥FS identical residues

Default values:
IS  - 50% of nr of sequences + 1
FS  - 85% of nr of sequences
CP  - 8 positions
BL1 - 15 positions
BL2 - 10 positions


All stretches of contiguous nonconserved positions >CP are rejected

'''


#for i in range(len(seq_objects)):
#    for x in seq_objects[i].sequence:
#        print(x)













