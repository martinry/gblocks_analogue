# -*- coding: utf-8 -*-

from collections import defaultdict

msa = defaultdict(lambda: defaultdict(int))

# Input parser

name = ''
acc_id = 0
sequence = ''

with open('test.pir') as data:
    for i, line in enumerate(data):
        line = line.strip().replace(' ', '')
 
        if(line.startswith('>')): # Add accession to dict
            name = line
            acc_id = i
            msa[line]['Name'] = name.replace('>', '')
            
          
        elif(i == (acc_id + 1)): # Add description to dict
            msa[name]['Description'] = (line)
          
        else:
            sequence += line # Combine fragmented sequences
            if('*' in line): # until *
                sequence = sequence.replace('*', '')
                msa[name]['Sequence'] = (sequence) # Add sequence to dictionary
                sequence = '' # Reset sequence var



#for i, k,v in enumerate(input_alignment.items()):
 #   t.extend
    #print(k)
    #residue = list(v[1])
    #print(residue)


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


























