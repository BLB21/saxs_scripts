#!/opt/anaconda3/bin/python

import sys, os

if len(sys.argv) < 2:
    sys.exit('useage: PrintSeq.py XXX.pdb')
else:
    if os.path.isfile(sys.argv[1]) and sys.argv[1][-3:] == 'pdb':
        pass
    else:
        sys.exit(str(sys.argv[1])+' is not a valid pdb file')

file = open(sys.argv[1], 'r')
data = file.readlines()

translator={'ALA': 'A',\
            'CYS': 'C',\
            'ASP': 'D',\
            'GLU': 'E',\
            'PHE': 'F',\
            'GLY': 'G',\
            'HIS': 'H',\
            'ILE': 'I',\
            'LYS': 'K',\
            'LEU': 'L',\
            'MET': 'M',\
            'ASN': 'N',\
            'PRO': 'P',\
            'GLN': 'Q',\
            'ARG': 'R',\
            'SER': 'S',\
            'THR': 'T',\
            'VAL': 'V',\
            'TRP': 'W',\
            'TYR': 'Y',\
            'DA': 'A',\
            'DT': 'T',\
            'DC': 'C',\
            'DG': 'G'}

sequence = []
is_seqres = False
for line in data:
    line = line.split()
    if line[0] == 'SEQRES':
        is_seqres = True
        for part in line:
            if part in translator.keys():
                sequence.append(translator[part])
    elif line[0] == 'ATOM' and line[2] == 'CA' and is_seqres == False:
        if line[3] in translator.keys():
            sequence.append(translator[line[3]])
    elif line[0] == 'ATOM' and line[2] == 'P' and is_seqres == False:
        if line[3] in translator.keys():
            sequence.append(translator[line[3]])
    else:
        pass

window_size = 50
window = 0

while [ 1 ]:
    window_start = window
    window_end = window+window_size
    if window_end > len(sequence):
        print(''.join(sequence[window_start:]))
        break
    else:
        print(''.join(sequence[window_start:window_end]))
        window = window + window_size

print(f'{len(sequence)} residues')                


