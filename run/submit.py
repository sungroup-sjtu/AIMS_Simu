#!/usr/bin/env python3
# coding=utf-8

'''
Usage: submit.py mols.txt 'remark for this computation'
'''

import sys
import requests

json_dict = {
    'id': 1,
    'user_id': 1,
    'detail': {
        'procedures': ['npt'],
        'combinations': [],
        'p': [int(1E5), int(1E5)]
    }
}

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit()

with open(sys.argv[1]) as f:
    lines = f.read().splitlines()

remark = sys.argv[2]

for line in lines:
    line.strip()
    if line == '' or line.startswith('#'):
        continue
    words = line.split()
    name = words[0]
    smiles = words[1]
    t_min = int(float(words[2]))
    t_max = int(float(words[3]))
    p_min = int(float(words[4]))
    p_max = int(float(words[5]))

    json_dict['detail']['combinations'].append({'smiles': [smiles],
                                                'names': [name],
                                                't': [t_min, t_max],
                                                'p': [p_min, p_max]
                                                })
    json_dict['remark'] = remark

print(json_dict)
r = requests.post('http://localhost:5050/api/submit', json=json_dict)
print(r.text)
