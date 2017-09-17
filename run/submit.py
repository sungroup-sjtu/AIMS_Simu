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
        'p': [1, 1000]  # bar
    }
}

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit()

with open(sys.argv[1]) as f:
    lines = f.read().splitlines()

json_dict['remark'] = sys.argv[2]

for line in lines:
    line.strip()
    if line == '' or line.startswith('#'):
        continue
    words = line.split()
    name = words[2]
    smiles = words[3]
    t_min = int(round(float(words[4])))
    t_max = int(round(float(words[5])))

    if t_min > t_max:
        print('!Not Good', line)
        continue

    json_dict['detail']['combinations'].append({'smiles': [smiles],
                                                'names': [name],
                                                't': [t_min, t_max],
                                                })

print(json_dict)
r = requests.post('http://localhost:5050/api/submit', json=json_dict)
print(r.text)
