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
    t_fus = words[4]
    t_vap = words[5]
    t_c = words[6]

    if t_vap == 'None':
        print('!ERROR: Tvap is None: %s' % line)
        continue

    t_vap = int(round(float(t_vap)))

    if t_fus == 'None':
        print('!WARNING: Tfus is None: %s' % line)
        t_min = 298
    else:
        t_fus = int(round(float(t_fus)))
        t_min = t_fus + 25

    if t_c == 'None':
        print('!WARNING: Tc is None: %s' % line)
        t_max = t_vap + 100
    else:
        t_c = int(round(float(t_c)))
        t_max = t_c - 100

    if t_min > t_max:
        print('!ERROR: t_min > t_max: %s' % line)
        continue

    json_dict['detail']['combinations'].append({'smiles': [smiles],
                                                'names': [name],
                                                't': [t_min, t_max],
                                                })

print(json_dict)
r = requests.post('http://localhost:5050/api/submit', json=json_dict)
print(r.text)
