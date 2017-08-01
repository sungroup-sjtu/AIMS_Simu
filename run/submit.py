#!/usr/bin/env python3
# coding=utf-8

import sys
import requests

json_dict = {
    'id': 1,
    'user_id': 1,
    'detail': {
        'procedures': ['npt'],
        'smiles_list': [[]],
        'mol_names': [],
        'states': [],
    }
}

with open(sys.argv[1]) as f:
    lines = f.read().splitlines()

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

    json_dict['detail']['smiles_list'][0].append(smiles)
    json_dict['detail']['mol_names'].append(name)
    json_dict['detail']['states'].append({'t_min': t_min, 't_max': t_max, 'p_min': p_min, 'p_max': p_max})

print(json_dict)
r = requests.post('http://localhost:5050/api/submit', json=json_dict)
print(r.text)
