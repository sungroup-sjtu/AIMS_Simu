#!/usr/bin/env python3
# coding=utf-8

import sys
import json

sys.path.append('..')
from app import create_app
from app.api.actions import ComputeAction
from mstools.utils import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to generate a Compute and corresponding Tasks in high-throughput simulation')
parser.add_argument('-i', '--input', type=str, help='Input txt file: name smiles mol_ratio')
parser.add_argument('-p', '--procedure', type=str, help='Procedure line: npt, nvt-slab, ppm')
parser.add_argument('-r', '--remark', type=str, help='remark of this compute')
parser.add_argument('-tp', '--temppresstyle', type=str, help='The way to define t_list and p_list')
parser.add_argument('-a', '--altertask', type=bool, help='Alter t_list and p_list of exist task', default=False)
opt = parser.parse_args()


procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()
app.test_request_context().push()


def submit(json_dict):
    computeAction = ComputeAction()
    try:
        compute_id = computeAction.init_from_json_dict(json_dict) # set smiles_repeat True to
    except Exception as e:
        return json.dumps({'success': False,
                           'reason' : str(e)
                           })
    else:
        return json.dumps({'success'   : True,
                           'compute_id': compute_id,
                           })


json_dict = {
    'id': 1,
    'user_id': 1,
    'detail': {
        'procedures': [procedure],
        'combinations': [],
    },
    'remark': opt.remark,
    'alter_task': opt.altertask
}

with open(opt.input) as f:
    lines = f.read().splitlines()

for line in lines:
    line.strip()
    if line == '' or line.startswith('#'):
        continue
    words = line.split()
    name = words[0].split('.')
    smiles = words[1].split('.')
    n_mol_ratio = list(map(int, words[2].split('.')))
    if opt.temppresstyle == 'XY':
        t_min = 300
        t_max = 450
        T_list = get_T_list_from_range(t_min, t_max, 4)
        P_list = [1]
        json_dict['detail']['combinations'].append({'smiles': smiles,
                                                    'names': name,
                                                    'n_mol_ratio': n_mol_ratio,
                                                    't_list': T_list,
                                                    'p_list': P_list,
                                                    })
    elif opt.temppresstyle == 'fixed':
        T_list = [37.78, 60.0, 98.89]
        for i in range(len(T_list)):
            T_list[i] += 273.15
        P_list = [1, 200, 600, 1800, 3000, 6000, 7000, 8000]
        json_dict['detail']['combinations'].append({'smiles': smiles,
                                                    'names': name,
                                                    'n_mol_ratio': n_mol_ratio,
                                                    't_list': T_list,
                                                    'p_list': P_list,
                                                    })
    elif opt.temppresstyle == 'assigned':
        T_list = list(map(int, words[3].split(',')))

        if words[4] == 'None':
            P_list = []
        else:
            P = words[4].split(',')
            P_list = [int(p) for p in P]

        json_dict['detail']['combinations'].append({'smiles': smiles,
                                                'names': name,
                                                'n_mol_ratio': n_mol_ratio,
                                                't_list': T_list,
                                                'p_list': P_list,
                                                })
    elif opt.temppresstyle == 'nvt-slab':
        T_list = list(map(int, words[3].split(',')))
        P_list = []

        json_dict['detail']['combinations'].append({'smiles': smiles,
                                                'names': name,
                                                'n_mol_ratio': n_mol_ratio,
                                                't_list': T_list,
                                                'p_list': P_list,
                                                })
    elif opt.temppresstyle == 'prior':
        json_dict['detail']['combinations'].append({'smiles': smiles,
                                                'names': name,
                                                'n_mol_ratio': n_mol_ratio,
                                                't_list': [],
                                                'p_list': [],
                                                })

if __name__ == '__main__':
    print(json_dict)
    r = submit(json_dict)
    print(r)
