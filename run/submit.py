#!/usr/bin/env python3
# coding=utf-8

'''
Usage: submit.py procedure mols.txt 'remark for this computation'
'''

import sys
import json

sys.path.append('..')

from app import create_app
from app.api.actions import ComputeAction

procedure = sys.argv[1]
app = create_app(procedure)
app.app_context().push()
app.test_request_context().push()


def submit(json_dict):
    computeAction = ComputeAction()
    try:
        compute_id = computeAction.init_from_json_dict(json_dict)
    except Exception as e:
        return json.dumps({'success': False,
                           'reason' : str(e)
                           })
    else:
        return json.dumps({'success'   : True,
                           'compute_id': compute_id,
                           })


if len(sys.argv) != 4:
    print(__doc__)
    sys.exit()

json_dict = {
    'id'     : 1,
    'user_id': 1,
    'detail' : {
        'procedures'  : [procedure],
        'combinations': [],
        'p'           : [1, 1000]  # bar. This option is useless
    },
    'remark' : sys.argv[3]
}

with open(sys.argv[2]) as f:
    lines = f.read().splitlines()

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
    t_vap = float(t_vap)

    if t_fus == 'None':
        print('!WARNING: Tfus is None: %s' % line)
        t_min = int(round(t_vap * 0.4 + 100))
    else:
        t_fus = float(t_fus)
        t_min = int(round(t_fus + 25))

    if t_c == 'None':
        print('!WARNING: Tc is None: %s' % line)
        t_max = int(round(t_vap * 1.2))
    else:
        t_c = float(t_c)
        t_max = int(round(t_c * 0.85))

    t_max = min(t_max, 650)

    if t_min >= t_max:
        print('!ERROR: t_min > t_max: %s' % line)
        continue

    json_dict['detail']['combinations'].append({'smiles': [smiles],
                                                'names' : [name],
                                                't'     : [t_min, t_max],
                                                })

if __name__ == '__main__':
    print(json_dict)
    r = submit(json_dict)
    print(r)
