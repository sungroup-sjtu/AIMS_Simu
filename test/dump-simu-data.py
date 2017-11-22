#!/usr/bin/env python3

import sys

sys.path.append('..')

from app.models import *
from app.models_cv import *
from app.models_yaws import *

content = ''
for f in sys.argv[1:]:
    with open(f) as fin:
        content += fin.read()

molecules = []

for line in content.splitlines():
    line = line.strip()
    if line == '' or line.startswith('#'):
        continue

    mol = YawsMolecule.query.filter(YawsMolecule.cas == line.split()[2]).first()
    molecules.append(mol)

print('molecule,property,t,p,value,stderr,ff')
for mol in molecules:
    smiles = mol.smiles

    task = Task.query.filter(Task.smiles_list == '["%s"]' % smiles).first()
    if task is None:
        continue
    if task.post_result is None:
        continue

    cv = Cv.query.filter(Cv.cas == mol.cas).first()

    for job in task.jobs:
        T = job.t
        P = job.p

        post_result = task.get_post_result(T, P)
        if post_result['score-density'] < 0.999 or post_result['score-e_inter'] < 0.999:
            continue

        if cv is None:
            cp_sim = None
        else:
            cv_intra = cv.get_post_result(T)
            cp_sim = cv_intra + post_result['Cp_inter'] + post_result['Cp_PV']

        print('%s,%s,%i,%i,%f,,TEAM_LS' % (smiles, 'density', T, P*100000, post_result['density']))
        print('%s,%s,%i,%i,%f,,TEAM_LS' % (smiles, 'expansion', T, P*100000, post_result['expansion']))
        print('%s,%s,%i,%i,%f,,TEAM_LS' % (smiles, 'compressibility', T, P*100000, post_result['compressibility']))
        if cp_sim != None:
            print('%s,%s,%i,%i,%f,,TEAM_LS' % (smiles, 'Cp', T, P*100000, cp_sim))
