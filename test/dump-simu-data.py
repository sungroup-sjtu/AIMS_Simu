#!/usr/bin/env python3

import sys
import pybel

sys.path.append('..')
sys.path.append('../../ms-tools')

from app.models import *
from app.models_cv import Cv
from mstools.formula import Formula

npt = init_simulation('npt')

molecules = []

print(
        '#SMILES T(K) P(bar) density(g/mL) einter(kJ/mol) expansion(/K) compressibility(/bar) cp(J/mol.K) raw_density raw_einter raw_expansion raw_compressibility')
tasks = Task.query.filter(Task.procedure == 'npt')

# for task in tasks.limit(2):
for task in tasks:
    smiles = json.loads(task.smiles_list)[0]
    m = pybel.readstring('smi', smiles)
    if set(Formula.read(m.formula).atomdict.keys()) != {'C', 'H'}:
        continue

    if task.post_result is None:
        continue

    cv = Cv.query.filter(Cv.smiles == smiles).first()

    post_result = json.loads(task.post_result)
    if post_result['density-poly4-score'] < 0.999 or post_result['e_inter-poly4-score'] < 0.999:
        continue

    for job in task.jobs:
        T = job.t
        P = job.p
        if not job.converged:
            continue

        post_result = task.get_post_data(T, P)
        density = post_result['density']
        e_inter = post_result['e_inter']
        expansion = post_result['expansion']
        compressibility = post_result['compressibility']

        raw_result = json.loads(job.result)
        begin = raw_result['converged_from']
        raw_density = raw_result['density'][0]
        raw_e_inter = raw_result['e_inter'][0] / json.loads(task.n_mol_list)[0]
        raw_expansion = raw_result['expansion'][0]
        raw_compressibility = raw_result['compressibility'][0]

        if cv is None:
            cp_sim = 0  ### TODO Be careful, cv data not found
        else:
            cv_intra = cv.get_post_data(T)
            cp_sim = cv_intra + post_result['cp_inter'] + post_result['cp_pv']

        print('%s %i %i %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e %.2e' % (
            smiles, T, P,
            density, e_inter, expansion, compressibility, cp_sim,
            raw_density, raw_e_inter, raw_expansion, raw_compressibility))
