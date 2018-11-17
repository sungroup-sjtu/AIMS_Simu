#!/usr/bin/env python3

import sys
import pybel

sys.path.append('..')
sys.path.append('../../ms-tools')

from app import create_app
from app.models import Task
from app.models_cv import Cv
from mstools.formula import Formula


def get_npt_result():
    print('#SMILES T(K) P(bar) density(g/mL) einter(kJ/mol) cp(J/mol.K) raw_density(g/mL) raw_einter(kJ/mol)')
    tasks = Task.query
    for task in tasks:
        smiles = task.get_smiles_list()[0]

        ### Ignore isomeric and bad molecules
        if task.remark in ['chiral', 'bad']:
            continue

        ### Only select alkanes
        if task.remark != 'alkane':
            continue

        ### Only select hydrocarbon with n_heavy > 3
        # m = pybel.readstring('smi', smiles)
        # if m.OBMol.NumHvyAtoms() <= 3:
        #     continue
        # if set(Formula.read(m.formula).atomdict.keys()) != {'C', 'H'}:
        #     continue

        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue

        for i, (t, p, raw_density) in enumerate(post_result['density']):
            raw_einter = post_result['einter'][i][2]

            post_data = task.get_post_data(t, p)
            density = post_data['density']
            einter = post_data['einter']
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']

            cv = Cv.query.filter(Cv.smiles == smiles).first()
            if cv is None:
                ### Be careful, cv data not found
                cp = 0
            else:
                cv_intra = cv.get_post_data(t)
                cp = cv_intra + cp_inter + cp_pv

            print('%s %i %i %.3e %.3e %.3e %.3e %.3e' % (smiles, t, p, density, einter, cp, raw_density, raw_einter))


def get_slab_result():
    print('#SMILES T(K) Tc(K) Dc(g/mL) dliq(g/mL) dgas(g/mL) st(mN/m) raw_dliq(g/mL) raw_dgas(g/mL) raw_st(mN/m)')
    tasks = Task.query
    for task in tasks:
        smiles = task.get_smiles_list()[0]

        ### Ignore isomeric and bad molecules
        if task.remark in ['chiral', 'bad']:
            continue

        ### Only select alkanes
        # if task.remark != 'alkane':
        #     continue

        ## Only select hydrocarbon with n_heavy > 3
        m = pybel.readstring('smi', smiles)
        if m.OBMol.NumHvyAtoms() <= 3:
            continue
        if set(Formula.read(m.formula).atomdict.keys()) != {'C', 'H'}:
            continue

        post_result = task.get_post_result()
        if post_result is None:
            continue

        for job in task.jobs:
            if not job.converged:
                continue
        for i, (t, raw_dliq) in enumerate(post_result['dliq']):
            raw_dgas = post_result['dgas'][i][1]
            raw_st = post_result['st'][i][1]

            post_data = task.get_post_data(t)
            tc = post_data['tc']
            dc = post_data['dc']

            dliq = post_data['dliq']
            dgas = post_data['dgas']
            st = post_data['st']
            print('%s %i %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e' % (
                smiles, t, tc, dc, dliq, dgas, st, raw_dliq, raw_dgas, raw_st))


if __name__ == '__main__':
    procedure = sys.argv[1]

    app = create_app(procedure)
    app.app_context().push()

    if procedure == 'npt':
        get_npt_result()
    elif procedure == 'nvt-slab':
        get_slab_result()
    else:
        raise Exception('Invalid procedure: ' + procedure)
