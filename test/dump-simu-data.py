#!/usr/bin/env python3

import sys
import math

sys.path.append('..')
sys.path.append('../../ms-tools')

from app import create_app
from app.models import Task
from app.models_cv import Cv
from mstools.formula import Formula
from mstools.utils import is_alkane


def get_npt_result():
    print('SMILES T(K) P(bar) density(g/mL) einter(kJ/mol) cp(J/mol.K) raw_density(g/mL) raw_einter(kJ/mol)')
    tasks = Task.query
    for task in tasks:
        smiles = task.get_smiles_list()[0]

        ### Only select molecules with n_heavy >= 4 and n_C >= 2
        m = task.get_mol_list()[0]
        f = Formula.read(m.formula)
        if m.OBMol.NumHvyAtoms() < 4:
            continue
        if f.atomdict.get('C') < 2:
            continue

        ### Only select alkanes
        # if not is_alkane(m):
        #     continue

        ### Only select hydrocarbon
        # if set(f.atomdict.keys()) != {'C', 'H'}:
        #     continue

        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue

        cv = Cv.query.filter(Cv.smiles == smiles).first()
        for i, (t, p, raw_density) in enumerate(post_result['density']):
            raw_einter = post_result['einter'][i][2]

            post_data = task.get_post_data(t, p)
            density = post_data['density']
            einter = post_data['einter']
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']

            if cv is None:
                ### Be careful, cv data not found
                cp = math.inf
            else:
                cv_intra = cv.get_post_data(t)
                cp = cv_intra + cp_inter + cp_pv

            print('%s %i %i %.3e %.3e %.3e %.3e %.3e' % (smiles, t, p, density, einter, cp, raw_density[0], raw_einter[0]))


def get_slab_result():
    f_critical = open('_tmp_critical.txt', 'w')
    f_critical.write('#SMILES Tc(K) Dc(g/mL)\n')
    print('SMILES T(K) Tc(K) Dc(g/mL) dliq(g/mL) dgas(g/mL) st(mN/m) raw_dliq(g/mL) raw_dgas(g/mL) raw_st(mN/m)')
    tasks = Task.query
    for task in tasks:
        smiles = task.get_smiles_list()[0]

        ### Only select molecules with n_heavy >= 4 and n_C >= 2
        m = task.get_mol_list()[0]
        f = Formula.read(m.formula)
        if m.OBMol.NumHvyAtoms() < 4:
            continue
        if f.atomdict.get('C') < 2:
            continue

        ### Only select alkanes
        # if not is_alkane(m):
        #     continue

        ### Only select hydrocarbon
        # if set(f.atomdict.keys()) != {'C', 'H'}:
        #     continue

        post_result = task.get_post_result()
        if post_result is None:
            continue

        post_data = task.get_post_data(0)
        tc = post_data['tc']
        dc = post_data['dc']
        f_critical.write('%s %.3e %.3e\n' % (smiles, tc, dc))

        for i, (t, raw_dliq) in enumerate(post_result['dliq']):
            raw_dgas = post_result['dgas'][i][1]
            raw_st = post_result['st'][i][1]

            post_data = task.get_post_data(t)
            dliq = post_data['dliq']
            dgas = post_data['dgas']
            st = post_data['st']
            print('%s %i %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e' % (
                smiles, t, tc, dc, dliq, dgas, st, raw_dliq[0], raw_dgas[0], raw_st[0]))


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
