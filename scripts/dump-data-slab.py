#!/usr/bin/env python3

import sys
import math

sys.path.append('..')
sys.path.append('../../ms-tools')

from app import create_app
from app.models import Task
from mstools.formula import Formula
from mstools.utils import is_alkane

app = create_app('nvt-slab')
app.app_context().push()

smiles_list = []
category_list = []
ml_list = []
t_list = []
dl_list = []
dg_list = []
st_list = []
tc_list = []
dc_list = []


def get_slab_result(limit=None):
    tasks = Task.query
    if limit is not None:
        tasks = tasks.limit(limit)
    for i, task in enumerate(tasks):
        sys.stdout.write(f'\r\t{i}')

        post_result = task.get_post_result()
        if post_result is None:
            continue

        smiles = task.get_smiles_list()[0]
        py_mol = task.get_mol_list()[0]
        f = Formula(py_mol.formula)
        n_C = f.atomdict['C']

        category = 'All'
        if set(f.atomdict.keys()) == {'C', 'H'}:
            category = 'CH'
        if is_alkane(py_mol):
            category = 'Ane'

        ### For ML, Only select molecules with n_heavy >= 4 and n_C >= 2
        ml = True
        if py_mol.OBMol.NumHvyAtoms() < 4:
            ml = False
        elif f.atomdict.get('C') < 2:
            ml = False

        post_data = task.get_post_data(0)
        tc = post_data['tc']
        dc = post_data['dc']

        for i, (t, raw_dliq) in enumerate(post_result['dliq']):
            raw_dgas = post_result['dgas'][i][1]
            raw_st = post_result['st'][i][1]

            smiles_list.append(smiles)
            category_list.append(category)
            ml_list.append(ml)
            t_list.append(t)
            dl_list.append(raw_dliq)
            dg_list.append(raw_dgas)
            st_list.append(raw_st)

            tc_list.append(tc)
            dc_list.append(dc)


def write_data(n_mol_per_file=int(1e9)):
    i_mol = 0
    last_smiles = ''
    for i, smiles in enumerate(smiles_list):
        if smiles != last_smiles:
            i_mol += 1
            if (i_mol - 1) % n_mol_per_file == 0:
                f_All = open(f'result-All-slab-{i_mol}.csv', 'w')
                print('molecule,t,p,property,value,uncertainty,ff', file=f_All)

            print('%s,%f,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, math.nan, math.nan, 'tc', tc_list[i], math.nan), file=f_All)
            print('%s,%f,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, math.nan, math.nan, 'dc', dc_list[i], math.nan), file=f_All)
            last_smiles = smiles

        print('%s,%i,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], math.nan, 'st-lg', *st_list[i]), file=f_All)
        print('%s,%i,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], math.nan, 'density-lg', *dl_list[i]), file=f_All)
        if dg_list[i][0] >= 0.01:
            print('%s,%i,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], math.nan, 'density-gl', *dg_list[i]), file=f_All)


def write_ml_data():
    f_ml_All = open('result-ML-All-slab.txt', 'w')
    f_ml_Ane = open('result-ML-Ane-slab.txt', 'w')
    f_ml_CH = open('result-ML-CH-slab.txt', 'w')
    f_ml_critical_All = open('result-ML-All-critical.txt', 'w')
    f_ml_critical_Ane = open('result-ML-Ane-critical.txt', 'w')
    f_ml_critical_CH = open('result-ML-CH-critical.txt', 'w')
    print('SMILES T dliq dgas st', file=f_ml_All)
    print('SMILES T dliq dgas st', file=f_ml_Ane)
    print('SMILES T dliq dgas st', file=f_ml_CH)
    print('SMILES tc dc', file=f_ml_critical_All)
    print('SMILES tc dc', file=f_ml_critical_Ane)
    print('SMILES tc dc', file=f_ml_critical_CH)

    smiles_last_All = ''
    smiles_last_Ane = ''
    smiles_last_CH = ''
    for i, smiles in enumerate(smiles_list):
        if not ml_list[i]:
            continue

        print('%s %i %.3e %.3e %.3e' % (smiles, t_list[i], dl_list[i][0], dg_list[i][0], st_list[i][0]), file=f_ml_All)

        if smiles != smiles_last_All:
            print('%s %.3e %.3e' % (smiles, tc_list[i], dc_list[i]), file=f_ml_critical_All)
            smiles_last_All = smiles

        if category_list[i] == 'Ane':
            print('%s %i %.3e %.3e %.3e' % (smiles, t_list[i], dl_list[i][0], dg_list[i][0], st_list[i][0]), file=f_ml_Ane)

            if smiles != smiles_last_Ane:
                print('%s %.3e %.3e' % (smiles, tc_list[i], dc_list[i]), file=f_ml_critical_Ane)
                smiles_last_Ane = smiles

        if category_list[i] in ['Ane', 'CH']:
            print('%s %i %.3e %.3e %.3e' % (smiles, t_list[i], dl_list[i][0], dg_list[i][0], st_list[i][0]), file=f_ml_CH)

            if smiles != smiles_last_CH:
                print('%s %.3e %.3e' % (smiles, tc_list[i], dc_list[i]), file=f_ml_critical_CH)
                smiles_last_CH = smiles


if __name__ == '__main__':
    get_slab_result(limit=None)
    write_data()
    write_ml_data()
