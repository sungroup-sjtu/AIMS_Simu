#!/usr/bin/env python3

import sys
import math

sys.path.append('..')
sys.path.append('../../ms-tools')

from app import create_app
from app.models import Task
from app.models_cv import Cv
from mstools.formula import Formula

app = create_app('npt')
app.app_context().push()

smiles_list = []
category_list = []
ml_list = []
t_list = []
p_list = []
den_list = []
ei_list = []
cp_list = []
hvap_list = []
comp_list = []
expa_list = []


def get_npt_result(limit=None):
    print('\nGet NPT results...\n')
    tasks = Task.query.filter(Task.remark == None)
    if limit is not None:
        tasks = tasks.limit(limit)
    for i, task in enumerate(tasks):
        sys.stdout.write(f'\r\t{i}')

        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue

        smiles = task.get_smiles_list()[0]
        py_mol = task.get_mol_list()[0]
        f = Formula(py_mol.formula)
        n_CON = f.atomdict.get('C', 0) + f.atomdict.get('O', 0) + f.atomdict.get('N', 0)

        category = 'All'
        if set(f.atomdict.keys()) == {'C', 'H'}:
            category = 'CH'

        ### For ML, Only select molecules with n_heavy >= 4
        ml = True
        if py_mol.OBMol.NumHvyAtoms() < 4:
            ml = False

        cv = Cv.query.filter(Cv.smiles == smiles).first()
        for i, (t, p, raw_density) in enumerate(post_result['density']):
            raw_einter = post_result['einter'][i][2]
            raw_compress = post_result['compress'][i][2]

            post_data = task.get_post_data(t, p)
            density = post_data['density']
            einter = post_data['einter']
            hvap = post_data['hvap']
            compress = post_data['compress']
            expansion = post_data['expansion']
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']

            smiles_list.append(smiles)
            category_list.append(category)
            ml_list.append(ml)
            t_list.append(t)
            p_list.append(p)
            den_list.append(raw_density)
            ei_list.append(raw_einter)
            ### Hvap = RT - Ei - n_CON/15 * RT
            hvap_list.append([(1 - n_CON / 15) * 8.314 * t / 1000 - raw_einter[0], raw_einter[1]])
            comp_list.append(raw_compress)
            expa_list.append(expansion)

            if None in [cv, cp_inter, cp_pv]:
                ### Be careful, cv data not found
                cp = None
            else:
                cv_intra = cv.get_post_data(t)
                cp = cv_intra + cp_inter + cp_pv

            cp_list.append(cp)


def write_data(n_mol_per_file=int(1e9)):
    i_mol = 0
    last_smiles = ''
    for i, smiles in enumerate(smiles_list):
        if smiles != last_smiles:
            i_mol += 1
            if (i_mol - 1) % n_mol_per_file == 0:
                f_All = open(f'result-All-npt-{i_mol}.csv', 'w')
                print('molecule,t,p,property,value,uncertainty,ff', file=f_All)
            last_smiles = smiles

        print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'density-l', *den_list[i]), file=f_All)
        print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'einter-l', *ei_list[i]), file=f_All)
        print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'compressibility-l', *comp_list[i]), file=f_All)
        if expa_list[i] is not None:
            print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'expansion-l', expa_list[i], math.nan), file=f_All)
        if cp_list[i] is not None:
            print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'cp-l', cp_list[i], math.nan), file=f_All)


def write_ml_data():
    f_ml_All = open('result-ML-All-npt.txt', 'w')
    f_ml_CH = open('result-ML-CH-npt.txt', 'w')
    f_ml_cp_All = open('result-ML-All-cp.txt', 'w')
    f_ml_cp_CH = open('result-ML-CH-cp.txt', 'w')
    f_ml_expan_All = open('result-ML-All-expan.txt', 'w')
    f_ml_expan_CH = open('result-ML-CH-expan.txt', 'w')
    print('SMILES T P density einter compress', file=f_ml_All)
    print('SMILES T P density einter compress', file=f_ml_CH)
    print('SMILES T P cp', file=f_ml_cp_All)
    print('SMILES T P cp', file=f_ml_cp_CH)
    print('SMILES T P expan', file=f_ml_expan_All)
    print('SMILES T P expan', file=f_ml_expan_CH)

    for i, smiles in enumerate(smiles_list):
        if not ml_list[i]:
            continue

        print('%s %i %i %.3e %.3e %.3e' % (smiles, t_list[i], p_list[i], den_list[i][0], ei_list[i][0], comp_list[i][0]), file=f_ml_All)
        if expa_list[i] is not None:
            print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], expa_list[i]), file=f_ml_expan_All)
        if cp_list[i] is not None:
            print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], cp_list[i]), file=f_ml_cp_All)

        if category_list[i] == 'CH':
            print('%s %i %i %.3e %.3e %.3e' % (smiles, t_list[i], p_list[i], den_list[i][0], ei_list[i][0], comp_list[i][0]), file=f_ml_CH)
            if expa_list[i] is not None:
                print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], expa_list[i]), file=f_ml_expan_CH)
            if cp_list[i] is not None:
                print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], cp_list[i]), file=f_ml_cp_CH)


if __name__ == '__main__':
    get_npt_result(limit=100)
    write_data(n_mol_per_file=1000)
    write_ml_data()
