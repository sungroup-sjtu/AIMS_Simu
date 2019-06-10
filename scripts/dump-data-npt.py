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
    tasks = Task.query
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
            # Hvap = RT - Ei - n_C/15 * RT
            hvap_list.append([(1 - n_C / 15) * 8.314 * t / 1000 - raw_einter[0], raw_einter[1]])
            comp_list.append(raw_compress)
            expa_list.append(expansion)

            if cv is None:
                ### Be careful, cv data not found
                cp = math.inf
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
        print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'expansion-l', expa_list[i], math.nan), file=f_All)

        if cp_list[i] != math.inf:
            print('%s,%i,%i,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], p_list[i], 'cp-l', cp_list[i], math.nan), file=f_All)

        if p_list[i] == 1:
            print('%s,%i,%f,%s,%.3e,%.1e,TEAM_MGI' % (smiles, t_list[i], math.nan, 'hvap-lg', *hvap_list[i]), file=f_All)


def write_ml_data():
    f_ml_All = open('result-ML-All-npt.txt', 'w')
    f_ml_Ane = open('result-ML-Ane-npt.txt', 'w')
    f_ml_CH = open('result-ML-CH-npt.txt', 'w')
    f_ml_hvap_All = open('result-ML-All-hvap.txt', 'w')
    f_ml_hvap_Ane = open('result-ML-Ane-hvap.txt', 'w')
    f_ml_hvap_CH = open('result-ML-CH-hvap.txt', 'w')
    f_ml_cp_All = open('result-ML-All-cp.txt', 'w')
    f_ml_cp_Ane = open('result-ML-Ane-cp.txt', 'w')
    f_ml_cp_CH = open('result-ML-CH-cp.txt', 'w')
    print('SMILES T P density einter compress expansion', file=f_ml_All)
    print('SMILES T P density einter compress expansion', file=f_ml_Ane)
    print('SMILES T P density einter compress expansion', file=f_ml_CH)
    print('SMILES T hvap', file=f_ml_hvap_All)
    print('SMILES T hvap', file=f_ml_hvap_Ane)
    print('SMILES T hvap', file=f_ml_hvap_CH)
    print('SMILES T P cp', file=f_ml_cp_All)
    print('SMILES T P cp', file=f_ml_cp_Ane)
    print('SMILES T P cp', file=f_ml_cp_CH)

    for i, smiles in enumerate(smiles_list):
        if not ml_list[i]:
            continue

        print('%s %i %i %.3e %.3e %.3e %.3e' % (
            smiles, t_list[i], p_list[i],
            den_list[i][0], ei_list[i][0], comp_list[i][0], expa_list[i]), file=f_ml_All)

        if cp_list[i] != math.inf:
            print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], cp_list[i]), file=f_ml_cp_All)

        if p_list[i] == 1:
            print('%s %i %.3e' % (smiles, t_list[i], hvap_list[i][0]), file=f_ml_hvap_All)

        if category_list[i] == 'Ane':
            print('%s %i %i %.3e %.3e %.3e %.3e' % (
                smiles, t_list[i], p_list[i],
                den_list[i][0], ei_list[i][0], comp_list[i][0], expa_list[i]), file=f_ml_Ane)

            if cp_list[i] != math.inf:
                print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], cp_list[i]), file=f_ml_cp_Ane)

            if p_list[i] == 1:
                print('%s %i %.3e' % (smiles, t_list[i], hvap_list[i][0]), file=f_ml_hvap_Ane)

        if category_list[i] in ['Ane', 'CH']:
            print('%s %i %i %.3e %.3e %.3e %.3e' % (
                smiles, t_list[i], p_list[i],
                den_list[i][0], ei_list[i][0], comp_list[i][0], expa_list[i]), file=f_ml_CH)

            if cp_list[i] != math.inf:
                print('%s %i %i %.3e' % (smiles, t_list[i], p_list[i], cp_list[i]), file=f_ml_cp_CH)

            if p_list[i] == 1:
                print('%s %i %.3e' % (smiles, t_list[i], hvap_list[i][0]), file=f_ml_hvap_CH)


if __name__ == '__main__':
    get_npt_result(limit=None)
    write_data(n_mol_per_file=1000)
    write_ml_data()
