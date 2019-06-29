#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *
from app.models_cv import Cv
from app.selection import *

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate input txt from experimental data')
parser.add_argument('-t', '--type', type=str, help='The type of the output(SIM, EXP)')
parser.add_argument('--database', type=str, help='The type of the experimental database(ILTHERMO, NIST)')
parser.add_argument('--procedure', type=str, help='The procedure of simulation(npt)')
parser.add_argument('--selection', type=bool, help='Use the selection function', default=True)
parser.add_argument('--training', type=str, help='The txt file for training sets smiles', default=None)

opt = parser.parse_args()
app = create_app('npt')
app.app_context().push()

if opt.type == 'EXP' and opt.database == 'NIST':
    f1 = open('Tvap.txt', 'w')
    f1.write('SMILES Tvap\n')
    f2 = open('Tm.txt', 'w')
    f2.write('SMILES Tm\n')
    f3 = open('Tc.txt', 'w')
    f3.write('SMILES Tc\n')
    molecules = NistMolecule.query.filter(NistMolecule.n_heavy > 5)
    for mol in molecules:
        print(mol.smiles)
        if not selection(mol.smiles, type='CH'):
            continue
        if mol.tb is not None:
            f1.write('%s %f\n' % (mol.smiles, mol.tb))
        if mol.tt is not None:
            f2.write('%s %f\n' % (mol.smiles, mol.tt))
        if mol.tc is not None:
            f3.write('%s %f\n' % (mol.smiles, mol.tc))
elif opt.type == 'SIM' and opt.procedure == 'npt':
    if opt.training is not None:
        f_train1 = open('result-ML-density-train.txt', 'w')
        f_train2 = open('result-ML-einter-train.txt', 'w')
        f_train3 = open('result-ML-compress-train.txt', 'w')
        f_train4 = open('result-ML-hvap-train.txt', 'w')
        f_train5 = open('result-ML-cp-train.txt', 'w')
        f_train1.write('SMILES T P density\n')
        f_train2.write('SMILES T P einter\n')
        f_train3.write('SMILES T P compress\n')
        f_train4.write('SMILES T P hvap\n')
        f_train5.write('SMILES T P cp\n')
        info = pd.read_csv(opt.training, sep='\s+', header=0)
        training_smiles_list = info['SMILES'].to_list()
    else:
        training_smiles_list = []
    f1 = open('result-ML-density.txt', 'w')
    f2 = open('result-ML-einter.txt', 'w')
    f3 = open('result-ML-compress.txt', 'w')
    f4 = open('result-ML-hvap.txt', 'w')
    f5 = open('result-ML-cp.txt', 'w')
    f1.write('SMILES T P density\n')
    f2.write('SMILES T P einter\n')
    f3.write('SMILES T P compress\n')
    f4.write('SMILES T P hvap\n')
    f5.write('SMILES T P cp\n')
    smiles_list = []
    t_list = []
    p_list = []
    den_list = []
    tasks = Task.query.filter(Task.procedure == opt.procedure)
    for task in tasks:
        if task.status in [Compute.Status.FAILED, Compute.Status.STARTED]:
            continue
        if not task_selection(task, select=opt.selection):
            continue
        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue
        smiles = task.get_smiles_list()[0]
        py_mol = task.get_mol_list()[0]
        if not 5 < get_heavy_atom_numbers(smiles) < 16:
            continue
        f = Formula(py_mol.formula)
        n_CON = f.atomdict.get('C', 0) + f.atomdict.get('O', 0) + f.atomdict.get('N', 0)
        cv = Cv.query.filter(Cv.smiles == smiles).first()
        for t in task.get_t_list():
            for p in task.get_p_list():
                post_data = task.get_post_data(t, p)
                density = post_data.get('density')
                einter = post_data.get('einter')
                expansion = post_data.get('expansion')
                if einter is not None:
                    hvap = (1 - n_CON / 15) * 8.314 * t / 1000 - einter
                else:
                    hvap = None
                compress = post_data.get('compressibility')
                if density is not None:
                    f1.write('%s %i %i %.3e\n' % (smiles, t, p, density))
                if einter is not None:
                    f2.write('%s %i %i %.3e\n' % (smiles, t, p, einter))
                if compress is not None:
                    f3.write('%s %i %i %.3e\n' % (smiles, t, p, compress))
                if hvap is not None:
                    f4.write('%s %i %i %.3e\n' % (smiles, t, p, hvap))

                if smiles in training_smiles_list:
                    if density is not None:
                        f_train1.write('%s %i %i %.3e\n' % (smiles, t, p, density))
                    if einter is not None:
                        f_train2.write('%s %i %i %.3e\n' % (smiles, t, p, einter))
                    if compress is not None:
                        f_train3.write('%s %i %i %.3e\n' % (smiles, t, p, compress))
                    if hvap is not None:
                        f_train4.write('%s %i %i %.3e\n' % (smiles, t, p, hvap))

                cp_inter = post_data.get('cp_inter')
                cp_pv = post_data.get('cp_pv')
                if None not in [cv, cp_inter, cp_pv]:
                    cv_intra = cv.get_post_cv(t)
                    cp = cv_intra + cp_inter + cp_pv
                    f5.write('%s %i %i %.3e\n' % (smiles, t, p, cp))
                    if smiles in training_smiles_list:
                        f_train5.write('%s %i %i %.3e\n' % (smiles, t, p, cp))
elif opt.type == 'SIM' and opt.procedure == 'nvt-slab':
    if opt.training is not None:
        f_train1 = open('result-ML-npt-train.txt', 'w')
        f_train2 = open('result-ML-cp-train.txt', 'w')
        f_train1.write('SMILES T P density einter compress hvap')
        f_train1.write('SMILES T P cp')
        info = pd.read_csv(opt.training, sep='\s+', header=0)
        training_smiles_list = info['SMILES'].to_list()
    else:
        training_smiles_list = []

