#!/usr/bin/env python3
# coding=utf-8

import sys, os
import pybel

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_cv import Cv
from config import Config

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.gauss import Cv as GaussCv
from mstools.jobmanager import Slurm
from mstools.utils import cd_or_create_and_cd
from mstools.analyzer.fitting import polyfit

# QUEUE = ('cpu', 8, 0, 8)
# GAUSS_BIN = '/share/apps/g16/g16'
QUEUE = ('fast', 6, 0, 6)
GAUSS_BIN = '/share/apps/g09/g09'

n_conformer = 1

slurm = Slurm(*QUEUE)
gauss = GaussCv(gauss_bin=GAUSS_BIN, jobmanager=slurm)

app = create_app('npt')
app.app_context().push()


class Mol:
    def __init__(self, name, smiles, formula=None):
        self.name = name
        self.smiles = smiles
        if formula is None:
            formula = pybel.readstring('smi', smiles).formula
        self.formula = formula
        self.T_list = []

    def __repr__(self):
        return '<Mol: %s %s>' % (self.name, self.smiles)

    def __eq__(self, m2):
        return self.smiles == m2.smiles


def get_mols():
    mols = []
    with open(sys.argv[2]) as f:
        lines = f.read().splitlines()
    for line in lines:
        if line.strip() == '' or line.startswith('#'):
            continue

        words = line.strip().split()
        name = words[0]  # str, could be 'None'
        smiles = words[1]

        mol = Mol(name, smiles)
        mol.T_list = [100, 200, 300, 400, 500, 600, 700]

        if mol not in mols:
            mols.append(mol)

    return mols


if __name__ == '__main__':
    CWD = os.getcwd()
    cmd = sys.argv[1]

    if cmd == 'prepare':
        fout = open('_cv_prepared.txt', 'w')
        n = 0
        mols = get_mols()
        for mol in mols:
            # if Cv.query.filter(Cv.smiles == mol.smiles).count() > 0:
                # continue
            n += 1
            print(mol.name + '_' + random_string(4), mol.smiles, n, mol.formula, sep='\t', file=fout)
        fout.close()

    if cmd == 'cv':
        mols = get_mols()
        cd_or_create_and_cd(os.path.join(Config.WORK_DIR, 'Cv'))
        for mol in mols:
            print(mol)
            if Cv.query.filter(Cv.smiles == mol.smiles).count() != 0:
                continue
            try:
                shutil.rmtree(mol.name)
            except:
                pass
            cd_or_create_and_cd(mol.name)
            gauss.set_system([mol.smiles], n_mol_list=[1])
            gauss.prepare(n_conformer=n_conformer, T_list=mol.T_list, jobname=mol.name)
            gauss.run()
            os.chdir('..')
        os.chdir(CWD)

    if cmd == 'get-cv':
        fout = open('_cv.log', 'w')
        fout2 = open('_enthalpy.log', 'w')
        mols = get_mols()
        print('Total mol = %i' % (len(mols)))
        count = 1
        for mol in mols:
            sys.stdout.write('\rpresent mol = %i' % (count))
            count += 1
            if Cv.query.filter(Cv.smiles == mol.smiles).count() != 0:
                continue
            try:
                os.chdir(os.path.join(Config.WORK_DIR, 'Cv', mol.name))
            except:
                print(mol, 'Error: Dir not exist', file=fout)
                continue

            gauss.logs = ['conf-%i.log' % i for i in range(n_conformer)]
            try:
                result = gauss.analyze()
            except Exception as e:
                print(mol, str(e), file=fout)
                print(mol, str(e), file=fout2)
                continue

            T_list = []
            Cv_list = []
            enthalpy_list = []
            for T, val_stderr in result['Cv-corrected'].items():
                T_list.append(T)
                Cv_list.append(val_stderr[0])
                enthalpy_list.append(result['enthalpy'][T][0])
            if T_list != [100, 200, 300, 400, 500, 600, 700]:
                print(mol, 'Some temperatures are failed', file=fout)
                print(mol, 'Some temperatures are failed', file=fout2)
            else:
                coef, score = polyfit(T_list, Cv_list, 4)
                print(mol, *coef, score, file=fout)
                coef, score = polyfit(T_list, enthalpy_list, 2)
                print(mol, *coef, score, file=fout2)

            os.chdir(CWD)
        fout.close()

    if cmd == 'save-db':
        Cv.load_from_log('_cv.log', '_enthalpy.log')
