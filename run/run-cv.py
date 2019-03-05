#!/usr/bin/env python3
# coding=utf-8

import sys
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

# QUEUE = ('cpu', 8, 0, 8)
# GAUSS_BIN = '/share/apps/g16/g16'
QUEUE = ('fast', 6, 0, 6)
GAUSS_BIN = '/share/apps/g09/g09'

n_conformer = 1

slurm = Slurm(*QUEUE)
gauss = GaussCv(gauss_bin=GAUSS_BIN, jobmanager=slurm)

app = create_app('npt')
app.app_context().push()


class Mol():
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
        name = words[2]  # str, could be 'None'
        smiles = words[3]

        mol = Mol(name, smiles)
        mol.T_list = [100, 200, 300, 400, 500, 600, 700]

        if mol not in mols:
            mols.append(mol)

    return mols


if __name__ == '__main__':
    CWD = os.getcwd()
    cmd = sys.argv[1]

    if cmd == 'print':
        n = 0
        mols = get_mols()
        for mol in mols:
            if Cv.query.filter(Cv.smiles == mol.smiles).count() > 0:
                continue
            n += 1
            print(n, mol.formula, mol.name + '_' + random_string(4), mol.smiles, sep='\t')

    if cmd == 'cv':
        mols = get_mols()
        cd_or_create_and_cd(os.path.join(Config.WORK_DIR, 'Cv'))
        for mol in mols:
            print(mol)
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
        mols = get_mols()
        from mstools.analyzer.fitting import polyfit

        for mol in mols:
            try:
                os.chdir(os.path.join(Config.WORK_DIR, 'Cv', mol.name))
            except:
                print(mol, 'Error: Dir not exist')
                continue

            gauss.logs = ['conf-%i.log' % i for i in range(n_conformer)]
            try:
                result = gauss.analyze()
            except Exception as e:
                print(mol, str(e))
                continue

            T_list = []
            Cv_list = []
            for T, val_stderr in result['Cv-corrected'].items():
                T_list.append(T)
                Cv_list.append(val_stderr[0])
            if T_list != [100, 200, 300, 400, 500, 600, 700]:
                print(mol, 'Some temperatures are failed')
            else:
                coef, score = polyfit(T_list, Cv_list, 4)
                print(mol, *coef, score)

            os.chdir(CWD)
