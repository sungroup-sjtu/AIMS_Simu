#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app.models import *
from app.models_cv import Cv
from config import Config

sys.path.append(Config.MS_TOOLS_DIR)
from mstools.simulation.gauss import Cv as GaussCv
from mstools.jobmanager import Torque
from mstools.utils import cd_or_create_and_cd

from collections import OrderedDict

QUEUE_DICT = OrderedDict([('cpu', 4)])
GAUSS_BIN = '/share/apps/g16/g16'
#QUEUE_DICT = OrderedDict([('fast', 6)])
#QUEUE_DICT = OrderedDict([('batch', 4)])
#GAUSS_BIN = '/share/apps/g09-b01/g09/g09'

n_conformer = 1

torque = Torque(queue_dict=QUEUE_DICT)
gauss_cv = GaussCv(gauss_bin=GAUSS_BIN, jobmanager=torque)


class Mol():
    def __init__(self, name, smiles):
        self.name = name
        self.smiles = smiles
        self.T_list = []

    def __repr__(self):
        return '<Mol: %s %s>' % (self.name, self.smiles)

    def __eq__(self, mol2):
        return self.name == mol2.name


def get_mols():
    mols = []
    with open(sys.argv[2]) as f:
        lines = f.read().splitlines()
    for line in lines:
        if line.strip() == '' or line.startswith('#'):
            continue
        words = line.strip().split()
        cas = words[1]
        smiles = words[2]
        mol = Mol(cas, smiles)
        mol.T_list = [100, 200, 300, 400, 500, 600, 700]
        mols.append(mol)
    return mols


if __name__ == '__main__':
    CWD = os.getcwd()
    cmd = sys.argv[1]

    if cmd == 'print':
        n = 0
        tasks = Task.query
        for task in tasks:
            cas = task.name[:-5]
            smiles = json.loads(task.smiles_list)[0]
            if Cv.query.filter(Cv.smiles == smiles).count() > 0:
                continue
            n += 1
            print(n, cas, smiles, sep='\t')

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
            gauss_cv.set_system([mol.smiles], n_mol_list=[1])
            gauss_cv.prepare(n_conformer=n_conformer, T_list=mol.T_list, jobname=mol.name)
            gauss_cv.run()
            os.chdir('..')
        os.chdir(CWD)

    if cmd == 'get-cv':
        mols = get_mols()
        from mstools.analyzer.fitting import polyfit

        cd_or_create_and_cd(os.path.join(Config.WORK_DIR, 'Cv'))
        for mol in mols:
            cd_or_create_and_cd(mol.name)
            gauss_cv.logs = ['conf-%i.log' % i for i in range(n_conformer)]
            try:
                result = gauss_cv.analyze()
            except Exception as e:
                print(mol, str(e))
            else:
                T_list = []
                Cv_list = []
                for T, val_stderr in result['Cv-corrected'].items():
                    T_list.append(T)
                    Cv_list.append(val_stderr[0])
                if T_list != [100, 200, 300, 400, 500, 600, 700]:
                    print(mol, 'Some temperatures are failed')
                else:
                    coef, score = polyfit(T_list, Cv_list, 4)
                    print(mol, *coef)

            os.chdir('..')
        os.chdir(CWD)
