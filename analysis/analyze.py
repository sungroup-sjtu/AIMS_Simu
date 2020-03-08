#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to reanalyze a specific compute')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

if procedure == 'npt':
    tasks = Task.query.filter(Task.procedure == procedure)
    smiles_tasks = []
    for task in tasks:
        smiles_tasks.append(json.loads(task.smiles_list)[0])
    print(len(smiles_tasks))

    nists = NistMolecule.query#.filter(NistMolecule.tb != None)
    smiles_nists = []
    for nist in nists:
        if nist.n_heavy <= 20:
            smiles_nists.append(nist.smiles)
    print(len(smiles_nists))
    for smiles in smiles_tasks:
        if smiles not in smiles_nists:
            print(smiles)








