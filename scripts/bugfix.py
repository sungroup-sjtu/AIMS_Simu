#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *
from app.models_cv import *

import argparse
parser = argparse.ArgumentParser(description='This is a code to fix the bug in the database')
parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
parser.add_argument('--smiles', type=bool, help='Check the smiles in database, make it canonical', default=True)
opt = parser.parse_args()

procedure = opt.procedure
if procedure == 'cv':
    procedure = 'npt'
app = create_app(procedure)
app.app_context().push()
if opt.procedure == 'cv' and opt.smiles:
    cvs = Cv.query
    count = 1
    for cv in cvs:
        sys.stdout.write('\rpresent task = %i' % count)
        count += 1
        if cv.smiles != get_canonical_smiles(cv.smiles):
            cv.smiles = get_canonical_smiles(cv.smiles)
            print('\ncv %i, %s: is not canonical, change it.' % (cv.id, cv.smiles))
    db.session.commit()
elif opt.smiles:
    tasks = Task.query.filter(Task.procedure == procedure)
    print('Total task = %i\n' % (tasks.count()))
    count = 1
    for task in tasks:
        sys.stdout.write('\rpresent task = %i' % count)
        count += 1
        smiles_list = json.loads(task.smiles_list)
        canonical_smiles_list = [get_canonical_smiles(smiles) for smiles in smiles_list]
        if task.smiles_list != json.dumps(canonical_smiles_list):
            task.smiles_list = json.dumps(canonical_smiles_list)
            print('\ntask %i, %s: is not canonical, change it.' % (task.id, task.smiles_list))
    db.session.commit()
