#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *

import argparse
import pybel
parser = argparse.ArgumentParser(description='This is a code to dump txt file')
parser.add_argument('-i', '--id', type=int, help='The specific compute id to dump', default=0)
parser.add_argument('-p', '--procedure', type=str, help='procedure: npt(ppm), or nvt-slab')
parser.add_argument('--nist', type=str, help='using the infomation in nist database')
opt = parser.parse_args()

procedure = opt.procedure
app = create_app(procedure)
app.app_context().push()

if opt.id!=0:
    compute = Compute.query.get(opt.id)
    tasks = compute.tasks
else:
    tasks = Task.query.filter(Task.procedure==procedure).filter(Task.remark==None)

f1 = open('input-p1.txt', 'w')
f2 = open('input-p2.txt', 'w')
f3 = open('input-t1.txt', 'w')
f1.write('#No Formula CAS SMILES Tfus Tvap Tc content_id\n')
count = 1
n = 1
print('Total task = %i' % (tasks.count()))
smiles_list = [] # part
for task in tasks:
    sys.stdout.write('\rpresent task = %i' % (count))
    count += 1
    # total  smiles_com = A.B.C.D.E
    smiles_com = '.'.join(json.loads(task.smiles_list))
    f3.write('%s\n' % (smiles_com))
    # part   smiles = A
    for i, smiles in enumerate(json.loads(task.smiles_list)):
        if smiles in smiles_list:
            continue
        if opt.nist:
            mol = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if mol != None:
                f1.write('%i %s %s %s %s %s %s %s\n' % (
                    n, mol.formula, mol.cas, smiles, str(mol.tt), str(mol.tb), str(mol.tc), mol.content_id))
            else:
                py_mol = pybel.readstring('smi', smiles)
                f1.write('%i %s %s %s %s %s %s %s\n' % (
                    n, py_mol.formula, task.name.split('_')[i], smiles, str(None), str(None), str(None), str(None)))
        else:
            py_mol = pybel.readstring('smi', smiles)
            f1.write('%i %s %s %s %s %s %s %s\n' % (
                n, py_mol.formula, task.name.split('_')[i], smiles, str(None), str(None), str(None), str(None)))
        smiles_list.append(smiles)
        f2.write(smiles + '\n')
        n += 1

db.session.commit()
