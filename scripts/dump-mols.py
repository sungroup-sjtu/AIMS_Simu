#!/usr/bin/env python3

'''
Dump molecules from NPT simulations
The IUPAC names are obtained from NIST database, or None if not exist.
'''

import sys
import csv

sys.path.append('..')

from app import create_app, db
from app.models import Task
from app.models_nist import NistMolecule

smiles_list = set()
npt = create_app('npt')
with npt.app_context():
    tasks = Task.query.filter(Task.remark == None)
    for task in tasks:
        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue
        smiles_list.add(task.get_smiles_list()[0])

slab = create_app('nvt-slab')
with slab.app_context():
    tasks = Task.query.filter(Task.remark == None)
    for task in tasks:
        post_result = task.get_post_result()
        if post_result is None:
            continue
        smiles_list.add(task.get_smiles_list()[0])

with npt.app_context():
    nist_smiles_iupac = dict(db.session.query(NistMolecule.smiles, NistMolecule.name).all())
    nist_smiles_cas = dict(db.session.query(NistMolecule.smiles, NistMolecule.cas).all())

fout = open('mols.csv', 'w')
fcsv = csv.writer(fout, quoting=csv.QUOTE_ALL)

fcsv.writerow(['smiles', 'iupac', 'cas', 'source', 'category'])
for smiles in smiles_list:
    iupac = nist_smiles_iupac.get(smiles, None)
    cas = nist_smiles_cas.get(smiles, None)
    source = None if iupac is None else 'NIST'

    fcsv.writerow([smiles, iupac or 'unknown', cas or 'unknown', source or 'unknown', sys.argv[1]])