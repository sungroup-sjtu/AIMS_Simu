#!/usr/bin/env python3

'''
Dump molecules from NPT simulations
The IUPAC names are obtained from NIST database, or YAWS if not exist.
'''

import sys

sys.path.append('..')

from app import create_app, db
from app.models import Task
from app.models_nist import NistMolecule
from app.models_yaws import YawsMolecule

smiles_list = set()
npt = create_app('npt')
with npt.app_context():
    for task in db.session.query(Task):
        post_result = task.get_post_result()
        if post_result is None or post_result['density-poly4'][-1] < 0.999:
            continue
        smiles_list.add(task.get_smiles_list()[0])

slab = create_app('nvt-slab')
with slab.app_context():
    for task in db.session.query(Task):
        post_result = task.get_post_result()
        if post_result is None:
            continue
        smiles_list.add(task.get_smiles_list()[0])

with npt.app_context():
    nist_smiles_iupac = dict(db.session.query(NistMolecule.smiles, NistMolecule.name).all())
    nist_smiles_cas = dict(db.session.query(NistMolecule.smiles, NistMolecule.cas).all())
    yaws_smiles_iupac = dict(db.session.query(YawsMolecule.smiles, YawsMolecule.iupac).all())
    yaws_smiles_cas = dict(db.session.query(YawsMolecule.smiles, YawsMolecule.cas).all())

print('smiles', 'iupac', 'cas', 'source', sep='\t')
for smiles in smiles_list:
    iupac = nist_smiles_iupac.get(smiles, None)
    cas = nist_smiles_cas.get(smiles, None)
    source = 'NIST'

    if iupac is None:
        iupac = yaws_smiles_iupac.get(smiles, None)
        cas = yaws_smiles_cas.get(smiles, None)
        source = 'Yaws'

        if iupac is None:
            source = None

    print(smiles, iupac or 'unknown', cas or 'unknown', source or 'unknown', sep='\t')
