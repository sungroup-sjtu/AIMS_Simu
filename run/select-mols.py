#!/usr/bin/env python3
# coding=utf-8

'''
./select-mols.py procedure mols.txt fraction
'''

import sys

sys.path.append('..')
from app import create_app
from app.models import *

import random
import pybel
from collections import OrderedDict


def get_N(N, n_heavy, frac=0.5):
    if n_heavy <= 5:
        return N

    Nmin = 15
    Nmax = 60

    if N <= Nmin:
        return N
    elif N >= Nmax:
        return int(N * frac)
    else:
        return Nmin + int((N - Nmin) / (Nmax - Nmin) * (Nmax * frac - Nmin))


class MolLine():
    def __init__(self, m, line):
        self.mol = m
        self.line = line
        self.smiles = m.write('can').strip()
        self.matched = False


app = create_app(sys.argv[1])
app.app_context().push()

with open(sys.argv[2]) as f:
    lines = f.read().splitlines()

smiles_list = []
n_heavy_mols = OrderedDict([(i, []) for i in range(1, 20)])

for line in lines:
    if line == '' or line.startswith('#'):
        continue
    words = line.strip().split()
    smiles = words[3]
    if smiles in smiles_list:
        continue

    smiles_list.append(smiles)
    m = pybel.readstring('smi', smiles)
    n_heavy = m.OBMol.NumHvyAtoms()
    n_heavy_mols[n_heavy].append(MolLine(m, line))

fout = open('out.txt', 'w')
print('n_heav n_mols n_need n_matc n_rema')
for k, molLine_list in n_heavy_mols.items():
    N = get_N(len(molLine_list), k, float(sys.argv[3]))

    n_matched = 0
    for molLine in molLine_list:
        if Task.query.filter(Task.smiles_list == json.dumps([molLine.smiles])).first() is not None:
            fout.write(molLine.line + '\n')
            molLine.matched = True
            n_matched += 1

    n_remain = N - n_matched
    print('%6i %6i %6i %6i %6i' % (k, len(molLine_list), N, n_matched, n_remain))
    n_remain = max(n_remain, 0)

    l = [molLine for molLine in molLine_list if not molLine.matched]
    random.shuffle(l)
    for molLine in l[:n_remain]:
        fout.write(molLine.line + '\n')
fout.close()
