#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app.models import *

import random
import pybel
from collections import OrderedDict


def get_N(n):
    if n <= 20:
        return n
    elif n >= 60:
        return n // 2
    else:
        return 20 + (n - 20) // 4


class MolLine():
    def __init__(self, m, line):
        self.mol = m
        self.line = line
        self.smiles = m.write('can').strip()
        self.matched = False


with open(sys.argv[1]) as f:
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

print('n_heav n_mols n_need n_matc')
for k, v in n_heavy_mols.items():
    N = get_N(len(v))

    n_matched = 0
    for molLine in v:
        if Task.query.filter(Task.smiles_list == json.dumps([molLine.smiles])).first() != None:
            n_matched += 1

    print('%6i %6i %6i %6i' % (k, len(v), N, n_matched))
