#!/usr/bin/env python3

import sys

sys.path.append('..')
sys.path.append('../../ms-tools')

from app.models import *
from app.models_nist import *

print('#SMILES T(K) P(bar) density(g/mL) e_inter(kJ/mol) cp(J/mol.K)')

for nist in NistMolecule.query.filter(NistMolecule.remark == 'alkane').filter(
        NistMolecule.n_heavy < 20).order_by(NistMolecule.n_heavy):
    task = Task.query.filter(Task.smiles_list == json.dumps([nist.smiles])).first()
    if task != None:
        continue

    t_list = [298]
    if nist.tt != None:
        t_list.append(nist.tt + 25)
    if nist.tb != None:
        t_list.append(nist.tb)
    if nist.tc != None:
        t_list.append(nist.tc * 0.8)

    spline_pvap = NistSpline.query.join(NistProperty).filter(NistSpline.molecule == nist).filter(
            NistProperty.name == 'pvap-lg').first()
    if spline_pvap == None:
        continue

    spline_dens = NistSpline.query.join(NistProperty).filter(NistSpline.molecule == nist).filter(
            NistProperty.name == 'density-lg').first()
    spline_hvap = NistSpline.query.join(NistProperty).filter(NistSpline.molecule == nist).filter(
            NistProperty.name == 'hvap-lg').first()
    spline_cp = NistSpline.query.join(NistProperty).filter(NistSpline.molecule == nist).filter(
            NistProperty.name == 'cp-lg').first()

    for t in t_list:
        p = spline_pvap.get_data(t)[0]  # kPa
        if p == None:
            continue

        t = int(round(t))  # K
        p = int(round(p / 100))  # bar

        if spline_dens != None:
            density = spline_dens.get_data(t)[0] or 0
        else:
            density = 0

        if spline_hvap != None:
            hvap = spline_hvap.get_data(t)[0] or 0
        else:
            hvap = 0
        if p <= 1 and hvap != 0:
            ei = 8.314 * t / 1000 - hvap  # kJ/mol
        else:
            ei = 0

        if spline_cp != None:
            cp = spline_cp.get_data(t)[0] or 0
        else:
            cp = 0

        print('%s %i %i %.3e %.3e %.3e' % (nist.smiles, t, p, density / 1000, ei, cp))
