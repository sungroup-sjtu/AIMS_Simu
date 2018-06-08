#!/usr/bin/env python3

import sys

sys.path.append('..')
sys.path.append('../../ms-tools')

from app.models import *
from app.models_nist import *

print('#SMILES T(K) P(bar) density(g/mL) u e_inter(kJ/mol) u cp(J/mol.K) u')

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
        p /= 100  # bar

        density, den_u, hvap, hvap_u, ei, ei_u, cp, cp_u = [None] * 8

        if spline_dens != None:
            density, den_u = spline_dens.get_data(t)

        if spline_hvap != None:
            hvap, hvap_u = spline_hvap.get_data(t)

        if spline_cp != None:
            cp, cp_u = spline_cp.get_data(t)

        if hvap != None:
            if p <= 1:
                ei = 8.314 * t / 1000 - hvap  # kJ/mol
            else:
                spline_dgl = NistSpline.query.join(NistProperty).filter(NistSpline.molecule == nist).filter(
                        NistProperty.name == 'density-gl').first()
                if spline_dgl != None:
                    dgas, _ = spline_dgl.get_data(t)
                    if dgas != None and density != None:
                        pVg_l = p * nist.weight * (1 / dgas - 1 / density) / 10  # kJ/mol
                        ei = pVg_l - hvap

        if density == None:
            density, den_u = 0, 0
        else:
            den_u = den_u / density

        if ei == None:
            ei, ei_u = 0, 0
        else:
            ei_u = hvap_u / hvap

        if cp == None:
            cp, cp_u = 0, 0
        else:
            cp_u = cp_u / cp

        t = int(round(t))  # K
        p = int(round(p))  # bar
        print('%s %i %i %.3e %.1e %.3e %.1e %.3e %.1e' % (nist.smiles, t, p,
                                                       density / 1000, den_u,
                                                       ei, ei_u,
                                                       cp, cp_u,
                                                       ))
