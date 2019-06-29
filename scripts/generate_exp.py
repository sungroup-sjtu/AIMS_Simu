#!/usr/bin/env python3
# coding=utf-8

import sys
sys.path.append('..')
from app import create_app
from app.models_nist import *
from mstools.smiles.smiles import *
from app.function import get_t_min_max

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate experimental data.')
parser.add_argument('-t', '--selectiontype', type=str, help='The type of the output(CH, LinearBranchAlkane)')

opt = parser.parse_args()
app = create_app('npt')
app.app_context().push()


prop_density = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
prop_hvap = NistProperty.query.filter(NistProperty.name == 'hvap-lg').first()
prop_cp = NistProperty.query.filter(NistProperty.name == 'cp-lg').first()
prop_st = NistProperty.query.filter(NistProperty.name == 'st-lg').first()
properties = [prop_density, prop_cp, prop_hvap, prop_st]
for property in properties:
    f = open('%s.txt' % property.name, 'w')
    f.write('SMILES T P Result Uncertainty\n')
    has_datas = NistHasData.query.filter(NistHasData.property == property).filter(NistHasData.has_exp == True)
    for has_data in has_datas:
        molecule = has_data.molecule
        t_min, t_max = get_t_min_max(Tvap=molecule.Tvap, Tfus=molecule.Tfus, Tc=molecule.Tc)
        if t_min is None or t_max is None:
            continue
        smiles = molecule.smiles
        if not 5 < get_heavy_atom_numbers(smiles) < 16:
            continue
        if not selection(smiles, type=opt.selectiontype):
            continue
        datas = NistData.query.filter(NistData.property == property).filter(NistData.molecule == molecule).\
            filter(NistData.t < t_max).filter(NistData.t > t_min)
        for data in datas:
            if data.t > 600:
                continue
            if property.name == 'density-lg':
                value = data.value / 1000
                uncertainty = data.uncertainty / 1000
            else:
                value = data.value
                uncertainty = data.uncertainty
            f.write('%s %.2f %.2f %.5e %.5e\n' % (smiles, data.t, data.p or 1, value, uncertainty))
