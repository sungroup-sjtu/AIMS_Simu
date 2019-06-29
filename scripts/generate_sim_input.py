#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app import create_app
from app.models_nist import *
from app.models_ilthermo import Property, Ion, Molecule, Spline, Data
from sqlalchemy import or_, and_
from app.selection import *

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate input txt from experimental data')
parser.add_argument('-t', '--type', type=str, help='The type of the output(SIMU, ML)')
parser.add_argument('-e', '--exp', type=str, help='The type of the experimental database(ILTHERMO, NIST)')
parser.add_argument('-s', '--smiles', type=str,
                    help='The smiles list(exp, CH). CH means you need CH-Tvap.txt and CH-Tc.txt files in current folder')
parser.add_argument('--temperature', type=int, help='Only consider the experimental data contain specific temperat',
                    default=None)
opt = parser.parse_args()
app = create_app('npt')
app.app_context().push()
# get input of ionic liquids with experimental data, used in force field validation
if opt.type == 'SIMU' and opt.exp == 'ILTHERMO' and opt.smiles == 'exp':
    def has_exp_data(mol, property, temperature=None, pressure=1):
        datas = Data.query.filter(Data.molecule == mol).filter(Data.property == property).filter(
            or_(Data.p == None, and_(Data.p > 98 * pressure, Data.p < 102 * pressure))).order_by(Data.t)
        if datas.count() > 0:
            if temperature is None or datas[0].t < temperature < datas[-1].t:
                return True
        return False


    def write_txt(property, txt):
        f = open(txt, 'w')
        mols = Molecule.query
        if property == prop_cp:
            cation_list = []
            anion_list = []
        for mol in mols:
            if mol.cation.force_field_support != 1 or mol.anion.force_field_support != 1:
                continue
            if has_exp_data(mol, property, temperature=opt.temperature):
                f.write('%i.%i %s.%s 1.1\n' % (mol.cation_id, mol.anion_id, mol.cation.smiles, mol.anion.smiles))
                if property == prop_cp:
                    if mol.cation not in cation_list:
                        cation_list.append(mol.cation)
                    if mol.anion not in anion_list:
                        anion_list.append(mol.anion)
        f.close()
        if property == prop_cp:
            f = open('cp_qm.txt', 'w')
            for cation in cation_list:
                f.write('%i %s 1\n' % (cation.id, cation.smiles))
            for anion in anion_list:
                f.write('%i %s 1\n' % (anion.id, anion.smiles))
            f.close()

    prop_density = Property.query.filter(Property.name == 'Density').first()
    prop_cp = Property.query.filter(Property.name == 'Heat capacity at constant pressure').first()
    prop_econ = Property.query.filter(Property.name == 'Electrical conductivity').first()
    prop_vis = Property.query.filter(Property.name == 'Viscosity').first()
    prop_diff = Property.query.filter(Property.name == 'Self-diffusion coefficient').first()
    write_txt(prop_density, 'den.txt')
    write_txt(prop_cp, 'cp.txt')
    write_txt(prop_econ, 'econ.txt')
    write_txt(prop_vis, 'vis.txt')
    write_txt(prop_diff, 'diff.txt')
# get input of organic liquids with experimental data, used in force field validation
elif opt.type == 'SIMU' and opt.exp == 'NIST' and opt.smiles == 'exp':
    def has_exp_data(mol, property, temperature=None, pressure=1):
        datas = NistData.query.filter(NistData.molecule == mol).filter(NistData.property == property).filter(
            or_(NistData.p == None, and_(NistData.p > 98 * pressure, NistData.p < 102 * pressure))).order_by(NistData.t)
        if datas.count() > 0:
            if temperature is None or datas[0].t < temperature < datas[-1].t:
                return True
        return False


    def write_txt(property, txt):
        f = open(txt, 'w')
        mols = NistMolecule.query.filter(NistMolecule.n_heavy < 21)
        count = 1
        for mol in mols:
            sys.stdout.write('\rpresent molecule = %i / %i\n' % (count, mols.count()))
            count += 1
            if mol.smiles is None:
                continue
            if not IsCH(mol.smiles):
                continue
            if has_exp_data(mol, property, temperature=opt.temperature):
                T_list = mol.get_sim_t_list()
                if T_list is None:
                    continue
                P_list = [1, 50, 100, 250, 500, 750, 1000]
                f.write('%s %s 1 %s %s\n' % (
                    mol.contend_id, mol.smiles, ','.join(list(map(str, T_list))), ','.join(list(map(str, P_list)))))
                f.flush()


    prop_vis = NistProperty.query.filter(NistProperty.name == 'viscosity-lg').first()
    write_txt(prop_vis, 'vis.txt')
# get input of organic liquids, you need to prepare name-Tvap.txt and name-Tc.txt files.
# The Tvap and Tc can be get from machine learning, using AIMS_ML
elif opt.type == 'SIMU' and opt.exp == 'NIST' and opt.smiles != 'exp':
    import pandas as pd
    tc_f = '%s-Tc.txt' % opt.smiles
    info = pd.read_csv(tc_f, sep='\s+', header=0)
    smiles_list = info['SMILES'].to_list()
    Tc_list = info['Result'].to_list()
    # print(smiles_list, Tc_list)
    tvap_f = '%s-Tvap.txt' % opt.smiles
    info = pd.read_csv(tvap_f, sep='\s+', header=0)
    if not smiles_list == info['SMILES'].to_list():
        raise Exception('SMILES in %s and %s must be the same' % (tc_f, tvap_f))
    Tvap_list = info['Result'].to_list()

    f = open('input.txt', 'w')
    for i, smiles in enumerate(smiles_list):
        mol = Molecule.query.filter(Molecule.smiles == smiles).first()
        if mol is None:
            T_list = get_sim_t_list(Tvap=Tvap_list[i], Tm=None, Tc=Tc_list[i])
        else:
            T_list = get_sim_t_list(Tvap=mol.Tvap or Tvap_list[i], Tm=mol.Tfus, Tc=mol.Tc or Tc_list[i])
        P_list = [1, 50, 100, 250, 500, 750, 1000]
        f.write('%s\t%s\t1\t%s\t%s\n' % (random_string(), smiles, ','.join(list(map(str, T_list))), ','.join(list(map(str, P_list)))))
