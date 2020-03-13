#!/usr/bin/env python3
# coding=utf-8

import sys

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_nist import *
from app.models_ilthermo import Property, Ion, Molecule, Spline, Data
from app.models_cv import Cv
from sqlalchemy import or_, and_
from app.selection import *
from config import ClassificationConfig

import argparse

parser = argparse.ArgumentParser(description='This is a code to generate input txt from experimental data')
parser.add_argument('-t', '--type', type=str, help='The type of the output(FromEXP, FromTXT)')
parser.add_argument('-e', '--expdb', type=str, help='The type of the experimental database(ILTHERMO, NIST)')
parser.add_argument('-s', '--smiles', type=str, help='The smiles list(exp, CH). CH means you need CH-Tvap.txt and '
                                                     'CH-Tc.txt files in current folder')
parser.add_argument('--temperature', type=int, help='Only consider the experimental data contain specific temperature',
                    default=None)
parser.add_argument('--tc', type=str, help='The txt file for smiles and Tc', default=None)
parser.add_argument('--tvap', type=str, help='The txt file for smiles and Tvap', default=None)
parser.add_argument('--selection', action='store_true', help='Selection molecules')
args = parser.parse_args()
app = create_app('npt')
app.app_context().push()


def has_exp_data(_Data, _mol, _property, temperature=None, pressure=1, least_data_points=5):
    datas = _Data.query.filter(_Data.molecule == _mol).filter(_Data.property == _property).filter(
        or_(_Data.p == None, and_(_Data.p > 98 * pressure, _Data.p < 102 * pressure))).order_by(_Data.t)
    if datas.count() >= least_data_points:
        if temperature is None or datas[0].t < temperature < datas[-1].t:
            return True
    return False


# get input of ionic liquids with experimental data, used in force field validation
if args.type == 'FromEXP' and args.expdb == 'ILTHERMO':
    def write_txt(_property, txt):
        file = open(txt, 'w')
        molecules = Molecule.query
        if _property == prop_cp:
            cation_list = []
            anion_list = []
        for mol in molecules:
            if mol.cation.force_field_support != 1 or mol.anion.force_field_support != 1:
                continue
            if has_exp_data(Data, mol, _property, temperature=args.temperature):
                file.write('%i.%i %s.%s 1.1\n' % (mol.cation_id, mol.anion_id, mol.cation.smiles, mol.anion.smiles))
                if _property == prop_cp:
                    if mol.cation not in cation_list:
                        cation_list.append(mol.cation)
                    if mol.anion not in anion_list:
                        anion_list.append(mol.anion)
        file.close()
        if _property == prop_cp:
            file = open('cp_qm.txt', 'w')
            for cation in cation_list:
                if Cv.query.filter(Cv.smiles == cation.smiles).first() is None:
                    file.write('%i %s 1\n' % (cation.id, cation.smiles))
            for anion in anion_list:
                if Cv.query.filter(Cv.smiles == anion.smiles).first() is None:
                    file.write('%i %s 1\n' % (anion.id, anion.smiles))
            file.close()

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
elif args.type == 'FromEXP' and args.expdb == 'NIST' and args.smiles == 'exp':
    def write_txt(property, txt):
        f = open(txt, 'w')
        mols = NistMolecule.query.filter(NistMolecule.n_heavy < 21)
        count = 1
        for mol in mols:
            sys.stdout.write('\rpresent molecule = %i / %i\n' % (count, mols.count()))
            count += 1
            if mol.smiles is None:
                continue
            if args.selection and not ClassificationConfig.classifyer.classify(mol.smiles):
                continue
            if has_exp_data(NistData, mol, property, temperature=args.temperature):
                T_list = mol.get_sim_t_list()
                if T_list is None:
                    continue
                P_list = [1, 50, 100, 250, 500, 750, 1000]
                f.write('%s %s 1 %s %s\n' % (
                    mol.content_id, mol.smiles, ','.join(list(map(str, T_list))), ','.join(list(map(str, P_list)))))
                f.flush()
    #prop_vis = NistProperty.query.filter(NistProperty.name == 'viscosity-lg').first()
    #write_txt(prop_vis, 'vis.txt')
    prop_den = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
    write_txt(prop_den, 'den.txt')
# get input of organic liquids, you need to prepare name-Tvap.txt and name-Tc.txt files.
# The Tvap and Tc can be get from machine learning, using AIMS_ML
elif args.type == 'FromTXT' and args.tc is not None and args.tvap is not None:
    import pandas as pd
    tc_f = '%s' % args.tc
    info = pd.read_csv(tc_f, sep='\s+', header=0)
    smiles_list = info['SMILES'].to_list()
    Tc_list = info['Result'].to_list()
    # print(smiles_list, Tc_list)
    tvap_f = '%s' % args.tvap
    info = pd.read_csv(tvap_f, sep='\s+', header=0)
    if not smiles_list == info['SMILES'].to_list():
        raise Exception('SMILES in %s and %s must be the same' % (tc_f, tvap_f))
    Tvap_list = info['Result'].to_list()
    smiles_list = [get_canonical_smiles(smiles) for smiles in smiles_list]

    f = open('input.txt', 'w')
    f_qm = open('qm.txt', 'w')
    for i, smiles in enumerate(smiles_list):
        mol = NistMolecule.query.filter(Molecule.smiles == smiles).first()
        if mol is None:
            T_list = get_sim_t_list(Tvap=Tvap_list[i], Tm=None, Tc=Tc_list[i])
        else:
            T_list = get_sim_t_list(Tvap=mol.Tvap or Tvap_list[i], Tm=mol.Tfus, Tc=mol.Tc or Tc_list[i])
        P_list = [1, 50, 100, 250, 500, 750, 1000]
        rs = random_string()
        f.write('%s\t%s\t1\t%s\t%s\n' % (rs, smiles, ','.join(list(map(str, T_list))), ','.join(list(map(str, P_list)))))
        if Cv.query.filter(Cv.smiles == smiles).first() is None:
            f_qm.write('%s\t%s\t\n' % (rs, smiles))
# get input of dft calculation
elif args.type == 'FromSIM':
    f_qm = open('qm.txt', 'w')
    tasks = Task.query.filter(Task.procedure == 'npt').filter(Task.status != -1)
    for task in tasks:
        smiles_list = json.loads(task.smiles_list)
        for smiles in smiles_list:
            if Cv.query.filter(Cv.smiles == smiles).first() is None:
                rs = random_string()
                f_qm.write('%s\t%s\t\n' % (rs, smiles))
