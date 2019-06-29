#!/usr/bin/env python3
# coding=utf-8

from app import create_app
from app.models_ilthermo import *
from sqlalchemy import or_, and_
import pybel, numpy as np, sys
from mstools.analyzer.fitting import *

procedure = 'npt'
app = create_app(procedure)
app.app_context().push()

def remove_repeat_data(t_list, v_list, u_list):
    t_list_new = list(set(t_list))
    v_list_new = []
    u_list_new = []
    for t in t_list_new:
        v_list_new.append([])
        u_list_new.append([])
    for i, t in enumerate(t_list):
        index = t_list_new.index(t)
        v_list_new[index].append(v_list[i])
        u_list_new[index].append(u_list[i])
    for i, v in enumerate(v_list_new):
        if np.std(v) / np.mean(v) > 0.1:
            print('v=',v, t_list_new[i])
        v_list_new[i] = np.mean(v)
    for i, u in enumerate(u_list_new):
        if np.std(u) / np.mean(u) > 0.1:
            print('u=',u, t_list_new[i])
        u_list_new[i] = np.mean(u)
    return t_list_new, v_list_new, u_list_new

def add_spline_info_ilthermo():
    molecule_property_paper = []
    datas = Data.query.filter(Data.phase == 'Liquid').order_by(Data.molecule_id, Data.property_id, Data.paper_id)
    for data in datas:
        m_p_p = [data.molecule, data.property, data.paper]
        if m_p_p not in molecule_property_paper:
            molecule_property_paper.append(m_p_p)

    for i, [mol, pro, paper] in enumerate(molecule_property_paper):
        sys.stdout.write('\r %i / %i' % (i, len(molecule_property_paper)))
        datas = Data.query.filter(Data.phase == 'Liquid').filter(
            and_(Data.molecule == mol, Data.property == pro, Data.paper == paper)).filter(
            or_(and_(Data.p < 120, Data.p > 80), Data.p == None)).order_by(Data.t)
        if datas.count() < 5:
            continue
        splines = Spline.query.filter(Spline.molecule == mol).filter(Spline.property == pro).filter(Spline.paper == paper)
        if splines.count() == 0:
            spline = Spline(molecule=mol, property=pro, paper=paper)
            spline.t_min = datas[0].t
            spline.t_max = datas[-1].t
            t_list = []
            v_list = []
            u_list = []
            for data in datas:
                t_list.append(data.t)
                v_list.append(data.value)
                u_list.append(data.stderr)
            if len(t_list) != len(set(t_list)):
                [t_list, v_list, u_list] = remove_repeat_data(t_list, v_list, u_list)
            if len(t_list) < 5:
                continue
            t_list, v_list, u_list = zip(*sorted(zip(t_list, v_list, u_list)))
            spl_v = interpolate.splrep(t_list, v_list)
            spl_u = interpolate.splrep(t_list, u_list)
            spline.coef_v = json.dumps([spl_v[0].tolist(), spl_v[1].tolist(), spl_v[2]])
            spline.coef_u = json.dumps([spl_u[0].tolist(), spl_u[1].tolist(), spl_u[2]])
            if pro.id in [10, 40, 50]:
                coef, score = VTFfit(t_list, v_list)
                spline.coef_VTF = json.dumps([coef, score])
            db.session.add(spline)
        elif splines.count() == 1:
            spline = splines.first()
            if pro.id not in [10, 40, 50] or spline.coef_VTF is not None:
                continue
            t_list = []
            v_list = []
            u_list = []
            for data in datas:
                if data.value != 0.:
                    t_list.append(data.t)
                    v_list.append(data.value)
                    u_list.append(data.stderr)
            if len(t_list) != len(set(t_list)):
                [t_list, v_list, u_list] = remove_repeat_data(t_list, v_list, u_list)
            if len(t_list) < 5:
                continue
            t_list, v_list, u_list = zip(*sorted(zip(t_list, v_list, u_list)))

            coef, score = VTFfit(t_list, v_list)
            spline.coef_VTF = json.dumps([coef.tolist(), score])
        db.session.commit()

def density_verify():
    datas = Data.query.filter(Data.phase == 'Liquid').filter(Data.property_id == 9)
    for data in datas:
        if data.value < 500:
            print(data.molecule_id, data.value)
# density_verify()

def ion_verify():
    '''
    ions = Ion.query
    for ion in ions:
        if 'R' in ion.name:
            print(ion.id, ion.smiles, ion.name)
        '''

    ions = Ion.query
    smiles_list = []
    id_list = []
    for ion in ions:
        if ion.smiles==None:
            continue
        if ion.smiles not in smiles_list:
            smiles_list.append(ion.smiles)
            id_list.append(ion.id)
            ion.duplicate = ion.id
        else:
            ion.duplicate = id_list[smiles_list.index(ion.smiles)]
    db.session.commit()

    mols = Molecule.query
    for mol in mols:
        if mol.cation.duplicate !=None:
            mol.cation_id = mol.cation.duplicate
        if mol.anion.duplicate != None:
            mol.anion_id = mol.anion.duplicate
    db.session.commit()

    def get_can_smiles(smiles):
        py_mol = pybel.readstring("smi", smiles)
        return py_mol.write('can', opt={'n': None}).strip()

    ions = Ion.query
    for ion in ions:
        if ion.smiles!=None:
            ion.smiles = get_can_smiles(ion.smiles)
    db.session.commit()

def molecule_verify():
    mols = Molecule.query
    c_a_id = []
    for mol in mols:
        if mol.cation.smiles == None or mol.anion.smiles ==None:
            print(mol.id, 'smiles not found')
        ca = [mol.cation_id, mol.anion_id]
        if ca not in c_a_id:
            c_a_id.append(ca)
        else:
            print(mol.id, 'repeated')

def ion_delete_repeat():
    ions = Ion.query
    for ion in ions:
        if ion.duplicate != None and int(ion.duplicate) != ion.id:
            db.session.delete(ion)
    db.session.commit()

# ion_verify()
# molecule_verify()
# ion_delete_repeat()
# db.create_all()
add_spline_info_ilthermo()




print('end')




