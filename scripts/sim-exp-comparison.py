#!/usr/bin/env python3
# coding=utf-8

import os, sys
sys.path.append('..')
import json
from sqlalchemy import or_
from app import create_app
from app.models import Task
from app.models_ilthermo import Data, Ion, Molecule, Property
from mstools.analyzer.fitting import polyfit
from app.models_qmcv import IonQM


procedure = sys.argv[1]
app = create_app(sys.argv[1])
app.app_context().push()
CWD = os.getcwd()


def GetUnsignedError(a, b):
    if a>b:
        return (a - b) / b
    elif a<b:
        return (b - a) / a
    else:
        return 0
def get_ion_id(smiles):
    ion = Ion.query.filter(Ion.smiles == smiles).first()
    return ion.duplicate

def get_data(cation_smiles, anion_smiles, p_name):
    prop = Property.query.filter(Property.name == p_name).first()
    cation_id = get_ion_id(cation_smiles)
    anion_id = get_ion_id(anion_smiles)
    mols = Molecule.query.filter(Molecule.cation_true_id == cation_id) \
        .filter(Molecule.anion_true_id == anion_id) \
        .all()
    data = []
    for mol in mols:
        datas = mol.datas.filter(Data.property == prop) \
            .filter(Data.phase == 'Liquid') \
            .filter(Data.value > 0.01) \
            .filter(or_(Data.p == None, Data.p<102)) \
            .all()
        data += datas
    data.sort(key=lambda x: x.t)
    return data

# density
def get_density_comparison(T='all', t_control=True):
    f = open('density_comparison_%s.xvg' % (str(T)), 'w')
    final_info = []
    for task in Task.query.filter(Task.procedure == procedure).filter(Task.post_result != None):
        cation_smiles, anion_smiles = json.loads(task.smiles_list)
        datas = get_data(cation_smiles, anion_smiles, 'Density')
        score = 0
        T_list = []
        den_list = []
        count_list = []
        if len(datas) > 5:
            for data in datas:
                if len(T_list) == 0:
                    T_list.append(data.t)
                    den_list.append(data.value / 1000)
                    count_list.append(1)
                else:
                    if data.t > T_list[-1]:
                        T_list.append(data.t)
                        den_list.append(data.value / 1000)
                        count_list.append(1)
                    elif data.t == T_list[-1]:
                        den_list[-1] += data.value / 1000
                        count_list[-1] += 1
            for i in range(len(T_list)):
                den_list[i] /= count_list[i]
            coeff, score = polyfit(T_list, den_list, 2)

        if score > 0.9:
            t_min = T_list[0]
            t_max = T_list[-1]
            t_p_den_stderr_list = json.loads(task.post_result).get('density')
            for info in t_p_den_stderr_list:
                t = info[0]
                p = 1
                if (t_control and t_min < t < t_max and p == 1) or (not t_control and p == 1):
                    if T =='all' or t==T:
                        den_sim = info[2][0]
                        den_exp = coeff[0] + coeff[1] * t + coeff[2] * t ** 2
                        final_info.append([den_sim, den_exp, GetUnsignedError(den_sim, den_exp), t, task.name, task.smiles_list])
    final_info.sort(key=lambda x: x[2])
    error = 0.
    for info in final_info:
        f.write('%.3f %.3f %.3f %i %s %s\n' % (info[0], info[1], info[2], info[3], info[4], info[5]))
        error += info[2]
    if len(final_info)!=0:
        error /= len(final_info)
        print('The average unsigned error for density for T=%s is %.3f' % (T, error))



# heat capacity
def get_Intra_CV(smiles, t):
    temp = IonQM.query.filter(IonQM.smiles==smiles).first()
    para = json.loads(temp.cv_parameter)
    return para.get('constant') + para.get('1') * t + para.get('2') * t**2 + para.get('3') * t**3 + para.get('4') * t**4
fp = open('p_mol.txt', 'w')
def get_cp_comparison(T='all', t_control=True):
    f = open('cp_comparison_%s.xvg' % (str(T)), 'w')
    final_info = []
    for task in Task.query.filter(Task.procedure == procedure).filter(Task.post_result != None):
        cation_smiles, anion_smiles = json.loads(task.smiles_list)
        if IonQM.query.filter(IonQM.smiles == cation_smiles).first() == None:
            fp.write('%s\n' % (cation_smiles))
            continue
        if IonQM.query.filter(IonQM.smiles == anion_smiles).first() == None:
            fp.write('%s\n' % (cation_smiles))
            continue
        datas = get_data(cation_smiles, anion_smiles, 'Heat capacity at constant pressure')
        score = 0
        T_list = []
        den_list = []
        count_list = []
        if len(datas) > 5:
            for data in datas:
                if len(T_list) == 0:
                    T_list.append(data.t)
                    den_list.append(data.value)
                    count_list.append(1)
                else:
                    if data.t > T_list[-1]:
                        T_list.append(data.t)
                        den_list.append(data.value)
                        count_list.append(1)
                    elif data.t == T_list[-1]:
                        den_list[-1] += data.value
                        count_list[-1] += 1
            for i in range(len(T_list)):
                den_list[i] /= count_list[i]
            coeff, score = polyfit(T_list, den_list, 2)

        if score > 0.9:
            t_min = T_list[0]
            t_max = T_list[-1]
            einter_t_poly3 = json.loads(task.post_result).get('einter-t-poly3')
            t_p_einter_stderr_list = json.loads(task.post_result).get('einter')
            for info in t_p_einter_stderr_list:
                t = info[0]
                p = 1
                if (t_control and t_min < t < t_max and p == 1) or (not t_control and p == 1):
                    if T == 'all' or t == T:
                        cp_sim = einter_t_poly3[0][1] + einter_t_poly3[0][2] * 2 * t + einter_t_poly3[0][
                            3] * 3 * t ** 2
                        cp_sim *= 1000
                        cp_qm = get_Intra_CV(cation_smiles, t) + get_Intra_CV(anion_smiles, t)
                        cp_cal = cp_sim + cp_qm
                        cp_exp = coeff[0] + coeff[1] * t + coeff[2] * t ** 2
                        final_info.append(
                            [cp_cal, cp_exp, GetUnsignedError(cp_cal, cp_exp), t, task.name, task.smiles_list])
    final_info.sort(key=lambda x: x[2])
    error = 0.
    for info in final_info:
        f.write('%.3f %.3f %.3f %i %s %s\n' % (info[0], info[1], info[2], info[3], info[4], info[5]))
        error += info[2]
    if len(final_info)!=0:
        error /= len(final_info)
        print('The average unsigned error for cp for T=%s is %.3f' % (T, error))

# econ
def get_econ_comparison(T='all', t_control=True):
    f = open('econ_comparison_%s.xvg' % (str(T)), 'w')
    final_info = []
    for task in Task.query.filter(Task.procedure == procedure).filter(Task.post_result != None):
        cation_smiles, anion_smiles = json.loads(task.smiles_list)
        datas = get_data(cation_smiles, anion_smiles, 'Electrical conductivity')
        score = 0
        T_list = []
        econ_list = []
        count_list = []
        if len(datas) > 5:
            for data in datas:
                if len(T_list) == 0:
                    T_list.append(data.t)
                    econ_list.append(data.value)
                    count_list.append(1)
                else:
                    if data.t > T_list[-1]:
                        T_list.append(data.t)
                        econ_list.append(data.value)
                        count_list.append(1)
                    elif data.t == T_list[-1]:
                        econ_list[-1] += data.value
                        count_list[-1] += 1
            for i in range(len(T_list)):
                econ_list[i] /= count_list[i]
            coeff, score = polyfit(T_list, econ_list, 2)

        if score > 0.9:
            t_min = T_list[0]
            t_max = T_list[-1]
            t_p_econ_stderr_list = json.loads(task.post_result).get('electrical conductivity from diffusion constant')
            for info in t_p_econ_stderr_list:
                t = info[0]
                p = 1
                if (t_control and t_min < t < t_max and p == 1) or (not t_control and p == 1):
                    if T =='all' or t==T:
                        econ_sim = info[2][0]
                        econ_exp = coeff[0] + coeff[1] * t + coeff[2] * t ** 2
                        if econ_sim > 1.0:
                            final_info.append([econ_sim, econ_exp, GetUnsignedError(econ_sim, econ_exp), t, task.name, task.smiles_list])
    final_info.sort(key=lambda x: x[2])
    error = 0.
    for info in final_info:
        f.write('%.3f %.3f %.3f %i %s %s\n' % (info[0], info[1], info[2], info[3], info[4], info[5]))
        error += info[2]
    if len(final_info)!=0:
        error /= len(final_info)
        print('The average unsigned error for econ for T=%s is %.3f' % (T, error))

#viscosity
def get_vis_comparison(T='all', t_control=True):
    from mstools.analyzer.fitting import fit_VTF
    import numpy as np
    f = open('vis_comparison_%s.xvg' % (str(T)), 'w')
    final_info = []
    for task in Task.query.filter(Task.procedure == procedure).filter(Task.post_result != None):
        cation_smiles, anion_smiles = json.loads(task.smiles_list)
        datas = get_data(cation_smiles, anion_smiles, 'Viscosity')
        score = 0
        T_list = []
        vis_list = []
        count_list = []
        if len(datas) > 5:
            for data in datas:
                if len(T_list) == 0:
                    T_list.append(data.t)
                    vis_list.append(data.value)
                    count_list.append(1)
                else:
                    if data.t > T_list[-1]:
                        T_list.append(data.t)
                        vis_list.append(data.value)
                        count_list.append(1)
                    elif data.t == T_list[-1]:
                        vis_list[-1] += data.value
                        count_list[-1] += 1
            for i in range(len(T_list)):
                vis_list[i] /= count_list[i]
            Y_fit = np.log(vis_list)
            coeff, score = fit_VTF(T_list, Y_fit)

        if score > 0.9:
            t_min = T_list[0]
            t_max = T_list[-1]
            t_p_vis_score_list = json.loads(task.post_result).get('t_p_vis_score_list')
            for info in t_p_vis_score_list:
                t = info[0]
                p = 1
                if (t_control and t_min < t < t_max and p == 1) or (not t_control and p == 1):
                    if T =='all' or t==T:
                        vis_sim = info[2][0]
                        vis_exp = np.exp(coeff[0] + coeff[2] / (t - coeff[1])) * 1000
                        final_info.append([vis_sim, vis_exp, GetUnsignedError(vis_sim, vis_exp), t, task.name, task.smiles_list, info[2][1]])
    final_info.sort(key=lambda x: x[2])
    error = 0.
    for info in final_info:
        f.write('%.3f %.3f %.3f %i %s %s %.3f\n' % (info[0], info[1], info[2], info[3], info[4], info[5], info[6]))
        error += info[2]
    if len(final_info)!=0:
        error /= len(final_info)
        print('The average unsigned error for viscosity for T=%s is %.3f' % (T, error))


T_list = [300, 350, 400, 450, 500]
if procedure == 'npt':
    get_density_comparison()
    for T in T_list:
        get_density_comparison(T=T)

    get_cp_comparison()
    for T in T_list:
        get_cp_comparison(T=T)

    get_econ_comparison()
    for T in T_list:
        get_econ_comparison(T=T)

if procedure == 'ppm':
    get_vis_comparison()
    get_vis_comparison()
    for T in T_list:
        get_vis_comparison(T=T)
