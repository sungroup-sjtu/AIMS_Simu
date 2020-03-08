#!/usr/bin/env python3
# coding=utf-8

import sys
import matplotlib
from sqlalchemy import or_

matplotlib.use('Agg')

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_ilthermo import Property, Ion, Molecule, Spline, Data
from app.models_cv import Cv
from app.models_nist import NistProperty, NistMolecule, NistData
from mstools.analyzer.plot import *
from app.selection import *
from mstools.analyzer.fitting import VTFfit, VTFval


# simulation part
def get_simulation_value_stderr(property, pressure, name, post_result, task=None, gk=False, diff_sum=False):
    if property == 'density':
        f = open('%s-%i-den-sim.txt' % (name, pressure), 'w')
    elif property == 'viscosity' and not gk:
        f = open('%s-%i-vis-sim.txt' % (name, pressure), 'w')
    elif property == 'viscosity' and gk:
        f = open('%s-%i-vis_gk-sim.txt' % (name, pressure), 'w')
    elif property == 'diffusion constant':
        f = open('%s-%i-diff-sim.txt' % (name, pressure), 'w')
    elif property == 'diffusion constant gk':
        f = open('%s-%i-diff_gk-sim.txt' % (name, pressure), 'w')
    elif property == 'electrical conductivity':
        f = open('%s-%i-econ_gk-sim.txt' % (name, pressure), 'w')
    elif property == 'Nernst-Einstein electrical conductivity':
        f = open('%s-%i-econ_NE-sim.txt' % (name, pressure), 'w')
    elif property == 'cp':
        if task is None:
            return False
        return get_simulation_cp(pressure, task)
    elif property == 'hvap':
        if task is None:
            return False
        return get_simulation_hvap(pressure, name.split('-'), task)
    elif property == 'pzz' and pressure is None:
        f = open('%s-pvap-sim.txt' % name, 'w')
    elif property == 'dliq':
        f = open('%s-dliq-sim.txt' % name, 'w')
    elif property == 'dgas':
        f = open('%s-dgas-sim.txt' % name, 'w')
    elif property == 'st':
        f = open('%s-st-sim.txt' % name, 'w')
    else:
        return False

    t_list = []
    value_list = []
    for info in post_result.get(property):
        if property == 'density':
            [t, p, [value, stderr]] = info
        elif property == 'viscosity' and not gk:
            [t, p, [value, stderr], score] = info
            t_list.append(t)
            value_list.append(value)
        elif property == 'viscosity' and gk:
            [t, p, value, score] = info
            stderr = None
        elif property == 'diffusion constant':
            if diff_sum:
                [t, p, diff_dict] = info
                keys = list(diff_dict.keys())
                keys.remove('System')
                value, stderr = 0., 0.
                for i in keys:
                    value += diff_dict[i][0]
                    stderr += diff_dict[i][1]
            else:
                [t, p, diff_dict] = info
                value, stderr = diff_dict.get('System')
            if t > 320:
                t_list.append(t)
                value_list.append(value)
        elif property == 'diffusion constant gk':
            if diff_sum:
                [t, p, diff_dict] = info
                keys = list(diff_dict.keys())
                keys.remove('System')
                value = 0.
                for i in keys:
                    value += diff_dict[i][0]
            else:
                [t, p, diff_dict] = info
                value, score = diff_dict.get('System')
            stderr = None
        elif property == 'electrical conductivity':
            [t, p, value, score] = info
            stderr = None
        elif property == 'Nernst-Einstein electrical conductivity':
            [t, p, [value, stderr]] = info
            if t > 320:
                t_list.append(t)
                value_list.append(value)
        elif property == 'pzz':
            [t, [value, stderr]] = info
            p = None
        elif property == 'dliq':
            [t, [value, stderr]] = info
            p = None
        elif property == 'dgas':
            [t, [value, stderr]] = info
            p = None
        elif property == 'st':
            [t, [value, stderr]] = info
            p = None
        else:
            return False

        if property in ['dliq', 'dgas']:
            f.write('%#.5e\t%#.5e\t%#.5e\n' % (value, t, stderr))
        elif p == pressure:
            if stderr is not None:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (t, value, stderr))
            else:
                f.write('%#.5e\t%#.5e\n' % (t, value))
    f.close()
    if len(t_list) > 5:
        coef, score = VTFfit(t_list, value_list)
        print('%s: Tg=%f' % (property, coef[1]))
        if property == 'viscosity' and not gk:
            f = open('%s-%i-vis-vtf.txt' % (name, pressure), 'w')
        elif property == 'diffusion constant':
            f = open('%s-%i-diff-vtf.txt' % (name, pressure), 'w')
        elif property == 'Nernst-Einstein electrical conductivity':
            f = open('%s-%i-econ_NE-vtf.txt' % (name, pressure), 'w')
        t_list = np.linspace(250, 500, 26)
        value_list = VTFval(t_list, coef)
        for t, value in zip(t_list, value_list):
            f.write('%#.5e\t%#.5e\n' % (t, value))
    return True


def get_simulation_cp(pressure, task):
    smiles_list = json.loads(task.smiles_list)
    name = '-'.join(smiles_list)
    f = open('%s-%i-cp-sim.txt' % (name, pressure), 'w')
    t_list = json.loads(task.t_list)
    for t in t_list:
        cp_qm = 0.
        for smiles in smiles_list:
            cv = Cv.query.filter(Cv.smiles == smiles).first()
            if cv is None:
                return False
            cp_qm += cv.get_post_cv(t)
        post_data = task.get_post_data(t, pressure)
        cp_inter = post_data.get('cp_inter')
        if cp_inter is None:
            return False
        cp_pv = post_data.get('cp_pv')
        if cp_pv is None:
            return False
        cp_md = cp_inter + cp_pv
        cp_sim = cp_md + cp_qm
        f.write('%#.5e\t%#.5e\t%#.5e\t%#.5e\n' % (t, cp_sim, cp_md, cp_qm))
    return True


def get_simulation_hvap(pressure, smiles_list, task):
    name = '-'.join(smiles_list)
    f = open('%s-%i-hvap-sim.txt' % (name, pressure), 'w')
    t_list = json.loads(task.t_list)
    for t in t_list:
        post_data = task.get_post_data(t, pressure)
        hvap_sim = post_data['hvap']
        if hvap_sim is None:
            return False
        f.write('%#.5e\t%#.5e\n' % (t, hvap_sim))
    return True


# experimental part
def get_exp(pressure, property, DATAS, name, type, simplify=False, dt=10):
    def get_simplify_data(_datas):
        _datas = _datas.order_by(Data.t)
        if _datas.count() < 10:
            return list(_datas)
        t_min, t_max = _datas[0].t, _datas[-1].t
        t = t_min
        data_list = []
        while t < t_max:
            best_data = _datas.filter(Data.t > t).filter(Data.t < t + dt).order_by(Data.stderr).first()
            if best_data is not None:
                data_list.append(best_data)
            t += dt
        return data_list

    if type == 'ilthermo':
        datas = DATAS.filter(Data.property == property).filter(Data.p > 98 * pressure).filter(
            Data.p < 102 * pressure).order_by(Data.t)
        if simplify:
            datas = get_simplify_data(datas)
        else:
            datas = list(datas)
        if len(datas) == 0:
            return False
        if property.name == 'Density':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'den'), 'w')
            for data in datas:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value / 1000, data.stderr / 1000))
        elif property.name == 'Heat capacity at constant pressure':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'cp'), 'w')
            for data in datas:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value, data.stderr))
        elif property.name == 'Viscosity':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'vis'), 'w')
            for data in datas:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value * 1000, data.stderr * 1000))
        elif property.name == 'Self-diffusion coefficient':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'diff'), 'w')
            for data in datas:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value * 1e-5, data.stderr * 1e-5))
        elif property.name == 'Electrical conductivity':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'econ'), 'w')
            for data in datas:
                f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value, data.stderr))
        else:
            return False
    elif type == 'nist':
        if pressure is None:
            datas = DATAS.filter(NistData.property == property).filter(NistData.p == None)
        elif pressure == 1:
            datas = DATAS.filter(NistData.property == property).filter(
                or_(NistData.p == None, and_(NistData.p > 98 * pressure, NistData.p < 102 * pressure))).order_by(
                NistData.t)
        else:
            datas = DATAS.filter(NistData.property == property).filter(
                and_(NistData.p > 98 * pressure, NistData.p < 102 * pressure)).order_by(NistData.t)
        if datas.count() == 0:
            return False
        if property.name == 'density-lg' and pressure is not None:
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'den'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value / 1000, data.uncertainty / 1000))
        elif property.name == 'cp-lg':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'cp'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value, data.uncertainty))
        elif property.name == 'hvap-lg':
            f = open('%s-%i-%s-exp.txt' % (name, pressure, 'hvap'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value, data.uncertainty))
        elif property.name == 'pvap-lg':
            f = open('%s-%s-exp.txt' % (name, 'pvap'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value / 100, data.uncertainty / 100))
        elif property.name == 'density-lg' and pressure is None:
            f = open('%s-%s-exp.txt' % (name, 'dliq'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.value / 1000, data.t, data.uncertainty / 1000))
        elif property.name == 'density-gl':
            f = open('%s-%s-exp.txt' % (name, 'dgas'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.value / 1000, data.t, data.uncertainty / 1000))
        elif property.name == 'st-lg':
            f = open('%s-%s-exp.txt' % (name, 'st'), 'w')
            for data in datas:
                if data.uncertainty < data.value:
                    f.write('%#.5e\t%#.5e\t%#.5e\n' % (data.t, data.value * 1000, data.uncertainty * 1000))
        else:
            return False
    return True


def compare_npt(exp_type, n_tasks=None, simplify_exp=False):
    tasks = Task.query.filter(Task.procedure == 'npt').filter(
        Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
    if n_tasks is not None:
        tasks = tasks.limit(n_tasks)
    for task in tasks:
        print(task)
        if not task_selection(task, select=opt.selection):
            continue

        # P_list = json.loads(task.p_list)
        P_list = [1]
        smiles_list = json.loads(task.smiles_list)
        post_result = json.loads(task.post_result)
        for P in P_list:
            # simulation part
            DEN_SIM = get_simulation_value_stderr(property='density', pressure=P, name='-'.join(smiles_list),
                                                  post_result=post_result)
            CP_SIM = get_simulation_value_stderr(property='cp', pressure=P, name='-'.join(smiles_list),
                                                 post_result=post_result, task=task)

            if exp_type == 'nist':
                HVAP_SIM = get_simulation_value_stderr(property='hvap', pressure=P, name='-'.join(smiles_list),
                                                       post_result=post_result, task=task)

            # experimental part
            if exp_type == 'ilthermo':
                prop_density = Property.query.filter(Property.name == 'Density').first()
                prop_cp = Property.query.filter(Property.name == 'Heat capacity at constant pressure').first()
                cation_smiles, anion_smiles = smiles_list
                cation_id = Ion.query.filter(Ion.smiles == cation_smiles).first().id
                anion_id = Ion.query.filter(Ion.smiles == anion_smiles).first().id
                molecules = Molecule.query.filter(Molecule.cation_id == cation_id).filter(Molecule.anion_id == anion_id)
                if molecules.count() == 1:
                    molecule = molecules.first()
                    datas = Data.query.filter(Data.molecule == molecule).filter(Data.phase == 'Liquid')
                    if get_exp(pressure=P, property=prop_density, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s-%s' % (cation_smiles, anion_smiles)) and DEN_SIM:
                        name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'den')
                        gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                                ylabel='density(g·cm^{-3})', font=40,
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'errorlines'],
                                title_list=['exp', 'sim'])

                    if get_exp(pressure=P, property=prop_cp, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s-%s' % (cation_smiles, anion_smiles)) and CP_SIM:
                        name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'cp')
                        gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                                ylabel='heat capacity(J·mol^{-1}·K^{-1})', font=40,
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'lines-3'],
                                title_list=['exp', ['sim', 'md', 'qm']])
            elif exp_type == 'nist':
                prop_density = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
                prop_cp = NistProperty.query.filter(NistProperty.name == 'cp-lg').first()
                prop_hvap = NistProperty.query.filter(NistProperty.name == 'hvap-lg').first()
                molecules = NistMolecule.query.filter(NistMolecule.smiles == smiles_list[0])
                if molecules.count() == 0:
                    molecules = NistMolecule.query.filter(
                        NistMolecule.smiles == get_canonical_smiles(smiles_list[0].replace('@', '')))
                if molecules.count() == 1:
                    molecule = molecules.first()
                    datas = NistData.query.filter(NistData.molecule == molecule)
                    if get_exp(pressure=P, property=prop_density, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s' % (smiles_list[0])) and DEN_SIM:
                        name = '%s-%i-%s' % (smiles_list[0], P, 'den')
                        gnuplot(output=name, title=name, xlabel='temperature(K)', ylabel='density(g·cm^{-3})', font=40,
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'errorlines'],
                                title_list=['exp', 'sim'])

                    if get_exp(pressure=P, property=prop_cp, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s' % (smiles_list[0])) and CP_SIM:
                        name = '%s-%i-%s' % (smiles_list[0], P, 'cp')
                        gnuplot(output=name, title=name, xlabel='temperature(K)', font=40,
                                ylabel='heat capacity(J·mol^{-1}·K^{-1})',
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'lines-3'],
                                title_list=['exp', ['sim', 'md', 'qm']])

                    if get_exp(pressure=P, property=prop_hvap, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s' % (smiles_list[0])) and HVAP_SIM:
                        name = '%s-%i-%s' % (smiles_list[0], P, 'hvap')
                        gnuplot(output=name, title=name, xlabel='temperature(K)', font=40,
                                ylabel='enthalpy of vaporization(kJ·mol^{-1})',
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'lines'],
                                title_list=['exp', 'sim'])


def compare_nvt_slab(exp_type, n_tasks=None, simplify_exp=False):
    tasks = Task.query.filter(Task.procedure == 'nvt-slab').filter(
        Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
    if n_tasks is not None:
        tasks = tasks.limit(n_tasks)
    for task in tasks:
        print(task)
        if not task_selection(task, select=opt.selection):
            continue

        smiles_list = json.loads(task.smiles_list)
        post_result = json.loads(task.post_result)
        # simulation part
        PVAP_SIM = get_simulation_value_stderr(property='pzz', pressure=None, name='-'.join(smiles_list),
                                               post_result=post_result)
        DLIQ_SIM = get_simulation_value_stderr(property='dliq', pressure=None, name='-'.join(smiles_list),
                                               post_result=post_result)
        DGAS_SIM = get_simulation_value_stderr(property='dgas', pressure=None, name='-'.join(smiles_list),
                                               post_result=post_result)
        ST_SIM = get_simulation_value_stderr(property='st', pressure=None, name='-'.join(smiles_list),
                                             post_result=post_result)
        if exp_type == 'nist':
            prop_pvap = NistProperty.query.filter(NistProperty.name == 'pvap-lg').first()
            prop_dliq = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
            prop_dgas = NistProperty.query.filter(NistProperty.name == 'density-gl').first()
            prop_st = NistProperty.query.filter(NistProperty.name == 'st-lg').first()
            molecules = NistMolecule.query.filter(NistMolecule.smiles == smiles_list[0])
            if molecules.count() == 1:
                molecule = molecules.first()
                datas = NistData.query.filter(NistData.molecule == molecule)
                if get_exp(pressure=None, property=prop_pvap, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s' % (smiles_list[0])) and PVAP_SIM:
                    name = '%s-%s' % (smiles_list[0], 'pvap')
                    gnuplot(output=name, title=name, xlabel='temperature(K)', ylabel='vapor pressure(bar)', font=40,
                            txt_list=[name + '-exp.txt', name + '-sim.txt'],
                            type_list=['errorbars', 'errorlines'],
                            title_list=['exp', 'sim'])
                if get_exp(pressure=None, property=prop_dliq, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s' % (smiles_list[0])) and DLIQ_SIM and \
                        get_exp(pressure=None, property=prop_dgas, DATAS=datas, type=exp_type, simplify=simplify_exp,
                                name='%s' % (smiles_list[0])) and DGAS_SIM:
                    name = smiles_list[0]
                    gnuplot(output=name + '-vle', title=name + '-vle', xlabel='density(g·cm^{-3})', font=40,
                            ylabel='temperature(K)',
                            txt_list=[name + '-dliq-exp.txt', name + '-dgas-exp.txt', name + '-dliq-sim.txt',
                                      name + '-dgas-sim.txt'],
                            type_list=['xerrorbars', 'xerrorbars', 'lines', 'lines'],
                            title_list=['dliq-exp', 'dgas-exp', 'dliq-sim', 'dgas-sim'])
                if get_exp(pressure=None, property=prop_st, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s' % (smiles_list[0])) and ST_SIM:
                    name = '%s-%s' % (smiles_list[0], 'st')
                    gnuplot(output=name, title=name, xlabel='temperature(K)', ylabel='surface tension', font=40,
                            txt_list=[name + '-exp.txt', name + '-sim.txt'],
                            type_list=['errorbars', 'errorlines'],
                            title_list=['exp', 'sim'])


def compare_ppm(exp_type, n_tasks=None, simplify_exp=False):
    tasks = Task.query.filter(Task.procedure == 'ppm').filter(
        Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
    if n_tasks is not None:
        tasks = tasks.limit(n_tasks)
    for task in tasks:
        print(task)
        if not task_selection(task, select=opt.selection):
            continue

        P_list = [1]
        # P_list = json.loads(task.p_list)
        smiles_list = json.loads(task.smiles_list)
        post_result = json.loads(task.post_result)
        for P in P_list:
            # simulation part
            VIS_SIM = get_simulation_value_stderr(property='viscosity', pressure=P, name='-'.join(smiles_list),
                                                  post_result=post_result)
            # experimental part
            if exp_type == 'nist':
                prop_vis = NistProperty.query.filter(NistProperty.name == 'viscosity-lg').first()
                molecules = NistMolecule.query.filter(NistMolecule.smiles == smiles_list[0])
                if molecules.count() == 0:
                    molecules = NistMolecule.query.filter(
                        NistMolecule.smiles == get_canonical_smiles(smiles_list[0].replace('@', '')))
                if molecules.count() == 1:
                    molecule = molecules.first()
                    datas = NistData.query.filter(NistData.molecule == molecule)
                    if get_exp(pressure=P, property=prop_vis, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s' % (smiles_list[0])) and VIS_SIM:
                        name = '%s-%i-%s' % (smiles_list[0], P, 'vis')
                        gnuplot(output=name, title=name, xlabel='temperature(K)', ylabel='viscosity(mPa·s)', font=40,
                                txt_list=[name + '-exp.txt', name + '-sim.txt'],
                                type_list=['errorbars', 'errorlines'],
                                title_list=['exp', 'sim-ppm'], y_log_scale=True, reciprical_x=True)

            elif exp_type == 'ilthermo':
                prop_vis = Property.query.filter(Property.name == 'Viscosity').first()
                cation_smiles, anion_smiles = smiles_list
                cation_id = Ion.query.filter(Ion.smiles == cation_smiles).first().id
                anion_id = Ion.query.filter(Ion.smiles == anion_smiles).first().id
                molecules = Molecule.query.filter(Molecule.cation_id == cation_id).filter(Molecule.anion_id == anion_id)
                if molecules.count() == 1:
                    molecule = molecules.first()
                    datas = Data.query.filter(Data.molecule == molecule).filter(Data.phase == 'Liquid')
                    if get_exp(pressure=P, property=prop_vis, DATAS=datas, type=exp_type, simplify=simplify_exp,
                               name='%s-%s' % (cation_smiles, anion_smiles)):
                        name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'vis')
                        txt_list = [name + '-exp.txt', name + '-sim.txt']
                        type_list = ['errorbars', 'errorlines']
                        title_list = ['exp', 'sim-ppm']
                        if os.path.exists(name + '-vtf.txt'):
                            txt_list += [name + '-vtf.txt']
                            type_list += ['lines']
                            title_list += ['sim-ppm-vtf']
                        gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                                ylabel='viscosity(mPa·s)', txt_list=txt_list, type_list=type_list, title_list=title_list,
                                y_log_scale=True, reciprical_x=True, font=40)


def compare_nvt_multi(exp_type, n_tasks=None, diff_gk=False, simplify_exp=False):
    if exp_type != 'ilthermo':
        return
    prop_econ = Property.query.filter(Property.name == 'Electrical conductivity').first()
    prop_vis = Property.query.filter(Property.name == 'Viscosity').first()
    prop_diff = Property.query.filter(Property.name == 'Self-diffusion coefficient').first()
    tasks = Task.query.filter(Task.procedure == 'nvt-multi').filter(
        Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
    if n_tasks is not None:
        tasks = tasks.limit(n_tasks)
    for task in tasks:
        print(task)
        if not task_selection(task, select=opt.selection):
            continue

        P_list = json.loads(task.p_list)
        cation_smiles, anion_smiles = json.loads(task.smiles_list)
        # simulation part
        post_result = json.loads(task.post_result)
        for P in P_list:
            get_simulation_value_stderr(property='viscosity', pressure=P, name='%s-%s' % (cation_smiles, anion_smiles),
                                        post_result=post_result, gk=True)
            get_simulation_value_stderr(property='diffusion constant', pressure=P,
                                        name='%s-%s' % (cation_smiles, anion_smiles),
                                        post_result=post_result, diff_sum=True)
            if diff_gk:
                get_simulation_value_stderr(property='diffusion constant gk', pressure=P,
                                            name='%s-%s' % (cation_smiles, anion_smiles),
                                            post_result=post_result, diff_sum=True)
            get_simulation_value_stderr(property='electrical conductivity', pressure=P,
                                        name='%s-%s' % (cation_smiles, anion_smiles),
                                        post_result=post_result)
            get_simulation_value_stderr(property='Nernst-Einstein electrical conductivity', pressure=P,
                                        name='%s-%s' % (cation_smiles, anion_smiles),
                                        post_result=post_result)

        cation_id = Ion.query.filter(Ion.smiles == cation_smiles).first().id
        anion_id = Ion.query.filter(Ion.smiles == anion_smiles).first().id
        # experimental part
        molecules = Molecule.query.filter(Molecule.cation_id == cation_id).filter(Molecule.anion_id == anion_id)
        if molecules.count() == 1:
            molecule = molecules.first()
            datas = Data.query.filter(Data.molecule == molecule).filter(Data.phase == 'Liquid')
            for P in P_list:
                if get_exp(pressure=P, property=prop_vis, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s-%s' % (cation_smiles, anion_smiles)):
                    name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'vis')
                    txt_list = [name + '-exp.txt', name + '_gk-sim.txt', name + '-sim.txt']
                    type_list = ['errorbars', 'lines', 'errorlines']
                    title_list = ['exp', 'sim-gk', 'sim-ppm']
                    if os.path.exists(name + '-vtf.txt'):
                        txt_list += [name + '-vtf.txt']
                        type_list += ['lines']
                        title_list += ['sim-ppm-vtf']
                    gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                            ylabel='viscosity(mPa·s)', txt_list=txt_list, type_list=type_list, title_list=title_list,
                            y_log_scale=True, reciprical_x=True, font=40)

                if get_exp(pressure=P, property=prop_diff, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s-%s' % (cation_smiles, anion_smiles)):
                    name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'diff')
                    txt_list = [name + '-exp.txt', name + '-sim.txt']
                    type_list = ['errorbars', 'errorlines']
                    title_list = ['exp', 'sim-Einstein']
                    if diff_gk:
                        txt_list += [name + '_gk-sim.txt']
                        type_list += ['lines']
                        title_list += ['sim-gk']
                    if os.path.exists(name + '-vtf.txt'):
                        txt_list += [name + '-vtf.txt']
                        type_list += ['lines']
                        title_list += ['sim-Einstein-vtf']
                    gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                            ylabel='diffusion constant(cm^2·s^{-1})', txt_list=txt_list, type_list=type_list,
                            title_list=title_list, y_log_scale=True, reciprical_x=True, xtics=0.0005, font=40)

                if get_exp(pressure=P, property=prop_econ, DATAS=datas, type=exp_type, simplify=simplify_exp,
                           name='%s-%s' % (cation_smiles, anion_smiles)):
                    name = '%s-%s-%i-%s' % (cation_smiles, anion_smiles, P, 'econ')
                    txt_list = [name + '-exp.txt', name + '_NE-sim.txt', name + '_gk-sim.txt']
                    type_list = ['errorbars', 'errorlines', 'lines']
                    title_list = ['exp', 'sim-Nernst-Einstein', 'sim-gk']
                    if os.path.exists(name + '_NE-vtf.txt'):
                        txt_list += [name + '_NE-vtf.txt']
                        type_list += ['lines']
                        title_list += ['sim-Nernst-Einstein-vtf']
                    gnuplot(output=name, title=name + '-%i-%i' % (cation_id, anion_id), xlabel='temperature(K)',
                            ylabel='electrical conductivity(S·m^{-1})', txt_list=txt_list, type_list=type_list,
                            title_list=title_list, y_log_scale=True, reciprical_x=True, xtics=0.0005, font=40)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='This is a code to compare simulation results with experimental \
        results. For one molecule and one property, a png file will generated show the detailed comparison between \
        simulation and experiment')
    parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
    parser.add_argument('-t', '--type', type=str, help='The type of the experimental database(ilthermo or nist)')
    parser.add_argument('--selection', help='select specific task', default=False, action='store_true')
    parser.add_argument('--simplify_exp', help='Only output the experimental data with lowest standard deviation',
                        default=False, action='store_true')
    opt = parser.parse_args()
    procedure = opt.procedure
    app = create_app(procedure)
    app.app_context().push()
    if procedure == 'npt':
        compare_npt(exp_type=opt.type, simplify_exp=opt.simplify_exp)
    elif procedure == 'nvt-slab':
        compare_nvt_slab(exp_type=opt.type, simplify_exp=opt.simplify_exp)
    elif procedure == 'ppm':
        compare_ppm(exp_type=opt.type, simplify_exp=opt.simplify_exp)
    elif procedure == 'nvt-multi':
        compare_nvt_multi(exp_type=opt.type, simplify_exp=opt.simplify_exp)
    else:
        print('Unknown procedure')
