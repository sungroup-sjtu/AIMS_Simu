#!/usr/bin/env python3
# coding=utf-8

import os, sys, time
import json
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO

sys.path.append('..')
from app import create_app
from app.models import Task, Compute, PbsJob
from app.models_nist import NistMolecule, NistProperty, NistSpline
from app.models_cv import Cv
from app.selection import task_selection

class StatAction():
    def __init__(self):
        self.nist_list = []
        self.task_list = []
        self.TP_list = {
            # 298: [],
            # 'Tm25': [],
            # 'Tvap': [],
            # 'Tcx8': [],
            # 'Tc85': [],
        }

    def update_mol_task_list(self, procedure):
        print('Get molecule list')
        tasks = Task.query.filter(Task.procedure == procedure).filter(Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
        for task in tasks:
            print(task)
            if task_selection(task, select=opt.selection):
                smiles = json.loads(task.smiles_list)[0]
                nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
                if nist is None:
                    continue

                if task.post_result is None:
                    continue
                if procedure == 'npt' and task.get_post_result()['density-poly4'][-1] < 0.999:
                    continue

                self.nist_list.append(nist)
                self.task_list.append(task)

                for _T in self.TP_list.keys():
                    T = self.get_T_nist(nist, _T)

                    if T is None:
                        P = None
                    else:
                        P = self.get_P_nist(nist, T)

                    self.TP_list[_T].append((T, P))
    '''
    def update_mol_slab_list(self):
        print('Get molecule list')

        tasks = Task.query
        for slab in tasks:
            print(slab)
            if task_selection(slab):
                smiles = json.loads(slab.smiles_list)[0]
                nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
                if nist is None:
                    continue

                if slab.post_result is None:
                    continue

                self.nist_list.append(nist)
                self.slab_list.append(slab)

                for _T in self.TP_list.keys():
                    T = self.get_T_nist(nist, _T)
                    P = None

                    self.TP_list[_T].append((T, P))
    '''

    def get_T_nist(self, nist, _T):
        if type(_T) == str and nist is None:
            return None

        if _T == 'Tm25':
            T = nist.tt
            if T is not None:
                T += 25
        elif _T == 'Tvap':
            T = nist.tb
        elif _T == 'Tcx8':
            T = nist.tc
            if T is not None:
                T *= 0.8
                T = min(T, 600)
        elif _T == 'Tc85':
            T = nist.tc
            if T is not None:
                T *= 0.85
                T = min(T, 650)
        else:
            T = _T

        if T is not None:
            T = int(round(T))

        return T

    def get_P_nist(self, nist: NistMolecule, T):
        if nist is None or nist.tb is None:
            return None

        if T < nist.tb + 1:
            return 1 # bar

        spline_pvap = nist.splines.filter(NistSpline.property_id == prop_pvap.id).first()
        if spline_pvap is None:
            P = None
        else:
            P = spline_pvap.get_data(T)[0]  # kPa
            if P is not None:
                P /= 100  # bar
        return P

    def get_property(self, property, property_name, _T):
        print('Get %s from nist' % (property.name))
        value_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))
            # get comparison temperature and pressure
            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue
            # get experimental value
            spline = nist.splines.filter(NistSpline.property_id == property.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue
            # unit transform
            if property == prop_density:
                value /= 1000
                uncertainty /= 1000
            elif property == prop_st:
                value *= 1000
                uncertainty *= 1000
            elif property == prop_pvap:
                value /= 100
                uncertainty /= 100
            # get simulation value
            task = self.task_list[i]
            post_data = task.get_post_data(T, P)
            if property_name == 'cp':
                cp_inter = post_data.get('cp_inter')
                cp_pv = post_data.get('cp_pv')
                if cp_inter is None or cp_pv is None:
                    continue
                cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
                if cv is None:
                    continue
                value_sim = cp_inter + cp_pv + cv.get_post_cv(T)
            else:
                value_sim = post_data.get(property_name)
            if value_sim is None:
                continue
            # correction for simulation value
            ### Correction for hvap
            if property_name == 'hvap':
                value_sim = value_sim - nist.n_nothx / 15 * (8.314 * T) / 1000
            value_list.append([nist.smiles, value, uncertainty, value_sim])
        return value_list

    def get_hl_nist(self, _T=298):
        print('Get liquid enthalpy from nist')
        hl_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.filter(NistSpline.property_id == prop_hl.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue
            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv is None:
                continue

            task = self.task_list[i]
            post_data = task.get_post_data(T, P)
            hl_sim = post_data['liquid enthalpy']
            if hl_sim is None:
                continue
            hl_sim += cv.get_post_enthalpy(T)

            hl_exp_list = []
            hl_sim_list = []
            tmin = max(spline.t_min, task.t_min)
            tmax = min(spline.t_max, task.t_max)
            t_list = np.linspace(tmin, tmax, 10)
            for t in t_list:
                v, u = spline.get_data(t)
                p = self.get_P_nist(nist, t)
                if p==None:
                    continue
                p_d = task.get_post_data(t, p)
                if p_d.get('liquid enthalpy')==None:
                    continue
                hl_exp_list.append(v)
                hl_sim_list.append(p_d.get('liquid enthalpy') + cv.get_post_enthalpy(T))
            # shift = np.array(hl_exp_list).mean() - np.array(hl_sim_list).mean()
            # hl_list.append([nist.smiles, value, uncertainty, hl_sim + shift])
            shift1 = hl_exp_list[0]
            shift2 = hl_sim_list[0]
            hl_list.append([nist.smiles, value-shift1, uncertainty, hl_sim - shift2])
            '''
            f = open(str(i),'w')
            for j, t in enumerate(t_list):
                f.write('%f %f %f\n' % (t, hl_exp_list[j], hl_sim_list[j]+shift))
            f.close()
            '''
        print('')
        return hl_list
    
    def get_sound_nist(self, _T=298):
        print('Get cSound from nist')
        sound_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.filter(NistSpline.property_id == prop_sound.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            cv = Cv.query.filter(Cv.smiles == nist.smiles).first()
            if cv is None:
                continue

            task = self.task_list[i]
            post_data = task.get_post_data(T, P)
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']
            dens = post_data['density']  # g/mL
            expan = post_data['expansion']  # /K
            compr = post_data['compress']  # /bar
            if None in [cp_inter, cp_pv, dens, expan, compr]:
                continue

            cp_sim = cp_inter + cp_pv + cv.get_post_data(T)
            cv_sim = cp_sim - T * nist.weight / dens * expan ** 2 / compr * 0.1  # J/mol.K
            sound_sim = (cp_sim / cv_sim / compr / dens) ** 0.5 * 10  # m/s

            sound_list.append([nist.smiles, value, uncertainty, sound_sim])

        print('')
        return sound_list

    def get_tc_nist(self):
        print('Get Tc from Nist')
        tc_list = []
        nist: NistMolecule
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            value, uncertainty = nist.tc, nist.tc_u
            if value is None:
                continue
            # if not nist.tc_has_exp:
            #     continue

            slab = self.task_list[i]
            tc = slab.get_post_data(100)['tc']

            tc_list.append([nist.smiles, value, uncertainty, tc])

        print('')
        return tc_list

    def get_dc_nist(self):
        print('Get Dc from Nist')
        dc_list = []
        nist: NistMolecule
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            value, uncertainty = nist.dc, nist.dc_u
            if value is None:
                continue
            # if not nist.dc_has_exp:
            #     continue

            slab = self.task_list[i]
            dc = slab.get_post_data(100)['dc']

            dc_list.append([nist.smiles, value / 1000, uncertainty / 1000, dc])

        print('')
        return dc_list

    def get_st_nist(self, _T=298):
        print('Get surface tension from Nist')
        st_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, _ = self.TP_list[_T][i]
            if None in (T,):
                continue

            # if nist.has_datas.filter(NistHasData.property_name == 'st-lg').filter(
            #         NistHasData.has_exp == True).count() == 0:

            spline = nist.splines.filter(NistSpline.property_id == prop_st.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            slab = self.task_list[i]
            st_sim = slab.get_post_data(T)['st']

            st_list.append([nist.smiles, value * 1000, uncertainty * 1000, st_sim])

        print('')
        return st_list

    def get_dgas_nist(self, _T=298):
        print('Get vapor density from Nist')
        dgas_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, _ = self.TP_list[_T][i]
            if None in (T,):
                continue

            # if nist.has_datas.filter(NistHasData.property_name == 'density-gl').filter(
            #         NistHasData.has_exp == True).count() == 0:
            #     continue

            spline = nist.splines.filter(NistSpline.property_id == prop_dgas.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            slab = self.task_list[i]
            try:
                # if T is larger than predicted Tc, dgas will not exist
                dgas_sim = slab.get_post_data(T)['dgas']
            except:
                continue

            dgas_list.append([nist.smiles, value / 1000, uncertainty / 1000, dgas_sim])

        print('')
        return dgas_list


def get_png_from_data(name, data_exp_sim_list, threthold=0.5):
    exp_list = []
    sim_list = []
    dev_list = []
    absdev_list = []
    u_list = []
    bad_list = []
    for smiles, ref, u, sim in data_exp_sim_list:
        if np.isnan(ref) or np.isnan(u) or np.isnan(sim):
            continue
        dev = (sim / ref - 1) * 100
        bad_list.append([smiles, '%.3g' % ref, '%.3g' % sim, round(dev, 1)])
        if threthold > 0 and abs(dev) > threthold * 100:  # incorrect expt. data
            continue
        dev_list.append(dev)
        absdev_list.append(abs(dev))
        exp_list.append(ref)
        sim_list.append(sim)
        u_list.append(u / ref * 100)

    bad_list.sort(key=lambda x: abs(x[-1]), reverse=True)

    p = BytesIO()
    fig = plt.figure(figsize=(13, 4))

    sp1 = fig.add_subplot(131)
    sp1.set_xlabel('Expt.')
    sp1.set_ylabel('Simu.')
    sp1.plot(exp_list, exp_list, '-')
    sp1.plot(exp_list, sim_list, '.', alpha=0.7)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    if len(exp_list)!=0:
        text = '%s\n%i molecules\nMDEV = %.1f %%\nMUD = %.1f %%' \
               % (name, len(exp_list), np.mean(dev_list), np.mean(absdev_list))
    else:
        text = '%s\n%i molecules\nMDEV = None %%\nMUD = None %%' \
               % (name, len(exp_list))
    sp1.text(0.05, 0.95, text, transform=sp1.transAxes, va='top', bbox=props)

    y, _ = np.histogram(absdev_list, bins=30, range=[0, 30])
    x = (_[1:] + _[:-1]) / 2
    sp2 = fig.add_subplot(132)
    sp2.set_xlabel('Unsigned deviation (%)')
    sp2.set_ylabel('Number of molecules')
    sp2.bar(x, y, color='C1', alpha=0.7)

    if len(exp_list) != 0:
        text = 'Unsigned deviation to expt. data\n%s\n%i molecules\nMUD = %.1f %%' \
           % (name, len(exp_list), np.mean(absdev_list))
    else:
        text = 'Unsigned deviation to expt. data\n%s\n%i molecules\nMUD = None %%' \
               % (name, len(exp_list))
    sp2.text(0.95, 0.95, text, transform=sp2.transAxes, va='top', ha='right', bbox=props)

    y, _ = np.histogram(u_list, bins=30, range=[0, 30])
    x = (_[1:] + _[:-1]) / 2
    sp3 = fig.add_subplot(133)
    sp3.set_xlabel('Uncertainty (%)')
    sp3.set_ylabel('Number of molecules')
    sp3.bar(x, y, alpha=0.7)

    if len(exp_list) != 0:
        text = 'Uncertainty of expt. data\n%s\n%i molecules\nMean uncertainty = %.1f %%' \
           % (name, len(exp_list), np.mean(u_list))
    else:
        text = 'Uncertainty of expt. data\n%s\n%i molecules\nMean uncertainty = None %%' \
               % (name, len(exp_list))
    sp3.text(0.95, 0.95, text, transform=sp3.transAxes, va='top', ha='right', bbox=props)

    plt.tight_layout()
    plt.savefig(p, format='png')
    plt.close()

    return p, bad_list


def write_plot(name, data):
    png, bad_list = get_png_from_data(name, data)
    name = name.replace(' ', '_')
    with open(f'_{name}.png', 'wb') as f:
        f.write(png.getvalue())
    with open(f'_{name}.txt', 'w') as f:
        for item in bad_list:
            f.write('\t'.join(map(str, item)) + '\n')


def compare_npt():

    action = StatAction()
    action.TP_list = {
        'Tm25': [],
        'Tvap': [],
        'Tc85': [],
    }
    action.update_mol_task_list(procedure='npt')
    print(len(action.nist_list))

    write_plot('density @ Tm+', action.get_property(property=prop_density, property_name='density', _T='Tm25'))
    write_plot('density @ Tvap', action.get_property(property=prop_density, property_name='density', _T='Tvap'))
    write_plot('density @ T0.85*', action.get_property(property=prop_density, property_name='density', _T='Tc85'))
    write_plot('Cp @ Tm+', action.get_property(property=prop_cp, property_name='cp', _T='Tm25'))
    write_plot('Cp @ Tvap', action.get_property(property=prop_cp, property_name='cp', _T='Tvap'))
    write_plot('Cp @ T0.85*', action.get_property(property=prop_cp, property_name='cp', _T='Tc85'))
    write_plot('Hvap @ Tm+', action.get_property(property=prop_hvap, property_name='hvap', _T='Tm25'))
    write_plot('Hvap @ Tvap', action.get_property(property=prop_hvap, property_name='hvap', _T='Tvap'))
    write_plot('Hvap @ T0.85*', action.get_property(property=prop_hvap, property_name='hvap', _T='Tc85'))
    # write_plot('Hliquid @ Tm+', action.get_hl_nist(_T='Tm25'))
    # write_plot('Hliquid @ Tvap', action.get_hl_nist(_T='Tvap'))
    # write_plot('Hliquid @ T0.85*', action.get_hl_nist(_T='Tc85'))


def compare_slab():

    action = StatAction()
    action.TP_list = {
        # 'Tm25': [],
        'Tvap': [],
        # 'Tc85': [],
    }
    action.update_mol_task_list(procedure='nvt-slab')
    print(len(action.nist_list))

    write_plot('surface tension @ Tvap', action.get_property(property=prop_st, property_name='st', _T='Tvap'))
    # write_plot('vapor pressure @ Tvap', action.get_property(property=prop_pvap, property_name='pzz', _T='Tvap'))
    write_plot('critical temperature', action.get_tc_nist())
    write_plot('critical density', action.get_dc_nist())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='This is a code to compare simulation results with experimental \
        results. For one property, a png file will generated show the detailed comparison between \
        simulation and experiment')
    parser.add_argument('-p', '--procedure', type=str, help='procedure of the compute: npt(ppm), or nvt-slab')
    parser.add_argument('--selection', type=bool, help='select specific task', default=False)
    opt = parser.parse_args()
    procedure = opt.procedure
    app = create_app(procedure)
    app.app_context().push()

    prop_pvap = NistProperty.query.filter(NistProperty.name == 'pvap-lg').first()
    prop_density = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
    prop_hvap = NistProperty.query.filter(NistProperty.name == 'hvap-lg').first()
    prop_cp = NistProperty.query.filter(NistProperty.name == 'cp-lg').first()
    prop_sound = NistProperty.query.filter(NistProperty.name == 'sound-lg').first()
    prop_st = NistProperty.query.filter(NistProperty.name == 'st-lg').first()
    prop_dgas = NistProperty.query.filter(NistProperty.name == 'density-gl').first()
    prop_hl = NistProperty.query.filter(NistProperty.name == 'h-lg').first()

    if procedure == 'npt':
        compare_npt()
    elif procedure == 'nvt-slab':
        compare_slab()
    else:
        print('Unknown procedure')
