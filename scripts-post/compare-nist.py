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


class StatAction():
    def __init__(self):
        self.nist_list = []
        self.task_list = []
        self.slab_list = []
        self.TP_list = {
            # 298: [],
            # 'Tm25': [],
            # 'Tvap': [],
            # 'Tcx8': [],
            # 'Tc85': [],
        }

    def update_mol_task_list(self, smiles_list):
        print('Get molecule list')

        for smiles in smiles_list:
            nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if nist is None:
                continue

            task = Task.query.filter(Task.smiles_list == json.dumps([smiles])).first()
            if task is None or task.post_result is None or task.get_post_result()['density-poly4'][-1] < 0.999:
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

    def update_mol_slab_list(self, smiles_list):
        print('Get molecule list')

        for smiles in smiles_list:
            nist = NistMolecule.query.filter(NistMolecule.smiles == smiles).first()
            if nist is None:
                continue

            slab = Task.query.filter(Task.smiles_list == json.dumps([smiles])).first()
            if slab is None or slab.post_result is None:
                continue

            self.nist_list.append(nist)
            self.slab_list.append(slab)

            for _T in self.TP_list.keys():
                T = self.get_T_nist(nist, _T)
                P = None

                self.TP_list[_T].append((T, P))

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
        if nist is None:
            return None

        spline_pvap = nist.splines.filter(NistSpline.property_id == prop_pvap.id).first()
        if spline_pvap is None:
            P = None
        else:
            P = spline_pvap.get_data(T)[0]  # kPa
            if P is not None:
                P /= 100  # bar
        return P

    def get_density_nist(self, _T=298):
        print('Get density from nist')
        dens_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.filter(NistSpline.property_id == prop_density.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            task = self.task_list[i]
            post_data = task.get_post_data(T, P)
            if post_data.get('error') is not None:
                continue
            dens_sim = post_data['density']

            dens_list.append([nist.smiles, value / 1000, uncertainty / 1000, dens_sim])

        print('')
        return dens_list

    def get_cp_nist(self, _T=298):
        print('Get Cp from nist')
        cp_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.filter(NistSpline.property_id == prop_cp.id).first()
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
            if post_data.get('error') is not None:
                continue
            cp_inter = post_data['cp_inter']
            cp_pv = post_data['cp_pv']
            cp_sim = cp_inter + cp_pv + cv.get_post_data(T)

            cp_list.append([nist.smiles, value, uncertainty, cp_sim])

        print('')
        return cp_list

    def get_hvap_nist(self, _T=298):
        print('Get Hvap from nist')
        hvap_list = []
        for i, nist in enumerate(self.nist_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, nist.smiles))

            T, P = self.TP_list[_T][i]
            if None in (T, P):
                continue

            spline = nist.splines.filter(NistSpline.property_id == prop_hvap.id).first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            task = self.task_list[i]
            post_data = task.get_post_data(T, P)
            if post_data.get('error') is not None:
                continue
            hvap_sim = post_data['hvap']

            ### Correction for hvap
            hvap_sim = hvap_sim - nist.n_nothx / 15 * (8.314 * T) / 1000

            hvap_list.append([nist.smiles, value, uncertainty, hvap_sim])

        print('')
        return hvap_list

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
            cp_sim = cp_inter + cp_pv + cv.get_post_data(T)

            dens = post_data['density']  # g/mL
            expan = post_data['expansion']  # /K
            compr = post_data['compress']  # /bar

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

            slab = self.slab_list[i]
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

            slab = self.slab_list[i]
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

            slab = self.slab_list[i]
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

            spline = nist.splines.join(NistProperty).filter(NistProperty.name == 'density-gl').first()
            if spline is None:
                continue
            value, uncertainty = spline.get_data(T)
            if value is None:
                continue

            slab = self.slab_list[i]
            dgas_sim = slab.get_post_data(T)['density_gas']

            dgas_list.append([nist.smiles, value / 1000, uncertainty / 1000, dgas_sim])

        print('')
        return dgas_list


def get_png_from_data(name, data_exp_sim_list):
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
        if abs(dev) > 50:  # incorrect expt. data
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
    text = '%s\n%i molecules\nMDEV = %.1f %%\nMUD = %.1f %%' \
           % (name, len(exp_list), np.mean(dev_list), np.mean(absdev_list))
    sp1.text(0.05, 0.95, text, transform=sp1.transAxes, va='top', bbox=props)

    y, _ = np.histogram(absdev_list, bins=30, range=[0, 30])
    x = (_[1:] + _[:-1]) / 2
    sp2 = fig.add_subplot(132)
    sp2.set_xlabel('Unsigned deviation (%)')
    sp2.set_ylabel('Number of molecules')
    sp2.bar(x, y, color='C1', alpha=0.7)

    text = 'Unsigned deviation to expt. data\n%s\n%i molecules\nMUD = %.1f %%' \
           % (name, len(exp_list), np.mean(absdev_list))
    sp2.text(0.95, 0.95, text, transform=sp2.transAxes, va='top', ha='right', bbox=props)

    y, _ = np.histogram(u_list, bins=30, range=[0, 30])
    x = (_[1:] + _[:-1]) / 2
    sp3 = fig.add_subplot(133)
    sp3.set_xlabel('Uncertainty (%)')
    sp3.set_ylabel('Number of molecules')
    sp3.bar(x, y, alpha=0.7)

    text = 'Uncertainty of expt. data\n%s\n%i molecules\nMean uncertainty = %.1f %%' \
           % (name, len(exp_list), np.mean(u_list))
    sp3.text(0.95, 0.95, text, transform=sp3.transAxes, va='top', ha='right', bbox=props)

    plt.tight_layout()
    plt.savefig(p, format='png')
    plt.close()

    return p, bad_list


def compare_npt():
    with open(sys.argv[2]) as f:
        smiles_list = f.read().splitlines()

    action = StatAction()
    action.TP_list = {
        # 'Tm25': [],
        'Tvap': [],
        # 'Tc85': [],
    }
    action.update_mol_task_list(smiles_list)
    print(len(action.nist_list))

    k_data_exp_sim_list = [
        # ('density @ Tm+', action.get_density_nist(_T='Tm25')),
        ('density @ Tvap', action.get_density_nist(_T='Tvap')),
        # ('density @ T0.85*', action.get_density_nist(_T='Tc85'),
        # ('Cp @ Tm+', action.get_cp_nist(_T='Tm25')),
        ('Cp @ Tvap', action.get_cp_nist(_T='Tvap')),
        # ('Cp @ T0.85*', action.get_cp_nist(_T='Tc85')),
        # ('Hvap @ Tm+', action.get_hvap_nist(_T='Tm25')),
        ('Hvap @ Tvap', action.get_hvap_nist(_T='Tvap')),
        # ('Hvap @ T0.85*', action.get_hvap_nist(_T='Tc85')),
    ]
    for k, v in k_data_exp_sim_list:
        png, bad_list = get_png_from_data(k, v)
        name = k.replace(' ', '_')
        with open(f'_{name}.png', 'wb') as f:
            f.write(png.getvalue())
        with open(f'_{name}.txt', 'w') as f:
            for item in bad_list:
                f.write('\t'.join(map(str, item)) + '\n')


def compare_slab():
    with open(sys.argv[2]) as f:
        smiles_list = f.read().splitlines()

    action = StatAction()
    action.TP_list = {
        # 'Tm25': [],
        'Tvap': [],
        # 'Tc85': [],
    }
    action.update_mol_slab_list(smiles_list)
    print(len(action.nist_list))

    k_data_exp_sim_list = [
        ('surface tension @ Tvap', action.get_st_nist(_T='Tvap')),
        ('critical temperature', action.get_tc_nist()),
        ('critical density', action.get_dc_nist()),
    ]
    for k, v in k_data_exp_sim_list:
        png, bad_list = get_png_from_data(k, v)
        name = k.replace(' ', '_')
        with open(f'_{name}.png', 'wb') as f:
            f.write(png.getvalue())
        with open(f'_{name}.txt', 'w') as f:
            for item in bad_list:
                f.write('\t'.join(map(str, item)) + '\n')


if __name__ == '__main__':
    procedure = sys.argv[1]
    app = create_app(procedure)
    app.app_context().push()

    prop_pvap = NistProperty.query.filter(NistProperty.name == 'pvap-lg').first()
    prop_density = NistProperty.query.filter(NistProperty.name == 'density-lg').first()
    prop_hvap = NistProperty.query.filter(NistProperty.name == 'hvap-lg').first()
    prop_cp = NistProperty.query.filter(NistProperty.name == 'cp-lg').first()
    prop_sound = NistProperty.query.filter(NistProperty.name == 'sound-lg').first()
    prop_st = NistProperty.query.filter(NistProperty.name == 'st-lg').first()
    prop_dgas = NistProperty.query.filter(NistProperty.name == 'density-gl').first()

    if procedure == 'npt':
        compare_npt()
    elif procedure == 'nvt-slab':
        compare_slab()
    else:
        print('Unknown procedure')
