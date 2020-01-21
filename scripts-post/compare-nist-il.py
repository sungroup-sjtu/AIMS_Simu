#!/usr/bin/env python3
# coding=utf-8

import os, sys, time
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO

sys.path.append('..')
from app import create_app
from app.models import *
from app.models_ilthermo import Property, Ion, Molecule, Spline, Data
from app.models_cv import Cv
from app.selection import task_selection
from mstools.smiles.smiles import *


class StatAction:
    def __init__(self):
        self.mol_list = []
        self.task_list = []
        self.Nslice = 5 # 6 temperature from t_min to t_max will be chosen in comparison

    def update_mol_task_list(self, procedure):
        print('Get molecule list')

        tasks = Task.query.filter(Task.procedure == procedure).filter(Task.status.in_((Compute.Status.ANALYZED, Compute.Status.DONE)))
        for task in tasks:
            print(task)
            if task_selection(task, select=opt.selection):
                smiles_list = json.loads(task.smiles_list)
                if get_charge(smiles_list[0]) < 0 or get_charge(smiles_list[1]) > 0:
                    raise Exception('the input of smiles must in cation.anion order')
                cation_id = Ion.query.filter(Ion.smiles == smiles_list[0]).first().id
                anion_id = Ion.query.filter(Ion.smiles == smiles_list[1]).first().id
                mols = Molecule.query.filter(Molecule.cation_id == cation_id).filter(Molecule.anion_id == anion_id)
                if mols.count() > 1:
                    raise Exception('molecule %s exists %i times in Molecule database' % ('-'.join(smiles_list), mols.count()))
                mol = mols.first()
                if mol is None:
                    continue

                if task.post_result is None or task.status == 1:
                    continue

                self.mol_list.append(mol)
                self.task_list.append(task)

    # T = '2/5' or 300K
    def get_property(self, property, property_name, T=None, P=1, VTF=True):
        print('Get %s from nist-ilthermo' % (property.name))
        value_list = []
        for i, mol in enumerate(self.mol_list):
            if i % 10 == 0:
                sys.stdout.write('\r\t%5i %s' % (i, mol.smiles()))

            # pressure check, only condider 1 atm pressure
            task = self.task_list[i]
            p_list = json.loads(task.p_list)
            if P not in p_list:
                continue

            # get t_min and t_max
            t_list = json.loads(task.t_list)
            t_min_list, t_max_list = [min(t_list)], [max(t_list)]
            splines = mol.splines.filter(Spline.property_id == property.id)
            if splines.count() == 0:
                continue
            for spline in splines:
                t_min_list.append(spline.t_min)
                t_max_list.append(spline.t_max)
            t_min = max(t_min_list)
            t_max = min(t_max_list)
            # get comparison temperature
            if type(T) == str:
                a, b = list(map(int, T.split('/')))
                T = a / b * (t_max - t_min) + t_min
            elif type(T) == int:
                if not t_min-1 < T < t_max+1:
                    continue

            # get simulation value
            post_data = task.get_post_data(T, P)
            if property_name == 'cp':
                cp_inter = post_data['cp_inter']
                cp_pv = post_data['cp_pv']
                if cp_inter is None or cp_pv is None:
                    continue
                cation_smiles = mol.cation.smiles
                anion_smiles = mol.anion.smiles
                cv_c = Cv.query.filter(Cv.smiles == cation_smiles).first()
                cv_a = Cv.query.filter(Cv.smiles == anion_smiles).first()
                if cv_c is None or cv_a is None:
                    continue
                value_sim = cp_inter + cp_pv + cv_c.get_post_cv(T) + cv_a.get_post_cv(T)
            else:
                value_sim = post_data.get(property_name)
            if value_sim is None:
                continue

            re_min = 100.
            # get experimental value
            for spline in splines:
                v, u = spline.get_data(T)
                if v is None or u is None:
                    continue
                if u < 0:
                    continue
                if property in [prop_vis, prop_diff, prop_econ] and VTF:
                    v = spline.get_VTF_data(T)
                # unit transform
                if property == prop_density:
                    v /= 1000
                    u /= 1000
                elif property == prop_vis:
                    v *= 1000
                    u *= 1000
                elif property == prop_diff:
                    value_sim *= 1e5
                relative_error = (value_sim - v) / v
                if relative_error < re_min:
                    re_min = relative_error
                    value = v
                    uncertainty = u

            value_list.append([mol.smiles(), value, uncertainty, value_sim])
        return value_list


def get_png_from_data(name, data_exp_sim_list, threthold=100, error_range=30):
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

    y, _ = np.histogram(absdev_list, bins=error_range, range=[0, error_range])
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

    y, _ = np.histogram(u_list, bins=error_range, range=[0, error_range])
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
    action.update_mol_task_list('npt')
    print(len(action.mol_list))

    write_plot('density-T=0', action.get_property(property=prop_density, property_name='density', T='0/2'))
    write_plot('density-T=1', action.get_property(property=prop_density, property_name='density', T='1/2'))
    write_plot('density-T=2', action.get_property(property=prop_density, property_name='density', T='2/2'))
    write_plot('density-T=300', action.get_property(property=prop_density, property_name='density', T=300))
    write_plot('density-T=350', action.get_property(property=prop_density, property_name='density', T=350))
    write_plot('density-T=400', action.get_property(property=prop_density, property_name='density', T=400))
    write_plot('density-T=450', action.get_property(property=prop_density, property_name='density', T=450))

    write_plot('cp-T=0', action.get_property(property=prop_cp, property_name='cp', T='0/2'))
    write_plot('cp-T=1', action.get_property(property=prop_cp, property_name='cp', T='1/2'))
    write_plot('cp-T=2', action.get_property(property=prop_cp, property_name='cp', T='2/2'))
    write_plot('cp-T=300', action.get_property(property=prop_cp, property_name='cp', T=300))
    write_plot('cp-T=350', action.get_property(property=prop_cp, property_name='cp', T=350))
    write_plot('cp-T=400', action.get_property(property=prop_cp, property_name='cp', T=400))
    write_plot('cp-T=450', action.get_property(property=prop_cp, property_name='cp', T=450))


def compare_ppm():
    action = StatAction()
    action.update_mol_task_list('ppm')
    print(len(action.mol_list))

    write_plot('vis_ppm-T=0', action.get_property(property=prop_vis, property_name='viscosity', T='0/2'))
    write_plot('vis_ppm-T=1', action.get_property(property=prop_vis, property_name='viscosity', T='1/2'))
    write_plot('vis_ppm-T=2', action.get_property(property=prop_vis, property_name='viscosity', T='2/2'))
    write_plot('vis_ppm-T=300', action.get_property(property=prop_vis, property_name='viscosity', T=300))
    write_plot('vis_ppm-T=350', action.get_property(property=prop_vis, property_name='viscosity', T=350))
    write_plot('vis_ppm-T=400', action.get_property(property=prop_vis, property_name='viscosity', T=400))
    write_plot('vis_ppm-T=450', action.get_property(property=prop_vis, property_name='viscosity', T=450))


def compare_nvt_multi():
    action = StatAction()
    action.update_mol_task_list('nvt-multi')
    print(len(action.mol_list))

    write_plot('econ_gk-T=0', action.get_property(property=prop_econ, property_name='electrical conductivity', T='0/2'))
    write_plot('econ_gk-T=1', action.get_property(property=prop_econ, property_name='electrical conductivity', T='1/2'))
    write_plot('econ_gk-T=2', action.get_property(property=prop_econ, property_name='electrical conductivity', T='2/2'))
    write_plot('econ_gk-T=300', action.get_property(property=prop_econ, property_name='electrical conductivity', T=300))
    write_plot('econ_gk-T=350', action.get_property(property=prop_econ, property_name='electrical conductivity', T=350))
    write_plot('econ_gk-T=400', action.get_property(property=prop_econ, property_name='electrical conductivity', T=400))
    write_plot('econ_gk-T=450', action.get_property(property=prop_econ, property_name='electrical conductivity', T=450))

    write_plot('econ_ne-T=0', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T='0/2'))
    write_plot('econ_ne-T=1', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T='1/2'))
    write_plot('econ_ne-T=2', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T='2/2'))
    write_plot('econ_ne-T=300', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T=300))
    write_plot('econ_ne-T=350', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T=350))
    write_plot('econ_ne-T=400', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T=400))
    write_plot('econ_ne-T=450', action.get_property(property=prop_econ, property_name='Nernst-Einstein electrical conductivity', T=450))

    write_plot('vis_gk-T=0', action.get_property(property=prop_vis, property_name='viscosity', T='0/2'))
    write_plot('vis_gk-T=1', action.get_property(property=prop_vis, property_name='viscosity', T='1/2'))
    write_plot('vis_gk-T=2', action.get_property(property=prop_vis, property_name='viscosity', T='2/2'))
    write_plot('vis_gk-T=300', action.get_property(property=prop_vis, property_name='viscosity', T=300))
    write_plot('vis_gk-T=350', action.get_property(property=prop_vis, property_name='viscosity', T=350))
    write_plot('vis_gk-T=400', action.get_property(property=prop_vis, property_name='viscosity', T=400))
    write_plot('vis_gk-T=450', action.get_property(property=prop_vis, property_name='viscosity', T=450))

    write_plot('diff-T=0', action.get_property(property=prop_diff, property_name='diffusion constant sum', T='0/2'))
    write_plot('diff-T=1', action.get_property(property=prop_diff, property_name='diffusion constant sum', T='1/2'))
    write_plot('diff-T=2', action.get_property(property=prop_diff, property_name='diffusion constant sum', T='2/2'))
    write_plot('diff-T=300', action.get_property(property=prop_diff, property_name='diffusion constant sum', T=300))
    write_plot('diff-T=350', action.get_property(property=prop_diff, property_name='diffusion constant sum', T=350))
    write_plot('diff-T=400', action.get_property(property=prop_diff, property_name='diffusion constant sum', T=400))
    write_plot('diff-T=450', action.get_property(property=prop_diff, property_name='diffusion constant sum', T=450))


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

    prop_density = Property.query.filter(Property.name == 'Density').first()
    prop_cp = Property.query.filter(Property.name == 'Heat capacity at constant pressure').first()
    prop_econ = Property.query.filter(Property.name == 'Electrical conductivity').first()
    prop_vis = Property.query.filter(Property.name == 'Viscosity').first()
    prop_diff = Property.query.filter(Property.name == 'Self-diffusion coefficient').first()

    if procedure == 'npt':
        compare_npt()
    elif procedure == 'ppm':
        compare_ppm()
    elif procedure == 'nvt-multi':
        compare_nvt_multi()
    else:
        print('Unknown procedure')
