#!/usr/bin/env python3
# coding=utf-8


import sys
import os
import math
from collections import OrderedDict
from typing import Dict

import pybel
import shutil
from scipy.optimize import root
import numpy as np
import copy
from lmfit import Parameters, Minimizer

sys.path.append('..')

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd
from mstools.simulation.gmx import GmxSimulation, Npt
from config import Config
from app import jobmanager

jobmanager.queue = 'fast'
jobmanager.nprocs = 6
kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
          'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
simulation = GmxSimulation(**kwargs)


class Target():
    def __int__(self):
        self.id = 0
        self.n_mol = 0
        self.smiles = None
        self.T = None
        self.P = None
        self.density = None
        self.wDensity = None
        self.hvap = None
        self.wHvap = None

    def __repr__(self):
        return '%s,%i' % (self.smiles, self.T)

    @property
    def dir(self):
        return os.path.join(Config.WORK_DIR, '%i-%i-%i' % (self.id, self.T, self.P))

    def n_mol(self, n_atoms=3000):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        return math.ceil(n_atoms / len(py_mol.atoms))

    def build(self, n_atoms=3000, ppf=None):
        os.chdir(self.dir)
        pdb = 'mol.pdb'
        mol2 = 'mol.mol2'
        py_mol = create_mol_from_smiles(self.smiles, pdb_out=pdb, mol2_out=mol2)
        n_mol = math.ceil(n_atoms / len(py_mol.atoms))
        mass = py_mol.molwt * n_mol
        length = (10 / 6.022 * mass / self.density) ** (1 / 3)  # assume cubic box

        print('Build coordinates using Packmol: %s molecules ...' % n_mol)
        simulation.packmol.build_box([pdb], [n_mol], 'init.pdb', length=length - 2, silent=True)

        print('Create box using DFF ...')
        simulation.dff.build_box_after_packmol([mol2], [n_mol], 'init.msd', mol_corr='init.pdb', length=length)
        if ppf is not None:
            shutil.copy(ppf, 'TEAM_LS.ppf')
        simulation.export(ppf='TEAM_LS.ppf', minimize=True)

        nprocs = simulation.jobmanager.nprocs
        commands = []
        simulation.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-eq.mdp', T=self.T, nsteps=int(4E5), nstxtcout=0)
        cmd = simulation.gmx.grompp(mdp='grompp-eq.mdp', tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-nvt.mdp', T=self.T, nsteps=int(2E5),
                                                 nstxout=1000, nstvout=1000, nstxtcout=0,
                                                 restart=True)
        cmd = simulation.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', tpr_out='nvt.tpr', cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        top_hvap = 'topol-hvap.top'
        simulation.gmx.generate_top_for_hvap('topol.top', top_hvap)
        cmd = simulation.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='hvap', nprocs=nprocs, rerun='nvt.trr', get_cmd=True)
        commands.append(cmd)

        simulation.jobmanager.generate_sh(os.getcwd(), commands, name='%i-%i' % (self.id, self.T))
        simulation.run()

    def get_pres_hvap_from_paras(self, ppf_file=None, d: OrderedDict = None) -> (float, float):
        paras = copy.copy(d)
        os.chdir(self.dir)
        if ppf_file != None:
            ppf = PPF(ppf_file)
        else:
            ppf = PPF('TEAM_LS.ppf')
        if not ppf.set_lj_para(paras):
            pres = simulation.gmx.get_property('nvt.edr', 'Pressure', begin=10)
            ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=10)
            hvap = 8.314 * self.T / 1000 - ei / self.n_mol
            return pres, hvap

        ppf.write('rerun.ppf')
        top = 'diff.top'

        shutil.copy('init.msd', 'diff.msd')
        simulation.dff.set_charge(['diff.msd'], 'rerun.ppf')
        simulation.dff.export_gmx('diff.msd', 'rerun.ppf', top_out=top)
        nprocs = simulation.jobmanager.nprocs

        # Press
        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), nstxtcout=0)
        simulation.gmx.grompp(top=top, tpr_out='diff.tpr', silent=True)
        simulation.gmx.mdrun(name='diff', nprocs=nprocs, rerun='nvt.trr', silent=True)

        pres = simulation.gmx.get_property('diff.edr', 'Pressure', begin=10)

        # Hvap
        top_hvap = 'diff-hvap.top'
        simulation.gmx.generate_top_for_hvap(top, top_hvap)
        nprocs = simulation.jobmanager.nprocs

        simulation.gmx.grompp(top=top_hvap, tpr_out='diff-hvap.tpr', silent=True)
        simulation.gmx.mdrun(name='diff-hvap', nprocs=nprocs, rerun='nvt.trr', silent=True)

        ei = simulation.gmx.get_property('diff-hvap.edr', 'Potential', begin=10)
        hvap = 8.314 * self.T / 1000 - ei / self.n_mol
        return pres, hvap

    def get_dPres_dHvap_from_paras(self, ppf_file, d: OrderedDict) -> ([float], [float]):
        paras = copy.copy(d)
        os.chdir(self.dir)
        dPres = []
        dHvap = []
        for k, v in paras.items():
            dP, dH = self.get_dPres_dHvap_for_para(ppf_file, paras, k)
            dPres.append(dP)
            dHvap.append(dH)
        return dPres, dHvap

    def get_dPres_dHvap_for_para(self, ppf_file, d: OrderedDict, k) -> (float, float):
        paras = copy.copy(d)
        os.chdir(self.dir)
        v = paras[k]
        pres_list = []
        hvap_list = []

        h = 0.02
        if k.endswith('r0'):
            h = 0.02
        elif k.endswith('e0'):
            h = 0.01

        for i in [-1, 1]:
            new_v = v + i * h
            paras[k] = new_v
            if ppf_file is not None:
                ppf = PPF(ppf_file)
            else:
                ppf = PPF('TEAM_LS.ppf')
            if not ppf.set_lj_para(paras):
                return 0, 0
            ppf.write('rerun.ppf')
            #print('    %s %10.5f %10.5f' %(k, v, new_v))
            top = 'diff%i.top' % i

            shutil.copy('init.msd', 'diff.msd')
            simulation.dff.set_charge(['diff.msd'], 'rerun.ppf')
            simulation.dff.export_gmx('diff.msd', 'rerun.ppf', top_out=top)
            nprocs = simulation.jobmanager.nprocs

            # Pres
            simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), nstxtcout=0)
            simulation.gmx.grompp(top=top, tpr_out='diff%i.tpr' % i, silent=True)
            simulation.gmx.mdrun(name='diff%i' % i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            pres = simulation.gmx.get_property('diff%i.edr' % i, 'Pressure', begin=10)
            pres_list.append(pres)

            # HVap
            top_hvap = 'diff%i-hvap.top'
            simulation.gmx.generate_top_for_hvap(top, top_hvap)

            simulation.gmx.grompp(top=top_hvap, tpr_out='diff%i-hvap.tpr' % i, silent=True)
            simulation.gmx.mdrun(name='diff%i-hvap' % i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            ei = simulation.gmx.get_property('diff%i-hvap.edr' % i, 'Potential', begin=10)
            hvap = 8.314 * self.T / 1000 - ei / self.n_mol
            hvap_list.append(hvap)

        dPres = (pres_list[1] - pres_list[0]) / 2 / h
        dHvap = (hvap_list[1] - hvap_list[0]) / 2 / h
        return dPres, dHvap

    @property
    def dir_npt(self):
        return os.path.join(Config.WORK_DIR, 'NPT-%i' % (self.id))

    def run_npt(self, ppf=None, n_atoms=3000):
        cd_or_create_and_cd(self.dir_npt)

        if not os.path.exists('init.msd'):
            pdb = 'mol.pdb'
            mol2 = 'mol.mol2'
            py_mol = create_mol_from_smiles(self.smiles, pdb_out=pdb, mol2_out=mol2)
            n_mol = math.ceil(n_atoms / len(py_mol.atoms))
            mass = py_mol.molwt * n_mol
            length = (10 / 6.022 * mass / (self.density-0.1)) ** (1 / 3)  # assume cubic box

            print('Build coordinates using Packmol: %s molecules ...' % n_mol)
            simulation.packmol.build_box([pdb], [n_mol], 'init.pdb', length=length - 2, silent=True)

            print('Create box using DFF ...')
            simulation.dff.build_box_after_packmol([mol2], [n_mol], 'init.msd', mol_corr='init.pdb', length=length)

        cd_or_create_and_cd(os.path.basename('%s' %ppf)[:-4])

        if not os.path.exists('conf.gro'):
            if ppf is not None:
                shutil.copy(ppf, 'TEAM_LS.ppf')
            simulation.msd = '../init.msd'
            simulation.export(ppf='TEAM_LS.ppf', minimize=True)

        cd_or_create_and_cd('%i-%i' %(self.T, self.P))

        npt = Npt(**kwargs)
        npt.prepare(model_dir='..', T=self.T, P=self.P, jobname=str(self))
        npt.run()


class PPF():
    def __init__(self, ppf):
        with open(ppf) as f:
            self.terms = f.read().splitlines()

    @property
    def adj_lj_paras(self) -> OrderedDict:
        terms = OrderedDict()
        for term in self.terms:
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                continue

            words = term.split(':')
            words = [w.strip() for w in words]
            if term.startswith('N12_6'):
                a_type = words[1]
                paras = words[2]

                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]
                if not r0.endswith('*'):
                    terms['%s_r0' % a_type] = float(r0)
                if not e0.endswith('*'):
                    terms['%s_e0' % a_type] = float(e0)
            elif term.startswith('BINC'):
                a_types = words[1]
                para = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]
                if not para.endswith('*'):
                    terms['%s_%s_bi' %(a1_type, a2_type)] = float(para)
        return terms

    def set_lj_para(self, new_paras: Dict):
        terms = [''] * len(self.terms)
        replaced = False
        for i, term in enumerate(self.terms):
            if not (term.startswith('N12_6') or term.startswith('BINC')):
                terms[i] = term
                continue

            if term.startswith('N12_6'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_type = words[1]
                paras = words[2]
                words = paras.split(',')
                words = [w.strip() for w in words]
                r0 = words[0]
                e0 = words[1]

                r0_key = '%s_r0' % a_type
                if r0_key in new_paras.keys():
                    r0 = new_paras[r0_key]
                    replaced = True
                e0_key = '%s_e0' % a_type
                if e0_key in new_paras.keys():
                    e0 = new_paras[e0_key]
                    replaced = True

                new_term = 'N12_6: %s: %s, %s:' % (a_type, str(r0), str(e0))
                terms[i] = new_term
            elif term.startswith('BINC'):
                words = term.split(':')
                words = [w.strip() for w in words]
                a_types = words[1]
                bi = words[2]

                words = a_types.split(',')
                words = [w.strip() for w in words]
                a1_type = words[0]
                a2_type = words[1]

                bi_key = '%s_%s_bi' %(a1_type, a2_type)
                if bi_key in new_paras.keys():
                    bi = new_paras[bi_key]
                    replaced = True

                new_term = 'BINC: %s, %s: %s:' %(a1_type, a2_type, str(bi))
                terms[i] = new_term

        self.terms = terms
        return replaced

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            for term in self.terms:
                f.write(term + '\n')



def build(ppf_file):
    for p in targets:
        print(p.dir)
        cd_or_create_and_cd(p.dir)
        p.build(ppf=ppf_file)


def optimize(ppf_file):
    def residual(params: Parameters):
        paras = OrderedDict()
        for k, v in params.items():
            paras[k] = v.value

        f = []
        for target in targets:
            pres, hvap = target.get_pres_hvap_from_paras(ppf_file, paras)
            f.append((pres - target.P) * target.wDensity)
            f.append((hvap - target.hvap) * target.wHvap)
        #print('CURRENT', params)
        #print('RESIDUE:', list(map(lambda x: round(x, 1), f)))
        #print('RSQ:', round(np.sum(list(map(lambda x: x ** 2, f))), 1))

        return f

    def dfunc(params: Parameters):
        paras = OrderedDict()
        for k, v in params.items():
            paras[k] = v.value

        df = []
        for target in targets:
            dPres, dHvap = target.get_dPres_dHvap_from_paras(ppf_file, paras)
            df.append([ i * target.wDensity for i in dPres ])
            df.append([ i * target.wHvap for i in dHvap ])
        print('\nDIRECTIVE:')
        for l in df:
            print(list(map(lambda x: round(x, 1), l)))

        return df

    def print_result(params:Parameters, iter:int, resid:[float]):
        print('\nITERATION:', iter)
        print('PARAMETERS:')
        for k, v in params.items():
            print('\t', k, round(v.value, 5))
        print('RESIDUE:', list(map(lambda x: round(x, 1), resid)))
        print('RSQ:', round(np.sum(list(map(lambda x: x ** 2, resid))), 1))
        print('')

    params = Parameters()
    for k,v in adj_lj_paras.items():
        if k.endswith('r0'):
            params.add(k, value=v, min=2, max=5)
        elif k.endswith('e0'):
            params.add(k, value=v, min=0.01, max=3)
        elif k.endswith('bi'):
            params.add(k, value=v, min=-1, max=1)

    params_init = copy.deepcopy(params)

    mini = Minimizer(residual, params, iter_cb=print_result)
    result = mini.leastsq(Dfun=dfunc)

    print('INIT PARAMETERS:')
    for k, v in params_init.items():
        print('\t', k, round(v.value, 5))
    res_init = residual(params_init)
    print('INIT RESIDUE:', list(map(lambda x: round(x, 1), res_init)))
    print('INIT RSQ:', round(np.sum(list(map(lambda x: x ** 2, res_init))), 1))

def run_npt(ppf_file):
    for target in targets:
        target.run_npt(ppf_file)

def plot(ppfs):
    ppfs = [ppf[:-4] for ppf in ppfs]
    props = dict()
    for target in targets:
        if not target.id in props.keys():
            props[target.id] = {'smiles': target.smiles, 'T':[], 'd_exp':[], 'h_exp': [],
                                'd_sim': OrderedDict(), 'h_sim': OrderedDict()}
        props[target.id]['T'].append(target.T)
        props[target.id]['d_exp'].append(target.density)
        props[target.id]['h_exp'].append(target.hvap)

        for ppf in ppfs:
            if ppf not in props[target.id]['d_sim'].keys():
                props[target.id]['d_sim'][ppf] = []
                props[target.id]['h_sim'][ppf] = []

            os.chdir(target.dir_npt)
            os.chdir(ppf)
            os.chdir('%i-%i' %(target.T, target.P))
            density = simulation.gmx.get_property('npt.edr', 'Density', begin=250)
            density /= 1000
            ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=250)
            hvap = 8.314 * target.T / 1000 - ei / target.n_mol

            props[target.id]['d_sim'][ppf].append(density)
            props[target.id]['h_sim'][ppf].append(hvap)

    os.chdir(sys.path[0])
    import pylab
    pylab.rcParams.update({'font.size': 16})
    for tid, prop in props.items():
        pylab.figure()
        pylab.plot(prop['T'], prop['d_exp'], '--')
        for ppf, points in prop['d_sim'].items():
            marker = 'x' if ppf.endswith('init') else 'o'
            pylab.plot(prop['T'], points, marker, label=ppf)
        y_mean = np.mean(prop['d_exp'])
        pylab.ylim(y_mean-0.2, y_mean+0.2)
        pylab.legend()
        pylab.title('Density %i %s (g/mL)' %(tid, prop['smiles']))
        pylab.savefig('density-%02i.png' %tid)

        pylab.figure()
        pylab.plot(prop['T'], prop['h_exp'], '--')
        for ppf, points in prop['h_sim'].items():
            marker = 'x' if ppf.endswith('init') else 'o'
            pylab.plot(prop['T'], points, marker, label=ppf)
        y_mean = np.mean(prop['h_exp'])
        pylab.ylim(y_mean-10, y_mean+10)
        pylab.legend()
        pylab.title('HVap %i %s (kJ/mol)' %(tid, prop['smiles']))
        pylab.savefig('hvap-%02i.png' %tid)


if __name__ == '__main__':
    def read_data(filename):
        targets = []
        with open(filename) as f:
            lines = f.read().splitlines()
        smiles_list = []
        for line in lines:
            if line.startswith('#') or line.strip() == '':
                continue
            words = line.split()

            smiles = words[0]
            if smiles not in smiles_list:
                smiles_list.append(smiles)

            target = Target()
            target.id = int(words[0])
            target.smiles = words[1]
            target.T = int(words[2])
            target.P = int(words[3])
            target.density = float(words[4])
            target.wDensity = float(words[5])
            target.hvap = float(words[6])
            target.wHvap = float(words[7])
            target.n_mol = target.n_mol()
            targets.append(target)
        return targets

    datafile = sys.argv[1]
    global targets
    targets = read_data(datafile)
    cmd = sys.argv[2]
    ppf_file = None
    if len(sys.argv) == 4:
        ppf_file = os.path.abspath(sys.argv[3])

    if cmd == 'build':
        build(ppf_file)
    if cmd == 'optimize':
        ppf = PPF(ppf_file)
        global adj_lj_paras
        adj_lj_paras = ppf.adj_lj_paras
        jobmanager.nprocs = 8
        optimize(ppf_file)
    if cmd == 'npt':
        run_npt(ppf_file)
    if cmd == 'plot':
        ppfs = sys.argv[3:]
        plot(ppfs)

