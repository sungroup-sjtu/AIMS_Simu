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

sys.path.append('..')

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd
from mstools.simulation.gmx import GmxSimulation
from config import Config
from app import jobmanager

jobmanager.queue = 'cpu'
jobmanager.nprocs = 8
kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
          'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
simulation = GmxSimulation(**kwargs)


class Target():
    def __int__(self):
        self.id = 0
        self.smiles = None
        self.T = None
        self.P = None
        self.density = None
        self.wDensity = None
        self.hvap = None
        self.wHvap = None

    def __repr__(self):
        return '%s,%i' %(self.smiles, self.T)

    @property
    def dir(self):
        return os.path.join(Config.WORK_DIR, '%i-%i-%i' % (self.id, self.T, self.P))

    def n_mol(self, n_atoms=3000):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        return math.ceil(n_atoms / len(py_mol.atoms))

    def convert_einter(self):
        n_mol = self.n_mol()
        self.einter = (8.314 * self.T - self.hvap) * n_mol
        self.wEinter = self.wHvap / n_mol

    def build(self, n_atoms=3000):
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
        simulation.export(minimize=True)

        commands = []
        simulation.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', T=self.T, nsteps=int(5E5), nstxtcout=0)
        cmd = simulation.gmx.grompp(tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='eq', nprocs=simulation.jobmanager.nprocs, get_cmd=True)
        commands.append(cmd)

        simulation.jobmanager.generate_sh(os.getcwd(), commands, name='eq-%i' % self.id)
        simulation.run()

    def run(self):
        nprocs = simulation.jobmanager.nprocs
        commands = []

        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5),
                                                 nstxout=1000, nstvout=1000, nstxtcout=0,
                                                 restart=True)
        cmd = simulation.gmx.grompp(gro='eq.gro', tpr_out='nvt.tpr', cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        top_hvap = 'topol-hvap.top'
        simulation.gmx.generate_top_for_hvap('topol.top', top_hvap)
        cmd = simulation.gmx.grompp(gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='hvap', nprocs=nprocs, rerun='nvt.trr', get_cmd=True)
        commands.append(cmd)

        simulation.jobmanager.generate_sh(os.getcwd(), commands, name='nvt-%i' % self.id)
        simulation.run()

    def get_pres(self):
        os.chdir(self.dir)
        return simulation.gmx.get_property('nvt.edr', 'Pressure', begin=100)

    def get_hvap(self):
        os.chdir(self.dir)
        ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=100)
        return 8.314 * self.T / 1000 - ei / self.n_mol

    def get_pres_hvap_from_paras(self, d:OrderedDict = None)->(float, float):
        paras = copy.copy(d)
        os.chdir(self.dir)
        ppf = PPF('TEAM_LS.ppf')
        if not ppf.set_lj_para(paras):
            return self.get_pres(), self.get_hvap()
        ppf.write('rerun.ppf')
        top = 'diff.top'
        simulation.dff.export_gmx('init.msd', 'rerun.ppf', top_out=top)
        nprocs = simulation.jobmanager.nprocs

        # Press
        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), nstxtcout=0)
        simulation.gmx.grompp(top=top, tpr_out='diff.tpr', silent=True)
        simulation.gmx.mdrun(name='diff', nprocs=nprocs, rerun='nvt.trr', silent=True)

        pres = simulation.gmx.get_property('diff.edr', 'Pressure', begin=100)

        # Hvap
        top_hvap = 'diff-hvap.top'
        simulation.gmx.generate_top_for_hvap(top, top_hvap)
        nprocs = simulation.jobmanager.nprocs

        simulation.gmx.grompp(top=top_hvap, tpr_out='diff-hvap.tpr', silent=True)
        simulation.gmx.mdrun(name='diff-hvap', nprocs=nprocs, rerun='nvt.trr', silent=True)

        ei = simulation.gmx.get_property('diff-hvap.edr', 'Potential', begin=100)
        hvap = 8.314 * self.T / 1000 - ei / self.n_mol
        return pres, hvap

    def get_dPres_dHvap_from_paras(self, d:OrderedDict)->([float], [float]):
        paras = copy.copy(d)
        dPres = []
        dHvap = []
        for k, v in paras.items():
            dP, dH = self.get_dPres_dHvap_for_para(paras, k)
            dPres.append(dP)
            dHvap.append(dH)
        return dPres, dHvap


    def get_dPres_dHvap_for_para(self, d:OrderedDict, k)->(float, float):
        paras = copy.copy(d)
        v = paras[k]
        pres_list = []
        hvap_list = []

        if k.endswith('r0'):
            h = 0.02
        else:
            h = 0.01

        for i in [-1, 1]:
            new_v = v + i * h
            paras[k] = new_v
            ppf = PPF('TEAM_LS.ppf')
            if not ppf.set_lj_para(paras):
                return 0, 0
            ppf.write('rerun.ppf')
            #print('    %s %10.5f %10.5f' %(para_k, v, new_v))
            top = 'diff%i.top' %i
            simulation.dff.export_gmx('init.msd', 'rerun.ppf', top_out=top)
            nprocs = simulation.jobmanager.nprocs

            # Pres
            simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), nstxtcout=0)
            simulation.gmx.grompp(top=top, tpr_out='diff%i.tpr' %i, silent=True)
            simulation.gmx.mdrun(name='diff%i' %i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            pres = simulation.gmx.get_property('diff%i.edr' %i, 'Pressure', begin=100)
            pres_list.append(pres)

            # HVap
            top_hvap = 'diff%i-hvap.top'
            simulation.gmx.generate_top_for_hvap(top, top_hvap)

            simulation.gmx.grompp(top=top_hvap, tpr_out='diff%i-hvap.tpr' %i, silent=True)
            simulation.gmx.mdrun(name='diff%i-hvap' %i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            ei = simulation.gmx.get_property('diff%i-hvap.edr' %i, 'Potential', begin=100)
            hvap = 8.314 * self.T / 1000 - ei / self.n_mol
            hvap_list.append(hvap)

        dPres = (pres_list[1] - pres_list[0]) / 2 / h
        dHvap = (hvap_list[1] - hvap_list[0]) / 2 / h
        return dPres, dHvap


class PPF():
    def __init__(self, ppf):
        with open(ppf) as f:
            self.terms = f.read().splitlines()

    @property
    def adj_lj_paras(self) -> OrderedDict:
        terms = OrderedDict()
        for term in self.terms:
            if not term.startswith('N12_6'):
                continue

            words = term.split(':')
            words = [w.strip() for w in words]
            a_type = words[1]
            paras = words[2]

            words = paras.split(',')
            words = [w.strip() for w in words]
            r0 = words[0]
            e0 = words[1]
            if not r0.endswith('*'):
                terms['%s-r0' % a_type] = float(r0)
            if not e0.endswith('*'):
                terms['%s-e0' % a_type] = float(e0)
        return terms

    def set_lj_para(self, new_paras: Dict):
        terms = [''] * len(self.terms)
        replaced = False
        for i, term in enumerate(self.terms):
            if not term.startswith('N12_6'):
                terms[i] = term
                continue

            words = term.split(':')
            words = [w.strip() for w in words]
            a_type = words[1]
            paras = words[2]
            words = paras.split(',')
            words = [w.strip() for w in words]
            r0 = words[0]
            e0 = words[1]

            r0_key = '%s-r0' % a_type
            if r0_key in new_paras.keys():
                r0 = new_paras[r0_key]
                replaced = True
            e0_key = '%s-e0' % a_type
            if e0_key in new_paras.keys():
                e0 = new_paras[e0_key]
                replaced = True

            new_term = 'N12_6: %s: %s, %s:' % (a_type, str(r0), str(e0))
            terms[i] = new_term

        self.terms = terms
        return replaced

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            for term in self.terms:
                f.write(term + '\n')


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


def build():
    for p in targets:
        print(p.dir)
        cd_or_create_and_cd(p.dir)
        if not os.path.exists('eq.gro'):
            p.build()


def run_nvt():
    for p in targets:
        os.chdir(p.dir)
        if not os.path.exists('eq.gro'):
            raise Exception('Should EQ first')
        p.run()


def optimize_():
    pres_list = []
    hvap_list = []
    rsq_total = 0
    dPres = []
    dHvap = []
    for p in targets:
        print(p)
        dPres.append([])
        dHvap.append([])
        os.chdir(p.dir)
        if not os.path.exists('nvt.gro'):
            raise Exception('Should run NVT first')

        pres, hvap, rsq = p.get_result()
        pres_list.append(pres)
        hvap_list.append(hvap)
        rsq_total += rsq
        for k, v in adj_lj_paras.items():
            print('  ' + k)
            dP, dH = p.diff(k, v)
            dPres[-1].append(dP)
            dHvap[-1].append(dH)

    print('\nRSQ = %.1f' %rsq_total)
    print('Parameter', end='')
    for j, target in enumerate(targets):
        print(' ', target, end='')
    print('')
    for (i, k) in enumerate(adj_lj_paras.keys()):
        print(k, end='')
        for j, target in enumerate(targets):
            print(' %.1f %.1f' %(dPres[j][i], dHvap[j][i]), end='')
        print('')


def func(x):
    paras = OrderedDict()
    for i, k in enumerate(adj_lj_paras.keys()):
        paras[k] = x[i]

    f = []
    df = []
    for target in targets:
        pres, hvap = target.get_pres_hvap_from_paras(paras)
        f.append(pres - target.P)
        f.append((hvap - target.hvap) * target.wHvap)
        dPres, dHvap = target.get_dPres_dHvap_from_paras(paras)
        df.append(dPres)
        df.append(dHvap)
    print('CURRENT', x)
    print('RESIDUE:', list(map(lambda x: round(x, 1), f)))
    print('RSQ:', round(np.sum(list(map(lambda x: x**2, f))), 1))
    print('DIRECTIVE:')
    for l in df:
        print(list(map(lambda x: round(x, 1), l)))

    return f, df

def print_result(x, f):
    print(x, f)

def optimize():
    print(adj_lj_paras)
    x0 = [v for k, v in adj_lj_paras.items()]
    sol = root(func, x0, jac=True, method='lm', options={'maxiter': 100})
    print(sol.x)
    print(sol.success)
    print(sol.message)

if __name__ == '__main__':
    datafile = sys.argv[1]
    cmd = sys.argv[2]
    global targets
    targets = read_data(datafile)

    if cmd == 'build':
        build()
    if cmd == 'run':
        run_nvt()
    if cmd == 'optimize':
        ppf = PPF(sys.argv[3])
        global adj_lj_paras
        adj_lj_paras = ppf.adj_lj_paras

        optimize()

