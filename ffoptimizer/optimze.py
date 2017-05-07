import sys
import os
import math
from collections import OrderedDict
from typing import Dict

import pybel
import shutil

sys.path.append('..')

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd
from mstools.simulation.gmx import GmxSimulation
from config import Config
from app import jobmanager

kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
          'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
simulation = GmxSimulation(**kwargs)


def f1d5p(f, h, x0):
    fm2, fm1, f1, f2 = [f(x0 + i * h) for i in [-2, -1, 1, 2]]
    fp = (-f2 + 8 * f1 - 8 * fm1 + fm2) / 12 / h
    return fp


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
        self.einter = None
        self.wEinter = None

    def __repr__(self):
        return self.smiles

    @property
    def dir(self):
        return os.path.join(Config.WORK_DIR, '%i-%i-%i' % (self.id, self.T, self.P))

    def n_mol(self, n_atoms=3000):
        py_mol = pybel.readstring('smi', self.smiles)
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
        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=1000, nsteps=int(4E5), nstxtcout=0)
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
                                                 nstxout=100, nstvout=100, nstxtcout=0)
        cmd = simulation.gmx.grompp(gro='eq.gro', tpr_out='nvt.tpr', get_cmd=True)
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

    def rsq(self):
        pressure = simulation.gmx.get_property('nvt.edr', 'Pressure', begin=20)
        ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=20)
        return (pressure - self.P) ** 2 * self.wDensity ** 2 + (ei - self.einter) ** 2 * self.wEinter

    def diff(self, ppf_file, para_k, para_v):
        P_list = []
        EI_list = []

        if para_k.endswith('r0'):
            h = 0.05
        else:
            h = 0.01

        for i in [-2, -1, 1, 2]:
            new_v = para_v + i * h
            ppf = PPF(ppf_file)
            ppf.set_lj_para({para_k: new_v})
            ppf.write('rerun.ppf')
            simulation.export(ppf='rerun.ppf', top_out='diff.top', minimize=False)
            nprocs = simulation.jobmanager.nprocs

            simulation.gmx.grompp(top='diff.top', tpr_out='diff.tpr')
            simulation.gmx.mdrun(name='diff', nprocs=nprocs, rerun='nvt.trr')

            top_hvap = 'diff-hvap.top'
            simulation.gmx.generate_top_for_hvap('diff.top', top_hvap)
            simulation.gmx.grompp(top=top_hvap, tpr_out='diff-hvap.tpr')
            simulation.gmx.mdrun(name='diff-hvap', nprocs=nprocs, rerun='nvt.trr')

            pressure = simulation.gmx.get_property('diff.edr', 'Pressure', begin=20)
            ei = simulation.gmx.get_property('diff-hvap.edr', 'Potential', begin=20)
            P_list.append(pressure)
            EI_list.append(ei)

        diff_P = (P_list[0] - 8 * P_list[1] + 8 * P_list[2] - P_list[3]) / 12 / h
        diff_EI = (EI_list[0] - 8 * EI_list[1] + 8 * EI_list[2] - EI_list[3]) / 12 / h
        return diff_P, diff_EI


class PPF():
    def __init__(self, ppf):
        with open(ppf) as f:
            self.terms = f.read().splitlines()

    @property
    def adj_lj_terms(self) -> OrderedDict:
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
            if r0.endswith('*'):
                terms['%s-r0' % a_type] = float(r0)
            if e0.endswith('*'):
                terms['%s-e0' % a_type] = float(e0)
        return terms

    def set_lj_para(self, new_paras: Dict):
        terms = [''] * len(self.terms)
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
            e0_key = '%s-e0' % a_type
            if e0_key in new_paras.keys():
                e0 = new_paras[e0_key]

            new_term = 'N12_6: %s: %s, %s:' % (a_type, str(r0), str(e0))
            terms[i] = new_term

        self.terms = terms

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            for term in self.terms:
                f.write(term + '\n')


def read_data(filename):
    targets = []
    with open(filename) as f:
        lines = f.read().splitlines()
    n = 0
    smiles_list = []
    for line in lines:
        if line.startswith('#') or line.strip() == '':
            continue
        words = line.split()

        smiles = words[0]
        if smiles not in smiles_list:
            smiles_list.append(smiles)
            n += 1

        target = Target()
        target.id = n
        target.smiles = words[0]
        target.T = int(words[1])
        target.P = int(words[2])
        target.density = float(words[3])
        target.wDensity = float(words[4])
        target.hvap = float(words[5])
        target.wHvap = float(words[6])
        target.convert_einter()
        targets.append(target)
    return targets


def build(datafile):
    targets = read_data(datafile)
    for p in targets:
        print(p.dir)
        cd_or_create_and_cd(p.dir)
        if not os.path.exists('eq.gro'):
            p.build()


def run_nvt(datafile):
    targets = read_data(datafile)
    for p in targets:
        os.chdir(p.dir)
        if not os.path.exists('eq.gro'):
            raise Exception('Should EQ first')
        p.run()


def optimize(datafile, ppf_file):
    targets = read_data(datafile)
    ppf = PPF(ppf_file)
    adj_lj_paras = ppf.adj_lj_terms

    rsq = 0
    directive = []
    for p in targets:
        directive.append([])
        os.chdir(p.dir)
        if not os.path.exists('nvt.gro'):
            raise Exception('Should run NVT first')

        rsq += p.rsq()
        for k, v in adj_lj_paras.items():
            directive.append(p.diff(ppf_file, k, v))

    print(rsq)
    print(directive)


if __name__ == '__main__':
    datafile = sys.argv[1]
    cmd = sys.argv[2]
    if cmd == 'build':
        build(datafile)
    if cmd == 'run':
        run_nvt(datafile)
    if cmd == 'optimize':
        ppf_file = os.path.abspath(sys.argv[3])
        optimize(datafile, ppf_file)
