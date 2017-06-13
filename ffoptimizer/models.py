import copy
import math
import os
import shutil
from collections import OrderedDict
from functools import partial

import pybel
import panedr
from sqlalchemy import Column, Integer, Text, Float, String

NotNullColumn = partial(Column, nullable=False)

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd, n_diff_lines
from mstools.simulation.gmx import GmxSimulation, Npt
from config import Config
from app import jobmanager

kwargs = {'packmol_bin': Config.PACKMOL_BIN, 'dff_root': Config.DFF_ROOT,
          'gmx_bin': Config.GMX_BIN, 'jobmanager': jobmanager}
simulation = GmxSimulation(**kwargs)

from ffoptimizer.ppf import PPF

from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()
metadata = Base.metadata


class Target(Base):
    __tablename__ = 'target'
    id = NotNullColumn(Integer, primary_key=True)
    name = NotNullColumn(String(200))
    smiles = NotNullColumn(Text)
    n_mol = Column(Integer, nullable=True)
    T = NotNullColumn(Integer)
    P = NotNullColumn(Integer)
    density = NotNullColumn(Float)
    hvap = NotNullColumn(Float)
    wDensity = NotNullColumn(Float)
    wHvap = NotNullColumn(Float)
    iteration = NotNullColumn(Integer, default=0)

    def __repr__(self):
        return '<Target: %s %s %i>' % (self.name, self.smiles, self.T)

    def calc_n_mol(self, n_atoms=3000, n_mol=0):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        self.n_mol = math.ceil(n_atoms / len(py_mol.atoms))
        if self.n_mol < n_mol:
            self.n_mol = n_mol

    @property
    def dir_base_npt(self):
        return os.path.join(Config.WORK_DIR, 'NPT-%s' % (self.name))

    @property
    def dir_child(self):
        return '%i-%i-%i' % (self.T, self.P, self.iteration)

    def run_npt(self, ppf_file=None, diff: OrderedDict = None):
        cd_or_create_and_cd(self.dir_base_npt)

        if not os.path.exists('init.msd'):
            pdb = 'mol.pdb'
            mol2 = 'mol.mol2'
            py_mol = create_mol_from_smiles(self.smiles, pdb_out=pdb, mol2_out=mol2)
            mass = py_mol.molwt * self.n_mol
            length = (10 / 6.022 * mass / (self.density - 0.1)) ** (1 / 3)  # assume cubic box

            print('Build coordinates using Packmol: %s molecules ...' % self.n_mol)
            simulation.packmol.build_box([pdb], [self.n_mol], 'init.pdb', length=length - 2, tolerance=1.7, silent=True)

            print('Create box using DFF ...')
            simulation.dff.build_box_after_packmol([mol2], [self.n_mol], 'init.msd', mol_corr='init.pdb', length=length)

        cd_or_create_and_cd(os.path.basename('%s' % ppf_file)[:-4])

        simulation.msd = '../init.msd'
        simulation.export(ppf=ppf_file, minimize=True)

        cd_or_create_and_cd(self.dir_child)

        npt = Npt(**kwargs)
        commands = npt.prepare(model_dir='..', T=self.T, P=self.P, jobname='NPT-%s-%i' % (self.name, self.T),
                               dt=0.002, nst_eq=int(2E5), nst_run=int(2E5), nst_trr=250, nst_xtc=250)

        if diff is not None:
            commands.append('export GMX_MAXCONSTRWARN=-1')
            for k in diff.keys():
                paras = copy.copy(diff)
                if k.endswith('r0'):
                    delta = 0.01
                elif k.endswith('e0'):
                    delta = 0.001
                elif k.endswith('bi'):
                    delta = 0.005
                else:
                    raise Exception('Unknown parameter: ' + k)

                msd_out = 'diff-%s.msd' % k
                ppf_out = 'diff-%s.ppf' % k
                top_out = 'diff-%s.top' % k
                top_out_hvap = 'diff-%s-hvap.top' % k

                paras[k] += delta
                ppf = PPF(ppf_file)
                ppf.set_lj_para(paras)
                ppf.write(ppf_out)

                shutil.copy('../../init.msd', msd_out)
                simulation.dff.set_charge([msd_out], ppf_out)
                simulation.dff.export_gmx(msd_out, ppf_out, gro_out='diff.gro', top_out=top_out)

                nprocs = simulation.jobmanager.nprocs

                cmd = simulation.gmx.grompp(mdp='grompp-npt.mdp', top=top_out, tpr_out='diff-%s.tpr' % k, get_cmd=True)
                commands.append(cmd)
                cmd = simulation.gmx.mdrun(name='diff-%s' % k, nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                commands.append(cmd)

                simulation.gmx.generate_top_for_hvap(top_out, top_out_hvap)

                cmd = simulation.gmx.grompp(mdp='grompp-npt.mdp', top=top_out_hvap, tpr_out='diff-%s-hvap.tpr' % k,
                                            get_cmd=True)
                commands.append(cmd)
                cmd = simulation.gmx.mdrun(name='diff-%s-hvap' % k, nprocs=nprocs, rerun='npt.trr', get_cmd=True)
                commands.append(cmd)

        commands.append('touch _finished_')
        npt.jobmanager.generate_sh(os.getcwd(), commands, name='NPT-%s-%i' % (self.name, self.T))
        npt.run()

    def get_npt_result(self, subdir) -> (float, float):
        os.chdir(self.dir_base_npt)
        os.chdir(subdir)
        os.chdir(self.dir_child)
        print(os.getcwd())
        density = simulation.gmx.get_property('npt.edr', 'Density')
        density /= 1000
        ei = simulation.gmx.get_property('hvap.edr', 'Potential')
        hvap = 8.314 * self.T / 1000 - ei / self.n_mol
        return density, hvap

    def get_dDens_dHvap_from_paras(self, ppf_file, d: OrderedDict):
        paras = copy.copy(d)
        dDens = []
        dHvap = []
        for k, v in paras.items():
            dD, dH = self.get_dDens_dHvap_from_para(ppf_file, paras, k)
            dDens.append(dD)
            dHvap.append(dH)
        return dDens, dHvap

    def get_dDens_dHvap_from_para(self, ppf_file, d: OrderedDict, k) -> (float, float):
        os.chdir(self.dir_base_npt)
        subdir = os.path.basename(ppf_file)[:-4]
        os.chdir(subdir)
        os.chdir(self.dir_child)

        if k.endswith('r0'):
            delta = 0.01
        elif k.endswith('e0'):
            delta = 0.001
        elif k.endswith('bi'):
            delta = 0.005
        else:
            raise Exception('Unknown parameter: ' + k)

        # energy and Hvap after diff
        df = panedr.edr_to_df('diff-%s.edr' % k)
        pene_series_diff = df.Potential

        df = panedr.edr_to_df('diff-%s-hvap.edr' % k)
        eint_series_diff = df.Potential
        hvap_series_diff = 8.314 * self.T / 1000 - eint_series_diff / self.n_mol

        # density, energy and Hvap
        df = panedr.edr_to_df('npt.edr')
        dens_series = df.Density.loc[pene_series_diff.index]
        pene_series = df.Potential.loc[pene_series_diff.index]

        df = panedr.edr_to_df('hvap.edr')
        hvap_series = df.Potential.loc[pene_series_diff.index]

        # calculate the derivative
        dPene_series = (pene_series_diff - pene_series) / delta
        dHvap_series = (hvap_series_diff - hvap_series) / delta

        dens_dPene = dens_series * dPene_series
        hvap_dPene = hvap_series * dPene_series

        dDdP = -1 / 8.314 / self.T * (dens_dPene.mean() - dens_series.mean() * dPene_series.mean())
        dHdP = dHvap_series.mean() - 1 / 8.314 / self.T * (hvap_dPene.mean() - hvap_series.mean() * dPene_series.mean())
        return dDdP, dHdP

    def npt_finished(self, ppf_file) -> bool:
        subdir = os.path.basename(ppf_file)[:-4]
        log_finished = os.path.join(self.dir_base_npt, subdir, self.dir_child, '_finished_')
        if os.path.exists(log_finished):
            return True

        return False

    def npt_started(self, ppf_file) -> bool:
        subdir = os.path.basename(ppf_file)[:-4]
        sh_job = os.path.join(self.dir_base_npt, subdir, self.dir_child, jobmanager.sh)
        if os.path.exists(sh_job):
            return True

        return False

    def clear_npt_result(self, ppf_file):
        subdir = os.path.basename(ppf_file)[:-4]
        dir = os.path.join(self.dir_base_npt, subdir, self.dir_child)
        log_finished = os.path.join(dir, '_finished_')
        sh_job = os.path.join(dir, jobmanager.sh)
        try:
            os.remove(log_finished)
            shutil.move(sh_job, sh_job + '.bak')
        except:
            pass


class Result(Base):
    __tablename__ = 'result'
    id = NotNullColumn(Integer, primary_key=True)
    ppf = NotNullColumn(String(200))
    parameter = Column(Text, nullable=True)
    residual = Column(Text, nullable=True)
    jacobian = Column(Text, nullable=True)
