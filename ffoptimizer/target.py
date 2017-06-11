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

from mstools.utils import create_mol_from_smiles, cd_or_create_and_cd
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
    cycle = NotNullColumn(Integer, default=0)

    def __repr__(self):
        return '<Target: %s %s %i>' % (self.name, self.smiles, self.T)

    @property
    def dir_nvt(self):
        base_dir = os.path.join(Config.WORK_DIR, 'NVT-%s' % self.name, str(self.T))
        if self.cycle == 0:
            return base_dir
        else:
            return os.path.join(base_dir, str(self.cycle))

    def calc_n_mol(self, n_atoms=3000, n_mol=100):
        py_mol = pybel.readstring('smi', self.smiles)
        py_mol.addh()
        self.n_mol = math.ceil(n_atoms / len(py_mol.atoms))
        if self.n_mol < n_mol:
            self.n_mol = n_mol

    def build_nvt(self, ppf=None):
        cd_or_create_and_cd(self.dir_nvt)
        if self.cycle == 0:
            pdb = 'mol.pdb'
            mol2 = 'mol.mol2'
            py_mol = create_mol_from_smiles(self.smiles, pdb_out=pdb, mol2_out=mol2)
            mass = py_mol.molwt * self.n_mol
            length = (10 / 6.022 * mass / self.density) ** (1 / 3)  # assume cubic box

            print('Build coordinates using Packmol: %s molecules ...' % self.n_mol)
            simulation.packmol.build_box([pdb], [self.n_mol], 'init.pdb', length=length - 2, tolerance=1.7, silent=True)

            print('Create box using DFF ...')
            simulation.dff.build_box_after_packmol([mol2], [self.n_mol], 'init.msd', mol_corr='init.pdb', length=length)
            if ppf is not None:
                shutil.copy(ppf, 'TEAM_LS.ppf')
                simulation.export(ppf='TEAM_LS.ppf', minimize=True)
            else:
                simulation.export(minimize=True)

        else:
            shutil.copy('../init.msd', 'init.msd')
            shutil.copy(ppf, 'TEAM_LS.ppf')
            simulation.export(ppf='TEAM_LS.ppf', minimize=False)
            shutil.copy('../nvt.gro', 'conf.gro')

        nprocs = simulation.jobmanager.nprocs
        commands = []
        simulation.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-eq.mdp', T=self.T,
                                                 nsteps=int(3E5), dt=0.001, nstxtcout=0)
        cmd = simulation.gmx.grompp(mdp='grompp-eq.mdp', tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = simulation.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-nvt.mdp', T=self.T,
                                                 nsteps=int(2E5), dt=0.002, restart=True,
                                                 nstxout=500, nstvout=500, nstxtcout=0)
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

        simulation.jobmanager.generate_sh(os.getcwd(), commands, name='%s-%i' % (self.name, self.T))
        simulation.run()

    def get_pres_hvap_from_paras(self, ppf_file=None, d: OrderedDict = None) -> (float, float):
        paras = copy.copy(d)
        os.chdir(self.dir_nvt)
        if ppf_file is not None:
            ppf = PPF(ppf_file)
        else:
            ppf = PPF('TEAM_LS.ppf')
        if not ppf.set_lj_para(paras):
            pres = simulation.gmx.get_property('nvt.edr', 'Pressure')
            ei = simulation.gmx.get_property('hvap.edr', 'Potential')
            hvap = 8.314 * self.T / 1000 - ei / self.n_mol
            return pres, hvap

        ppf.write('rerun.ppf')
        top = 'diff.top'

        shutil.copy('init.msd', 'diff.msd')
        simulation.dff.set_charge(['diff.msd'], 'rerun.ppf')
        simulation.dff.export_gmx('diff.msd', 'rerun.ppf', top_out=top)
        nprocs = simulation.jobmanager.nprocs

        # Press
        simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), dt=0.002, nstxtcout=0)
        simulation.gmx.grompp(top=top, tpr_out='diff.tpr', silent=True)
        simulation.gmx.mdrun(name='diff', nprocs=nprocs, rerun='nvt.trr', silent=True)

        pres = simulation.gmx.get_property('diff.edr', 'Pressure')

        # Hvap
        top_hvap = 'diff-hvap.top'
        simulation.gmx.generate_top_for_hvap(top, top_hvap)
        nprocs = simulation.jobmanager.nprocs

        simulation.gmx.grompp(top=top_hvap, tpr_out='diff-hvap.tpr', silent=True)
        simulation.gmx.mdrun(name='diff-hvap', nprocs=nprocs, rerun='nvt.trr', silent=True)

        ei = simulation.gmx.get_property('diff-hvap.edr', 'Potential')
        hvap = 8.314 * self.T / 1000 - ei / self.n_mol
        return pres, hvap

    def get_dPres_dHvap_from_paras(self, ppf_file, d: OrderedDict) -> ([float], [float]):
        paras = copy.copy(d)
        dPres = []
        dHvap = []
        for k, v in paras.items():
            dP, dH = self.get_dPres_dHvap_for_para(ppf_file, paras, k)
            dPres.append(dP)
            dHvap.append(dH)
        return dPres, dHvap

    def get_dPres_dHvap_for_para(self, ppf_file, d: OrderedDict, k) -> (float, float):
        paras = copy.copy(d)
        os.chdir(self.dir_nvt)
        v = paras[k]
        pres_list = []
        hvap_list = []

        if k.endswith('r0'):
            delta = 0.02
        elif k.endswith('e0'):
            delta = 0.01
        elif k.endswith('bi'):
            delta = 0.02
        else:
            raise Exception('Unknown parameter: ' + k)

        for i in [-1, 1]:
            new_v = v + i * delta
            paras[k] = new_v
            if ppf_file is not None:
                ppf = PPF(ppf_file)
            else:
                ppf = PPF('TEAM_LS.ppf')
            if not ppf.set_lj_para(paras):
                return 0, 0
            ppf.write('rerun.ppf')
            # print('    %s %10.5f %10.5f' %(k, v, new_v))
            top = 'diff%i.top' % i

            shutil.copy('init.msd', 'diff.msd')
            simulation.dff.set_charge(['diff.msd'], 'rerun.ppf')
            simulation.dff.export_gmx('diff.msd', 'rerun.ppf', top_out=top)
            nprocs = simulation.jobmanager.nprocs

            # Pres
            simulation.gmx.prepare_mdp_from_template('t_nvt.mdp', T=self.T, nsteps=int(2E5), dt=0.002, nstxtcout=0)
            simulation.gmx.grompp(top=top, tpr_out='diff%i.tpr' % i, silent=True)
            simulation.gmx.mdrun(name='diff%i' % i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            pres = simulation.gmx.get_property('diff%i.edr' % i, 'Pressure')
            pres_list.append(pres)

            # HVap
            top_hvap = 'diff%i-hvap.top' % i
            simulation.gmx.generate_top_for_hvap(top, top_hvap)

            simulation.gmx.grompp(top=top_hvap, tpr_out='diff%i-hvap.tpr' % i, silent=True)
            simulation.gmx.mdrun(name='diff%i-hvap' % i, nprocs=nprocs, rerun='nvt.trr', silent=True)

            ei = simulation.gmx.get_property('diff%i-hvap.edr' % i, 'Potential')
            hvap = 8.314 * self.T / 1000 - ei / self.n_mol
            hvap_list.append(hvap)

        dPres = (pres_list[1] - pres_list[0]) / 2 / delta
        dHvap = (hvap_list[1] - hvap_list[0]) / 2 / delta
        return dPres, dHvap

    @property
    def dir_base_npt(self):
        return os.path.join(Config.WORK_DIR, 'NPT-%s' % (self.name))

    def run_npt(self, ppf=None):
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

        cd_or_create_and_cd(os.path.basename('%s' % ppf)[:-4])

        simulation.msd = '../init.msd'
        simulation.export(ppf=ppf, minimize=True)

        cd_or_create_and_cd('%i-%i' % (self.T, self.P))

        npt = Npt(**kwargs)
        npt.prepare(model_dir='..', T=self.T, P=self.P, jobname='NPT-%s-%i' % (self.name, self.T),
                    dt=0.002, nst_eq=int(2E5), nst_run=int(1E5), nst_trr=200)
        npt.run()

    def get_npt_result(self, subdir) -> (float, float):
        os.chdir(self.dir_base_npt)
        os.chdir(subdir)
        os.chdir('%i-%i' % (self.T, self.P))
        print(os.getcwd())
        density = simulation.gmx.get_property('npt.edr', 'Density', begin=250)
        density /= 1000
        ei = simulation.gmx.get_property('hvap.edr', 'Potential', begin=250)
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
        paras = copy.copy(d)

        os.chdir(self.dir_base_npt)
        subdir = os.path.basename(ppf_file)[:-4]
        os.chdir(subdir)
        os.chdir('%i-%i' % (self.T, self.P))

        v = paras[k]
        pene_series_list = []
        hvap_series_list = []

        if k.endswith('r0'):
            delta = 0.02
        elif k.endswith('e0'):
            delta = 0.01
        elif k.endswith('bi'):
            delta = 0.02
        else:
            raise Exception('Unknown parameter: ' + k)

        for i in [-1, 1]:
            new_v = v + i * delta
            paras[k] = new_v
            if ppf_file is not None:
                ppf = PPF(ppf_file)
            else:
                ppf = PPF('TEAM_LS.ppf')
            if not ppf.set_lj_para(paras):
                return 0, 0
            ppf.write('diff.ppf')
            # print('    %s %10.5f %10.5f' %(k, v, new_v))
            top = 'diff.top'

            shutil.copy('init.msd', 'diff.msd')
            simulation.dff.set_charge(['diff.msd'], 'diff.ppf')
            simulation.dff.export_gmx('diff.msd', 'diff.ppf', top_out=top)
            nprocs = simulation.jobmanager.nprocs

            # Pres
            simulation.gmx.grompp(mdp='grompp-npt.mdp', top=top, tpr_out='diff.tpr', silent=True)
            simulation.gmx.mdrun(name='diff', nprocs=nprocs, rerun='npt.trr', silent=True)

            df = panedr.edr_to_df('diff.edr')
            pene_series = df.Potential
            pene_series_list.append(pene_series)

            # HVap
            top_hvap = 'diff-hvap.top'
            simulation.gmx.generate_top_for_hvap(top, top_hvap)

            simulation.gmx.grompp(mdp='grompp-npt.mdp', top=top_hvap, tpr_out='diff-hvap.tpr', silent=True)
            simulation.gmx.mdrun(name='diff-hvap', nprocs=nprocs, rerun='npt.trr', silent=True)

            df = panedr.edr_to_df('diff.edr')
            eint_series = df.Potential
            hvap_series = 8.314 * self.T / 1000 - eint_series / self.n_mol
            hvap_series_list.append(hvap_series)

        dPene_series = (pene_series_list[1] - pene_series_list[0]) / 2 / delta
        dHvap_series = (hvap_series_list[1] - hvap_series_list[0]) / 2 / delta

        df = panedr.edr_to_df('npt.edr')
        dens_series = df.Density
        dens_series = dens_series.loc[dens_series.index in dPene_series]

        df = panedr.edr_to_df('hvap.edr')
        hvap_series = df.Potential
        hvap_series = hvap_series.loc[hvap_series.index in dPene_series]

        dens_dPene = dens_series * dPene_series
        hvap_dPene = hvap_series * dPene_series

        dDdP = -1 / 8.314 / self.T * (dens_dPene.mean() - dens_series.mean() * dPene_series.mean())
        dHdP = dHvap_series.mean() - 1 / 8.314 / self.T * (hvap_dPene.mean() - hvap_series.mean() * dPene_series.mean())
        return dDdP, dHdP
