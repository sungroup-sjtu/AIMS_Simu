import os
import shutil

from .gmx import GmxSimulation
from ...unit import Unit
from ...utils import check_convergence


class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.requirement = []

    def build(self, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        self.export(minimize=minimize)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, nproc=1, jobname=None):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        commands = []
        self.gmx.prepare_mdp_from_template('t_npt.mdp', T=T, P=P / Unit.bar, nsteps=int(1E6))
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out=self.procedure, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name=self.procedure, nprocs=nproc, get_cmd=True)
        commands.append(cmd)
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self, dirs=None):
        if dirs is None:
            dirs = ['.']
        import numpy as np
        import pandas as pd
        import panedr

        # TODO check convergence, utilizing previous cycles
        temp_series = pd.Series()
        press_series = pd.Series()
        pe_series = pd.Series()
        density_series = pd.Series()
        for dir in dirs:
            df = panedr.edr_to_df(os.path.join(dir, self.procedure + '.edr'))
            temp_series = temp_series.append(df.Temperature)
            press_series = press_series.append(df.Pressure)
            pe_series = pe_series.append(df.Potential)
            density_series = density_series.append(df.Density)

        converged, when = check_convergence(density_series)

        return converged, {
            'temperature': np.mean(temp_series[when:]),
            'pressure': np.mean(press_series[when:]),
            'potential': np.mean(pe_series[when:]),
            'density': np.mean(density_series[when:])
        }
