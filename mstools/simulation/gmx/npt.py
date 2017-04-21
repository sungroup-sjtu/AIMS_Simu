import os
import shutil

from .gmx import GmxSimulation
from ...unit import Unit
from ...utils import check_converged


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

    def analyze(self):
        # TODO check convergence, utilizing previous cycles
        import numpy as np
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        temperature_series = df.Temperature[500:]
        pressure_series = df.Pressure[500:]
        potential_series = df.Potential[500:]
        density_series = df.Density[500:]

        converged = check_converged(potential_series)
        if converged:
            converged = check_converged(density_series)

        return converged, {
            'temperature': np.mean(temperature_series),
            'pressure': np.mean(pressure_series),
            'potential': np.mean(potential_series),
            'density': np.mean(density_series)
        }
