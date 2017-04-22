import os
import shutil

from .npt import Npt
from ...utils import check_convergence


class NptHvap(Npt):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt-hvap'
        self.requirement = []

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, nproc=1, jobname=None):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        top_hvap = 'topol-hvap.top'
        self.gmx.generate_top_for_hvap(top, top_hvap)

        commands = []
        self.gmx.prepare_mdp_from_template('t_npt.mdp', T=T, P=0, nsteps=int(1E6), nstxtcout=1000)
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out='npt.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='npt', nprocs=nproc, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.grompp(gro=gro, top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='hvap', nprocs=nproc, rerun='npt.xtc', get_cmd=True)
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
        density_series = pd.Series()
        inter_series = pd.Series()
        for dir in dirs:
            df = panedr.edr_to_df(os.path.join(dir, 'npt.edr'))
            temp_series = temp_series.append(df.Temperature)
            density_series = density_series.append(df.Density)

            df = panedr.edr_to_df(os.path.join(dir, 'hvap.edr'))
            lj_series = inter_series.append(df['LJ (SR)'])
            tail_series = inter_series.append(df['Disper. corr.'])
            coul_series = inter_series.append(df['Coulomb (SR)'])
            pme_series = inter_series.append(df['Coul. recip.'])
            inter_series = lj_series + tail_series + coul_series + pme_series

        converged, when = check_convergence(inter_series)

        return converged, {
            'temperature': np.mean(temp_series[when:]),
            'density': np.mean(density_series[when:]),
            'inter': np.mean(inter_series[when:])
        }
