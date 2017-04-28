import os
import shutil

from .gmx import GmxSimulation
from ...analyzer.series import is_converged
from ...unit import Unit


class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.requirement = []
        self.logs = ['npt.log', 'hvap.log', 'cp.log']

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
        # NPT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_npt_sd.mdp', mdp_out='grompp-eq.mdp', T=T, P=P / Unit.bar,
                                           nsteps=int(5E5))
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro, top=top, tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nproc, get_cmd=True)
        commands.append(cmd)

        # NPT production with Velocity Rescaling thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P / Unit.bar,
                                           nsteps=int(1E6), nstxtcout=1000, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-npt.mdp', gro='eq.gro', top=top, tpr_out='npt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='npt', nprocs=nproc, get_cmd=True)
        commands.append(cmd)

        # Enthalpy of vaporization
        top_hvap = 'topol-hvap.top'
        self.gmx.generate_top_for_hvap(top, top_hvap)
        cmd = self.gmx.grompp(mdp='grompp-npt.mdp', gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='hvap', nprocs=nproc, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        # Heat capacity using 2-Phase Thermodynamics
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-cp.mdp', T=T, P=P / Unit.bar,
                                           nsteps=int(2E4), nstvout=4, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-cp.mdp', gro='npt.gro', top=top, tpr_out='cp.tpr',
                              cpt='npt.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='cp', nprocs=nproc, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.dos(trr='cp.trr', tpr='cp.tpr', T=T, get_cmd=True)
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
        inter_series = pd.Series()
        for dir in dirs:
            df = panedr.edr_to_df(os.path.join(dir, 'npt.edr'))
            temp_series = temp_series.append(df.Temperature)
            press_series = press_series.append(df.Pressure)
            pe_series = pe_series.append(df.Potential)
            density_series = density_series.append(df.Density)

            df = panedr.edr_to_df(os.path.join(dir, 'hvap.edr'))
            inter_series = inter_series.append(df.Potential)

        converged, when = is_converged(density_series)
        if converged:
            Cp = 0
            try:
                with open('dos.log') as f:
                    lines = f.readlines()
                for line in lines:
                    if line.startswith('Heat capacity'):
                        Cp = float(line.split()[2])
            except:
                pass

            return converged, {'temperature': np.mean(temp_series[when:]),
                               'pressure': np.mean(press_series[when:]),
                               'potential': np.mean(pe_series[when:]),
                               'density': np.mean(density_series[when:]),
                               'inter': np.mean(inter_series[when:]),
                               'cp': Cp
                               }
        else:
            return converged, None
