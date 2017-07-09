import os, shutil
import math
from typing import Dict
from .gmx import GmxSimulation
from ...unit import Unit


class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.requirement = []
        self.logs = []
        for i in range(5):
            self.logs.append('cv%i.log' % i)
            self.logs.append('dos%i.log' % i)

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, gro='conf.gro', top='topol.top', T=None, P=None, jobname=None,
                prior_job_dir=None, prior_job_result: Dict = None, nst_eq=int(1E5), nst_run=int(5E4), **kwargs):
        # Copy topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')

        # Slice structures from prior NPT trajectory, named conf0.gro, conf1.gro ...
        trr = os.path.join(prior_job_dir, 'npt.trr')
        tpr = os.path.join(prior_job_dir, 'npt.tpr')

        if prior_job_result is not None:
            simulation_length = prior_job_result['simulation_length']
            converged_from = prior_job_result['converged_from']
        else:
            simulation_length = self.gmx.get_length_of_traj(trr)
            converged_from = 0
        dt = math.floor((simulation_length - converged_from) / 40) * 10  # the dt should be 10 ps at least
        begin = simulation_length - 4 * dt
        self.gmx.slice_gro_from_traj(trr, tpr, 'conf.gro', begin, simulation_length, dt)

        # Scale gro box for NVT simulation
        if prior_job_result is not None:
            box = prior_job_result['box']
        else:
            box = self.gmx.get_box(os.path.join(prior_job_dir, 'npt.edr'), converged_from)
        for i in range(5):
            gro_nvt = 'conf%i.gro' % i
            self.gmx.scale_box(gro_nvt, gro_nvt, box)

        nprocs = self.jobmanager.nprocs
        commands = []

        # Heat capacity using 2-Phase Thermodynamics
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-eq.mdp', T=T,
                                           nsteps=nst_eq, tcoupl='nose-hoover', restart=True)
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-nvt.mdp', T=T,
                                           nsteps=nst_run, nstvout=4, nstenergy=1, tcoupl='nose-hoover', restart=True)
        for i in range(5):
            gro_conf = 'conf%i.gro' % i
            gro_eq = 'eq%i.gro' % i
            name_eq = 'eq%i' % i
            name_nvt = 'nvt%i' % i
            tpr_eq = name_eq + '.tpr'
            tpr_nvt = name_nvt + '.tpr'
            cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro_conf, top=top, tpr_out=tpr_eq, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_eq, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro=gro_eq, top=top, tpr_out=tpr_nvt,
                                  cpt=name_eq + '.cpt', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_nvt, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            cmd = self.gmx.dos(trr=name_nvt + '.trr', tpr=tpr, T=T, log_out='dos%i.log' % i, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, dirs=None):
        cv_list = []
        for i in range(5):
            with open('dos%i.log' % i) as f_dos:
                lines = f_dos.readlines()
            for line in lines:
                if line.startswith('Heat capacity'):
                    cv_list.append(float(line.split()[2]))
                    break
            else:
                raise Exception('Heat capacity not found')

        import numpy as np
        return {'cv': [np.mean(cv_list), np.std(cv_list, ddof=1) / math.sqrt(5)]}
