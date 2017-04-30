import os, shutil
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

    def build(self, minimize=False):
        pass

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, jobname=None,
                prior_job_dir=None, prior_job_result: Dict = None):
        # Copy topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(model_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(model_dir, f), '.')

        # Slice structures from prior NPT trajectory, named conf0.gro, conf1.gro ...
        trr = os.path.join(prior_job_dir, 'npt.trr')
        tpr = os.path.join(prior_job_dir, 'npt.tpr')
        self.gmx.slice_gro_from_traj(trr, tpr, 'conf.gro', 60, 100, 10)

        # Scale gro box for NVT simulation
        box = prior_job_result['box']
        for i in range(5):
            gro_nvt = 'conf%i.gro' % i
            self.gmx.scale_box(gro_nvt, gro_nvt, box)

        nprocs = self.jobmanager.nprocs
        commands = []
        # Heat capacity using 2-Phase Thermodynamics
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-cv.mdp', T=T, P=P / Unit.bar,
                                           nsteps=int(4E4), nstvout=4, restart=True)
        for i in range(5):
            gro_nvt = 'conf%i.gro' % i
            name = 'cv%i' % i
            tpr = name + '.tpr'
            cmd = self.gmx.grompp(mdp='grompp-cv.mdp', gro=gro_nvt, top=top, tpr_out=tpr, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.dos(trr=name + '.trr', tpr=tpr, T=T, log_out='dos%i.log' % i, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self, dirs=None):
        cv = 0
        for i in range(5):
            with open('dos%i.log' % i) as f_dos:
                lines = f_dos.readlines()
            for line in lines:
                if line.startswith('Heat capacity'):
                    cv += float(line.split()[2])
            raise Exception('Heat capacity not found')
        cv /= 5
        return True, {'cv': cv}
