import os, shutil
from .gmx import GmxSimulation
from ...unit import Unit
import math


class Npt(GmxSimulation):
    requirement = []

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
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        density = df.Density[500:]
        potential = df.Potential[500:]
        return {
            'density': density,
            'potential': potential
        }
