import os
import shutil

from .gmx import GmxSimulation


class NvtSlab(GmxSimulation):
    requirement = []

    def build(self, ppf=None, minimize=False):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         size=[self.length, self.length, self.length * 5])
        self.export(ppf=ppf, minimize=minimize)

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=None, P=None, jobname=None, **kwargs):
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        nprocs = self.jobmanager.nprocs
        commands = []
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', T=T, nsteps=int(1E6))
        cmd = self.gmx.grompp(gro=gro, top=top, tpr_out=self.procedure, get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name=self.procedure, nprocs=nprocs, get_cmd=True)
        commands.append(cmd)
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)

    def analyze(self):
        import panedr
        df = panedr.edr_to_df(self.procedure + '.edr')
        surface_tension = df['#Surf*SurfTen']
        return {
            'surface_tension': surface_tension
        }
